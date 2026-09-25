"""The read endpoints: the gene record, and the catalog that names it."""

import contextlib
from typing import Optional
from urllib.error import HTTPError, URLError

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import RedirectResponse

from . import catalog, suggest, tracks
from .catalog import CatalogUnavailable
from .genbank_parser import GeneNotFound, extract_gene
from .record_cache import fetch as fetch_genbank_record
from .schemas import (
    CatalogPage,
    GeneResponse,
    ProteinDetail,
    SearchResponse,
    SuggestResponse,
    Track,
)
from .impact_explanations import NOT_COVERED, read_impact_explanations, stored_url

router = APIRouter()


@contextlib.contextmanager
def _catalog_errors():
    """Report an unreadable catalog as unavailable, never as an empty one.

    There is nothing to read through to here. `record_cache` can lose its
    database and still answer from NCBI; the catalog has no second source, and
    answering an unreadable catalog with an empty one would tell a client that
    a protein does not exist when the truth is that nobody could look. The
    client keeps its own copy, and that is what keeps search working meanwhile.
    """
    try:
        yield
    except CatalogUnavailable as exc:
        raise HTTPException(status_code=503, detail=str(exc))


@router.get("/gene/{id}/{gene}/impact-explanations", response_model=None)
def impact_explanations(id: str, gene: str):
    """A saved attribution payload; no Atlas credentials or inference at runtime.

    Where storage holds it, the client is redirected there and fetches the
    bytes itself: cftr's payload is 3.98 MB, and carrying it would put it
    through 512 MB of memory and a tenth of a CPU for no gain. The redirect
    keeps the endpoint's contract -- a client that follows it, as Dio does by
    default, still receives the JSON it always received.
    """
    found = stored_url(id, gene)
    if found is NOT_COVERED:
        raise HTTPException(
            status_code=404,
            detail="AVI explanations are not included for this gene.",
        )
    if found is not None:
        # 307 rather than 302: the method must survive the hop, and the target
        # is where the bytes are today rather than where they will always be.
        return RedirectResponse(found, status_code=307)
    return read_impact_explanations(id, gene)


def _upstream_message(exc: HTTPError) -> str:
    """Best-effort read of the error body NCBI sent with an HTTPError."""
    try:
        body = exc.read()
    except (OSError, ValueError, AttributeError):
        return exc.reason or str(exc)
    if isinstance(body, bytes):
        body = body.decode("utf-8", errors="replace")
    return body.strip() or (exc.reason or str(exc))


def _no_record(id: str) -> str:
    return "No GenBank record found for id {!r}.".format(id)


@contextlib.contextmanager
def _upstream_errors(id: str):
    """Translate urllib failures from an Entrez call into HTTPExceptions.

    A context manager rather than inline handling so the mapping is stated once
    and stays reusable as endpoints are added. The clauses stay
    most-specific-first because HTTPError subclasses URLError, which subclasses
    OSError -- reordering them would silently swallow the specific cases.

    Any upstream 400 is read as "not found": ``id`` is the only part of the
    request that can be invalid.
    """
    try:
        yield
    except HTTPError as exc:
        if exc.code == 400:
            raise HTTPException(status_code=404, detail=_no_record(id))
        raise HTTPException(
            status_code=502,
            detail="NCBI returned {}: {}".format(exc.code, _upstream_message(exc)),
        )
    except URLError as exc:
        raise HTTPException(
            status_code=502,
            detail="Could not reach NCBI: {}".format(exc.reason),
        )
    except OSError as exc:
        raise HTTPException(
            status_code=502,
            detail="Could not reach NCBI: {}".format(exc),
        )


@router.get("/gene/{id}/{gene}", response_model=GeneResponse)
def read_gene(id: str, gene: str) -> dict:
    """Fetch one GenBank record and return just ``gene`` from it, structured.

    Deliberately a sync ``def``: ``Entrez.efetch`` is blocking urllib, so
    FastAPI runs this in a threadpool. Making it ``async def`` would block the
    event loop for the whole round trip to NCBI.

    ``id`` and ``gene`` are both free parameters: the record to read and the
    gene to keep are the client's choice, not this service's.

    The fetch is cache-aware. A cache miss, or no database at all, reads
    through to NCBI exactly as before.
    """
    with _upstream_errors(id):
        try:
            record = fetch_genbank_record(id)
        except ValueError:
            # Not a GenBank record at all -- NCBI answers an unusable id with
            # plain text and an HTTP 200.
            raise HTTPException(status_code=404, detail=_no_record(id))

    try:
        return extract_gene(record, gene)
    except GeneNotFound:
        raise HTTPException(
            status_code=404,
            detail="No gene {!r} in record {!r}.".format(gene, id),
        )


@router.get("/catalog", response_model=CatalogPage)
def read_catalog(
    limit: int = Query(200, ge=1, le=500),
    cursor: Optional[str] = None,
) -> dict:
    """The whole catalog, a page at a time.

    The client hydrates its on-device copy from this and then searches that
    copy, so the page size is about one request rather than about a screen.
    """
    with _catalog_errors():
        proteins, next_cursor = catalog.page(limit=limit, cursor=cursor)
    return {"proteins": proteins, "next": next_cursor}


@router.get("/catalog/search", response_model=SearchResponse)
def search_catalog(
    q: str = "",
    limit: int = Query(20, ge=1, le=100),
) -> dict:
    """Catalog rows matching ``q``, ranked the way the app ranks them.

    ``candidates`` -- proteins the catalog does not hold yet -- stays empty
    until the resolver lands. It is a separate list rather than a flag on a
    protein so the screen can keep saying only what it knows: nothing about a
    candidate has been resolved, and nothing about it has been checked.
    """
    with _catalog_errors():
        proteins = catalog.search(q, limit=limit)
    return {"proteins": proteins, "candidates": []}


@router.get("/proteins/suggest", response_model=SuggestResponse)
def suggest_proteins(
    q: str = "",
    limit: int = Query(12, ge=1, le=25),
) -> dict:
    """Every reviewed human protein ``q`` might mean, best first.

    Not the catalog: the twenty are among the answers, tagged ``listed``, but
    so is everything else UniProt has reviewed, each saying whether the app
    can build it. The client asks this as the reader types, so it is one index
    range and, only when no prefix matches, one near-miss lookup.
    """
    with _catalog_errors():
        return suggest.suggest(q, limit=limit)


@router.get("/protein/{slug}", response_model=ProteinDetail)
def read_protein(slug: str) -> dict:
    """One protein, whole: what `targets.py` and `protein_catalog.dart` both held."""
    with _catalog_errors():
        found = catalog.detail(slug)
    if found is None:
        raise HTTPException(
            status_code=404,
            detail="No protein {!r} in the catalog.".format(slug),
        )
    return found


@router.get("/protein/{slug}/tracks", response_model=dict)
def read_protein_tracks(slug: str) -> dict:
    """Where each of this protein's tracks is, and what state it is in.

    Separate from `/protein/{slug}` because the two change on different clocks:
    the row is written once when the protein resolves, while the tracks flip
    from pending to ready over the minutes after it. A walk already open
    re-reads this without re-reading the protein.
    """
    with _catalog_errors():
        known = catalog.exists(slug)
    if not known:
        raise HTTPException(
            status_code=404,
            detail="No protein {!r} in the catalog.".format(slug),
        )
    return {
        kind: Track(**track).model_dump()
        for kind, track in tracks.for_slug(slug).items()
    }
