"""The read endpoint: GET /gene/{id}/{gene}."""

import contextlib
from urllib.error import HTTPError, URLError

from fastapi import APIRouter, HTTPException

from .entrez_client import fetch_genbank_record
from .genbank_parser import GeneNotFound, extract_gene
from .schemas import GeneResponse

router = APIRouter()


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
