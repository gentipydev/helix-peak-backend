"""The single endpoint: GET /fetch/{id}."""

from urllib.error import HTTPError, URLError

from fastapi import APIRouter, HTTPException

from .entrez_client import fetch_genbank

router = APIRouter()

# NCBI does not signal "no such record" with an empty body. A well-formed but
# unknown accession comes back as HTTP 400, while an id it cannot even parse
# comes back as HTTP 200 whose body is literally
# "Error: F a i l e d  t o  u n d e r s t a n d  i d : ..." -- so a valid
# response has to be recognised positively, by its GenBank header.
_GENBANK_HEADER = "LOCUS"


def _upstream_message(exc: HTTPError) -> str:
    """Best-effort read of the error body NCBI sent with an HTTPError."""
    try:
        body = exc.read()
    except (OSError, ValueError, AttributeError):
        return exc.reason or str(exc)
    if isinstance(body, bytes):
        body = body.decode("utf-8", errors="replace")
    return body.strip() or (exc.reason or str(exc))


@router.get("/fetch/{id}")
def fetch(id: str) -> dict:
    """Fetch one GenBank record from NCBI and return its raw text.

    Deliberately a sync ``def``: ``Entrez.efetch`` is blocking urllib, so
    FastAPI runs this in a threadpool. Making it ``async def`` would block the
    event loop for the whole round trip to NCBI.
    """
    try:
        content = fetch_genbank(id)
    except HTTPError as exc:
        # db/rettype/retmode are hardcoded, so `id` is the only thing that can
        # make the request invalid -- an upstream 400 means "no such record".
        if exc.code == 400:
            raise HTTPException(
                status_code=404,
                detail="No GenBank record found for id {!r}.".format(id),
            )
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
        # socket.timeout and friends; also the fallback for transport errors
        # that Biopython re-raises without a URLError wrapper.
        raise HTTPException(
            status_code=502,
            detail="Could not reach NCBI: {}".format(exc),
        )

    if not content.lstrip().startswith(_GENBANK_HEADER):
        raise HTTPException(
            status_code=404,
            detail="No GenBank record found for id {!r}.".format(id),
        )

    return {"id": id, "content": content}
