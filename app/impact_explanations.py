"""Serve generated evidence by exact accession/gene identity.

Two paths to the same payload. The first asks `protein_track` where the bytes
are and sends the client straight to Supabase storage, so a 4 MB attribution
file never passes through this service. The second reads a local directory, and
exists for a checkout with no database -- development, and the test suite, which
runs without a socket.

The validation that used to happen on every read now happens once, on the way
in: `tool/upload_tracks.py` refuses a payload whose scorer, schema version or
`impact_sha256` disagree with the AVI track it explains, so a row that says
`ready` is a row that was checked. That is the same order `record_cache` uses
when it parses before it writes.
"""
import json
import logging
from typing import Optional

from fastapi import HTTPException

from . import db
from .config import settings
from .tracks import public_url

logger = logging.getLogger(__name__)

_UNAVAILABLE = "AVI explanations are temporarily unavailable."

# Left join, not inner: a catalogued gene with no explanations has to be
# distinguishable from a gene the catalog has never heard of. The first is a
# settled 404; the second is a question this database cannot answer.
_SELECT = """
select t.state, t.bucket, t.object_path
from protein p
left join protein_track t
       on t.slug = p.slug and t.kind = 'impact_explanations'
where p.accession = %s and p.gene = %s
"""

# The catalog holds this gene and has no explanations for it. Not a fallback
# case: falling back would let a missing mount answer a question the database
# already answered.
NOT_COVERED = object()


def stored_url(accession: str, gene: str):
    """Where storage holds this gene's explanations.

    Returns a URL to redirect to, ``NOT_COVERED`` when the catalog knows this
    gene and has nothing for it, or None when the question is unanswerable
    here -- no database, a broken one, or a gene that is not in the catalog --
    and the local directory should be tried instead.
    """
    if db.pool is None:
        return None
    try:
        with db.pool.connection() as conn:
            row = conn.execute(_SELECT, (accession, gene)).fetchone()
    except Exception:
        logger.exception("Track lookup failed for %s/%s", accession, gene)
        return None
    if row is None:
        return None
    state, bucket, object_path = row
    if state != "ready":
        return NOT_COVERED
    url = public_url(bucket, object_path)
    if url is None:
        # The row says the object is there and this service cannot say where.
        logger.error("SUPABASE_URL is unset; stored explanations are unreachable")
        raise HTTPException(503, _UNAVAILABLE)
    return url


def read_impact_explanations(accession: str, gene: str) -> dict:
    """The saved payload, read off the local directory."""
    # A directory that is not there is a deployment fault, not an answer about
    # this gene. `Path.glob` yields nothing for a missing directory rather than
    # raising, so without this check an unmounted volume reports "not included
    # for this gene" for every gene -- which is the one thing the 404 must not
    # be able to mean, because three genes ship expecting a payload.
    directory = settings.impact_explanations_dir
    if not directory.is_dir():
        logger.error(
            "No stored track and IMPACT_EXPLANATIONS_DIR is not a directory: %s. "
            "Every request will report explanations as unavailable.",
            directory,
        )
        raise HTTPException(503, _UNAVAILABLE)

    # Request values never become paths. The header identifies each generated
    # bundle, including the transcript, assembly and original impact digest.
    for path in sorted(directory.glob("*.json")):
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError):
            raise HTTPException(503, _UNAVAILABLE)
        if not isinstance(data, dict):
            raise HTTPException(503, _UNAVAILABLE)
        if data.get("accession") != accession or data.get("gene") != gene:
            continue
        if data.get("schema_version") != 1 or data.get("scorer") != "AVI_SCORE_FEATURE_IMPORTANCE":
            raise HTTPException(503, _UNAVAILABLE)
        return data
    raise HTTPException(404, "AVI explanations are not included for this gene.")
