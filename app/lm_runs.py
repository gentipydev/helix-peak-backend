"""Where a recorded training run of the de novo peptide model is.

The app plays one, stage by stage. Like a track (`tracks.py`), this names the
bytes and never carries them: the row says where the run is in storage and the
client fetches it from there. `pipeline/denovo/publish_trace.py` writes the row
(`lm_run`, migration 0013).
"""

import logging
from typing import Optional

from . import db
from .catalog import CatalogUnavailable
from .tracks import public_url

logger = logging.getLogger(__name__)

_SELECT = """
select name, bucket, object_path, bytes, sha256, format, provenance
from lm_run
where name = %s
"""


def get(name: str) -> Optional[dict]:
    """The run called ``name``, or None where there is none.

    A database that is not there, or cannot be read, is not an answer: that is
    ``CatalogUnavailable`` (503), never a 404 saying the run does not exist.
    Nor is a row whose bytes this service cannot address.
    """
    if db.pool is None:
        raise CatalogUnavailable("No database is configured for training runs.")
    try:
        with db.pool.connection() as conn:
            row = conn.execute(_SELECT, (name,)).fetchone()
    except Exception as exc:
        logger.exception("Training run read failed for %s", name)
        raise CatalogUnavailable("The training runs could not be read.") from exc
    if row is None:
        return None
    (run, bucket, object_path, size, sha256, fmt, provenance) = row
    url = public_url(bucket, object_path)
    if url is None:
        logger.error("SUPABASE_URL is unset; training run %s is unreachable", run)
        raise CatalogUnavailable("This run is stored but the service cannot address it.")
    return {
        "name": run,
        "url": url,
        "bytes": size,
        "sha256": sha256,
        "format": fmt,
        "provenance": provenance or {},
    }
