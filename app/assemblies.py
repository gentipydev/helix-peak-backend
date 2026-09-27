"""Where each track for an assembly is, and what state it is in.

An assembly (migrations/0006_assemblies.sql) is a molecule made of more than
one gene's chains -- the hemoglobin tetramer -- which the catalog cannot hold:
a catalog row is one gene's product. Its rows are in tables of their own,
which `/catalog`, `/catalog/search` and `/proteins/suggest` never read, so the
catalog stays twenty. This serves an assembly's tracks the way `tracks.py`
serves a protein's: it names bytes, and never carries them.
"""

import logging
from typing import Dict

from . import db
from .catalog import CatalogUnavailable
from .tracks import _absent, _row_to_track

logger = logging.getLogger(__name__)

# What an assembly is baked into: its states' morph, the models beside it
# named in the row's provenance.
KINDS = ("morph",)

_EXISTS = "select 1 from assembly where slug = %s"

_SELECT = """
select kind, state, reason, bucket, object_path, bytes, sha256,
       content_encoding, format, provenance
from assembly_track
where slug = %s
"""


def exists(slug: str) -> bool:
    """Whether the assembly table holds this slug.

    Unreadable -- no database, or a database the table has not been migrated
    into -- is ``CatalogUnavailable``, a 503: our own source, not an answer.
    """
    if db.pool is None:
        raise CatalogUnavailable("No database is configured for the catalog.")
    try:
        with db.pool.connection() as conn:
            row = conn.execute(_EXISTS, (slug,)).fetchone()
    except Exception as exc:
        logger.exception("Assembly existence check failed for %s", slug)
        raise CatalogUnavailable("The assemblies could not be read.") from exc
    return row is not None


def tracks_for(slug: str) -> Dict[str, dict]:
    """Every track for one assembly, keyed by kind, with every ``KINDS`` present."""
    found = {kind: _absent() for kind in KINDS}
    if db.pool is None:
        return found
    try:
        with db.pool.connection() as conn:
            rows = conn.execute(_SELECT, (slug,)).fetchall()
    except Exception:
        logger.exception("Assembly track read failed for %s", slug)
        return found
    for row in rows:
        found[row[0]] = _row_to_track(row)
    return found
