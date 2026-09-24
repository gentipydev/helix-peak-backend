"""Read-through cache for whole GenBank records.

Keyed by accession, not by (accession, gene). The expensive step is the NCBI
round trip for a record; ``extract_gene`` is pure and costs ~0.1 ms, so one
cached record serves every gene in it -- NG_007114 answers INS, TH and
INS-IGF2 from a single row.

What is stored is the raw text NCBI returned, not a parsed payload. Reparsing
is ~2 ms against a ~0.5-3 s fetch, and raw text keeps the cache independent of
both the parser and the response schema: a fix to ``genbank_parser`` changes
what clients see on the next request without invalidating anything.

Every cache operation is best-effort. A database that is missing, asleep or
broken must cost latency and nothing else, so failures here are logged and the
request continues to NCBI.
"""

import logging

from Bio.SeqRecord import SeqRecord

from . import db
from .config import settings
from .entrez_client import fetch_genbank_text, parse_genbank_text

logger = logging.getLogger(__name__)

SCHEMA = """
create table if not exists genbank_record (
    accession  text primary key,
    record     text not null,
    fetched_at timestamptz not null default now()
)
"""

_SELECT = """
select record from genbank_record
where accession = %s and fetched_at > now() - make_interval(days => %s)
"""

_UPSERT = """
insert into genbank_record (accession, record, fetched_at)
values (%s, %s, now())
on conflict (accession) do update
set record = excluded.record, fetched_at = excluded.fetched_at
"""


def create_schema() -> None:
    """Create the cache table if it is not there yet.

    Runs at startup rather than as a separate migration step: one table, no
    history to migrate, and a deploy that cannot half-apply. A failure is
    logged and swallowed -- without the table every read simply misses.
    """
    if db.pool is None:
        return
    try:
        with db.pool.connection() as conn:
            conn.execute(SCHEMA)
    except Exception:
        logger.exception("Could not create the cache table; running without a cache")


def _get(accession: str):
    """The cached text for ``accession``, or None on a miss or any failure."""
    if db.pool is None:
        return None
    try:
        with db.pool.connection() as conn:
            row = conn.execute(_SELECT, (accession, settings.cache_ttl_days)).fetchone()
    except Exception:
        logger.exception("Cache read failed for %s", accession)
        return None
    return row[0] if row else None


def _put(accession: str, text: str) -> None:
    if db.pool is None:
        return
    try:
        with db.pool.connection() as conn:
            conn.execute(_UPSERT, (accession, text))
    except Exception:
        logger.exception("Cache write failed for %s", accession)


def _forget(accession: str) -> None:
    if db.pool is None:
        return
    try:
        with db.pool.connection() as conn:
            conn.execute("delete from genbank_record where accession = %s", (accession,))
    except Exception:
        logger.exception("Could not evict %s", accession)


def fetch(accession: str) -> SeqRecord:
    """One GenBank record, from the cache when possible.

    Raises ValueError when NCBI answers with something that is not a GenBank
    record, which the router reads as a 404. Nothing unparseable is ever
    stored: the parse happens before the write, so a bad upstream response
    cannot be served back from the cache later.
    """
    cached = _get(accession)
    if cached is not None:
        try:
            return parse_genbank_text(cached)
        except ValueError:
            # A row written by an older, looser version of this code, or a
            # truncated write. Drop it and read through.
            logger.warning("Evicting unparseable cached record %s", accession)
            _forget(accession)

    text = fetch_genbank_text(accession)
    record = parse_genbank_text(text)
    _put(accession, text)
    return record
