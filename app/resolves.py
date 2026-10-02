"""Proteins the catalog does not list, asked for by a reader: Phase 6.

The service does not resolve anything itself. A request is a row in
`resolve_request`; the resolver on Modal (`pipeline/resolver/`) works it, and
writes the protein row and its record track when it is done. What this module
does is decide what an ask means, write the row when it is a new one, and wake
the resolver so the reader does not wait for its schedule.

What a protein can be, as `ResolveResponse.state` says it:

- ``ready``: there is a protein row for this gene, listed or resolved. ``slug``
  is where it opens, which for the twenty is the slug they shipped with.
- ``pending``: a request is queued or running.
- ``refused``: the resolver declined, and ``reason`` is why. Final, for the
  resolver version that said it; asking again changes nothing.
- ``failed``: it broke more often than it is retried. Asking again queues a new
  request.
- ``unavailable``: the index cannot build this protein, and ``reason`` is why
  (no MANE Select transcript, or one for another isoform).
- ``buildable``: nobody has asked yet, and nothing stands in the way. Only a
  read says this; an ask turns it into ``pending``.

A gene the index does not hold at all is not a state but a 404, and an
unreadable database is a 503, never an answer: the catalog's rule.

The one thing that costs money is a request, so a day's are capped
(``RESOLVES_PER_DAY``). Two readers asking for one gene make one request, which
the database holds to (`resolve_request_inflight_idx`).
"""

import logging
import urllib.request
from typing import Optional, Tuple

from .catalog import CatalogUnavailable, _require_pool
from .config import settings

logger = logging.getLogger(__name__)

HUMAN = 9606


class DailyCapReached(Exception):
    """Today's requests are used up."""


_PROTEIN = """
select slug from protein where lower(gene) = lower(%s) and taxon_id = 9606
"""

# A buildable entry first: a gene can name several entries, at most one of
# them buildable (`protein_index_buildable_gene`).
_INDEX = """
select uniprot, gene, buildable, unavailable_reason from protein_index
where lower(gene) = lower(%s) and gene <> ''
order by buildable desc, uniprot
limit 1
"""

_LATEST = """
select state, reason from resolve_request where lower(gene) = lower(%s)
order by requested_at desc, id desc
limit 1
"""

_TODAY = """
select count(*) from resolve_request where requested_at > now() - interval '1 day'
"""

_QUEUE = """
insert into resolve_request (gene, uniprot, slug) values (%s, %s, %s)
on conflict (lower(gene)) where state in ('queued', 'running') do nothing
returning id
"""


def _said(state: str, slug: Optional[str], reason: Optional[str] = None) -> dict:
    return {"slug": slug, "state": state, "reason": reason}


def _current(conn, gene: str) -> Tuple[Optional[dict], Optional[tuple]]:
    """What this gene is now, and its index entry when a request could be made.

    The first is None for a gene the index does not hold; the second is set
    only where nothing stands in the way of a new request.
    """
    found = conn.execute(_PROTEIN, (gene,)).fetchone()
    if found is not None:
        return _said("ready", found[0]), None
    entry = conn.execute(_INDEX, (gene,)).fetchone()
    if entry is None:
        return None, None
    uniprot, symbol, buildable, unavailable = entry
    slug = symbol.lower()
    if not buildable:
        return _said("unavailable", None, unavailable), None
    latest = conn.execute(_LATEST, (gene,)).fetchone()
    if latest is not None:
        state, reason = latest
        if state in ("queued", "running"):
            return _said("pending", slug), None
        if state == "refused":
            return _said("refused", None, reason), None
        if state == "failed":
            return _said("failed", slug, reason), (uniprot, symbol, slug)
    return None, (uniprot, symbol, slug)


def current(gene: str, taxon: int = HUMAN) -> Optional[dict]:
    """What a gene's protein is now. None where the index does not hold it."""
    if taxon != HUMAN:
        return None
    pool = _require_pool()
    try:
        with pool.connection() as conn:
            said, askable = _current(conn, gene)
    except Exception as exc:
        logger.exception("Resolve state read failed for %r", gene)
        raise CatalogUnavailable("The catalog could not be read.") from exc
    if said is None and askable is not None:
        # Buildable, and never asked for.
        return _said("buildable", askable[2])
    return said


def request(gene: str, taxon: int = HUMAN) -> Tuple[Optional[dict], bool]:
    """Ask for a gene's protein. Returns what it is now, and whether this ask
    queued a new request (and so should wake the resolver).

    Raises DailyCapReached where a new request would be one too many today.
    """
    if taxon != HUMAN:
        return None, False
    pool = _require_pool()
    try:
        with pool.connection() as conn:
            said, askable = _current(conn, gene)
            if askable is None:
                return said, False
            if conn.execute(_TODAY).fetchone()[0] >= settings.resolves_per_day:
                raise DailyCapReached()
            uniprot, symbol, slug = askable
            queued = conn.execute(_QUEUE, (symbol, uniprot, slug)).fetchone()
    except DailyCapReached:
        raise
    except Exception as exc:
        logger.exception("Resolve request failed for %r", gene)
        raise CatalogUnavailable("The request could not be recorded.") from exc
    # Not queued means another reader's request got there in the same moment:
    # it is pending either way, and only the one that queued it wakes anyone.
    return _said("pending", slug), queued is not None


def wake() -> None:
    """Tell the resolver on Modal a request is waiting. Best effort.

    Unconfigured, or unanswered, this does nothing a reader would notice except
    wait: the resolver's own schedule takes every queued request within five
    minutes. So a failure is logged and never raised.
    """
    url = settings.modal_wake_url
    if not url:
        return
    headers = {}
    if settings.modal_key and settings.modal_secret:
        headers = {"Modal-Key": settings.modal_key, "Modal-Secret": settings.modal_secret}
    call = urllib.request.Request(url, data=b"", method="POST", headers=headers)
    try:
        with urllib.request.urlopen(call, timeout=10):
            return
    except Exception:
        logger.warning("Could not wake the resolver at %s; its schedule will take the request",
                       url, exc_info=True)
