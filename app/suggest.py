"""Search suggestions over every reviewed human protein.

The catalog is the twenty the app lists; this is everything a reader might
type instead. A suggestion says what the app can do with the protein today:

- ``listed``: one of the twenty, opened from the list;
- ``ready``: built on demand earlier, and opened at once;
- ``buildable``: not built yet, and the pipeline can start from its row;
- ``unavailable``: not buildable, and ``reason`` says why.

``building`` is true while a build is under way: a request queued or running,
or a constraint bake not yet over. A protein is ``ready`` once its row is
written, minutes before ESM-2 has scored it, and a walk opened then keeps the
track it read as pending; so the app watches the build instead of opening it.
``stopped`` is true for one whose scoring the reader who asked stopped: it
opens unscored, and asking for it again scores it.

Ranking is in tiers. An exact match of a whole term comes first -- an
accession, a symbol, a synonym or a name, never one word of a name -- then a
prefix of the symbol, of a synonym, of the protein's name, of another of its
names, and of a word of any of them. Within a tier the listed proteins come
first, then the built ones, then the closest match: a shorter matching term
is nearer to what was typed. Only when nothing matches by prefix is a near miss
looked for, so a typo such as "insuln" still finds insulin, and an exact
accession is not padded out with accessions that merely look like it.

The proteins built on demand are listed here too (`/proteins/built`), as the
same suggestions: every install's, since a built protein is shared, newest
first. One still being built is left out until it is done, since the app
watches that one rather than opens it; one whose scoring was stopped is kept,
and says so.

Reads only, and like the catalog it has no second source: an unreadable index is
reported as unavailable, never as an empty one.
"""

import logging
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple

from . import db
from .catalog import CatalogUnavailable
from .protein_index import bounds, normalize

logger = logging.getLogger(__name__)

# A needle this short matches symbols and synonyms only. "a" as a prefix of
# every word of every name would rank most of the index to answer one key.
_SHORT_NEEDLE = 2

# Near misses are looked for only past this length: two letters have too few
# trigrams to say anything is near them.
_FUZZY_NEEDLE = 3

_ROW = """
    i.uniprot, i.gene, i.name, i.length, i.buildable, i.unavailable_reason,
    p.slug, p.display, p.catalog_order
"""

# The term range is every tier but the near miss: each prefix is a plain range
# on a "C"-collated column, so the btree serves it prepared or not. A term
# scores its tier times a thousand plus its length, so a protein's best score
# is its best tier and, within it, its closest term; no term is that long.
_RANKED = """
with hit as (
    select t.uniprot, t.gene,
           min(case when t.term = %(needle)s and t.kind < 5 then 0 else t.kind end * 1000
               + length(t.term)) as score
    from protein_index_term t
    where t.term >= %(lo)s and t.term < %(hi)s
      and (not %(short)s or t.kind in (1, 2) or (t.term = %(needle)s and t.kind < 5))
    group by t.uniprot, t.gene
)
select {row}, h.score / 1000 as tier
from hit h
join protein_index i on i.uniprot = h.uniprot and i.gene = h.gene
left join protein p on p.gene = i.gene and p.taxon_id = 9606 and i.gene <> ''
order by h.score / 1000, p.catalog_order nulls last, (p.slug is null), h.score %% 1000,
         i.annotation_score desc, i.existence, i.gene, i.uniprot
limit %(limit)s
""".format(row=_ROW)

# `%%` is the trigram similarity operator, escaped for the driver. It uses the
# trigram index and pg_trgm's default threshold of 0.3.
#
# Ranked whole before the limit, never cut first: "insuln" is exactly as near
# the word "insulin" in a hundred names as it is to insulin's own name, and a
# limit taken over those ties kept whichever twelve the plan met first. A word
# gives way to a whole term at equal nearness, and then the order is the
# prefix tiers' own.
_FUZZY = """
with near as (
    select t.uniprot, t.gene,
           max(similarity(t.term, %(needle)s) - case when t.kind = 5 then 0.01 else 0 end) as sim
    from protein_index_term t
    where t.term %% %(needle)s
    group by t.uniprot, t.gene
)
select {row}, 6 as tier
from near n
join protein_index i on i.uniprot = n.uniprot and i.gene = n.gene
left join protein p on p.gene = i.gene and p.taxon_id = 9606 and i.gene <> ''
order by n.sim desc, p.catalog_order nulls last, (p.slug is null),
         i.annotation_score desc, i.existence, i.gene, i.uniprot
limit %(limit)s
""".format(row=_ROW)

_RELEASE = "select uniprot, mane from protein_index_release"

# Which of these genes have a build under way (a request not yet resolved, or a
# protein resolved on demand whose latest ESM-2 bake is not over), and which a
# reader stopped while ESM-2 scored, with nothing on the track since.
_BUILDS = """
select lower(gene), 'building' from resolve_request
where state in ('queued', 'running') and lower(gene) = any(%(genes)s)
union
select lower(p.gene), case when j.state = 'stopped' then 'stopped' else 'building' end
from protein p
join lateral (
    select b.state from bake_job b where b.slug = p.slug and b.kind = 'constraint'
    order by b.requested_at desc, b.id desc limit 1
) j on true
join protein_track t on t.slug = p.slug and t.kind = 'constraint'
where p.catalog_order is null and lower(p.gene) = any(%(genes)s)
  and (j.state in ('queued', 'running') or (j.state = 'stopped' and t.state = 'absent'))
"""

# Every protein built on demand, newest first: a row with a null
# `catalog_order`, and the index row it was built from. Joined on the accession
# as well as the gene, because a gene can have more than one index row (SIRPB1
# has a buildable one and one that is not), and a built protein is listed once.
#
# Paged by (resolved_at, slug) rather than by offset, as the catalog is by slug,
# so a protein built between two requests cannot shift a page and hide a row.
_BUILT = """
select {row}, p.resolved_at
from protein p
join protein_index i on i.uniprot = p.uniprot and i.gene = p.gene
where p.catalog_order is null and p.taxon_id = 9606 and i.gene <> ''
  and (%(at)s::timestamptz is null or p.resolved_at < %(at)s
       or (p.resolved_at = %(at)s and p.slug > %(slug)s))
order by p.resolved_at desc, p.slug
limit %(limit)s
""".format(row=_ROW)

# A cursor is the last row's time to the microsecond, in UTC, and its slug:
# "2026-10-09T08:45:20.800935Z|ttr". Opaque to the client, which passes it back.
_CURSOR_TIME = "%Y-%m-%dT%H:%M:%S.%fZ"


class BadCursor(ValueError):
    """A ``before`` this service did not give out."""


def _suggestion(row, building: bool = False, stopped: bool = False) -> dict:
    (uniprot, gene, name, length, buildable, reason,
     slug, display, catalog_order) = row[:9]
    if catalog_order is not None:
        status, reason = "listed", None
    elif slug is not None:
        status, reason = "ready", None
    elif buildable:
        # The slug a build would give it. New proteins take lower(gene); only
        # the twenty kept the slugs they shipped with.
        status, slug, reason = "buildable", gene.lower(), None
    else:
        status, slug = "unavailable", None
    return {
        "uniprot": uniprot,
        "gene": gene or None,
        "name": name,
        "display": display,
        "length": length,
        "slug": slug,
        "status": status,
        "reason": reason,
        "building": building,
        "stopped": stopped,
    }


def suggest(query: str, limit: int = 12) -> Dict:
    """Up to ``limit`` proteins for what has been typed so far, best first."""
    needle = normalize(query)
    if not needle:
        return {"q": query, "release": None, "suggestions": []}

    pool = db.pool
    if pool is None:
        raise CatalogUnavailable("No database is configured for the protein index.")

    lo, hi = bounds(needle)
    params = {
        "needle": needle,
        "lo": lo,
        "hi": hi,
        "short": len(needle) <= _SHORT_NEEDLE,
        "limit": limit,
    }
    try:
        with pool.connection() as conn:
            rows = conn.execute(_RANKED, params).fetchall()
            if not rows and len(needle) >= _FUZZY_NEEDLE:
                rows = conn.execute(_FUZZY, params).fetchall()
            release = conn.execute(_RELEASE).fetchone()
            genes = sorted({row[1].lower() for row in rows if row[1]})
            said = []
            if genes:
                said = conn.execute(_BUILDS, {"genes": genes}).fetchall()
            building = {gene for gene, what in said if what == "building"}
            stopped = {gene for gene, what in said if what == "stopped"} - building
    except Exception as exc:
        logger.exception("Protein index search failed for %r", query)
        raise CatalogUnavailable("The protein index could not be searched.") from exc

    label: Optional[str] = None
    if release is not None:
        label = "UniProt {} · MANE {}".format(release[0], release[1])
    suggestions: List[dict] = [
        _suggestion(row, bool(row[1]) and row[1].lower() in building,
                    bool(row[1]) and row[1].lower() in stopped)
        for row in rows
    ]
    return {"q": query, "release": label, "suggestions": suggestions}


def _cursor(resolved_at: datetime, slug: str) -> str:
    return "{}|{}".format(resolved_at.astimezone(timezone.utc).strftime(_CURSOR_TIME), slug)


def _parse_cursor(cursor: str) -> Tuple[datetime, str]:
    at, bar, slug = cursor.partition("|")
    try:
        when = datetime.strptime(at, _CURSOR_TIME).replace(tzinfo=timezone.utc)
    except ValueError:
        when = None
    if when is None or not bar or not slug:
        raise BadCursor("Not a cursor this service gave out: {!r}.".format(cursor))
    return when, slug


def built(limit: int = 50, before: Optional[str] = None) -> Tuple[List[dict], Optional[str]]:
    """Up to ``limit`` proteins built on demand, newest first, and the cursor
    for the next page: null on the last.

    A protein whose build is under way is left out after the page is cut, so a
    page can come back short with more to follow; ``next`` says which.
    """
    at, after = _parse_cursor(before) if before else (None, None)

    pool = db.pool
    if pool is None:
        raise CatalogUnavailable("No database is configured for the protein index.")

    params = {"at": at, "slug": after, "limit": limit + 1}
    try:
        with pool.connection() as conn:
            rows = conn.execute(_BUILT, params).fetchall()
            more = len(rows) > limit
            rows = rows[:limit]
            genes = sorted({row[1].lower() for row in rows})
            said = []
            if genes:
                said = conn.execute(_BUILDS, {"genes": genes}).fetchall()
    except Exception as exc:
        logger.exception("Built proteins read failed")
        raise CatalogUnavailable("The proteins built on demand could not be read.") from exc

    building = {gene for gene, what in said if what == "building"}
    stopped = {gene for gene, what in said if what == "stopped"}
    proteins = [
        _suggestion(row, stopped=row[1].lower() in stopped)
        for row in rows
        if row[1].lower() not in building
    ]
    next_cursor = _cursor(rows[-1][9], rows[-1][6]) if more else None
    return proteins, next_cursor
