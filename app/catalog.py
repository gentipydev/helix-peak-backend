"""The protein catalog: what `protein_catalog.dart` held as 614 lines of const.

Reads only. Rows are written by the seeder (`tool/seed_catalog.py`) and, later,
by the resolver.

Unlike `record_cache`, a failure here has nothing to read through to. The cache
can lose its database and still answer from NCBI; a catalog that cannot reach
its rows has no second source, and inventing one on this side would put a
second copy of the truth where it can disagree with the first. So these raise
``CatalogUnavailable`` and the router reports it as a 503. The client holds its
own on-device copy, which is what keeps search working while this is down.
"""

import logging
from typing import Dict, List, Optional, Tuple

from . import db
from .config import settings
from . import tracks

logger = logging.getLogger(__name__)


class CatalogUnavailable(Exception):
    """The catalog rows could not be read."""


# The five names `matching()` compares whole, in `protein_catalog.dart`:
# display, gene, slug, uniprot, accession. Kept exactly as they are on the Dart
# side so the twenty rank identically to the const list they replace.
_WHOLE_NAME_KINDS = ["display", "gene", "slug", "uniprot", "accession"]

_SUMMARY_COLUMNS = """
    p.slug, p.display, p.gene, p.uniprot, p.accession, p.summary,
    p.residues, p.exons, p.chains, p.bridges, p.catalog_order,
    p.chain_name, p.structure
"""

_DETAIL_COLUMNS = _SUMMARY_COLUMNS + """,
    p.mature_peptides, p.transcript_id, p.protein_id,
    p.regions, p.disulfides, p.provenance, p.resolver_version
"""

# Paged by slug, and ordered by slug because of it: a keyset cursor is only
# stable when the ordering is the column being compared. The reading order
# travels as `catalog_order` on each row and is restored by the client, which
# holds the whole catalog anyway. It used to travel as the order of a const
# list; it has to travel as something.
#
# The catalog is the list, and the list is the curated rows only. A protein
# built on demand has a row too, with a null `catalog_order`, but it is found
# through `/proteins/suggest` and never joins the list someone chose.
_PAGE = """
select {columns}
from protein p
where p.catalog_order is not null
  and (%s::text is null or p.slug > %s)
order by p.slug
limit %s
""".format(columns=_SUMMARY_COLUMNS)

_DETAIL = "select {} from protein p where p.slug = %s".format(_DETAIL_COLUMNS)

# The four tiers of `_closeness`, in SQL. 0 exact, 1 prefix (on a whole name or
# on any word of the display name), 2 contained in a whole name, 3 contained in
# the summary. Synonyms -- which the const catalog had no notion of -- join at
# tier 2, so they add recall without reordering anything that ranked before.
_SEARCH = """
with matched as (
    select a.slug,
           min(case
                 when a.kind = any(%(whole)s) and lower(a.alias) = %(needle)s then 0
                 when (a.kind = any(%(whole)s) or a.kind = 'word')
                      and lower(a.alias) like %(prefix)s then 1
                 when (a.kind = any(%(whole)s) or a.kind = 'synonym')
                      and lower(a.alias) like %(infix)s then 2
               end) as tier
    from protein_alias a
    group by a.slug
)
select {columns},
       coalesce(m.tier,
                case when position(%(needle)s in lower(p.summary)) > 0 then 3 end) as tier
from protein p
left join matched m on m.slug = p.slug
where p.catalog_order is not null
  and (m.tier is not null or position(%(needle)s in lower(p.summary)) > 0)
order by tier, p.catalog_order nulls last, p.display, p.slug
limit %(limit)s
""".format(columns=_SUMMARY_COLUMNS)


def _summary(row) -> dict:
    (slug, display, gene, uniprot, accession, summary,
     residues, exons, chains, bridges, catalog_order, chain_name, structure) = row[:13]
    structure = structure or {}
    return {
        "slug": slug,
        "display": display,
        "gene": gene,
        "uniprot": uniprot,
        "accession": accession,
        "summary": summary,
        "facts": {
            "residues": residues,
            "exons": exons,
            "chains": chains,
            "bridges": bridges,
        },
        "catalog_order": catalog_order,
        "chain": chain_name,
        "chains": structure.get("chains", []),
        "structure": structure.get("chrome"),
        "tracks": {},
    }


def _attach_states(summaries: List[dict]) -> List[dict]:
    """One track-state query for the whole page, not one per protein."""
    if not summaries:
        return summaries
    states = tracks.states_for([s["slug"] for s in summaries])
    for summary in summaries:
        summary["tracks"] = states.get(summary["slug"], {})
    return summaries


def _require_pool():
    if db.pool is None:
        raise CatalogUnavailable("No database is configured for the catalog.")
    return db.pool


def page(limit: int = 200, cursor: Optional[str] = None) -> Tuple[List[dict], Optional[str]]:
    """One page of the catalog, and the cursor for the next.

    Paged by slug rather than by offset so a protein resolved between two
    requests cannot shift the page under the client and hide a row.
    """
    pool = _require_pool()
    try:
        with pool.connection() as conn:
            rows = conn.execute(_PAGE, (cursor, cursor, limit + 1)).fetchall()
    except Exception as exc:
        logger.exception("Catalog page read failed")
        raise CatalogUnavailable("The catalog could not be read.") from exc

    more = len(rows) > limit
    rows = rows[:limit]
    summaries = _attach_states([_summary(row) for row in rows])
    return summaries, (summaries[-1]["slug"] if more and summaries else None)


def search(query: str, limit: int = 20) -> List[dict]:
    """The catalog rows matching ``query``, ranked as `matching()` ranks them."""
    needle = query.strip().lower()
    if not needle:
        found, _ = page(limit=limit)
        return found

    pool = _require_pool()
    # LIKE metacharacters in a user's query are matched literally: someone
    # typing "p53%" is looking for a protein, not writing a pattern.
    escaped = needle.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")
    params = {
        "whole": _WHOLE_NAME_KINDS,
        "needle": needle,
        "prefix": escaped + "%",
        "infix": "%" + escaped + "%",
        "limit": limit,
    }
    try:
        with pool.connection() as conn:
            rows = conn.execute(_SEARCH, params).fetchall()
    except Exception as exc:
        logger.exception("Catalog search failed for %r", query)
        raise CatalogUnavailable("The catalog could not be searched.") from exc
    return _attach_states([_summary(row) for row in rows])


def exists(slug: str) -> bool:
    """Whether the catalog holds this slug, without reading the row.

    `/protein/{slug}/tracks` is polled while a protein's bakes land, and it
    asks about tracks; reading the whole protein row to decide whether to 404
    would put the heaviest read on the most repeated request.
    """
    pool = _require_pool()
    try:
        with pool.connection() as conn:
            row = conn.execute("select 1 from protein where slug = %s", (slug,)).fetchone()
    except Exception as exc:
        logger.exception("Catalog existence check failed for %s", slug)
        raise CatalogUnavailable("The catalog could not be read.") from exc
    return row is not None


def detail(slug: str) -> Optional[Dict]:
    """One protein, whole, or None when the catalog does not hold it."""
    pool = _require_pool()
    try:
        with pool.connection() as conn:
            row = conn.execute(_DETAIL, (slug,)).fetchone()
    except Exception as exc:
        logger.exception("Catalog detail read failed for %s", slug)
        raise CatalogUnavailable("The catalog could not be read.") from exc
    if row is None:
        return None

    found = _summary(row)
    (mature_peptides, transcript_id, protein_id,
     regions, disulfides, provenance, resolver_version) = row[13:]

    found.update({
        "mature_peptides": mature_peptides,
        "transcript_id": transcript_id,
        "protein_id": protein_id,
        "regions": regions or [],
        "disulfides": disulfides or [],
        "provenance": provenance or {},
        "resolver_version": resolver_version,
    })
    found["tracks"] = tracks.states_for([slug]).get(slug, {})
    return found
