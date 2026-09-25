"""Search suggestions over every reviewed human protein.

The catalog is the twenty the app lists; this is everything a reader might
type instead. A suggestion says what the app can do with the protein today:

- ``listed``: one of the twenty, opened from the list;
- ``ready``: built on demand earlier, and opened at once;
- ``buildable``: not built yet, and the pipeline can start from its row;
- ``unavailable``: not buildable, and ``reason`` says why.

Ranking is in tiers. An exact match of a whole term comes first -- an
accession, a symbol, a synonym or a name, never one word of a name -- then a
prefix of the symbol, of a synonym, of the protein's name, of another of its
names, and of a word of any of them. Within a tier the listed proteins come
first, then the built ones, then the closest match: a shorter matching term
is nearer to what was typed. Only when nothing matches by prefix is a near miss
looked for, so a typo such as "insuln" still finds insulin, and an exact
accession is not padded out with accessions that merely look like it.

Reads only, and like the catalog it has no second source: an unreadable index is
reported as unavailable, never as an empty one.
"""

import logging
from typing import Dict, List, Optional

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
_FUZZY = """
with near as (
    select t.uniprot, t.gene, max(similarity(t.term, %(needle)s)) as sim
    from protein_index_term t
    where t.term %% %(needle)s
    group by t.uniprot, t.gene
    order by sim desc
    limit %(limit)s
)
select {row}, 6 as tier
from near n
join protein_index i on i.uniprot = n.uniprot and i.gene = n.gene
left join protein p on p.gene = i.gene and p.taxon_id = 9606 and i.gene <> ''
order by n.sim desc, p.catalog_order nulls last, i.annotation_score desc, i.gene
""".format(row=_ROW)

_RELEASE = "select uniprot, mane from protein_index_release"


def _suggestion(row) -> dict:
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
    except Exception as exc:
        logger.exception("Protein index search failed for %r", query)
        raise CatalogUnavailable("The protein index could not be searched.") from exc

    label: Optional[str] = None
    if release is not None:
        label = "UniProt {} · MANE {}".format(release[0], release[1])
    suggestions: List[dict] = [_suggestion(row) for row in rows]
    return {"q": query, "release": label, "suggestions": suggestions}
