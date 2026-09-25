"""Seed the Supabase catalog from the two tables the twenty are curated in.

Run by hand. `targets.py` holds the coordinates, regions and disulfides;
`curated/catalog.json` holds the prose, the facts and the chain tints -- what
`protein_catalog.dart` held until Phase 3 of HANDOFF-ONDEMAND.md moved it here.
Neither is the whole row, which is the reason this exists: after the seed, one
table holds both, and `check_assets.py --against <base-url>` checks that the
service still agrees with both, field for field.

    DATABASE_URL=... python3 pipeline/seed_catalog.py
    python3 pipeline/seed_catalog.py --dry-run     # print the rows, touch nothing
    DATABASE_URL=... python3 pipeline/seed_catalog.py --check
                                                   # diff the rows with the live ones

A real run upserts the protein rows and replaces their aliases. It only inserts
track rows that are not there, so it never undoes an upload.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.paths import CURATED  # noqa: E402
from pipeline.targets import TARGETS, Target, partition  # noqa: E402

RESOLVER_VERSION = 0  # 0 means "not resolved": these rows were written by hand.


# --- reading the curated rows -----------------------------------------------


def curated_rows() -> dict[str, dict]:
    """Every curated row, keyed by slug, in the shape the builders below take.

    `order` is the row's place in the list; `chrome` is the fold page's prose,
    or None where a protein has no structure.
    """
    rows: dict[str, dict] = {}
    for protein in json.loads(CURATED.read_text())["proteins"]:
        rows[protein["slug"]] = {
            "order": protein["catalog_order"],
            "display": protein["display"],
            "summary": protein["summary"],
            "chain": protein.get("chain"),
            "facts": protein["facts"],
            "chains": protein["chains"],
            "chrome": protein.get("structure"),
            "impact_explanations": "impact_explanations" in protein.get("tracks", []),
        }
    return rows


# --- building the rows ------------------------------------------------------


def _aliases(target: Target, curated: dict) -> list[tuple[str, str]]:
    """Every string `matching()` can rank against, as (alias, kind).

    The five whole names it compares, each word of the display name, and the
    slug this protein shipped with -- which is not `lower(gene)` for most of
    them, and is the key every existing asset path is built from.
    """
    display = curated["display"]
    found = [
        (display, "display"),
        (target.gene, "gene"),
        (target.slug, "slug"),
        (target.uniprot, "uniprot"),
        (target.source.accession, "accession"),
    ]
    for word in re.split(r"[\s()\-]+", display):
        if word and word.lower() != display.lower():
            found.append((word, "word"))
    if target.slug != target.gene.lower():
        found.append((target.slug, "legacy_slug"))
    # One row per (alias, kind); the primary key would reject a repeat.
    return sorted({(alias, kind) for alias, kind in found if alias})


def protein_row(target: Target, curated: dict) -> dict:
    source = target.source
    structure = None
    if curated["chrome"] is not None:
        structure = {"chrome": curated["chrome"], "chains": curated["chains"]}

    return {
        "slug": target.slug,
        "gene": target.gene,
        "uniprot": target.uniprot,
        "accession": source.accession,
        "slice_start": source.seq_start,
        "slice_end": source.seq_stop,
        # A slice is always fetched plus-strand and re-based to 1; the record's
        # own strand is a fact about the gene, not about how it was cut.
        "slice_strand": 1 if source.seq_start is not None else None,
        "transcript_id": source.transcript_id,
        "protein_id": source.protein_id,
        "taxon_id": 9606,
        "display": curated["display"],
        "summary": curated["summary"],
        "chain_name": curated["chain"],
        "mature_peptides": target.mature_peptides,
        "residues": curated["facts"]["residues"],
        "exons": curated["facts"]["exons"],
        "chains": curated["facts"]["chains"],
        "bridges": curated["facts"]["bridges"],
        "regions": partition(target),
        "disulfides": [list(pair) for pair in target.disulfides],
        "structure": structure,
        "provenance": {
            # Written by a person, checked against the live entries by hand, and
            # recorded in docs/protein-verification.md. Nothing here was derived.
            "source": "targets.py + protein_catalog.dart",
            "prose": "hand",
            "uniprot_variants": [list(v) for v in target.uniprot_variants],
            "cleaved": target.cleaved,
        },
        "resolver_version": RESOLVER_VERSION,
        "catalog_order": curated["order"],
    }


def track_rows(target: Target, curated: dict) -> list[dict]:
    """The four booleans the app shipped with, as track states.

    Every one of the twenty is baked, so these all start `absent` with the
    object still in the app bundle -- the state flips to `ready` in the phase
    that uploads the bytes, not in this one. Seeding them `ready` now would
    name objects that are not there yet.
    """
    baked = {
        "constraint": target.scored,
        "impact": target.impact_scored,
        "clinvar": target.clinvar_available,
        "impact_explanations": curated["impact_explanations"],
        "structure": curated["chrome"] is not None,
    }
    rows = [{
        "slug": target.slug,
        "kind": "record",
        "state": "absent",
        "reason": None,
        "format": "json",
        "provenance": {"source": "NCBI Entrez", "accession": target.source.accession},
    }]
    for kind, is_baked in baked.items():
        rows.append({
            "slug": target.slug,
            "kind": kind,
            "state": "absent",
            "reason": None if is_baked else "Not included for this protein.",
            "format": "glb" if kind == "structure" else "json",
            "provenance": {},
        })
    return rows


# --- writing ----------------------------------------------------------------

_UPSERT_PROTEIN = """
insert into protein (
    slug, gene, uniprot, accession, slice_start, slice_end, slice_strand,
    transcript_id, protein_id, taxon_id, display, summary, chain_name,
    mature_peptides, residues, exons, chains, bridges,
    regions, disulfides, structure, provenance, resolver_version, catalog_order
) values (
    %(slug)s, %(gene)s, %(uniprot)s, %(accession)s, %(slice_start)s, %(slice_end)s,
    %(slice_strand)s, %(transcript_id)s, %(protein_id)s, %(taxon_id)s, %(display)s,
    %(summary)s, %(chain_name)s, %(mature_peptides)s, %(residues)s, %(exons)s,
    %(chains)s, %(bridges)s, %(regions)s, %(disulfides)s, %(structure)s,
    %(provenance)s, %(resolver_version)s, %(catalog_order)s
)
on conflict (slug) do update set
    gene = excluded.gene, uniprot = excluded.uniprot,
    accession = excluded.accession, slice_start = excluded.slice_start,
    slice_end = excluded.slice_end, slice_strand = excluded.slice_strand,
    transcript_id = excluded.transcript_id, protein_id = excluded.protein_id,
    display = excluded.display, summary = excluded.summary,
    chain_name = excluded.chain_name, mature_peptides = excluded.mature_peptides,
    residues = excluded.residues, exons = excluded.exons,
    chains = excluded.chains, bridges = excluded.bridges,
    regions = excluded.regions, disulfides = excluded.disulfides,
    structure = excluded.structure, provenance = excluded.provenance,
    resolver_version = excluded.resolver_version,
    catalog_order = excluded.catalog_order,
    resolved_at = now()
"""

_UPSERT_ALIAS = """
insert into protein_alias (slug, alias, kind) values (%s, %s, %s)
on conflict (slug, alias, kind) do nothing
"""

# Track rows are written only where there is not one already: this seeds the
# shape, and a later phase that has actually uploaded bytes owns the state.
_INSERT_TRACK = """
insert into protein_track (slug, kind, state, reason, format, provenance)
values (%(slug)s, %(kind)s, %(state)s, %(reason)s, %(format)s, %(provenance)s)
on conflict (slug, kind) do nothing
"""


def write(rows: list[tuple[dict, list[tuple[str, str]], list[dict]]], url: str) -> None:
    import psycopg
    from psycopg.types.json import Jsonb

    jsonb = ("regions", "disulfides", "structure", "provenance")
    with psycopg.connect(url) as conn:
        with conn.cursor() as cur:
            for protein, aliases, tracks in rows:
                payload = dict(protein)
                for key in jsonb:
                    payload[key] = Jsonb(payload[key]) if payload[key] is not None else None
                cur.execute(_UPSERT_PROTEIN, payload)
                cur.execute("delete from protein_alias where slug = %s", (protein["slug"],))
                for alias, kind in aliases:
                    cur.execute(_UPSERT_ALIAS, (protein["slug"], alias, kind))
                for track in tracks:
                    cur.execute(_INSERT_TRACK, {**track, "provenance": Jsonb(track["provenance"])})
        conn.commit()


_PROTEIN_COLUMNS = (
    "slug", "gene", "uniprot", "accession", "slice_start", "slice_end", "slice_strand",
    "transcript_id", "protein_id", "taxon_id", "display", "summary", "chain_name",
    "mature_peptides", "residues", "exons", "chains", "bridges", "regions", "disulfides",
    "structure", "provenance", "resolver_version", "catalog_order",
)


def check(rows: list[tuple[dict, list[tuple[str, str]], list[dict]]], url: str) -> int:
    """What a real run would change in the live rows, with nothing written.

    The proof that the curated file is the table the live rows were seeded from:
    once the prose left the Dart source, a seed from the file has to reproduce
    the database column for column and alias for alias, or the move lost
    something. `resolved_at` is the one column not compared; a run stamps it.
    """
    import psycopg

    differences: list[str] = []
    with psycopg.connect(url, prepare_threshold=None) as conn:
        for protein, aliases, _ in rows:
            slug = protein["slug"]
            live = conn.execute(
                f"select {', '.join(_PROTEIN_COLUMNS)} from protein where slug = %s",
                (slug,),
            ).fetchone()
            if live is None:
                differences.append(f"{slug}: no live row")
                continue
            for column, value in zip(_PROTEIN_COLUMNS, live):
                if protein[column] != value:
                    differences.append(
                        f"{slug}.{column}: live {value!r}, seed {protein[column]!r}")
            served = {
                (alias, kind) for alias, kind in conn.execute(
                    "select alias, kind from protein_alias where slug = %s", (slug,)
                ).fetchall()
            }
            if served != set(aliases):
                differences.append(
                    f"{slug}: aliases the seed adds {sorted(set(aliases) - served)}, "
                    f"drops {sorted(served - set(aliases))}")
    for line in differences:
        print(f"  {line}", file=sys.stderr)
    print(f"{len(rows)} proteins checked against the live rows, "
          f"{len(differences)} difference(s)", file=sys.stderr)
    return 1 if differences else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="print the rows as JSON and touch no database")
    parser.add_argument("--target", action="append", default=None,
                        help="one slug; repeatable. Default is every row.")
    parser.add_argument("--check", action="store_true",
                        help="diff the rows a run would write against the live ones; "
                             "write nothing")
    args = parser.parse_args()

    curated = curated_rows()
    wanted = set(args.target) if args.target else None

    missing = [t.slug for t in TARGETS if t.slug not in curated]
    if missing:
        raise SystemExit(f"no curated/catalog.json row for: {sorted(missing)}")
    extra = sorted(set(curated) - {t.slug for t in TARGETS})
    if extra:
        raise SystemExit(f"curated/catalog.json rows with no targets.py row: {extra}")

    rows = []
    for target in TARGETS:
        if wanted and target.slug not in wanted:
            continue
        row = curated[target.slug]
        rows.append((protein_row(target, row), _aliases(target, row),
                     track_rows(target, row)))

    if args.dry_run:
        print(json.dumps(
            [{"protein": p, "aliases": [list(a) for a in al], "tracks": t}
             for p, al, t in rows],
            indent=2, sort_keys=True,
        ))
        print(f"\n{len(rows)} proteins, "
              f"{sum(len(a) for _, a, _ in rows)} aliases, "
              f"{sum(len(t) for _, _, t in rows)} track rows", file=sys.stderr)
        return 0

    url = os.environ.get("DATABASE_URL")
    if not url:
        raise SystemExit("DATABASE_URL is not set. Use --dry-run to see the rows.")
    if args.check:
        return check(rows, url)
    write(rows, url)
    print(f"seeded {len(rows)} proteins", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
