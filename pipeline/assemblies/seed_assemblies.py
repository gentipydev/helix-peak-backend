"""Seed the `assembly` table from `assemblies.py`, and hold the catalog to twenty.

    python3 pipeline/assemblies/seed_assemblies.py --dry-run    # print the rows
    set -a && . ./.env && set +a
    .venv/Scripts/python pipeline/assemblies/seed_assemblies.py --check   # diff, and assert
    .venv/Scripts/python pipeline/assemblies/seed_assemblies.py           # write, and assert

Writes `assembly` and `assembly_track` (migrations/0006_assemblies.sql) and
nothing else: never `protein`, `protein_alias` or `protein_track`, which is
what keeps an assembly out of `/catalog`. A track row is inserted `absent`
only where there is none, so a seed never undoes an upload.

After a write, and in `--check`, it asserts both of the things an assembly
must not change: `seed_catalog.py --check` finds no difference in the
twenty's live rows, and `GET /catalog` on the service returns exactly twenty.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import seed_catalog  # noqa: E402
from pipeline.assemblies.assemblies import ASSEMBLIES, Assembly  # noqa: E402
from pipeline.fetch_tracks import SERVICE  # noqa: E402
from pipeline.targets import TARGETS  # noqa: E402

KINDS = ("morph",)


def assembly_row(assembly: Assembly) -> dict:
    return {
        "slug": assembly.slug,
        "display": assembly.display,
        "subunits": [{"node": u.node, "gene": u.gene, "uniprot": u.uniprot,
                      "pdb_chain": u.pdb_chain} for u in assembly.subunits],
        "states": [{"name": s.name, "pdb": s.pdb, "ligand": s.ligand, "biomt": s.biomt}
                   for s in assembly.states],
        "provenance": {"source": "RCSB PDB",
                       "entries": [s.pdb for s in assembly.states],
                       "superposed_on": list(assembly.superpose_on),
                       "seeded_by": "pipeline/assemblies/seed_assemblies.py"},
    }


def track_rows(assembly: Assembly) -> list[dict]:
    return [{"slug": assembly.slug, "kind": kind, "state": "absent", "reason": None,
             "format": "json", "provenance": {}} for kind in KINDS]


# The only tables written: an assembly's own.
_UPSERT = """
insert into assembly (slug, display, subunits, states, provenance)
values (%(slug)s, %(display)s, %(subunits)s, %(states)s, %(provenance)s)
on conflict (slug) do update set
    display = excluded.display, subunits = excluded.subunits,
    states = excluded.states, provenance = excluded.provenance, updated_at = now()
"""

_INSERT_TRACK = """
insert into assembly_track (slug, kind, state, reason, format, provenance)
values (%(slug)s, %(kind)s, %(state)s, %(reason)s, %(format)s, %(provenance)s)
on conflict (slug, kind) do nothing
"""

_COLUMNS = ("slug", "display", "subunits", "states", "provenance")


def write(url: str) -> None:
    import psycopg
    from psycopg.types.json import Jsonb

    with psycopg.connect(url) as conn:
        with conn.cursor() as cur:
            for assembly in ASSEMBLIES:
                row = assembly_row(assembly)
                cur.execute(_UPSERT, {**row, **{k: Jsonb(row[k]) for k in
                                                ("subunits", "states", "provenance")}})
                for track in track_rows(assembly):
                    cur.execute(_INSERT_TRACK, {**track, "provenance": Jsonb(track["provenance"])})
        conn.commit()


def check(url: str) -> list[str]:
    """What a write would change in the live `assembly` rows."""
    import psycopg

    differences: list[str] = []
    with psycopg.connect(url, prepare_threshold=None) as conn:
        exists = conn.execute("select to_regclass('public.assembly')").fetchone()[0]
        if exists is None:
            return ["the assembly table is not there: apply migrations/0006_assemblies.sql"]
        for assembly in ASSEMBLIES:
            row = assembly_row(assembly)
            live = conn.execute(f"select {', '.join(_COLUMNS)} from assembly where slug = %s",
                                (assembly.slug,)).fetchone()
            if live is None:
                differences.append(f"{assembly.slug}: no live row")
                continue
            for column, value in zip(_COLUMNS, live):
                if row[column] != value:
                    differences.append(f"{assembly.slug}.{column}: live {value!r}, seed {row[column]!r}")
    return differences


def catalog_count(service: str) -> int:
    """How many proteins `GET /catalog` serves, every page of it."""
    count, cursor = 0, None
    while True:
        url = f"{service.rstrip('/')}/catalog" + (f"?cursor={cursor}" if cursor else "")
        with urllib.request.urlopen(url, timeout=120) as response:
            page = json.loads(response.read())
        count += len(page.get("proteins", []))
        cursor = page.get("next")
        if not cursor:
            return count


def catalog_holds(url: str, service: str) -> list[str]:
    """The two things an assembly must not change, asserted."""
    found = []
    rows = [(seed_catalog.protein_row(t, seed_catalog.curated_rows()[t.slug]),
             seed_catalog._aliases(t, seed_catalog.curated_rows()[t.slug]),
             seed_catalog.track_rows(t, seed_catalog.curated_rows()[t.slug]))
            for t in TARGETS]
    if seed_catalog.check(rows, url) != 0:
        found.append("seed_catalog.py --check finds the twenty's live rows changed")
    served = catalog_count(service)
    if served != len(TARGETS) or len(TARGETS) != 20:
        found.append(f"GET /catalog serves {served} proteins, not twenty")
    return found


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--service", default=SERVICE)
    args = parser.parse_args()

    if args.dry_run:
        print(json.dumps([{"assembly": assembly_row(a), "tracks": track_rows(a)}
                          for a in ASSEMBLIES], indent=2))
        print(f"\n{len(ASSEMBLIES)} assembly row(s)", file=sys.stderr)
        return 0
    url = os.environ.get("DATABASE_URL")
    if not url:
        raise SystemExit("DATABASE_URL is not set. Use --dry-run to see the rows.")
    problems = check(url) if args.check else []
    if not args.check:
        write(url)
        print(f"seeded {len(ASSEMBLIES)} assembly row(s)", file=sys.stderr)
    problems += catalog_holds(url, args.service)
    for problem in problems:
        print(f"  {problem}", file=sys.stderr)
    print(f"{len(problems)} problem(s)", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
