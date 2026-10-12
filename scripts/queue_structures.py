"""Queue the fold of proteins built on demand to be made again: an experimental
entry where one passes, with the fold animation, and AlphaFold DB's model only
where a whole search finds none.

    .venv/bin/python scripts/queue_structures.py --dry-run --all       # what it would queue; reads only
    .venv/bin/python scripts/queue_structures.py TTR PRL               # queue them
    .venv/bin/python scripts/queue_structures.py --all                 # every protein built on demand
    .venv/bin/python scripts/queue_structures.py --refresh TTR         # search again, past what stands

From the repository root: the rows are reached through `.env`'s DATABASE_URL.
Each protein is one transaction (`store.queue_structure`): its structure bake
queued, as the worker queues one for a new build, and a structure track with
no model yet said to be on its way (`pending`). A model a protein draws stays
drawn, and its fold page as it is, until the worker publishes the one that
replaces it, the tracks and the row's words together (`store.publish_fold`);
where nothing better can be drawn, it keeps the one it has. Its build is never
held for it, as a new build's is for its first fold.

A protein whose experimental fold stands for its record, regions and rule is
left as it is (`kept`), and one with a bake under way already (`waiting`).
`--refresh` asks for a search that starts afresh: past a search that found
nothing not long ago, past one cut short, and past a fold that stands. One of
the twenty is refused: it is baked by hand.

The Mac's worker takes them at its next cycle, one protein at a time, before
any scoring. Every fold it publishes is kept as a generation, and
`scripts/fold_rollback.py` puts back the one before.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from app.config import settings  # noqa: E402
from pipeline.resolver import store  # noqa: E402
from pipeline.structure.pdb_policy import POLICY_VERSION  # noqa: E402


def built(conn, genes: list) -> list:
    """The slugs of these genes' proteins built on demand, or of all of them;
    a gene with none is said, and left out."""
    rows = conn.execute(
        "select slug, gene from protein where catalog_order is null and resolver_version > 0 "
        "order by resolved_at").fetchall()
    if not genes:
        return [slug for slug, _ in rows]
    by_gene = {gene.lower(): slug for slug, gene in rows}
    found = []
    for gene in genes:
        if gene.lower() in by_gene:
            found.append(by_gene[gene.lower()])
        else:
            print(f"{gene}: no protein built on demand", file=sys.stderr)
    return found


def drawn(conn, slug: str) -> str:
    """What a protein's fold page draws now, in a few words."""
    track = conn.execute(
        "select state, reason, provenance ->> 'source', provenance ->> 'entry' "
        "from protein_track where slug = %s and kind = 'structure'", (slug,)).fetchone()
    folding = conn.execute(
        "select state from protein_track where slug = %s and kind = 'folding'", (slug,)).fetchone()
    if track is None:
        return "no model asked for"
    state, reason, source, entry = track
    if state != "ready":
        return f"structure {state}" + (f" ({reason})" if reason else "")
    animated = folding is not None and folding[0] == "ready"
    return f"{source or 'a model'} {entry or ''}".strip() + (", animated" if animated else "")


def main(arguments: list) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("genes", nargs="*")
    parser.add_argument("--all", action="store_true", help="every protein built on demand")
    parser.add_argument("--refresh", action="store_true",
                        help="search afresh, past what earlier searches decided")
    parser.add_argument("--dry-run", action="store_true", help="say what would be queued")
    args = parser.parse_args(arguments)
    if bool(args.genes) == args.all:
        parser.error("name the genes, or say --all")

    with store.connect(settings.database_url) as conn:
        if not store.folds_ready(conn):
            print("The database cannot record a fold yet: apply "
                  "migrations/0012_fold_generation.sql first.", file=sys.stderr)
            return 1
        if args.dry_run:
            # Every statement after this one may read and nothing more.
            conn.execute("set session characteristics as transaction read only")
        slugs = built(conn, args.genes)
        refused: list = []
        for slug in slugs:
            now = drawn(conn, slug)
            if args.dry_run:
                prior = store.fold_prior(conn, slug)
                record = store.ready_track(conn, slug, "record")
                waiting = conn.execute(
                    "select 1 from bake_job where slug = %s and kind = 'structure' "
                    "and state in ('queued', 'running')", (slug,)).fetchone() is not None
                if waiting:
                    would = "nothing: a bake is under way"
                elif record is None:
                    would = "nothing: it has no ready record"
                elif not args.refresh and store.pinned(prior, record["sha256"],
                                                        store.regions_of(conn, slug),
                                                        POLICY_VERSION):
                    would = "nothing: its experimental fold stands"
                else:
                    would = "a structure bake" + (" (afresh)" if args.refresh else "")
                    if not args.refresh and prior["none"] is not None:
                        would += (f", which will not search the PDB again before "
                                  f"{prior['none']['decided_at'][:10]} + 28 days")
                print(f"{slug}: draws {now}; would queue {would}")
                continue
            try:
                done = store.queue_structure(conn, slug, refresh=args.refresh,
                                             policy=POLICY_VERSION)
            except RuntimeError as exc:
                print(f"{slug}: draws {now}; not queued: {exc}", file=sys.stderr)
                refused.append(slug)
                continue
            print(f"{slug}: draws {now}; {done}")
    # A gene named with no protein built on demand is a mistake in the asking.
    return 1 if refused or (args.genes and len(slugs) < len(args.genes)) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
