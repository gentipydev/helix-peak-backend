"""Queue AVI and ClinVar for proteins built before the worker made them.

    .venv/bin/python scripts/queue_evidence.py --dry-run TTR PRL   # what it would queue; reads only
    .venv/bin/python scripts/queue_evidence.py TTR PRL             # queue them
    .venv/bin/python scripts/queue_evidence.py --dry-run --all     # every protein built on demand

From the repository root: the rows are reached through `.env`'s DATABASE_URL.
Each protein is one transaction (`store.queue_evidence`): its tracks pending,
and a bake queued, for each of the two kinds not ready already, as the worker
queues them for a new build. One of the twenty is refused: it is baked by hand.
The Mac's worker takes them at its next cycle, AVI before ClinVar
(`store.WAITS_FOR`).

A protein opened while its evidence is pending is drawn without it until the
app restarts (`TrackClient` keeps a protein's rows), so this is best run when
nobody is walking the proteins it names.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from app.config import settings  # noqa: E402
from pipeline.resolver import store  # noqa: E402


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


def tracks(conn, slug: str) -> dict:
    return dict(conn.execute(
        "select kind, state from protein_track where slug = %s and kind = any(%s)",
        (slug, list(store.EVIDENCE))).fetchall())


def main(arguments: list) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("genes", nargs="*")
    parser.add_argument("--all", action="store_true", help="every protein built on demand")
    parser.add_argument("--dry-run", action="store_true", help="say what would be queued")
    args = parser.parse_args(arguments)
    if bool(args.genes) == args.all:
        parser.error("name the genes, or say --all")

    with store.connect(settings.database_url) as conn:
        if args.dry_run:
            # Every statement after this one may read and nothing more.
            conn.execute("set session characteristics as transaction read only")
        slugs = built(conn, args.genes)
        for slug in slugs:
            before = tracks(conn, slug)
            if args.dry_run:
                waiting = [kind for kind in store.EVIDENCE if before.get(kind) != "ready"]
                print(f"{slug}: {before or 'no evidence rows'}; would queue {waiting or 'nothing'}")
                continue
            queued = store.queue_evidence(conn, slug)
            print(f"{slug}: {before or 'no evidence rows'}; queued {queued or 'nothing'}, "
                  f"now {tracks(conn, slug)}")
    # A gene named with no protein built on demand is a mistake in the asking.
    return 1 if args.genes and len(slugs) < len(args.genes) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
