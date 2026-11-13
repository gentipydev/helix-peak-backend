"""Make the locus track, the zoom's way in, for proteins built before the
resolver made it.

    .venv/bin/python scripts/queue_locus.py --dry-run B2M PFN1   # what it would write; writes nothing
    .venv/bin/python scripts/queue_locus.py B2M PFN1             # make and write them
    .venv/bin/python scripts/queue_locus.py --dry-run --all      # every protein built on demand

From the repository root: the rows are reached through `.env`'s DATABASE_URL,
and storage through its SUPABASE_URL and SUPABASE_SERVICE_KEY. A protein's
locus is a pure function of its index row and the pinned snapshot of the
locus bake's inputs (`pipeline/locus/inputs.py`, built by `snapshot.py`), so
nothing is queued: each is made at once, exactly as the resolver makes one for
a protein it builds (`worker.make_locus`), the same bytes. Its object is stored
as the uploader names every track, `tracks/locus/<slug>.<sha12>.json`, and its
row written ready in a transaction of its own (`store.write_locus`), or
refused, with the reader's sentence, for a gene no zoom can be drawn for. Only
a protein resolved on demand is touched (`catalog_order` null): a curated one
keeps the locus baked by hand. A protein whose index row is gone is said, and
left as it is.

A protein opened while this runs has its tracks as they were when opened, and
offers the zoom from its next open.
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from app.config import settings  # noqa: E402
from pipeline.locus import bake_locus  # noqa: E402
from pipeline.locus import inputs as locus_inputs  # noqa: E402
from pipeline.resolver import store, worker  # noqa: E402
from pipeline.resolver.resolve import target_of  # noqa: E402


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


def state(conn, slug: str) -> str:
    row = conn.execute("select state from protein_track where slug = %s and kind = 'locus'",
                       (slug,)).fetchone()
    return row[0] if row else "absent"


def main(arguments: list) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("genes", nargs="*")
    parser.add_argument("--all", action="store_true", help="every protein built on demand")
    parser.add_argument("--dry-run", action="store_true", help="say what would be written")
    args = parser.parse_args(arguments)
    if bool(args.genes) == args.all:
        parser.error("name the genes, or say --all")

    snapshot = locus_inputs.load()
    print(f"locus inputs {snapshot.sha256[:12]} ({snapshot.created})", file=sys.stderr)
    storage = None
    if not args.dry_run:
        url, key = os.environ.get("SUPABASE_URL"), os.environ.get("SUPABASE_SERVICE_KEY")
        if not (url and key):
            from dotenv import dotenv_values

            found = dotenv_values(os.path.join(os.path.dirname(__file__), "..", ".env"))
            url, key = found.get("SUPABASE_URL"), found.get("SUPABASE_SERVICE_KEY")
        if not (url and key):
            raise SystemExit("SUPABASE_URL and SUPABASE_SERVICE_KEY must be set.")
        storage = store.TrackStorage(url, key)

    with store.connect(settings.database_url) as conn:
        # A dry run reads in one transaction that may do nothing more: `set
        # transaction read only` holds through any pooler, where a session's
        # setting may not.
        reading = conn.transaction() if args.dry_run else contextlib.nullcontext()
        with reading:
            if args.dry_run:
                conn.execute("set transaction read only")
            release = store.index_release(conn) or {}
            slugs = built(conn, args.genes)
            missing = _each(conn, storage, slugs, snapshot, release, args.dry_run)
    # A gene named with no protein built on demand is a mistake in the asking.
    asked_for = args.genes and len(slugs) < len(args.genes)
    return 1 if asked_for or missing else 0


def _each(conn, storage, slugs: list, snapshot, release: dict, dry_run: bool) -> int:
    """Make each protein's locus, or say what would be made; how many could not be."""
    missing = 0
    for slug in slugs:
        protein = store.protein(conn, slug)
        target = target_of(protein)
        row = store.index_row(conn, protein["uniprot"], protein["gene"])
        before = state(conn, slug)
        if row is None:
            print(f"{slug}: the index no longer holds {protein['uniprot']} "
                  f"({protein['gene']}); left {before}")
            missing += 1
            continue
        if dry_run:
            try:
                track = bake_locus.bake(
                    dataclasses.asdict(row), snapshot, slug=slug,
                    uniprot_source={"release": release.get("uniprot")},
                    mane_release=release.get("mane"))
                payload = bake_locus.encode(track)
                path = track["path"]
                print(f"{slug}: {before}; would write ready: {track['locus']} on "
                      f"{track['sequence']}, {path['tissue']} / {path['cell_type']} "
                      f"({path['tissue_from']}, {path['cell_from']}), "
                      f"{len(payload):,} bytes")
            except bake_locus.LocusRefused as exc:
                print(f"{slug}: {before}; would write refused: {exc}")
            continue
        made, refused, said = worker.make_locus(conn, storage, row, target)
        if made is None and refused is None:
            print(f"{slug}: {before}; nothing written, {said}")
            missing += 1
            continue
        with conn.transaction():
            store.write_locus(conn, slug, made, refused)
        print(f"{slug}: {before} -> {state(conn, slug)}; {said}")
    return missing


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
