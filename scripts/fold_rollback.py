"""Put back the fold a protein built on demand drew before, whole.

    .venv/bin/python scripts/fold_rollback.py TTR --list          # its folds, newest first; reads only
    .venv/bin/python scripts/fold_rollback.py TTR --dry-run       # what would be put back; reads only
    .venv/bin/python scripts/fold_rollback.py TTR                 # the fold before the one drawn now
    .venv/bin/python scripts/fold_rollback.py TTR --to 7          # that generation, back or forward

From the repository root: the rows are reached through `.env`'s DATABASE_URL.
Every fold the worker publishes is a generation (`fold_generation`, 0012), and
so is the one a protein drew before its first: an AlphaFold model, or the
sentence why it had none. Putting one back is one transaction
(`store.roll_back_fold`): the row's words, the structure track and the folding
track as that generation wrote them (an AlphaFold model has no folding track,
so it reads absent again), and that generation the one drawn. No object is
deleted when a fold is published, so every generation's objects are still
where its rows point.

Refused while a structure bake for the protein is queued or running, which
would publish over what is put back: let it finish, or stop the worker, first.
An app reads the fold put back the next time the protein is opened.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from app.config import settings  # noqa: E402
from pipeline.resolver import store  # noqa: E402


def said(generation: dict) -> str:
    """One generation, on one line."""
    track = generation["structure_track"] or {}
    provenance = track.get("provenance") or {}
    what = (f"{provenance.get('source') or generation['source']} {provenance.get('entry') or ''}"
            .strip() if generation["state"] == "ready" else f"refused: {track.get('reason')}")
    animated = ", animated" if generation["folding_track"] else ""
    return (f"{'*' if generation['active'] else ' '} {generation['id']:>4}  "
            f"{generation['created_at'].astimezone():%Y-%m-%d %H:%M}  {generation['source']:<9} "
            f"{what}{animated}")


def main(arguments: list) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("gene")
    parser.add_argument("--to", type=int, default=None, help="the generation to put back")
    parser.add_argument("--list", action="store_true", help="its folds, newest first")
    parser.add_argument("--dry-run", action="store_true", help="say what would be put back")
    args = parser.parse_args(arguments)

    with store.connect(settings.database_url) as conn:
        if not store.folds_ready(conn):
            print("The database records no fold yet: migrations/0012_fold_generation.sql is "
                  "not applied.", file=sys.stderr)
            return 1
        found = conn.execute(
            "select slug from protein where lower(gene) = lower(%s) and catalog_order is null "
            "and resolver_version > 0", (args.gene,)).fetchone()
        if found is None:
            print(f"{args.gene}: no protein built on demand", file=sys.stderr)
            return 1
        slug = found[0]
        every = store.generations(conn, slug)
        if args.list or args.dry_run:
            for generation in every:
                print(said(generation))
            if not every:
                print(f"{slug}: no fold published since 0012")
            if args.list:
                return 0
            active = next((g for g in every if g["active"]), None)
            if args.to is not None:
                chosen = next((g for g in every if g["id"] == args.to), None)
            else:
                chosen = next((g for g in every if active is None or g["id"] < active["id"]),
                              None)
            if chosen is None or chosen is active:
                print(f"{slug}: nothing to put back")
                return 1
            print(f"{slug}: would put back generation {chosen['id']} ({chosen['source']})")
            return 0
        try:
            done = store.roll_back_fold(conn, slug, args.to)
        except RuntimeError as exc:
            print(f"{slug}: {exc}", file=sys.stderr)
            return 1
        print(f"{slug}: generation {done['to']} ({done['source']}, {done['state']}) put back"
              + (f" over {done['from']}" if done["from"] is not None else ""))
        for generation in store.generations(conn, slug):
            print(said(generation))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
