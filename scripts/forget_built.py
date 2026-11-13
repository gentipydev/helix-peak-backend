"""Delete proteins built on demand, whole: their rows, their requests and their objects.

    .venv/bin/python scripts/forget_built.py --all --dry-run    # what would go; reads only
    .venv/bin/python scripts/forget_built.py --all              # every protein built on demand
    .venv/bin/python scripts/forget_built.py TTR PRL            # these, by gene

From the repository root: the rows are reached through `.env`'s DATABASE_URL, and storage
through its SUPABASE_URL and SUPABASE_SERVICE_KEY (the environment's, where it sets them). A
protein built on demand is a `protein` row with a null `catalog_order`. A curated protein is
never touched: naming one is refused before anything is written.

What goes, for each:
- its `protein` row, and with it (`on delete cascade`) its aliases, track rows, bakes, fold
  selections and every fold generation;
- every `resolve_request` for its gene, so the gene reads as never asked for;
- every object of its in storage: those its rows name (each track's, a structure's `.glb`,
  every generation's model, folding track and audit), and any other named
  `<folder>/<slug>.<sha12>.<ext>` in the `tracks` and `models` buckets, which a re-bake
  leaves behind.

Refused while a request for one of them, or a bake of one, is queued or running: let it end,
or stop it, first. The rows go in one transaction, and only then the objects, so nothing is
deleted while a row still points at it. A protein deleted can be asked for again; whether it is
offered to be built is the index's to say (`app/protein_index.py`).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.request
from typing import Callable, Dict, Iterable, List, Optional, Set, Tuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

BUCKETS = ("tracks", "models")

# What each protein holds in each table the delete reaches, in the order it is said.
_TABLES = (
    ("protein_alias", "aliases"),
    ("protein_track", "tracks"),
    ("bake_job", "bakes"),
    ("fold_selection", "fold selections"),
    ("fold_generation", "fold generations"),
)

Object = Tuple[str, str]  # (bucket, path)


class Refused(Exception):
    """Nothing was written: the sentence says why."""


def owned(slug: str, path: str) -> bool:
    """Whether `path` is an object of `slug`'s: `<folder>/<slug>.<sha12>.<ext>`, exactly.

    The name is the protein's own, never a longer one that starts with it
    (`b2m` is not `b2mx`), and sits one folder deep, as every bake writes it."""
    return re.fullmatch(r"[a-z0-9_-]+/" + re.escape(slug) + r"\.[0-9a-f]{12}\.[a-z]+", path) is not None


def _named(holder: Optional[dict], key: str = "object_path") -> List[Object]:
    """The objects one row of a track, as stored or as a generation keeps it, names: its
    own and, for a structure, the `.glb` it was compiled from, in the same bucket."""
    if not holder or not holder.get("bucket"):
        return []
    found = []
    if holder.get(key):
        found.append((holder["bucket"], holder[key]))
    glb = ((holder.get("provenance") or {}).get("glb") or {}).get("path")
    if glb:
        found.append((holder["bucket"], glb))
    return found


def objects_of(slug: str, tracks: Iterable[dict], generations: Iterable[dict]) -> Set[Object]:
    """Every object `slug`'s rows name. Raises Refused on one that is not its own, which
    would be another protein's object."""
    found = set()  # type: Set[Object]
    for track in tracks:
        found.update(_named(track))
    for generation in generations:
        found.update(_named(generation.get("structure_track")))
        found.update(_named(generation.get("folding_track")))
        found.update(_named(generation.get("audit"), key="path"))
    strangers = sorted("{}/{}".format(*o) for o in found if not owned(slug, o[1]))
    if strangers:
        raise Refused("{}'s rows name objects that are not its own: {}".format(
            slug, ", ".join(strangers)))
    return found


class Storage:
    """Lists and deletes objects in Supabase storage, with the service key."""

    def __init__(self, url: str, key: str):
        from pipeline import upload_tracks

        self._uploads = upload_tracks.Storage(url, key)
        self._folders = {}  # type: Dict[str, List[str]]

    def _list(self, bucket: str, prefix: str, search: str = "") -> List[dict]:
        found, offset = [], 0  # type: List[dict], int
        while True:
            body = json.dumps({"prefix": prefix, "search": search, "limit": 100, "offset": offset,
                               "sortBy": {"column": "name", "order": "asc"}}).encode()
            request = urllib.request.Request(
                "{}/object/list/{}".format(self._uploads.base, bucket), data=body, method="POST",
                headers={**self._uploads.headers, "Content-Type": "application/json"})
            with urllib.request.urlopen(request, timeout=120) as response:
                page = json.loads(response.read())
            found.extend(page)
            if len(page) < 100:
                return found
            offset += 100

    def folders(self, bucket: str) -> List[str]:
        if bucket not in self._folders:
            self._folders[bucket] = [o["name"] for o in self._list(bucket, "") if o.get("id") is None]
        return self._folders[bucket]

    def objects(self, bucket: str, folder: str, start: str) -> List[str]:
        """The names in `folder` of `bucket` that start with `start`, as paths."""
        return ["{}/{}".format(folder, o["name"])
                for o in self._list(bucket, folder + "/", start) if o.get("id") is not None]

    def delete(self, bucket: str, paths: List[str]) -> None:
        self._uploads.delete(bucket, paths)


def stored(storage, slug: str) -> Set[Object]:
    """Every object in storage named as `slug`'s, whatever row names it or none."""
    found = set()  # type: Set[Object]
    for bucket in BUCKETS:
        for folder in storage.folders(bucket):
            found.update((bucket, path) for path in storage.objects(bucket, folder, slug + ".")
                         if owned(slug, path))
    return found


def chosen(conn, genes: List[str], every: bool) -> List[Tuple[str, str]]:
    """(slug, gene) of each protein to delete. Raises Refused on a gene named that has no
    protein, or whose protein is curated, before anything is written."""
    if every:
        return [tuple(row) for row in conn.execute(
            "select slug, gene from protein where catalog_order is null order by slug")]
    rows = conn.execute(
        "select slug, gene, catalog_order from protein where lower(gene) = any(%s)",
        ([g.lower() for g in genes],)).fetchall()
    by_gene = {gene.lower(): (slug, gene, order) for slug, gene, order in rows}
    unknown = [g for g in genes if g.lower() not in by_gene]
    curated = [g for g in genes if g.lower() in by_gene and by_gene[g.lower()][2] is not None]
    if unknown or curated:
        raise Refused("; ".join(
            ["no protein for {}".format(", ".join(unknown))] * bool(unknown)
            + ["curated, never deleted here: {}".format(", ".join(curated))] * bool(curated)))
    return sorted({(slug, gene) for slug, gene, _order in by_gene.values()})


def counts(conn, slug: str, gene: str) -> Dict[str, int]:
    said = {name: conn.execute("select count(*) from {} where slug = %s".format(table),
                               (slug,)).fetchone()[0]
            for table, name in _TABLES}
    said["requests"] = conn.execute(
        "select count(*) from resolve_request where lower(gene) = lower(%s)", (gene,)).fetchone()[0]
    return said


def named(conn, slug: str) -> Set[Object]:
    """Every object `slug`'s rows name now."""
    tracks = [{"bucket": bucket, "object_path": path, "provenance": provenance}
              for bucket, path, provenance in conn.execute(
                  "select bucket, object_path, provenance from protein_track where slug = %s",
                  (slug,))]
    generations = [{"structure_track": s, "folding_track": f, "audit": a}
                   for s, f, a in conn.execute(
                       "select structure_track, folding_track, audit from fold_generation "
                       "where slug = %s", (slug,))]
    return objects_of(slug, tracks, generations)


_BUSY = """
select 'request for ' || gene || ' ' || state from resolve_request
where lower(gene) = any(%(genes)s) and state in ('queued', 'running')
union all
select kind || ' bake of ' || slug || ' ' || state from bake_job
where slug = any(%(slugs)s) and state in ('queued', 'running')
"""


def forget(conn, storage, genes: List[str], every: bool, dry_run: bool,
           say: Callable[[str], None] = print) -> int:
    """Delete what `genes` (or `every`) names, saying each step. Returns how many proteins
    were deleted; a dry run deletes none."""
    proteins = chosen(conn, genes, every)
    if not proteins:
        say("No protein built on demand.")
        return 0
    objects = {}  # type: Dict[str, Set[Object]]
    for slug, gene in proteins:
        objects[slug] = named(conn, slug) | stored(storage, slug)
        held = counts(conn, slug, gene)
        say("{} ({}): {}; {} objects".format(
            slug, gene, ", ".join("{} {}".format(n, name) for name, n in held.items()),
            len(objects[slug])))
        for bucket, path in sorted(objects[slug]):
            say("    {}/{}".format(bucket, path))
    if dry_run:
        say("Dry run: nothing deleted ({} proteins, {} objects).".format(
            len(proteins), sum(len(o) for o in objects.values())))
        return 0

    slugs = [slug for slug, _gene in proteins]
    lowered = [gene.lower() for _slug, gene in proteins]
    with conn.transaction():
        conn.execute("select slug from protein where slug = any(%s) and catalog_order is null "
                     "for update", (slugs,))
        busy = [row[0] for row in conn.execute(_BUSY, {"genes": lowered, "slugs": slugs})]
        if busy:
            raise Refused("under way: {}".format("; ".join(busy)))
        for slug in slugs:
            # Read again under the lock: a bake that ended since names objects too.
            objects[slug] |= named(conn, slug)
        conn.execute("delete from resolve_request where lower(gene) = any(%s)", (lowered,))
        gone = conn.execute("delete from protein where slug = any(%s) and catalog_order is null",
                            (slugs,)).rowcount
        if gone != len(slugs):
            raise Refused("{} of {} protein rows deleted; nothing kept".format(gone, len(slugs)))
    say("Deleted {} proteins' rows.".format(gone))

    # Only now that no row points at them; and whatever storage holds under their names,
    # listed again.
    for slug in slugs:
        objects[slug] |= stored(storage, slug)
    for bucket in BUCKETS:
        paths = sorted(p for every_one in objects.values() for b, p in every_one if b == bucket)
        for start in range(0, len(paths), 100):
            storage.delete(bucket, paths[start:start + 100])
    left = sorted("{}/{}".format(*o) for slug in slugs for o in stored(storage, slug))
    if left:
        say("Still in storage: {}".format(", ".join(left)))
        return -1
    say("Deleted {} objects.".format(sum(len(o) for o in objects.values())))
    return gone


def _credentials() -> Tuple[str, str, str]:
    from dotenv import dotenv_values

    dotenv = dotenv_values(os.path.join(os.path.dirname(__file__), "..", ".env"))
    said = [os.environ.get(key) or dotenv.get(key)
            for key in ("DATABASE_URL", "SUPABASE_URL", "SUPABASE_SERVICE_KEY")]
    missing = [key for key, value in zip(
        ("DATABASE_URL", "SUPABASE_URL", "SUPABASE_SERVICE_KEY"), said) if not value]
    if missing:
        raise SystemExit("{} must be set, in the environment or .env".format(", ".join(missing)))
    return said[0], said[1], said[2]


def main(arguments: list) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("genes", nargs="*", help="the genes whose built proteins go")
    parser.add_argument("--all", action="store_true", help="every protein built on demand")
    parser.add_argument("--dry-run", action="store_true", help="say what would go; read only")
    args = parser.parse_args(arguments)
    if bool(args.genes) == args.all:
        parser.error("name genes, or --all, not both")

    from pipeline.resolver import store

    database, url, key = _credentials()
    with store.connect(database) as conn:
        try:
            done = forget(conn, Storage(url, key), args.genes, args.all, args.dry_run)
        except Refused as exc:
            print("Refused, nothing deleted: {}".format(exc), file=sys.stderr)
            return 1
    return 1 if done < 0 else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
