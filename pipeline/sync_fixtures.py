"""Write curated rows into the app's test fixtures: their stored tracks and their `/catalog` row.

    .venv/bin/python pipeline/sync_fixtures.py --target phospholamban            # from storage
    .venv/bin/python pipeline/sync_fixtures.py --target phospholamban --from local
    .venv/bin/python pipeline/sync_fixtures.py --target phospholamban --check     # write nothing

The app's tests read a curated protein from three places (its CLAUDE.md):

| fixture | track |
|---|---|
| `test/fixtures/mock/gene_<gene>.json` | record |
| `test/fixtures/constraint/<slug>_esm_constraint.json` | constraint |
| `test/fixtures/impact/<slug>_avi.json` | impact |
| `test/fixtures/clinvar/<slug>_clinvar.json` | clinvar |
| `test/fixtures/impact_explanations/<slug>.json` | impact_explanations, where it has one |
| `test/fixtures/models/<slug>.glb` | structure (the `.glb` its scene is compiled from) |
| `test/features/zoom/fixtures/<slug>_locus.json` | locus |
| `test/shared/folding/fixtures/<slug>_folding.json` | folding |

and `test/fixtures/catalog.json`, a `/catalog` response, whose rows carry the
track kinds the app's tests have always read there (not `folding` or `locus`,
which they take from the files above).

`--from live` (the default) writes storage's bytes, each held to the digest
its row names, and the row `/catalog` serves. `--from local` writes what the
bakes left in `pipeline/data/` and the row the service will serve once it is
seeded and uploaded (`seed_catalog.protein_row`, shaped as `/catalog` shapes
it): the app's tests can then run on a row before anything is live, and the
live sync afterwards must leave every file as it found it.

It adds and never changes what the app has committed: a row or a file already
in the app's `HEAD` must come out the same, or nothing is written. A row it
wrote and nobody committed yet may be written again. `--check` compares and
writes nothing; it exits 1 on any difference.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

BACKEND = Path(__file__).resolve().parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import fetch_tracks  # noqa: E402
from pipeline.paths import CLIENT, DATA  # noqa: E402
from pipeline.targets import BY_SLUG, Target  # noqa: E402

CATALOG = "test/fixtures/catalog.json"


def fixture_paths(target: Target) -> Dict[str, str]:
    """Each track kind's fixture, relative to the app's checkout."""
    return {
        "record": f"test/fixtures/mock/gene_{target.gene.lower()}.json",
        "constraint": f"test/fixtures/constraint/{target.slug}_esm_constraint.json",
        "impact": f"test/fixtures/impact/{target.slug}_avi.json",
        "clinvar": f"test/fixtures/clinvar/{target.slug}_clinvar.json",
        "impact_explanations": f"test/fixtures/impact_explanations/{target.slug}.json",
        "structure": f"test/fixtures/models/{target.slug}.glb",
        "locus": f"test/features/zoom/fixtures/{target.slug}_locus.json",
        "folding": f"test/shared/folding/fixtures/{target.slug}_folding.json",
    }


def encode(catalog: dict) -> bytes:
    """The catalog fixture's own formatting."""
    return (json.dumps(catalog, indent=2, separators=(",", ": "), ensure_ascii=False)
            + "\n").encode()


def fixture_kinds(catalog: dict) -> List[str]:
    """The track kinds the fixture's rows carry, in their order."""
    kinds: List[str] = []
    for row in catalog["proteins"]:
        for kind in row.get("tracks") or {}:
            if kind not in kinds:
                kinds.append(kind)
    return kinds


def local_row(target: Target, kinds: List[str]) -> dict:
    """The row `/catalog` will serve once [target] is seeded and its tracks are
    uploaded, from the tables and the bakes here."""
    from pipeline import seed_catalog
    from pipeline.structure.experimental import served_row

    curated = seed_catalog.curated_rows()[target.slug]
    row = served_row(seed_catalog.protein_row(target, curated))
    baked = {
        "record": True,
        "constraint": target.scored,
        "impact": target.impact_scored,
        "clinvar": target.clinvar_available,
        "structure": target.structure is not None,
        "impact_explanations": curated["impact_explanations"],
    }
    row["tracks"] = {kind: "ready" if baked.get(kind) and (DATA / fetch_tracks.asset_path(
        kind, target)).exists() else "absent" for kind in kinds}
    return row


def live_row(service: str, target: Target, kinds: List[str]) -> dict:
    """The row `/catalog` serves, its tracks those the fixture's rows carry."""
    cursor: Optional[str] = None
    while True:
        url = f"{service.rstrip('/')}/catalog" + (f"?cursor={cursor}" if cursor else "")
        with urllib.request.urlopen(url, timeout=120) as answer:
            page = json.loads(answer.read())
        for row in page["proteins"]:
            if row["slug"] == target.slug:
                tracks = row.get("tracks") or {}
                row["tracks"] = {kind: tracks.get(kind, "absent") for kind in kinds}
                return row
        cursor = page.get("next")
        if not cursor:
            raise SystemExit(f"{target.slug}: /catalog does not list it; seed it first")


def payloads(target: Target, source: str, service: str) -> Dict[str, bytes]:
    """Each kind's bytes: storage's, held to their rows' digests, or the bakes'."""
    found: Dict[str, bytes] = {}
    if source == "local":
        for kind in fetch_tracks.KINDS:
            path = DATA / fetch_tracks.asset_path(kind, target)
            if path.exists():
                found[kind] = path.read_bytes()
        return found
    for kind, (url, digest) in fetch_tracks.stored_objects(service, target).items():
        with urllib.request.urlopen(url, timeout=300) as answer:
            blob = answer.read()
        if hashlib.sha256(blob).hexdigest() != digest:
            raise SystemExit(f"{target.slug} {kind}: the object is not the one its row names")
        found[kind] = blob
    return found


class App:
    """The app's checkout: its working tree, and what its `HEAD` has committed."""

    def __init__(self, root: Path):
        self.root = root

    def committed(self, path: str) -> Optional[bytes]:
        ran = subprocess.run(["git", "-C", str(self.root), "show", f"HEAD:{path}"],
                             capture_output=True)
        return ran.stdout if ran.returncode == 0 else None

    def read(self, path: str) -> Optional[bytes]:
        full = self.root / path
        return full.read_bytes() if full.exists() else None


def planned(app: App, targets: List[Target], source: str,
            service: str) -> Tuple[Dict[str, bytes], List[str]]:
    """What the fixtures would hold, path by path, and why it may not be written."""
    refused: List[str] = []
    writes: Dict[str, bytes] = {}
    catalog = json.loads(app.read(CATALOG))
    committed = {row["slug"]: row for row in json.loads(app.committed(CATALOG))["proteins"]}
    kinds = fixture_kinds({"proteins": list(committed.values())})
    rows = {row["slug"]: row for row in catalog["proteins"]}
    for target in targets:
        if target.slug not in BY_SLUG or target.structure is None:
            refused.append(f"{target.slug}: not a curated row with a fold")
            continue
        row = (local_row(target, kinds) if source == "local"
               else live_row(service, target, kinds))
        if target.slug in committed and committed[target.slug] != row:
            refused.append(f"{target.slug}: the app has committed another row for it")
        rows[target.slug] = row
        held = payloads(target, source, service)
        for kind, path in fixture_paths(target).items():
            if kind not in held:
                continue
            before = app.committed(path)
            if before is not None and before != held[kind]:
                refused.append(f"{path}: the app has committed other bytes")
            writes[path] = held[kind]
    order = [row["slug"] for row in catalog["proteins"]]
    order += [slug for slug in rows if slug not in order]
    writes[CATALOG] = encode({**catalog, "proteins": [rows[slug] for slug in order]})
    return writes, refused


def main(arguments: List[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--target", action="append", required=True)
    parser.add_argument("--from", dest="source", choices=("live", "local"), default="live")
    parser.add_argument("--app", type=Path, default=CLIENT)
    parser.add_argument("--service", default=fetch_tracks.SERVICE)
    parser.add_argument("--check", action="store_true", help="compare; write nothing")
    args = parser.parse_args(arguments)
    unknown = sorted(set(args.target) - set(BY_SLUG))
    if unknown:
        parser.error(f"no such target: {unknown}")
    app = App(args.app)
    if app.read(CATALOG) is None or app.committed(CATALOG) is None:
        raise SystemExit(f"{args.app} is not the app's checkout")
    writes, refused = planned(app, [BY_SLUG[s] for s in args.target], args.source, args.service)
    for problem in refused:
        print(f"  refused: {problem}", file=sys.stderr)
    if refused:
        print("nothing written", file=sys.stderr)
        return 1
    differ = [path for path, blob in writes.items() if app.read(path) != blob]
    if args.check:
        for path in differ:
            print(f"  differs: {path}", file=sys.stderr)
        print(f"{len(writes)} fixture(s) compared, {len(differ)} differ", file=sys.stderr)
        return 1 if differ else 0
    for path in differ:
        full = app.root / path
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_bytes(writes[path])
    print(f"{len(writes)} fixture(s) for {', '.join(args.target)} from {args.source}: "
          f"{len(differ)} written, {len(writes) - len(differ)} already so", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
