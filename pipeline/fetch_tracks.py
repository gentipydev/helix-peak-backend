"""Bring the stored tracks back to disk, for the tools that read files.

Storage is where the twenty's tracks live since Phase 3 of HANDOFF-ONDEMAND.md:
the app fetches them from there and nothing ships them. Some tools still want
files -- `check_assets.py` reads every family, the ClinVar bake reads the record
and the AVI map, the constraint and impact tests hold stored tracks to their
gates -- so this writes storage's copy into `pipeline/data/`, in the `assets/...`
layout `targets.py` names, verified against the sha256 each row carries.

    python3 pipeline/fetch_tracks.py                       # all twenty, every family
    python3 pipeline/fetch_tracks.py --target insulin --kind constraint
    python3 pipeline/fetch_tracks.py --verify-client ../helix-peek

Reads the public service and the public buckets only; it needs no key. A file
already on disk with the right digest is not fetched again.

`--verify-client` is the other direction and writes nothing: every file under a
client checkout's six asset directories has to equal the stored object its row
names, byte for byte. It is the gate Phase 4 runs before the client stops
bundling those files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.paths import CLIENT, DATA  # noqa: E402
from pipeline.targets import TARGETS, Target  # noqa: E402

SERVICE = "https://helix-peak-backend.onrender.com"
KINDS = ("record", "constraint", "impact", "clinvar", "impact_explanations", "structure")


def asset_path(kind: str, target: Target) -> str:
    """Where a family's file sits, relative to the data directory or a client.

    The structure family names the `.fsceneb` the app draws; what the bake
    makes, and what `check_assets.py` reads, is the `.glb` it was compiled from,
    which the row names in its provenance. So that is the file kept here.
    """
    if kind == "record":
        return target.mock_asset
    if kind == "constraint":
        return target.constraint_asset
    if kind == "impact":
        return target.impact_asset
    if kind == "clinvar":
        return f"assets/clinvar/{target.slug}_clinvar.json"
    if kind == "impact_explanations":
        return f"assets/impact_explanations/{target.slug}.json"
    if kind == "structure":
        return target.structure_asset
    raise ValueError(f"unknown kind {kind!r}")


def _get(url: str, timeout: int) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return response.read()


def stored_objects(service: str, target: Target) -> dict[str, tuple[str, str]]:
    """Each ready family's (url, sha256), as the service names them.

    For structure, the `.glb` beside the container: its object path is in the
    row's provenance and it sits in the same bucket, so its URL is the
    container's with the object path swapped.
    """
    rows = json.loads(_get(f"{service.rstrip('/')}/protein/{target.slug}/tracks", 60))
    found: dict[str, tuple[str, str]] = {}
    for kind in KINDS:
        row = rows.get(kind) or {}
        if row.get("state") != "ready" or not row.get("url"):
            continue
        if kind == "structure":
            glb = (row.get("provenance") or {}).get("glb") or {}
            if not glb.get("path") or not glb.get("sha256"):
                continue
            prefix = row["url"][: -len(row["url"].split("/models/", 1)[1])]
            found[kind] = (prefix + glb["path"], glb["sha256"])
        else:
            found[kind] = (row["url"], row["sha256"])
    return found


def fetch(service: str, targets: list[Target], kinds: tuple[str, ...]) -> int:
    fetched = kept = failed = 0
    total = 0
    for target in targets:
        try:
            objects = stored_objects(service, target)
        except (urllib.error.URLError, OSError, ValueError) as exc:
            print(f"  {target.slug}: could not read its tracks: {exc}", file=sys.stderr)
            failed += 1
            continue
        for kind in kinds:
            if kind not in objects:
                continue
            url, digest = objects[kind]
            path = DATA / asset_path(kind, target)
            if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() == digest:
                kept += 1
                total += path.stat().st_size
                continue
            try:
                blob = _get(url, 300)
            except (urllib.error.URLError, OSError) as exc:
                print(f"  {target.slug} {kind}: {exc}", file=sys.stderr)
                failed += 1
                continue
            if hashlib.sha256(blob).hexdigest() != digest:
                print(f"  {target.slug} {kind}: the object does not match its row's sha256",
                      file=sys.stderr)
                failed += 1
                continue
            path.parent.mkdir(parents=True, exist_ok=True)
            part = path.with_name(path.name + ".part")
            part.write_bytes(blob)
            part.replace(path)
            fetched += 1
            total += len(blob)
    print(f"{fetched} fetched, {kept} already current, {failed} failed; "
          f"{total / 1e6:.2f} MB under {DATA}", file=sys.stderr)
    return 1 if failed else 0


def verify_client(service: str, client: Path) -> int:
    """Every client asset file against the stored object its row names."""
    mismatched: list[str] = []
    orphans: list[str] = []
    checked = 0
    stored: dict[str, str] = {}
    for target in TARGETS:
        for kind, (_, digest) in stored_objects(service, target).items():
            stored[asset_path(kind, target)] = digest
    # A client checkout keeps them under `assets/`; its test fixtures, from
    # Phase 4 on, keep the same six directories at their own root.
    root = client / "assets" if (client / "assets").is_dir() else client
    families = ("mock", "constraint", "impact", "clinvar", "models", "impact_explanations")
    for family in families:
        directory = root / family
        if not directory.is_dir():
            continue
        for path in sorted(directory.iterdir()):
            if not path.is_file() or path.name.startswith("."):
                continue
            relative = "assets/" + path.relative_to(root).as_posix()
            want = stored.get(relative)
            if want is None:
                orphans.append(relative)
                continue
            checked += 1
            if hashlib.sha256(path.read_bytes()).hexdigest() != want:
                mismatched.append(relative)
    for relative in mismatched:
        print(f"  differs from storage: {relative}", file=sys.stderr)
    for relative in orphans:
        print(f"  no stored object names: {relative}", file=sys.stderr)
    missing = sorted(p for p in stored if not (root / p[len("assets/"):]).exists())
    for relative in missing:
        print(f"  stored, but not in the client: {relative}", file=sys.stderr)
    print(f"{checked} client files checked against storage: {len(mismatched)} differ, "
          f"{len(orphans)} unnamed, {len(missing)} stored but absent", file=sys.stderr)
    return 1 if mismatched or orphans else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--service", default=SERVICE, help=f"default {SERVICE}")
    parser.add_argument("--target", action="append", default=None,
                        help="one slug; repeatable. Default is all twenty.")
    parser.add_argument("--kind", action="append", choices=KINDS, default=None,
                        help="one family; repeatable. Default is every family.")
    parser.add_argument("--verify-client", nargs="?", const=str(CLIENT), default=None,
                        metavar="CLIENT", help=f"compare a client checkout instead "
                                               f"(default {CLIENT})")
    args = parser.parse_args()

    if args.verify_client is not None:
        return verify_client(args.service, Path(args.verify_client))
    wanted = set(args.target or [])
    unknown = wanted - {t.slug for t in TARGETS}
    if unknown:
        raise SystemExit(f"no such target: {sorted(unknown)}")
    targets = [t for t in TARGETS if not wanted or t.slug in wanted]
    return fetch(args.service, targets, tuple(args.kind or KINDS))


if __name__ == "__main__":
    raise SystemExit(main())
