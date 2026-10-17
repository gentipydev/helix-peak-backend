"""Hold a curated row's baked fold to the app's own check, as a fold made on demand is.

    pipeline/structure/venv/bin/python pipeline/structure/check_bundle.py --target phospholamban

For a row cut from a working copy (`pdb.source_of`; the twenty are not), after
`bake.py` and `bake_folding.py`:

- The model is the one the draft's own bake made from the same working copy,
  which the app's check passed then: its `.glb` has the digest the working
  copy's sidecar (`structures/working/<slug>.json`) records.
- The app's check passes it again, now with the curated row the app will read,
  the fold this checkout baked and the scene compiled from the model: the
  bundle `experimental.bundle` lays out (`row.json`, `model.glb`,
  `model.fsceneb`, `folding.json`), held to `test/fold_bundle_check.dart` in
  the app's checkout by `experimental.app_check`. Every ring and rod binds, the
  fold ends on the stored mesh float for float, and nothing jumps.

The bundle is left in `structure/output/<slug>/bundle/`, which git ignores.
Exits 1 with the reason where either fails.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.folding.bake_folding import folding_asset  # noqa: E402
from pipeline.paths import CLIENT, DATA  # noqa: E402
from pipeline.structure import pdb  # noqa: E402
from pipeline.targets import BY_SLUG  # noqa: E402

OUTPUT = HERE / "output"


def served(target) -> dict:
    """The row as the service will serve it once seeded, as far as the app's
    check reads one (`experimental.served_row`)."""
    from pipeline import seed_catalog
    from pipeline.structure.experimental import served_row

    return served_row(seed_catalog.protein_row(target, seed_catalog.curated_rows()[target.slug]))


def check(target, app: Path, dart: str, flutter: str) -> dict:
    """The bundle's verdict, with what the check measured. Raises where the
    model is not the draft's or the check does not pass it."""
    from pipeline.structure import alphafold, experimental

    copy = pdb.working_copy_of(target)
    if copy is None:
        raise ValueError(f"{target.slug} is not cut from a working copy")
    said = json.loads(copy.with_suffix(".json").read_text(encoding="utf-8"))
    glb = (DATA / target.structure_asset).read_bytes()
    digest = hashlib.sha256(glb).hexdigest()
    if digest != said["sha256"]["glb"]:
        raise ValueError(f"{target.slug}: the model baked here ({digest[:12]}) is not the one "
                         f"the draft's check passed ({said['sha256']['glb'][:12]})")
    folder = OUTPUT / target.slug / "bundle"
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    (folder / "row.json").write_text(json.dumps(served(target), ensure_ascii=False),
                                     encoding="utf-8")
    (folder / "model.glb").write_bytes(glb)
    alphafold.compile_scene(folder / "model.glb", folder / "model.fsceneb", app, dart)
    (folder / "folding.json").write_bytes((DATA / folding_asset(target)).read_bytes())
    return experimental.app_check(folder, app, flutter)


def main(arguments: list[str]) -> int:
    from pipeline.resolver.local_worker import default_dart, default_flutter

    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--target", action="append", required=True)
    parser.add_argument("--app", type=Path, default=CLIENT)
    parser.add_argument("--dart", default=None)
    parser.add_argument("--flutter", default=None)
    args = parser.parse_args(arguments)
    unknown = sorted(set(args.target) - set(BY_SLUG))
    if unknown:
        parser.error(f"no such target: {unknown}")
    dart = args.dart or str(default_dart(Path.home()))
    flutter = args.flutter or str(default_flutter(Path.home()))
    failed = 0
    for slug in args.target:
        try:
            measured = check(BY_SLUG[slug], args.app, dart, flutter)
        except Exception as exc:  # noqa: BLE001 -- said, then the next
            print(f"{slug}: {type(exc).__name__}: {exc}", file=sys.stderr)
            failed += 1
            continue
        print(f"{slug}: the model is the draft's, and the app's check passes it "
              f"({measured.get('rings')} rings, {measured.get('bridges')} bridges, "
              f"bound in {measured.get('bind_ms')} ms)")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
