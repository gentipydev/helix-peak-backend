"""`alphafold.build` for one protein, as a process of the structure bake's environment.

    pipeline/structure/venv/bin/python -u -m pipeline.resolver.model_local \\
        job.json out/ --app ../helix-peek --dart /path/to/dart

Run from the repository root by `local_worker.py`, the resolver's worker on a
Mac. That worker's environment has psycopg and Biopython and no numpy; the
structure bake's has numpy, scipy and trimesh, and is the one the twenty's
models were baked in. So the worker runs beside it, as it runs beside the
scorer's (`score_local.py`), and hands over what a model is made for: the
`alphafold.Job` as JSON (`worker.model_job`).

How it ends is the whole protocol:

- exit 0: `out/` holds `model.glb`, the scene compiled from it as
  `model.fsceneb`, and `described.json`: the fold page's words (`chrome`), its
  `chains`, and the bake's `provenance`.
- exit 3: no model is drawn for this protein (`alphafold.Refused`: AlphaFold
  DB has none, or none of this sequence, or none sure enough, or none that
  fits) and the reader's sentence is in `out/refusal.txt`. The worker raises
  it as `worker.Unmodelled`, which `worker.structure_next` tells as a refusal.
- anything else: it broke (the API not answering, PyMOL, the importer), stderr
  says how, and the worker tries again.

A bake that breaks three times is refused for good, so the worker first asks
whether a model could be made here at all, and claims nothing while not:

    pipeline/structure/venv/bin/python -u -m pipeline.resolver.model_local \\
        --check --app ../helix-peek --dart /path/to/dart

exits 0 where the bake's packages import, PyMOL runs, and the app's checkout
compiles a scene with the flutter_scene every stored scene was compiled by.
Otherwise it says on stderr which is missing.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

REFUSED = 3
_CHECK_TIMEOUT = 60


def check(app: Path, dart: str) -> int:
    """Whether a model could be made here, without asking AlphaFold DB for one."""
    import scipy  # noqa: F401
    import trimesh

    from pipeline.structure import alphafold, bake

    unready = alphafold.importer_unready(app)
    if unready:
        print(unready, file=sys.stderr)
        return 1
    try:
        pymol = subprocess.run([bake.PYMOL, "-cq", "-d", "quit"], capture_output=True,
                               text=True, timeout=_CHECK_TIMEOUT)
    except (OSError, subprocess.TimeoutExpired) as exc:
        print(f"PyMOL ({bake.PYMOL}) could not be run: {exc}", file=sys.stderr)
        return 1
    if pymol.returncode != 0:
        print(f"PyMOL ({bake.PYMOL}) exited with status {pymol.returncode}: "
              f"{' '.join(pymol.stderr.split())[-200:]}", file=sys.stderr)
        return 1
    # The importer, on the smallest model there is: it has no way to be asked
    # whether it would run, short of running it.
    with tempfile.TemporaryDirectory() as folder:
        glb = Path(folder) / "check.glb"
        scene = trimesh.Scene()
        scene.add_geometry(trimesh.creation.box(), geom_name="check", node_name="check")
        scene.export(glb)
        try:
            alphafold.compile_scene(glb, Path(folder) / "check.fsceneb", app, dart)
        except (OSError, RuntimeError, subprocess.TimeoutExpired) as exc:
            print(f"The scene importer could not be run from {app}: "
                  f"{' '.join(str(exc).split())[-300:]}", file=sys.stderr)
            return 1
    return 0


def main(arguments: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--app", type=Path, required=True,
                        help="the app's checkout, for its scene importer")
    parser.add_argument("--dart", required=True)
    parser.add_argument("files", nargs="*", metavar="job.json out/")
    args = parser.parse_args(arguments)
    if args.check:
        return check(args.app, args.dart)
    if len(args.files) != 2:
        print(__doc__, file=sys.stderr)
        return 2

    from pipeline.structure import alphafold

    job_file, out = Path(args.files[0]), Path(args.files[1])
    said = json.loads(job_file.read_text(encoding="utf-8"))
    job = alphafold.Job(
        slug=said["slug"], accession=said["accession"], display=said["display"],
        protein=said["protein"], kept=tuple(tuple(span) for span in said["kept"]),
        disulfides=tuple(tuple(pair) for pair in said["disulfides"]), allowed=said["allowed"])
    out.mkdir(parents=True, exist_ok=True)
    # Only the bake's own `Refused` is a refusal. Anything else raised is this
    # process breaking, and is left to exit 1.
    try:
        built = alphafold.build(job, out / "work", args.app, args.dart)
    except alphafold.Refused as exc:
        (out / "refusal.txt").write_text(str(exc), encoding="utf-8")
        return REFUSED
    (out / "model.glb").write_bytes(built.glb)
    (out / "model.fsceneb").write_bytes(built.scene)
    (out / "described.json").write_text(json.dumps(
        {"chrome": built.chrome, "chains": built.chains, "provenance": built.provenance},
        ensure_ascii=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
