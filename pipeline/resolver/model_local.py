"""A protein's fold, made in the structure bake's environment: an experimental
entry where one passes, and AlphaFold DB's model where a whole search finds none.

    pipeline/structure/venv/bin/python -u -m pipeline.resolver.model_local \\
        job.json out/ --app ../helix-peek --dart /path/to/dart --flutter /path/to/flutter

Run from the repository root by `local_worker.py`, the resolver's worker on a
Mac. That worker's environment has psycopg and Biopython and none of the bake's
tools; the structure bake's has numpy, scipy and trimesh, and is the one the
twenty's models were baked in. So the worker runs beside it, as it runs beside
the scorer's (`score_local.py`), and hands over what a fold is made for in the
job file (`worker.model_job`): AlphaFold's `Job`, the protein's `Target`
(`experimental.target_json`), the row the app's check reads, the name of the
stored record's file beside the job's, and what earlier searches for this
protein decided (`prior`, `store.fold_prior`).

First the PDB (`structure/experimental.py`): every experimental entry of the
protein ranked by the rule, examined, baked and held to the app's own check,
until one passes. Then, only where that search ran whole and nothing passed,
AlphaFold DB's model (`structure/alphafold.py`), unchanged. A search that found
nothing for these inputs under this rule, and has not expired, is not made
again (`prior.none`); one cut short at its bounds goes on from where it stopped
(`prior.progress`). `prior.refresh` asks for neither.

How it ends is the whole protocol, and `out/outcome.json` says it:

- exit 0: a model. `out/` holds `model.glb`, the scene compiled from it as
  `model.fsceneb`, and `described.json`: the fold page's words (`chrome`), its
  `chains` and the bake's `provenance`. `outcome.json`'s `made` says whose:
  `pdb`, with the `folding` track (`folding.json`, its provenance in
  `described.json`) and the bake's `audit.json` beside them, or `alphafold`.
- exit 3: no model from either source. The reader's sentence is in
  `out/refusal.txt`. The worker raises it as `worker.Unmodelled`, which
  `worker.structure_next` tells as a refusal.
- exit 4: the search reached its bounds before it decided. How far it got is
  in `outcome.json`, and the worker asks again next cycle, without spending
  one of the protein's tries.
- anything else: it broke (a service not answering, PyMOL, the importer, the
  app's check), stderr says how, and the worker tries again.

Every ending but the last says, in `outcome.json`'s `selection`, what the
search decided (`fold_report.selection`), or names the earlier one it relied
on (`{"cached": id}`).

A bake that breaks three times is refused for good, so the worker first asks
whether a fold could be made here at all, and claims nothing while not:

    pipeline/structure/venv/bin/python -u -m pipeline.resolver.model_local \\
        --check --app ../helix-peek --dart /path/to/dart --flutter /path/to/flutter

exits 0 where the bake's packages import, PyMOL runs and is the version the
twenty were drawn with, the app's checkout compiles a scene with the
flutter_scene every stored scene was compiled by, and Flutter and the app's
check are there to hold a fold to. Otherwise it says on stderr which is missing.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

REFUSED = 3
UNFINISHED = 4
_CHECK_TIMEOUT = 60


def check(app: Path, dart: str, flutter: str) -> int:
    """Whether a fold could be made here, without asking any service for one."""
    import scipy  # noqa: F401
    import trimesh

    from pipeline.structure import alphafold, bake, experimental

    unready = alphafold.importer_unready(app)
    if unready:
        print(unready, file=sys.stderr)
        return 1
    with tempfile.TemporaryDirectory() as folder:
        # Asked in a script, as the bake asks: on the command line PyMOL
        # rewrites what it is given.
        script = Path(folder) / "version.pml"
        script.write_text("python\nfrom pymol import cmd\n"
                          "print('PYMOL_VERSION ' + cmd.get_version()[0])\npython end\n")
        try:
            pymol = subprocess.run([bake.PYMOL, "-cq", str(script)], capture_output=True,
                                   text=True, timeout=_CHECK_TIMEOUT)
        except (OSError, subprocess.TimeoutExpired) as exc:
            print(f"PyMOL ({bake.PYMOL}) could not be run: {exc}", file=sys.stderr)
            return 1
    if pymol.returncode != 0:
        print(f"PyMOL ({bake.PYMOL}) exited with status {pymol.returncode}: "
              f"{' '.join(pymol.stderr.split())[-200:]}", file=sys.stderr)
        return 1
    version = re.search(r"^PYMOL_VERSION (\S+)$", pymol.stdout, re.M)
    if version is None or version.group(1) != experimental.PYMOL:
        print(f"PyMOL ({bake.PYMOL}) is {version.group(1) if version else 'of no version it says'}"
              f"; every model is drawn with {experimental.PYMOL}.", file=sys.stderr)
        return 1
    if not (app / experimental.BUNDLE_CHECK).exists():
        print(f"The app's checkout at {app} has no {experimental.BUNDLE_CHECK} to hold a fold "
              f"to.", file=sys.stderr)
        return 1
    if not (Path(flutter).exists() and os.access(flutter, os.X_OK)):
        print(f"Flutter ({flutter}) is not there to run the app's check with.", file=sys.stderr)
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


def no_model(entities, alphafold_says: str) -> str:
    """The fold page's sentence for a protein drawn from neither source: what
    the PDB had, and then AlphaFold DB's own reason."""
    first = ("The PDB holds no experimental structure of it." if entities == 0 else
             "None of its experimental structures in the PDB passes the checks for drawing.")
    return f"{first} {alphafold_says}"


def _write(path: Path, said) -> None:
    path.write_text(json.dumps(said, ensure_ascii=False), encoding="utf-8")


def make(job: dict, record: bytes, out: Path, app: Path, dart: str, flutter: str,
         **services) -> int:
    """Make one protein's fold into `out`, and return the exit status that says
    how it ended. `services` stand in for the network in a test
    (`experimental.make`'s, and `alphafold.build`'s `get`)."""
    from pipeline.structure import alphafold, experimental, fold_report, pdb_policy

    out.mkdir(parents=True, exist_ok=True)
    target = experimental.target_of_json(job["target"])
    prior = job.get("prior") or {}
    afresh = bool(prior.get("refresh"))
    alphafold_get = services.pop("alphafold_get", None)
    inputs, release = experimental.inputs_for(
        target, record, fetch_entry=services.pop("fetch_entry", None),
        pfam=services.pop("pfam", None))
    digest = fold_report.inputs_digest(inputs)

    def holds(said) -> bool:
        return (not afresh and bool(said) and said.get("inputs") == digest
                and said.get("policy") == pdb_policy.POLICY_VERSION)

    none = prior.get("none")
    if holds(none):
        # Searched whole for these inputs under this rule, not long ago, and
        # nothing passed.
        selection = {"cached": none["id"]}
        entities = none.get("entities")
    else:
        progress = prior.get("progress")
        resume = (progress["search"], progress["refused"]) if holds(progress) else None
        made = experimental.make(target, record, job["row"], out / "work", app, dart, flutter,
                                 inputs=inputs, release=release, resume=resume, **services)
        choice = made.choice
        selection = fold_report.selection(made.inputs, choice, made.uniprot_release)
        if choice.outcome == "incomplete":
            _write(out / "outcome.json", {"made": None, "unfinished": True,
                                          "selection": selection})
            return UNFINISHED
        built = made.built
        if built is not None:
            (out / "model.glb").write_bytes(built.glb)
            (out / "model.fsceneb").write_bytes(built.scene)
            (out / "folding.json").write_bytes(built.folding)
            _write(out / "audit.json", {**built.audit, "selection": selection,
                                        "measured": built.measured})
            _write(out / "described.json", {
                "chrome": built.chrome, "chains": built.chains, "provenance": built.provenance,
                "folding_provenance": built.folding_provenance})
            _write(out / "outcome.json", {"made": "pdb", "selection": selection,
                                          "measured": built.measured})
            return 0
        entities = len(choice.discovery.ids)

    # The whole search found nothing that passes: AlphaFold DB's model, as before.
    model_job = alphafold.Job(
        slug=job["slug"], accession=job["accession"], display=job["display"],
        protein=job["protein"], kept=tuple(tuple(span) for span in job["kept"]),
        disulfides=tuple(tuple(pair) for pair in job["disulfides"]), allowed=job["allowed"])
    try:
        built = alphafold.build(model_job, out / "alphafold", app, dart,
                                **({"get": alphafold_get} if alphafold_get else {}))
    except alphafold.Refused as exc:
        (out / "refusal.txt").write_text(no_model(entities, str(exc)), encoding="utf-8")
        _write(out / "outcome.json", {"made": None, "selection": selection})
        return REFUSED
    (out / "model.glb").write_bytes(built.glb)
    (out / "model.fsceneb").write_bytes(built.scene)
    _write(out / "described.json",
           {"chrome": built.chrome, "chains": built.chains, "provenance": built.provenance})
    _write(out / "outcome.json", {"made": "alphafold", "selection": selection})
    return 0


def main(arguments: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--app", type=Path, required=True,
                        help="the app's checkout, for its scene importer and its check")
    parser.add_argument("--dart", required=True)
    parser.add_argument("--flutter", required=True)
    parser.add_argument("files", nargs="*", metavar="job.json out/")
    args = parser.parse_args(arguments)
    if args.check:
        return check(args.app, args.dart, args.flutter)
    if len(args.files) != 2:
        print(__doc__, file=sys.stderr)
        return 2

    job_file, out = Path(args.files[0]), Path(args.files[1])
    job = json.loads(job_file.read_text(encoding="utf-8"))
    record = (job_file.parent / job["record"]).read_bytes()
    # Only a refusal or a search cut short is said with an exit status of its
    # own. Anything else raised is this process breaking, and is left to exit 1.
    return make(job, record, out, args.app, args.dart, args.flutter)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
