"""Bake one curated row whole: every track, each in its own environment, then the checks.

    NCBI_EMAIL=you@example.com .venv/bin/python pipeline/bake_curated.py --target phospholamban
    .venv/bin/python pipeline/bake_curated.py --target phospholamban --dry-run      # the commands
    .venv/bin/python pipeline/bake_curated.py --target phospholamban --from impact  # resume

In order, each a process of its own, stopping at the first that fails:

| step | baker | environment |
|---|---|---|
| record | `mock/build_gene_record.py` | `.venv`, `NCBI_EMAIL` |
| structure | `structure/bake.py` | `pipeline/structure/venv`, PyMOL |
| folding | `folding/bake_folding.py` | `pipeline/structure/venv`, PyMOL |
| bundle | `structure/check_bundle.py`: the app's own check of the fold | `pipeline/structure/venv`, the app's checkout |
| constraint | `constraint/score_protein.py` (ESM-2) | `pipeline/.esm-venv`, never installed into |
| impact | `impact/bake_impact.py` (AVI) | `pipeline/impact/venv`, `ALPHAGENOME_API_KEY`, `uv` |
| clinvar | `clinvar/bake_clinvar.py` | `.venv`, `NCBI_EMAIL` |
| locus | `locus/bake_locus.py` | `.venv` |
| check | `check_assets.py`, `folding/check_folding.py`, `locus/check_locus.py` | `.venv` |

Every step is given `--target` and nothing else, so no other row is touched.
`bundle` runs for a row cut from a working copy (`structure/pdb.source_of`):
the twenty were never held to the app's bundle check. A row unscored, with no
impact track or no ClinVar snapshot (`targets.py`) skips that step. After the
record, a row cut from a working copy has its translation held to the one its
fold was chosen for (the working copy's sidecar).

`ALPHAGENOME_API_KEY` is read from the environment, else the backend's `.env`,
else `~/.env`, and handed to the impact step alone. A row whose record the
service already serves is a row that is live: it is baked again only with
`--rebake`. The bakes write under `pipeline/data/` (`HELIXPEEK_DATA`), as by
hand; nothing is uploaded or seeded.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional

BACKEND = Path(__file__).resolve().parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

SERVICE = "https://helix-peak-backend.onrender.com"
STEPS = ("record", "structure", "folding", "bundle", "constraint", "impact", "clinvar",
         "locus", "check")
VENV = BACKEND / ".venv" / "bin" / "python"
STRUCTURE = BACKEND / "pipeline" / "structure" / "venv" / "bin" / "python"
ESM = BACKEND / "pipeline" / ".esm-venv" / "bin" / "python"
IMPACT = BACKEND / "pipeline" / "impact" / "venv" / "bin" / "python"
SECRET = "ALPHAGENOME_API_KEY"


def _dotenv(path: Path) -> Dict[str, str]:
    """`KEY=value` lines, as a shell sources them; nothing else is read."""
    found: Dict[str, str] = {}
    if not path.exists():
        return found
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        found[key.removeprefix("export ").strip()] = value.strip().strip("\"'")
    return found


def api_key(environ: Dict[str, str], home: Path) -> Optional[str]:
    """The AlphaGenome key: the environment's, else the backend's `.env`, else `~/.env`."""
    return (environ.get(SECRET) or _dotenv(BACKEND / ".env").get(SECRET)
            or _dotenv(home / ".env").get(SECRET))


def plan(target, environ: Dict[str, str], home: Path, working_copy: bool) -> List[dict]:
    """Every step for [target], in order: its name, its command and its
    environment. A step the row does not have is left out."""
    from pipeline.resolver.local_worker import default_uv

    slug = target.slug
    base = {key: value for key, value in environ.items() if key != SECRET}
    ncbi = {**base, "NCBI_EMAIL": environ.get("NCBI_EMAIL", "")}
    uv = default_uv(home)
    impact = {**base, SECRET: api_key(environ, home) or "",
              "PATH": f"{uv.parent}{os.pathsep}{base.get('PATH', '')}"}

    def step(name: str, python: Path, script: str, env: Dict[str, str], *extra: str) -> dict:
        return {"step": name, "env": env,
                "command": [str(python), "-u", str(BACKEND / "pipeline" / script),
                            "--target", slug, *extra]}

    steps = [step("record", VENV, "mock/build_gene_record.py", ncbi),
             step("structure", STRUCTURE, "structure/bake.py", base),
             step("folding", STRUCTURE, "folding/bake_folding.py", base)]
    if working_copy:
        steps.append(step("bundle", STRUCTURE, "structure/check_bundle.py", base))
    if target.scored:
        steps.append(step("constraint", ESM, "constraint/score_protein.py", base))
    if target.impact_scored:
        steps.append(step("impact", IMPACT, "impact/bake_impact.py", impact))
    if target.clinvar_available:
        steps.append(step("clinvar", VENV, "clinvar/bake_clinvar.py", ncbi))
    steps.append(step("locus", VENV, "locus/bake_locus.py", base))
    steps += [step("check", VENV, script, base)
              for script in ("check_assets.py", "folding/check_folding.py",
                             "locus/check_locus.py")]
    return steps


def chosen(steps: List[dict], only: List[str], start: Optional[str]) -> List[dict]:
    """The steps asked for: `--only` these, or from `--from` on, or all."""
    if only:
        return [s for s in steps if s["step"] in only]
    if start:
        names = [s["step"] for s in steps]
        if start not in names:
            raise SystemExit(f"this row has no {start!r} step ({', '.join(dict.fromkeys(names))})")
        return steps[names.index(start):]
    return steps


def stored_record(service: str, slug: str) -> Optional[bool]:
    """Whether the service serves this row's record: True, False (no such
    protein yet), or None where it cannot be asked."""
    try:
        with urllib.request.urlopen(f"{service}/protein/{slug}/tracks", timeout=90) as answer:
            rows = json.loads(answer.read())
    except urllib.error.HTTPError as exc:
        return False if exc.code == 404 else None
    except (urllib.error.URLError, OSError, ValueError):
        return None
    return (rows.get("record") or {}).get("state") == "ready"


def translation_held(target) -> Optional[str]:
    """Why the record does not translate to what the row's fold was chosen
    for, or None where it does (or the row has no working copy)."""
    import hashlib

    from pipeline.paths import DATA
    from pipeline.structure import pdb

    copy = pdb.working_copy_of(target)
    if copy is None:
        return None
    said = json.loads(copy.with_suffix(".json").read_text(encoding="utf-8"))
    translation = json.loads((DATA / target.mock_asset).read_bytes())["protein"]["translation"]
    if hashlib.sha256(translation.encode()).hexdigest() != said["sha256"]["translation"]:
        return (f"{target.slug}: the record translates to another protein than the one its "
                f"fold was chosen for ({copy.with_suffix('.json').name})")
    return None


def main(arguments: List[str]) -> int:
    from pipeline.structure import pdb
    from pipeline.targets import BY_SLUG

    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--target", required=True, help="the row's slug")
    parser.add_argument("--only", action="append", default=[], choices=STEPS)
    parser.add_argument("--from", dest="start", default=None, choices=STEPS)
    parser.add_argument("--dry-run", action="store_true", help="print the commands; run none")
    parser.add_argument("--rebake", action="store_true",
                        help="bake a row whose record the service already serves")
    parser.add_argument("--service", default=SERVICE)
    args = parser.parse_args(arguments)
    target = BY_SLUG.get(args.target)
    if target is None:
        parser.error(f"no such target: {args.target!r}")
    home = Path.home()
    steps = chosen(plan(target, dict(os.environ), home, pdb.working_copy_of(target) is not None),
                   args.only, args.start)

    if args.dry_run:
        for s in steps:
            secret = " (with the AlphaGenome key)" if s["env"].get(SECRET) else ""
            print(f"[{s['step']}]{secret} {' '.join(s['command'])}")
        return 0
    if any(s["step"] in ("record", "clinvar") for s in steps) and not os.environ.get("NCBI_EMAIL"):
        raise SystemExit("Set NCBI_EMAIL: the record and ClinVar are fetched from NCBI.")
    if any(s["step"] == "impact" for s in steps) and not api_key(dict(os.environ), home):
        raise SystemExit(f"No {SECRET} in the environment, the backend's .env or ~/.env.")
    live = stored_record(args.service, target.slug)
    if live is not False and not args.rebake:
        raise SystemExit(f"{target.slug}: " + ("the service serves its record already: it is "
                                               "live. Pass --rebake to bake it again."
                                               if live else "the service could not be asked "
                                               "whether it is live. Pass --rebake to bake anyway."))

    started = time.monotonic()
    for s in steps:
        print(f"\n== {target.slug}: {s['step']} ==  {' '.join(s['command'][2:])}", flush=True)
        began = time.monotonic()
        ran = subprocess.run(s["command"], cwd=str(BACKEND), env=s["env"], stdin=subprocess.DEVNULL)
        took = time.monotonic() - began
        if ran.returncode != 0:
            print(f"\n{target.slug}: {s['step']} failed (exit {ran.returncode}) after "
                  f"{took:.0f} s; resume with --from {s['step']}", file=sys.stderr)
            return 1
        if s["step"] == "record":
            wrong = translation_held(target)
            if wrong:
                print(f"\n{wrong}", file=sys.stderr)
                return 1
        print(f"-- {s['step']} done in {took:.0f} s", flush=True)
    print(f"\n{target.slug}: {len(steps)} step(s) done in {time.monotonic() - started:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
