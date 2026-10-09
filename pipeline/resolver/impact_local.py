"""`bake_impact.bake` for one protein, as a process of the AVI bake's environment.

    pipeline/impact/venv/bin/python -u -m pipeline.resolver.impact_local \\
        target.pkl record.json out/ --state <the worker's state>/impact

Run from the repository root by `local_worker.py`, the resolver's worker on a
Mac. That worker's environment has psycopg and Biopython and no AlphaGenome
client; the AVI bake's, `pipeline/impact/venv`, has the client and numpy, and
is the one the twenty's tracks were baked in. So the worker runs beside it, as
it runs beside the scorer's (`score_local.py`), and hands over what a bake is
made from: the `Target`, pickled as `score_local.py` takes it, and the
record's bytes.

How it ends is the whole protocol:

- exit 0: `out/impact.json` holds the track's bytes.
- exit 3: one of the bake's own gates said no (`BakeError`: the record's
  exons do not pair with GENCODE's MANE Select ones, its bases are not
  GRCh38's where GENCODE puts the gene, too few were scored, or the scores
  say the map is wrong) and the gate's words are in `out/refusal.txt`. The
  worker raises them as `worker.Unplaced`, which `worker.impact_next` tells as
  a refusal.
- anything else: it broke (the Atlas or GENCODE's lookup not answering, the
  quota spent past the bake's own patience), stderr says how, and the worker
  tries again. What was pulled is checkpointed under `--state`, so the next
  try starts where this one stopped.

The bake is called unchanged. What is set from here is where it reads and
writes (`DATA`, `GENCODE_DIR`, `CHECKPOINT_DIR`, never `pipeline/impact/`), and
three of its functions are wrapped to tell its verdicts from its failures:

- `pull`'s `BakeError` is the Atlas not answering, not a verdict, and is
  raised as `Unreachable`.
- `gencode`'s is GENCODE's lookup failing, except "expected exactly one MANE
  Select transcript", which is a fact about the gene.
- `check_biology` is skipped for a gene with no intron, the user's choice of
  2026-10-09: exons outscoring intron interiors is a smoke alarm a gene of one
  exon cannot sound, and the sequence gate still holds every base to GRCh38.
  The track then says null for the intron and splice medians, which nothing
  reads.

A bake that breaks three times is refused for good, so the worker first asks
whether one could be made here at all, and claims nothing while not:

    pipeline/impact/venv/bin/python -u -m pipeline.resolver.impact_local --check

exits 0 where the client imports, `ALPHAGENOME_API_KEY` is set, the skill's
GENCODE lookup starts under `uv`, and the Atlas answers for one base with
that key (one request). Otherwise it says on stderr which is missing.
"""

from __future__ import annotations

import argparse
import os
import pickle
import statistics
import subprocess
import sys
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

REFUSED = 3
_CHECK_TIMEOUT = 60
# The base the check asks the Atlas about: one of TTR's, which it scores.
_PROBE = ("chr18", 31_591_767)
_ONE_MANE = "expected exactly one MANE Select transcript"


class Unreachable(RuntimeError):
    """The Atlas or GENCODE's lookup did not answer. Worth another try, and
    never a verdict on the protein."""


def interior(bake_impact, record: dict) -> int:
    """How many drawn intron bases lie further than `check_biology`'s edge (8
    bases, or half a short intron) from either end of their intron: the bases
    its alarm compares the exons with."""
    total = 0
    for low, high in bake_impact.intron_spans(record):
        span = high - low + 1
        total += max(0, span - 2 * min(8, span // 2))
    return total


def guard(bake_impact) -> None:
    """Wrap the three functions the module docstring names, once."""
    if getattr(bake_impact, "_guarded", False):
        return
    pull, gencode, biology = bake_impact.pull, bake_impact.gencode, bake_impact.check_biology

    def pulling(client, chromosome, low, high):
        try:
            return pull(client, chromosome, low, high)
        except bake_impact.BakeError as exc:
            raise Unreachable(str(exc)) from exc

    def looking_up(gene, skill, refresh=False):
        try:
            return gencode(gene, skill, refresh=refresh)
        except bake_impact.BakeError as exc:
            if _ONE_MANE in str(exc):
                raise
            raise Unreachable(str(exc)) from exc

    def checking(record, positions, gene):
        if interior(bake_impact, record):
            return biology(record, positions, gene)
        exonic = [max(positions[local]) for exon in record["exons"]
                  for local in range(exon["start"], exon["end"] + 1) if local in positions]
        if not exonic:
            return biology(record, positions, gene)
        return {"exon_median_phred": round(statistics.median(exonic), 2),
                "splice_median_phred": None, "intron_median_phred": None}

    bake_impact.pull, bake_impact.gencode, bake_impact.check_biology = \
        pulling, looking_up, checking
    bake_impact._guarded = True


def bake(target, record: bytes, data: Path, state: Path, client) -> bytes:
    """The AVI track's bytes for one protein: `bake_impact.bake`, unchanged,
    reading and writing in `data`, keeping GENCODE's answers and its
    checkpoints in `state`. A gate's verdict is raised as the bake's own
    `BakeError`; anything else is the bake breaking."""
    from pipeline.impact import bake_impact

    guard(bake_impact)
    bake_impact.DATA = data
    bake_impact.GENCODE_DIR = state / "gencode"
    bake_impact.CHECKPOINT_DIR = state / "checkpoints"
    path = data / target.mock_asset
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(record)
    try:
        bake_impact.bake(target, client, bake_impact.DEFAULT_SKILL, resume=True,
                         refresh_gencode=False)
    except bake_impact.BakeError:
        # A verdict: what was pulled for it will never be needed again.
        (bake_impact.CHECKPOINT_DIR / f"{target.slug}.jsonl").unlink(missing_ok=True)
        raise
    return (data / target.impact_asset).read_bytes()


def check() -> int:
    """Whether an AVI track could be made here, at the cost of one request."""
    from alphagenome.atlas import atlas
    from alphagenome.data import genome

    from pipeline.impact import bake_impact

    key = os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        print("ALPHAGENOME_API_KEY is not set in the worker's .env.", file=sys.stderr)
        return 1
    script = bake_impact.DEFAULT_SKILL / "scripts" / "alphagenome_atlas_avi.py"
    if not script.exists():
        print(f"There is no AlphaGenome skill at {bake_impact.DEFAULT_SKILL}.", file=sys.stderr)
        return 1
    try:
        lookup = subprocess.run(["uv", "run", str(script), "gtf", "--help"],
                                cwd=bake_impact.DEFAULT_SKILL, capture_output=True, text=True,
                                timeout=_CHECK_TIMEOUT)
    except (OSError, subprocess.TimeoutExpired) as exc:
        print(f"GENCODE's lookup could not be run under uv: {exc}", file=sys.stderr)
        return 1
    if lookup.returncode != 0:
        print(f"GENCODE's lookup exited with status {lookup.returncode}: "
              f"{' '.join(lookup.stderr.split())[-200:]}", file=sys.stderr)
        return 1
    chromosome, position = _PROBE
    try:
        scores = atlas.create(key).query_interval(
            genome.Interval(chromosome=chromosome, start=position - 1, end=position),
            requested_scorers=[bake_impact.SCORER])
    except Exception as exc:  # noqa: BLE001 -- grpc raises its own family
        print(f"The Atlas did not answer: {' '.join(str(exc).split())[-300:]}", file=sys.stderr)
        return 1
    if bake_impact.SCORER not in scores:
        print(f"The Atlas answered with no {bake_impact.SCORER}.", file=sys.stderr)
        return 1
    return 0


def main(arguments: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--state", type=Path,
                        help="what the bake keeps between tries: GENCODE's answers, checkpoints")
    parser.add_argument("files", nargs="*", metavar="target.pkl record.json out/")
    args = parser.parse_args(arguments)
    if args.check:
        return check()
    if len(args.files) != 3 or args.state is None:
        print(__doc__, file=sys.stderr)
        return 2

    target_file, record_file, out = (Path(name) for name in args.files)
    out.mkdir(parents=True, exist_ok=True)
    # Named before `pipeline.paths` is first imported, and set on the bake
    # besides: nothing is written into `pipeline/data/`.
    os.environ["HELIXPEEK_DATA"] = str(out / "data")
    key = os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        print("ALPHAGENOME_API_KEY is not set in the worker's .env.", file=sys.stderr)
        return 1
    target = pickle.loads(target_file.read_bytes())

    from alphagenome.atlas import atlas

    from pipeline.impact import bake_impact

    try:
        payload = bake(target, record_file.read_bytes(), out / "data", args.state,
                       atlas.create(key))
    except bake_impact.BakeError as exc:
        (out / "refusal.txt").write_text(str(exc), encoding="utf-8")
        return REFUSED
    (out / "impact.json").write_bytes(payload)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
