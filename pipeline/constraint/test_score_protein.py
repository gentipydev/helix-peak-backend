"""Fast scientific-contract checks; no weights or network required."""

import contextlib
import dataclasses
import io
import json
import math
from pathlib import Path
import sys
import unittest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.constraint.score_protein import (  # noqa: E402
    AMINO_ACIDS, alignment, gate, normalize_conservation,
)
from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import BY_SLUG, TARGETS, Region, Target, partition  # noqa: E402

# The tracks live in storage now. `fetch_tracks.py` brings them back into the
# data directory, and until it has, there is nothing on disk to hold to the gate.
HAVE_TRACKS = any((DATA / t.constraint_asset).exists() for t in TARGETS if t.scored)


def _positions(sequence: str, entropies=None, favoured=None):
    """A scored track, with the cysteines at the top as a real one has them.

    At each position the residue `favoured(i)` names scores above every other
    one. By default that is the residue that is there, which is what a track
    filed under the right residues looks like.
    """
    positions = []
    for i, aa in enumerate(sequence):
        best = favoured(i) if favoured else aa
        positions.append({
            "index": i,
            "wildtype": aa,
            "entropy": entropies[i] if entropies else (0.01 if aa == "C" else 1 + i / 1000),
            "conservation": 0.99 if aa == "C" else 0.2,
            "substitutions": {
                other: 0.0 if other == aa else (1.0 if other == best else -1.0)
                for other in AMINO_ACIDS
            },
        })
    return positions


INSULIN = (
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED"
    "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)


class ConstraintChecks(unittest.TestCase):
    def test_confidence_increases_when_entropy_decreases(self):
        self.assertEqual(normalize_conservation([0.0, 0.5, 2.0]), [1.0, 0.75, 0.0])

    def test_rejects_undefined_or_invalid_normalization(self):
        for values in ([], [1.0, 1.0], [0.0, math.nan], [0.0, math.inf], [-1.0, 0.0]):
            with self.subTest(values=values), self.assertRaises(ValueError):
                normalize_conservation(values)

    def test_insulin_cysteine_numbering_matches_the_table(self):
        target = BY_SLUG["insulin"]
        self.assertEqual(len(INSULIN), target.aa)
        bonded = sorted({n for pair in target.disulfides for n in pair})
        self.assertEqual([i + 1 for i, aa in enumerate(INSULIN) if aa == "C"], bonded)

    def test_every_target_partitions_its_precursor_with_no_gaps(self):
        for target in TARGETS:
            with self.subTest(target=target.slug):
                regions = partition(target)
                self.assertEqual(regions[0]["start"], 1)
                self.assertEqual(regions[-1]["end"], target.aa)
                for before, after in zip(regions, regions[1:]):
                    self.assertEqual(after["start"], before["end"] + 1)
                for number in {n for pair in target.disulfides for n in pair}:
                    self.assertTrue(any(r["start"] <= number <= r["end"] for r in regions))

    def test_an_unnamed_end_of_a_cut_precursor_is_an_extension(self):
        """Glucagon's trailing RK, named for the chain, read "removed with
        Proglucagon"."""
        last = partition(BY_SLUG["glucagon"])[-1]
        self.assertEqual(
            (last["label"], last["start"], last["end"], last["kept"]),
            ("C-terminal extension", 179, 180, False),
        )
        # Either end of a cut precursor. A folded chain's ends are more of that
        # chain, and a cut precursor with nothing named is named for its chain.
        cut = dataclasses.replace(
            BY_SLUG["glucagon"], aa=10, chain_label="Pro-X",
            regions=(Region("Chain", "", 5, 8),),
        )
        self.assertEqual(
            [r["label"] for r in partition(cut)],
            ["N-terminal extension", "Chain", "C-terminal extension"],
        )
        self.assertEqual(
            [r["label"] for r in partition(dataclasses.replace(cut, cleaved=False))],
            ["Pro-X", "Chain", "Pro-X"],
        )
        self.assertEqual([r["label"] for r in partition(dataclasses.replace(cut, regions=()))], ["Pro-X"])

    def test_gate_rejects_a_track_shifted_by_one_residue(self):
        """What an offset bug looks like: every residue scored as its neighbour."""
        target = BY_SLUG["insulin"]
        with contextlib.redirect_stdout(io.StringIO()):
            gate(target, INSULIN, _positions(INSULIN))
            for step in (-1, 1):
                shifted = _positions(
                    INSULIN, favoured=lambda i: INSULIN[min(max(i + step, 0), len(INSULIN) - 1)]
                )
                with self.subTest(step=step), self.assertRaisesRegex(ValueError, "gate FAILED"):
                    gate(target, INSULIN, shifted)

    def test_gate_publishes_bridges_less_constrained_than_the_rest(self):
        """SOD1's bridge ranks 46 and 51 of 154, under its copper and zinc
        histidines. That is a finding, and an aligned track publishes it."""
        target = BY_SLUG["sod1"]
        sequence = "".join(
            "C" if i + 1 in {58, 147} else "AG"[i % 2] for i in range(target.aa)
        )
        entropies = [1 + i / 1000 for i in range(target.aa)]
        entropies[146], entropies[57] = 1.0445, 1.0485
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            gate(target, sequence, _positions(sequence, entropies))
        self.assertIn("rank=46/154", output.getvalue())
        self.assertIn("rank=51/154", output.getvalue())

    def test_gate_rejects_a_bridge_on_a_residue_that_is_not_a_cysteine(self):
        target = BY_SLUG["sod1"]
        sequence = "AG" * (target.aa // 2)
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError, "not a cysteine"):
                gate(target, sequence, _positions(sequence))

    def test_alignment_is_not_measured_inside_a_run_of_one_residue(self):
        """In a polyQ tract the residue before and after each Q is Q, so a
        position there says nothing about which residue a score belongs to.
        Counted, the forty would read as forty misfiled scores."""
        sequence = "MKT" + "Q" * 40 + "LVR"
        self.assertEqual(alignment(sequence, _positions(sequence)), (1.0, 1.0))
        with self.assertRaisesRegex(ValueError, "cannot be measured"):
            alignment("QQQQ", _positions("QQQQ"))

    def test_gate_rejects_a_track_that_is_not_the_sequence(self):
        target = BY_SLUG["insulin"]
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError, "positions for"):
                gate(target, INSULIN, [])

    @unittest.skipUnless(HAVE_TRACKS, "no constraint tracks under pipeline/data; "
                                      "run pipeline/fetch_tracks.py")
    def test_every_track_on_disk_passes_the_gate(self):
        """The gate a track would be baked under today, over every track in the
        data directory. A row not baked yet is skipped: `check_assets.py` is what
        wants the file to exist."""
        checked = 0
        for target in TARGETS:
            path = DATA / target.constraint_asset
            if not target.scored or not path.exists():
                continue
            asset = json.loads(path.read_text())
            with self.subTest(target=target.slug), contextlib.redirect_stdout(io.StringIO()):
                gate(target, asset["sequence"], asset["positions"])
            checked += 1
        self.assertGreater(checked, 0)


if __name__ == "__main__":
    unittest.main()
