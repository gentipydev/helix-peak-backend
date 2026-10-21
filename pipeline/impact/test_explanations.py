import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.impact.bake_explanations import finite  # noqa: E402
from pipeline.impact.bake_impact import BakeError  # noqa: E402
from pipeline.impact.check_explanations import validate  # noqa: E402
from pipeline.impact.track_winners import reduce  # noqa: E402
from pipeline.paths import DATA  # noqa: E402

# Insulin's two tracks, from storage by way of `fetch_tracks.py`.
_AVI = DATA / "assets/impact/insulin_avi.json"
_EXPLANATIONS = DATA / "assets/impact_explanations/insulin.json"


class _Frame:
    """The slice of a pandas frame the reducer reads: `len` and `.iloc[i].get`."""

    def __init__(self, rows):
        self.iloc = rows

    def __len__(self):
        return len(self.iloc)


class _Scores:
    def __init__(self, x, tracks, genes=()):
        self.X = np.array(x, dtype=np.float32)
        self.var = _Frame(list(tracks))
        self.obs = _Frame(list(genes))


def _track(title, biosample, **extra):
    return {"Assay title": title, "biosample_name": biosample, "ontology_curie": "UBERON:1", **extra}


class ExplanationsTest(unittest.TestCase):
    def test_rejects_nonfinite_model_data(self):
        for value in (float("nan"), float("inf"), -float("inf")):
            with self.assertRaises(BakeError):
                finite([0.1, value])
        with self.assertRaises(BakeError):
            finite([])

    def test_assay_winner_is_the_largest_absolute_effect_with_its_sign(self):
        replies = {"RNA_SEQ": _Scores(
            [[0.2, -0.9], [0.5, 0.1]],
            [_track("total RNA-seq", "liver"), _track("total RNA-seq", "pancreas")],
            [{"gene_name": "INS"}, {"gene_name": "TH"}],
        )}
        *track, value = reduce(replies)["MAX_ABS_RNA_SEQ"]
        self.assertEqual(track, ["RNA_SEQ", "total RNA-seq", "pancreas", "UBERON:1", "INS", ""])
        self.assertAlmostEqual(value, -0.9, places=6)

    def test_chip_names_its_factor_or_mark_and_splice_sites_their_kind(self):
        replies = {
            "CHIP_TF": _Scores([[0.3]], [_track("TF ChIP-seq", "HepG2", transcription_factor="CTCF")]),
            "CHIP_HISTONE": _Scores([[0.3]], [_track("Histone ChIP-seq", "HepG2", histone_mark="H3K27ac")]),
            "SPLICE_SITES": _Scores([[0.1, 0.4]], [{"name": "donor"}, {"name": "acceptor"}],
                                    [{"gene_name": "INS"}]),
        }
        out = reduce(replies)
        self.assertEqual(out["MAX_ABS_CHIP_TF"][5], "CTCF")
        self.assertEqual(out["MAX_ABS_CHIP_HISTONE"][5], "H3K27ac")
        self.assertEqual(out["MERGED_SPLICING"][:6], ["SPLICE_SITES", "", "", "", "INS", "acceptor"])

    def test_splicing_weighs_junctions_at_a_fifth(self):
        replies = {
            "SPLICE_SITE_USAGE": _Scores([[0.15]], [_track("total RNA-seq", "liver")]),
            "SPLICE_JUNCTIONS": _Scores([[0.6]], [_track("polyA plus RNA-seq", "pancreas")]),
        }
        self.assertEqual(reduce(replies)["MERGED_SPLICING"][0], "SPLICE_SITE_USAGE")
        replies["SPLICE_JUNCTIONS"] = _Scores([[0.9]], [_track("polyA plus RNA-seq", "pancreas")])
        self.assertEqual(reduce(replies)["MERGED_SPLICING"][0], "SPLICE_JUNCTIONS")

    @unittest.skipUnless(_AVI.exists() and _EXPLANATIONS.exists(),
                         "no insulin tracks under pipeline/data; run pipeline/fetch_tracks.py")
    def test_rejects_stale_or_misassigned_evidence(self):
        raw = _AVI.read_bytes()
        data = json.loads(_EXPLANATIONS.read_text())
        validate(raw, data)
        stale = copy.deepcopy(data)
        stale["impact_sha256"] = hashlib.sha256(b"wrong").hexdigest()
        with self.assertRaisesRegex(AssertionError, "digest"):
            validate(raw, stale)
        swapped = copy.deepcopy(data)
        swapped["positions"]["5294"][0][0] = -1
        with self.assertRaisesRegex(AssertionError, "allele score"):
            validate(raw, swapped)
        indel = copy.deepcopy(data)
        indel["positions"]["5294"][0][2][data["features"].index("IS_DELETION")] = 0.1
        with self.assertRaisesRegex(AssertionError, "indel"):
            validate(raw, indel)
        orphan = copy.deepcopy(data)
        orphan["positions"]["5294"][0][3].append([0, len(data["tracks"]), 0.1])
        with self.assertRaisesRegex(AssertionError, "winner"):
            validate(raw, orphan)


if __name__ == "__main__":
    unittest.main()
