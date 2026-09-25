import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.impact.bake_explanations import compact  # noqa: E402
from pipeline.impact.bake_impact import BakeError  # noqa: E402
from pipeline.impact.check_explanations import validate  # noqa: E402
from pipeline.paths import DATA  # noqa: E402

# Insulin's two tracks, from storage by way of `fetch_tracks.py`.
_AVI = DATA / "assets/impact/insulin_avi.json"
_EXPLANATIONS = DATA / "assets/impact_explanations/insulin.json"


class ExplanationsTest(unittest.TestCase):
    def test_selects_by_absolute_magnitude_without_losing_negative_sign(self):
        self.assertEqual(compact([0.2, -0.9, 0.5, -0.6]), [[1, -0.9], [3, -0.6], [2, 0.5]])
        self.assertEqual(compact([0.0, 0.0]), [])

    def test_rejects_nonfinite_model_data(self):
        for value in (float("nan"), float("inf"), -float("inf")):
            with self.assertRaises(BakeError):
                compact([value])

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


if __name__ == "__main__":
    unittest.main()
