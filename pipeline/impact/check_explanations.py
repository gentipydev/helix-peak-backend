"""Offline cross-asset gates for the AVI contribution payloads (schema 2)."""
import hashlib
import json
import math

INDELS = ("IS_INSERTION", "IS_DELETION")


def validate(impact_bytes: bytes, data: dict) -> None:
    impact = json.loads(impact_bytes)
    assert data["schema_version"] == 2, "unsupported schema"
    assert data["scorer"] == "AVI_SCORE_FEATURE_IMPORTANCE", "wrong scorer"
    assert data["scope"] == "across_atlas_genes_and_biosamples", "wrong scope"
    assert data["units"] == "raw_score_attribution", "wrong units"
    assert data["selection"] == "all", "wrong selection"
    assert data["impact_sha256"] == hashlib.sha256(impact_bytes).hexdigest(), "stale impact digest"
    for field in ("gene", "uniprot", "accession", "assembly", "annotation", "chromosome",
                  "transcript", "start", "sequence", "complemented", "orientation", "runs", "alt_order"):
        assert data[field] == impact[field], f"mismatched {field}"
    features = data["features"]
    assert len(features) == 18 and len(features) == len(set(features)), "not the 18 features"
    indels = [features.index(name) for name in INDELS]
    tracks = data["tracks"]
    assert all(len(t) == 6 and all(type(s) is str for s in t) for t in tracks), "invalid track"
    assert len({tuple(t) for t in tracks}) == len(tracks), "duplicate track"
    offsets = data["raw_minus_sum"]
    baseline = data["baseline"]
    assert baseline is None or math.isclose(baseline, offsets["median"]), "baseline is not the median offset"
    assert data["positions"].keys() == impact["positions"].keys(), "incomplete or extra coverage"
    for position, alternatives in data["positions"].items():
        assert len(alternatives) == 3, f"wrong alternative count at {position}"
        for expected, row in zip(impact["positions"][position], alternatives):
            assert len(row) == 4 and row[0] == expected, f"stale allele score at {position}"
            _, raw, values, won = row
            assert math.isfinite(raw), "invalid raw score"
            assert len(values) == 18 and all(math.isfinite(v) for v in values), "invalid contribution"
            assert all(values[i] == 0 for i in indels), f"indel attribution at {position}"
            offset = raw - math.fsum(values)
            assert offsets["min"] - 1e-5 <= offset <= offsets["max"] + 1e-5, f"offset out of range at {position}"
            seen = set()
            for feature, track, value in won:
                assert type(feature) is int and 0 <= feature < 18 and feature not in seen, "invalid winner feature"
                assert type(track) is int and 0 <= track < len(tracks), "invalid winner track"
                assert math.isfinite(value), "invalid winner value"
                seen.add(feature)
