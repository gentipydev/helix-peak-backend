"""Offline cross-asset gates for the compact AVI contribution payloads."""
import hashlib
import json
import math


def validate(impact_bytes: bytes, data: dict) -> None:
    impact = json.loads(impact_bytes)
    assert data["schema_version"] == 1, "unsupported schema"
    assert data["scorer"] == "AVI_SCORE_FEATURE_IMPORTANCE", "wrong scorer"
    assert data["scope"] == "across_atlas_genes_and_biosamples", "wrong scope"
    assert data["units"] == "raw_score_attribution", "wrong units"
    assert data["selection"] == "top_3_absolute_signed", "wrong selection"
    assert data["impact_sha256"] == hashlib.sha256(impact_bytes).hexdigest(), "stale impact digest"
    for field in ("gene", "uniprot", "accession", "assembly", "annotation", "chromosome",
                  "transcript", "start", "sequence", "complemented", "orientation", "runs", "alt_order"):
        assert data[field] == impact[field], f"mismatched {field}"
    features = data["features"]
    assert features and len(features) == len(set(features)), "duplicate or empty features"
    assert data["positions"].keys() == impact["positions"].keys(), "incomplete or extra coverage"
    for position, alternatives in data["positions"].items():
        assert len(alternatives) == 3, f"wrong alternative count at {position}"
        for expected, row in zip(impact["positions"][position], alternatives):
            assert len(row) == 2 and row[0] == expected, f"stale allele score at {position}"
            values = row[1]
            assert len(values) <= 3, "too many contributions"
            seen = set()
            previous = math.inf
            for index, value in values:
                assert type(index) is int and 0 <= index < len(features) and index not in seen, "invalid feature"
                assert math.isfinite(value) and 0 < abs(value) <= previous, "invalid contribution"
                seen.add(index)
                previous = abs(value)
