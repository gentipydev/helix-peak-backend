import json

import pytest

from app.config import settings


@pytest.fixture
def evidence_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(settings, "impact_explanations_dir", tmp_path)
    return tmp_path


def test_returns_exact_saved_contract(client, evidence_dir, efetch_raising):
    efetch_raising(AssertionError("This route must not call NCBI"))
    data = {"schema_version": 1, "gene": "INS", "accession": "NG_007114",
            "scorer": "AVI_SCORE_FEATURE_IMPORTANCE", "positions": {"4999": [[1.2, [[0, -0.4]]]]}}
    (evidence_dir / "insulin.json").write_text(json.dumps(data))
    response = client.get("/gene/NG_007114/INS/impact-explanations")
    assert response.status_code == 200
    assert response.json() == data
    assert client.get("/gene/NG_007114/INS-IGF2/impact-explanations").status_code == 404
    assert client.get("/gene/NG_000000/INS/impact-explanations").status_code == 404


def test_unincluded_and_malformed_are_distinct(client, evidence_dir):
    assert client.get("/gene/NG_007114/INS/impact-explanations").status_code == 404
    (evidence_dir / "insulin.json").write_text("not json")
    assert client.get("/gene/NG_007114/INS/impact-explanations").status_code == 503


def test_local_generated_bundle_matches_mobile(client):
    path = settings.impact_explanations_dir / "insulin.json"
    if not path.exists():
        pytest.skip("Companion mobile checkout not present")
    data = json.loads(path.read_text())
    response = client.get(f"/gene/{data['accession']}/{data['gene']}/impact-explanations")
    assert response.status_code == 200
    assert response.json() == data
