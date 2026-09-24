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


def test_missing_directory_is_a_fault_not_a_gene_answer(client, tmp_path, monkeypatch):
    """An unmounted volume must not answer "not included for this gene".

    The deployed image copies only `app/`, so the default directory -- which is
    built from `__file__` and expects the companion checkout beside it --
    resolves to a path that is not there. `Path.glob` yields nothing for a
    missing directory rather than raising, so every request used to report a
    404 while insulin, hemoglobin and cftr all ship expecting a payload.
    """
    monkeypatch.setattr(settings, "impact_explanations_dir", tmp_path / "not-mounted")
    response = client.get("/gene/NG_007114/INS/impact-explanations")
    assert response.status_code == 503


class _TrackPool:
    """Answers the track lookup with one (state, bucket, object_path) row."""

    def __init__(self, row=None, failing=False):
        self.row = row
        self.failing = failing

    def connection(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def execute(self, sql, params=()):
        if self.failing:
            raise OSError("connection closed")
        assert "protein_track" in sql
        self.params = params
        return type("C", (), {"fetchone": lambda _self: self.row})()


def test_a_stored_payload_redirects_instead_of_being_carried(client, monkeypatch):
    """cftr's payload is 3.98 MB. It does not belong in 512 MB of memory."""
    from app import db
    from app.config import settings as live

    monkeypatch.setattr(live, "supabase_url", "https://ref.supabase.co")
    monkeypatch.setattr(db, "pool",
                        _TrackPool(("ready", "tracks", "impact_explanations/cftr.abc.json")))
    response = client.get("/gene/NG_016465/CFTR/impact-explanations",
                          follow_redirects=False)
    assert response.status_code == 307
    assert response.headers["location"] == (
        "https://ref.supabase.co/storage/v1/object/public/tracks/"
        "impact_explanations/cftr.abc.json"
    )


def test_the_lookup_is_by_accession_and_gene(client, monkeypatch):
    from app import db
    from app.config import settings as live

    monkeypatch.setattr(live, "supabase_url", "https://ref.supabase.co")
    pool = _TrackPool(("ready", "tracks", "impact_explanations/cftr.abc.json"))
    monkeypatch.setattr(db, "pool", pool)
    client.get("/gene/NG_016465/CFTR/impact-explanations", follow_redirects=False)
    assert pool.params == ("NG_016465", "CFTR")


def test_no_stored_row_falls_back_to_the_directory(client, evidence_dir, monkeypatch):
    """A checkout with no database still serves what it has on disk."""
    from app import db

    monkeypatch.setattr(db, "pool", _TrackPool(row=None))
    data = {"schema_version": 1, "gene": "INS", "accession": "NG_007114",
            "scorer": "AVI_SCORE_FEATURE_IMPORTANCE", "positions": {}}
    (evidence_dir / "insulin.json").write_text(json.dumps(data))
    response = client.get("/gene/NG_007114/INS/impact-explanations")
    assert response.status_code == 200
    assert response.json() == data


def test_a_broken_database_falls_back_rather_than_failing(client, evidence_dir,
                                                          monkeypatch):
    from app import db

    monkeypatch.setattr(db, "pool", _TrackPool(failing=True))
    data = {"schema_version": 1, "gene": "INS", "accession": "NG_007114",
            "scorer": "AVI_SCORE_FEATURE_IMPORTANCE", "positions": {}}
    (evidence_dir / "insulin.json").write_text(json.dumps(data))
    assert client.get("/gene/NG_007114/INS/impact-explanations").status_code == 200


def test_a_catalogued_gene_with_no_explanations_is_a_404(client, monkeypatch):
    """Not a fallback case, and not a 503.

    The deployed image has no assets directory, so falling through for a gene
    the catalog has already answered about would turn "not included for this
    gene" into "temporarily unavailable" for all seventeen of them.
    """
    from app import db

    monkeypatch.setattr(db, "pool", _TrackPool(("absent", None, None)))
    response = client.get("/gene/NG_017013/TP53/impact-explanations")
    assert response.status_code == 404
    assert "not included" in response.json()["detail"]


def test_a_stored_track_this_service_cannot_address_is_a_fault(client, monkeypatch):
    """Redirecting to a URL it cannot build would send the client nowhere."""
    from app import db
    from app.config import settings as live

    monkeypatch.setattr(live, "supabase_url", None)
    monkeypatch.setattr(db, "pool",
                        _TrackPool(("ready", "tracks", "impact_explanations/cftr.abc.json")))
    assert client.get("/gene/NG_016465/CFTR/impact-explanations").status_code == 503
