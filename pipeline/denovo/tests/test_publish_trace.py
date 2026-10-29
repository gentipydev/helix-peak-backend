from __future__ import annotations

import hashlib
import io
import json
import sys
import urllib.request
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[3]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.denovo import publish_trace
from pipeline.denovo.check_trace import validate

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "toy_trace.json"


@pytest.fixture
def no_network(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("a dry run must not call out")

    monkeypatch.setattr(urllib.request, "urlopen", forbidden)


def test_the_row_names_the_bytes_by_their_digest():
    payload = FIXTURE.read_bytes()
    row = publish_trace.row_for(payload, validate(payload), commit="abc", made="2026-10-30")
    digest = hashlib.sha256(payload).hexdigest()
    assert row["name"] == "toy"
    assert row["bucket"] == "tracks"
    assert row["object_path"] == "peptide-lm/toy.%s.json" % digest[:12]
    assert row["sha256"] == digest and row["bytes"] == len(payload)
    assert row["format"] == "json"
    assert row["provenance"]["commit"] == "abc"
    assert row["provenance"]["made"] == "2026-10-30"
    assert row["provenance"]["steps"] == 20
    assert row["provenance"]["engine"] == "pipeline/denovo/train.py train_step"


def test_a_dry_run_plans_and_calls_nobody(no_network, capsys):
    assert publish_trace.main(["--trace", str(FIXTURE), "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "tracks/peptide-lm/toy." in out
    assert "checkpoints 1, " in out


def test_a_dry_run_against_a_service_says_whether_it_serves_these_bytes(monkeypatch, capsys):
    digest = hashlib.sha256(FIXTURE.read_bytes()).hexdigest()

    def answer(request, timeout):
        assert request.full_url == "https://service.test/peptide-lm/runs/toy"
        return io.BytesIO(json.dumps({"sha256": digest}).encode())

    monkeypatch.setattr(urllib.request, "urlopen", answer)
    assert publish_trace.main(
        ["--trace", str(FIXTURE), "--dry-run", "--against", "https://service.test"]) == 0
    assert "already serves these bytes" in capsys.readouterr().out


def test_a_run_that_fails_its_check_is_never_published(tmp_path, no_network, capsys):
    bad = tmp_path / "trace.json"
    bad.write_bytes(b'{"schema": 9}')
    assert publish_trace.main(["--trace", str(bad)]) == 1
    assert "REFUSED" in capsys.readouterr().err


def test_publishing_needs_the_storage_and_database_credentials(monkeypatch, no_network):
    for name in ("SUPABASE_URL", "SUPABASE_SERVICE_KEY", "DATABASE_URL"):
        monkeypatch.delenv(name, raising=False)
    with pytest.raises(SystemExit, match="must be set"):
        publish_trace.main(["--trace", str(FIXTURE)])


def test_against_is_only_for_a_dry_run():
    with pytest.raises(SystemExit):
        publish_trace.main(["--trace", str(FIXTURE), "--against", "https://service.test"])
