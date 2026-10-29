"""Tests for `/peptide-lm/runs/{name}`.

No database: a fake pool answers the one statement by its shape, so what the
route promises -- the row with its storage URL, a missing run as a 404, and an
unreadable or unaddressable one as a 503 -- is testable without a socket.
"""

import pytest

from app import db
from app.config import settings

ROW = ("curated", "tracks", "peptide-lm/curated.6f3031599652.json", 2989968,
       "6f30315996523252df4c108ac2fd952682daa2deece6b1ad271eabc9af87595a", "json",
       {"steps": 200, "checkpoints": [1, 8, 14, 35, 55, 100, 200]})


class FakeCursor:
    def __init__(self, rows):
        self._rows = rows

    def fetchone(self):
        return self._rows[0] if self._rows else None


class FakePool:
    def __init__(self, rows=(), failing=False):
        self.rows = list(rows)
        self.failing = failing
        self.statements = []

    def connection(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def execute(self, sql, params=None):
        if self.failing:
            raise OSError("connection closed")
        flat = " ".join(sql.split())
        self.statements.append((flat, params))
        if flat.startswith("select name, bucket, object_path") and "from lm_run" in flat:
            return FakeCursor([row for row in self.rows if row[0] == params[0]])
        return FakeCursor([])


@pytest.fixture
def pool(monkeypatch):
    def _install(**kwargs):
        fake = FakePool(**kwargs)
        monkeypatch.setattr(db, "pool", fake)
        return fake

    return _install


@pytest.fixture
def storage(monkeypatch):
    monkeypatch.setattr(settings, "supabase_url", "https://ref.supabase.co")


def test_a_run_is_named_where_storage_holds_it(client, pool, storage, efetch_raising):
    efetch_raising(AssertionError("a training run never asks NCBI"))
    pool(rows=[ROW])
    response = client.get("/peptide-lm/runs/curated")
    assert response.status_code == 200
    assert response.json() == {
        "name": "curated",
        "url": "https://ref.supabase.co/storage/v1/object/public/tracks/"
               "peptide-lm/curated.6f3031599652.json",
        "bytes": 2989968,
        "sha256": ROW[4],
        "format": "json",
        "provenance": {"steps": 200, "checkpoints": [1, 8, 14, 35, 55, 100, 200]},
    }


def test_a_run_nobody_published_is_not_found(client, pool, storage):
    fake = pool(rows=[ROW])
    response = client.get("/peptide-lm/runs/corpus")
    assert response.status_code == 404
    assert "corpus" in response.json()["detail"]
    assert fake.statements[0][1] == ("corpus",)


def test_no_database_is_unavailable_not_missing(client, storage, monkeypatch):
    monkeypatch.setattr(db, "pool", None)
    assert client.get("/peptide-lm/runs/curated").status_code == 503


def test_a_database_that_fails_is_unavailable(client, pool, storage):
    pool(rows=[ROW], failing=True)
    response = client.get("/peptide-lm/runs/curated")
    assert response.status_code == 503
    assert "could not be read" in response.json()["detail"]


def test_a_run_the_service_cannot_address_is_unavailable(client, pool, monkeypatch):
    monkeypatch.setattr(settings, "supabase_url", None)
    pool(rows=[ROW])
    response = client.get("/peptide-lm/runs/curated")
    assert response.status_code == 503
    assert "cannot address" in response.json()["detail"]
