"""`/assembly/{slug}/tracks`, and the wall between assemblies and the catalog.

A fake pool answers by SQL shape, as `test_catalog.py`'s does, so nothing here
dials a database.
"""

import pathlib

import pytest

from app import assemblies, db
from app.config import settings

APP = pathlib.Path(__file__).resolve().parents[1] / "app"


class FakeCursor:
    def __init__(self, rows):
        self._rows = rows

    def fetchone(self):
        return self._rows[0] if self._rows else None

    def fetchall(self):
        return self._rows


class FakePool:
    """Knows [slugs] as assemblies, with [track_rows] for any of them."""

    def __init__(self, slugs=(), track_rows=(), failing=False):
        self.slugs = set(slugs)
        self.track_rows = list(track_rows)
        self.failing = failing

    def connection(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def execute(self, sql, params=()):
        if self.failing:
            raise OSError("relation \"assembly\" does not exist")
        flat = " ".join(sql.split())
        if flat.startswith("select 1 from assembly"):
            return FakeCursor([(1,)] if params[0] in self.slugs else [])
        if "from assembly_track" in flat:
            return FakeCursor(self.track_rows)
        raise AssertionError("an assembly read only its own tables: " + flat)


@pytest.fixture
def pool(monkeypatch):
    def _install(**kwargs):
        fake = FakePool(**kwargs)
        monkeypatch.setattr(db, "pool", fake)
        return fake

    return _install


def test_an_unknown_assembly_is_not_found(client, pool):
    pool(slugs=["hemoglobin-a"])
    assert client.get("/assembly/insulin/tracks").status_code == 404


def test_an_assembly_with_no_rows_has_its_morph_absent(client, pool):
    pool(slugs=["hemoglobin-a"])
    body = client.get("/assembly/hemoglobin-a/tracks").json()
    assert list(body) == ["morph"]
    assert body["morph"]["state"] == "absent"
    assert body["morph"]["url"] is None


def test_a_ready_morph_is_a_public_url(client, pool, monkeypatch):
    monkeypatch.setattr(settings, "supabase_url", "https://example.supabase.co")
    pool(slugs=["hemoglobin-a"], track_rows=[(
        "morph", "ready", None, "tracks", "morph/hemoglobin-a.a629e4ad525c.json",
        171064, "a629e4ad525c" + "0" * 52, None, "json", {"schema_version": 1},
    )])
    morph = client.get("/assembly/hemoglobin-a/tracks").json()["morph"]
    assert morph["state"] == "ready"
    assert morph["url"] == ("https://example.supabase.co/storage/v1/object/public/"
                            "tracks/morph/hemoglobin-a.a629e4ad525c.json")
    assert morph["provenance"] == {"schema_version": 1}


def test_unreadable_assemblies_are_a_503_not_a_404(client, pool):
    pool(failing=True)
    assert client.get("/assembly/hemoglobin-a/tracks").status_code == 503


def test_no_database_is_a_503(client, monkeypatch):
    monkeypatch.setattr(db, "pool", None)
    assert client.get("/assembly/hemoglobin-a/tracks").status_code == 503


def test_nothing_that_serves_the_catalog_reads_an_assembly():
    # The wall that keeps /catalog at twenty: the catalog, its search and the
    # suggestions read the protein tables, and never an assembly's.
    for module in ("catalog.py", "suggest.py", "tracks.py", "protein_index.py"):
        assert "assembly" not in (APP / module).read_text(encoding="utf-8"), module
    # And an assembly reads only its own.
    source = (APP / "assemblies.py").read_text(encoding="utf-8")
    for table in ("from protein ", "from protein_track", "from protein_alias"):
        assert table not in source
    assert assemblies.KINDS == ("morph",)
