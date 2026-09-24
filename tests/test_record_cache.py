"""Tests for the read-through record cache.

No database is involved: a fake pool stands in for psycopg so the policy --
what gets stored, what is served, and what happens when the cache misbehaves
-- is testable without a socket. The SQL itself is exercised against the real
database only by /health/db and by running the service.
"""

import pytest

from app import db, record_cache


class FakeCursor:
    def __init__(self, row):
        self._row = row

    def fetchone(self):
        return self._row


class FakePool:
    """Records every statement, and answers selects from ``rows``."""

    def __init__(self, rows=None, failing=False):
        self.rows = rows if rows is not None else {}
        self.failing = failing
        self.statements = []

    def connection(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def execute(self, sql, params=()):
        if self.failing:
            raise OSError("connection closed")
        self.statements.append((sql, params))
        if sql.lstrip().startswith("select record"):
            # psycopg hands back a row tuple, not the bare column.
            record = self.rows.get(params[0])
            return FakeCursor(None if record is None else (record,))
        return FakeCursor(None)


@pytest.fixture
def fake_pool(monkeypatch):
    def _install(**kwargs):
        pool = FakePool(**kwargs)
        monkeypatch.setattr(db, "pool", pool)
        return pool

    return _install


def test_a_hit_never_calls_ncbi(fake_pool, genbank_text, efetch_raising):
    fake_pool(rows={"NG_007114": genbank_text})
    efetch_raising(AssertionError("NCBI was called on a cache hit"))

    assert record_cache.fetch("NG_007114").id.startswith("NG_007114")


def test_a_miss_fetches_and_stores(fake_pool, genbank_text, efetch_returning):
    pool = fake_pool()
    efetch_returning(genbank_text)

    record_cache.fetch("NG_007114")
    writes = [p for sql, p in pool.statements if sql.lstrip().startswith("insert")]

    assert writes == [("NG_007114", genbank_text)]


def test_unparseable_upstream_is_never_stored(fake_pool, efetch_returning):
    """NCBI answers a bad id with plain text and HTTP 200; that must not be cached."""
    pool = fake_pool()
    efetch_returning("Error: CEFetchPApplication::proxy_stream()")

    with pytest.raises(ValueError):
        record_cache.fetch("NOT_A_REAL_ID")

    assert [sql for sql, _ in pool.statements if sql.lstrip().startswith("insert")] == []


def test_a_poisoned_row_is_evicted_and_read_through(
    fake_pool, genbank_text, efetch_returning
):
    pool = fake_pool(rows={"NG_007114": "not a genbank record"})
    efetch_returning(genbank_text)

    record = record_cache.fetch("NG_007114")
    deletes = [sql for sql, _ in pool.statements if sql.lstrip().startswith("delete")]

    assert record.id.startswith("NG_007114")
    assert len(deletes) == 1


def test_a_broken_database_still_serves_from_ncbi(
    fake_pool, genbank_text, efetch_returning
):
    """The whole point of the design: no cache means slower, not broken."""
    fake_pool(failing=True)
    efetch_returning(genbank_text)

    assert record_cache.fetch("NG_007114").id.startswith("NG_007114")


def test_no_database_configured_reads_through(
    monkeypatch, genbank_text, efetch_returning
):
    monkeypatch.setattr(db, "pool", None)
    efetch_returning(genbank_text)

    assert record_cache.fetch("NG_007114").id.startswith("NG_007114")


def test_the_ttl_is_part_of_the_lookup(fake_pool, genbank_text, efetch_returning):
    """A row older than the TTL must not count as a hit, so it is in the query."""
    pool = fake_pool()
    efetch_returning(genbank_text)

    record_cache.fetch("NG_007114")
    select = next(s for s in pool.statements if s[0].lstrip().startswith("select record"))

    assert "make_interval" in select[0]
    assert select[1][1] == 30


def test_the_cache_is_keyed_by_record_not_by_gene(
    fake_pool, genbank_text, efetch_raising, client
):
    """Two genes from one record must both be served by the single cached row."""
    fake_pool(rows={"NG_007114": genbank_text})
    efetch_raising(AssertionError("NCBI was called on a cache hit"))

    assert client.get("/gene/NG_007114/INS").status_code == 200
    assert client.get("/gene/NG_007114/TH").status_code == 200
