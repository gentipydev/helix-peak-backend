"""Tests for `/proteins/suggest`.

No database, the same trade `test_catalog.py` makes: a fake pool answers each
statement by its shape, so what the route promises -- statuses, the short
and near-miss rules, the parameters it binds, and an unreadable index reported
as unavailable -- is testable without a socket. The SQL itself is exercised
against Supabase by `scripts/check_suggest.py`.
"""

import pytest

from app import db


def _row(uniprot="P01308", gene="INS", name="Insulin", length=110, buildable=True,
         reason=None, slug="insulin", display="Insulin", catalog_order=0, tier=0):
    return (uniprot, gene, name, length, buildable, reason, slug, display,
            catalog_order, tier)


LISTED = _row()
READY = _row(uniprot="P61769", gene="B2M", name="Beta-2-microglobulin", length=119,
             slug="b2m", display="Beta-2-microglobulin", catalog_order=None, tier=3)
BUILDABLE = _row(uniprot="P06213", gene="INSR", name="Insulin receptor", length=1382,
                 slug=None, display=None, catalog_order=None, tier=3)
UNAVAILABLE = _row(uniprot="Q13625", gene="TP53BP2",
                   name="Apoptosis-stimulating of p53 protein 2", length=1128,
                   buildable=False,
                   reason="MANE Select encodes isoform Q13625-3, and UniProt numbers "
                          "its features on the canonical sequence.",
                   slug=None, display=None, catalog_order=None, tier=3)
NEAR = _row(uniprot="P14735", gene="IDE", name="Insulin-degrading enzyme", length=1019,
            slug=None, display=None, catalog_order=None, tier=4)


class FakeCursor:
    def __init__(self, rows):
        self._rows = rows

    def fetchone(self):
        return self._rows[0] if self._rows else None

    def fetchall(self):
        return list(self._rows)


class FakePool:
    def __init__(self, ranked=(), near=(), release=("2026_03", "v1.5"), failing=False):
        self.ranked = list(ranked)
        self.near = list(near)
        self.release = release
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
        if flat.startswith("with hit as"):
            return FakeCursor(self.ranked)
        if flat.startswith("with near as"):
            return FakeCursor(self.near)
        if "from protein_index_release" in flat:
            return FakeCursor([self.release] if self.release else [])
        return FakeCursor([])


@pytest.fixture
def fake_pool(monkeypatch):
    def _install(**kwargs):
        pool = FakePool(**kwargs)
        monkeypatch.setattr(db, "pool", pool)
        return pool

    return _install


@pytest.fixture
def no_pool(monkeypatch):
    monkeypatch.setattr(db, "pool", None)


# --- unavailable, never empty -----------------------------------------------


def test_no_database_reports_unavailable(client, no_pool):
    assert client.get("/proteins/suggest?q=ins").status_code == 503


def test_a_broken_database_reports_unavailable(client, fake_pool):
    fake_pool(failing=True)
    assert client.get("/proteins/suggest?q=ins").status_code == 503


def test_nothing_typed_asks_nothing(client, no_pool):
    """Punctuation alone normalises to nothing, and nothing is not a query."""
    for typed in ("", "   ", "%-_"):
        body = client.get("/proteins/suggest", params={"q": typed}).json()
        assert body["suggestions"] == []


# --- statuses ---------------------------------------------------------------


def test_each_status_says_what_the_app_can_do(client, fake_pool):
    fake_pool(ranked=[LISTED, READY, BUILDABLE, UNAVAILABLE])
    body = client.get("/proteins/suggest?q=ins").json()
    assert body["release"] == "UniProt 2026_03 · MANE v1.5"
    listed, ready, buildable, unavailable = body["suggestions"]

    assert (listed["status"], listed["slug"], listed["display"]) == ("listed", "insulin", "Insulin")
    assert (ready["status"], ready["slug"]) == ("ready", "b2m")
    # A buildable protein carries the slug a build would give it.
    assert (buildable["status"], buildable["slug"], buildable["reason"]) == (
        "buildable", "insr", None)
    assert unavailable["status"] == "unavailable"
    assert unavailable["slug"] is None
    assert "Q13625-3" in unavailable["reason"]


def test_a_protein_with_no_gene_says_so_with_null(client, fake_pool):
    fake_pool(ranked=[_row(uniprot="Q8N8Q1", gene="", name="Uncharacterized protein",
                           buildable=False, reason="UniProt names no gene.", slug=None,
                           display=None, catalog_order=None, tier=3)])
    (only,) = client.get("/proteins/suggest?q=unchar").json()["suggestions"]
    assert only["gene"] is None
    assert only["status"] == "unavailable"


# --- what is asked ----------------------------------------------------------


def test_the_needle_is_normalised_and_bound_as_a_range(client, fake_pool):
    pool = fake_pool(ranked=[LISTED] * 12)
    client.get("/proteins/suggest", params={"q": "  Insulin-Deg%"})
    (flat, params) = pool.statements[0]
    assert flat.startswith("with hit as")
    assert params["needle"] == "insulin deg"
    assert (params["lo"], params["hi"]) == ("insulin deg", "insulin deh")
    assert params["short"] is False
    assert params["limit"] == 12


def test_a_short_needle_matches_symbols_and_synonyms_only(client, fake_pool):
    pool = fake_pool(ranked=[])
    client.get("/proteins/suggest?q=hb")
    assert pool.statements[0][1]["short"] is True
    # Two letters are too few for a near miss to mean anything.
    assert not any(flat.startswith("with near as") for flat, _ in pool.statements)


def test_a_typo_that_matches_nothing_is_answered_by_near_misses(client, fake_pool):
    pool = fake_pool(ranked=[], near=[LISTED, NEAR])
    body = client.get("/proteins/suggest?q=insuln").json()
    assert [s["gene"] for s in body["suggestions"]] == ["INS", "IDE"]
    assert any(flat.startswith("with near as") for flat, _ in pool.statements)


def test_any_prefix_match_asks_for_no_near_misses(client, fake_pool):
    """One exact accession is the answer, not the first of twelve lookalikes."""
    pool = fake_pool(ranked=[LISTED], near=[NEAR])
    body = client.get("/proteins/suggest?q=P01308").json()
    assert [s["gene"] for s in body["suggestions"]] == ["INS"]
    assert not any(flat.startswith("with near as") for flat, _ in pool.statements)


def test_an_exact_match_is_a_whole_term_never_a_word(client, fake_pool):
    pool = fake_pool(ranked=[])
    client.get("/proteins/suggest?q=ins")
    flat = pool.statements[0][0]
    assert "t.term = %(needle)s and t.kind < 5 then 0" in flat


def test_the_limit_is_bounded(client, fake_pool):
    fake_pool(ranked=[])
    assert client.get("/proteins/suggest?q=ins&limit=0").status_code == 422
    assert client.get("/proteins/suggest?q=ins&limit=26").status_code == 422
    assert client.get("/proteins/suggest?q=ins&limit=25").status_code == 200


def test_an_index_with_no_release_row_still_answers(client, fake_pool):
    fake_pool(ranked=[LISTED], release=None)
    body = client.get("/proteins/suggest?q=ins").json()
    assert body["release"] is None
    assert len(body["suggestions"]) == 1
