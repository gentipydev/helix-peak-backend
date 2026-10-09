"""Tests for `/proteins/built`.

No database, as in `test_suggest.py`: a fake pool answers each statement by its
shape, and pages its rows the way the keyset does, so what the route promises --
the row, the order it keeps, which builds it leaves out, the cursor, and an
unreadable index reported as unavailable -- is testable without a socket. The
SQL itself is exercised against Postgres by `pipeline/resolver/test_worker_pg.py`.
"""

from datetime import datetime, timedelta, timezone

import pytest

from app import db

T0 = datetime(2026, 10, 9, 8, 45, 20, 800935, tzinfo=timezone.utc)


def _row(uniprot="P02766", gene="TTR", name="Transthyretin", length=147, slug="ttr",
         display="Transthyretin", resolved_at=T0):
    return (uniprot, gene, name, length, True, None, slug, display, None, resolved_at)


# Newest first, as the query orders them.
TTR = _row()
PRL = _row(uniprot="P01236", gene="PRL", name="Prolactin", length=227, slug="prl",
           display="Prolactin", resolved_at=T0 - timedelta(minutes=10))
HERC2 = _row(uniprot="O95714", gene="HERC2", name="E3 ubiquitin-protein ligase HERC2",
             length=4834, slug="herc2", display="E3 ubiquitin-protein ligase HERC2",
             resolved_at=T0 - timedelta(hours=19))


class FakeCursor:
    def __init__(self, rows):
        self._rows = rows

    def fetchall(self):
        return list(self._rows)


class FakePool:
    def __init__(self, rows=(), failing=False, building=(), stopped=()):
        self.rows = list(rows)
        self.failing = failing
        # Genes, lower-cased, with a build under way, and whose scoring was stopped.
        self.building = list(building)
        self.stopped = list(stopped)
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
        if "from protein p join protein_index i" in flat:
            at, slug = params["at"], params["slug"]
            after = [row for row in self.rows
                     if at is None or row[9] < at or (row[9] == at and row[6] > slug)]
            return FakeCursor(after[:params["limit"]])
        if flat.startswith("select lower(gene), 'building' from resolve_request"):
            return FakeCursor(
                [(gene, "building") for gene in self.building if gene in params["genes"]]
                + [(gene, "stopped") for gene in self.stopped if gene in params["genes"]])
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


def _pages(pool):
    """What each page query bound, in order."""
    return [params for flat, params in pool.statements
            if "from protein p join protein_index i" in flat]


def _built(client, **params):
    response = client.get("/proteins/built", params=params)
    assert response.status_code == 200
    return response.json()


# --- unavailable, never empty -----------------------------------------------


def test_no_database_reports_unavailable(client, no_pool):
    assert client.get("/proteins/built").status_code == 503


def test_a_broken_database_reports_unavailable(client, fake_pool):
    fake_pool(failing=True)
    assert client.get("/proteins/built").status_code == 503


# --- the rows ---------------------------------------------------------------


def test_each_row_is_a_ready_suggestion_that_opens_where_it_was_built(client, fake_pool):
    fake_pool(rows=[TTR])
    body = _built(client)
    assert body == {
        "proteins": [{
            "uniprot": "P02766", "gene": "TTR", "name": "Transthyretin",
            "display": "Transthyretin", "length": 147, "slug": "ttr", "status": "ready",
            "reason": None, "building": False, "stopped": False,
        }],
        "next": None,
    }


def test_newest_first_and_never_a_curated_row(client, fake_pool):
    pool = fake_pool(rows=[TTR, PRL, HERC2])
    body = _built(client)
    assert [s["slug"] for s in body["proteins"]] == ["ttr", "prl", "herc2"]
    (flat, params) = pool.statements[0]
    assert "where p.catalog_order is null" in flat
    assert flat.endswith("order by p.resolved_at desc, p.slug limit %(limit)s")
    # Its own index row, not every row its gene has.
    assert "on i.uniprot = p.uniprot and i.gene = p.gene" in flat


def test_a_build_under_way_is_left_out_and_a_stopped_one_kept(client, fake_pool):
    pool = fake_pool(rows=[TTR, PRL, HERC2], building=["prl"], stopped=["herc2"])
    body = _built(client)
    assert [(s["gene"], s["status"], s["building"], s["stopped"])
            for s in body["proteins"]] == [
        ("TTR", "ready", False, False), ("HERC2", "ready", False, True),
    ]
    (asked,) = [params for flat, params in pool.statements
                if flat.startswith("select lower(gene), 'building' from resolve_request")]
    assert asked == {"genes": ["herc2", "prl", "ttr"]}


def test_none_built_asks_nothing_about_builds(client, fake_pool):
    pool = fake_pool(rows=[])
    assert _built(client) == {"proteins": [], "next": None}
    assert not [flat for flat, _ in pool.statements if "resolve_request" in flat]


# --- pages ------------------------------------------------------------------


def test_next_leads_to_the_rest_and_is_null_on_the_last_page(client, fake_pool):
    pool = fake_pool(rows=[TTR, PRL, HERC2])
    first = _built(client, limit=2)
    assert [s["slug"] for s in first["proteins"]] == ["ttr", "prl"]
    assert first["next"] == "2026-10-09T08:35:20.800935Z|prl"
    # One more than the page is asked for, to know whether there is more.
    assert _pages(pool)[0] == {"at": None, "slug": None, "limit": 3}

    second = _built(client, limit=2, before=first["next"])
    assert [s["slug"] for s in second["proteins"]] == ["herc2"]
    assert second["next"] is None
    assert _pages(pool)[-1] == {"at": T0 - timedelta(minutes=10), "slug": "prl", "limit": 3}


def test_a_page_exactly_full_has_no_next(client, fake_pool):
    fake_pool(rows=[TTR, PRL])
    assert _built(client, limit=2)["next"] is None


def test_a_page_cut_short_by_a_build_under_way_still_has_a_next(client, fake_pool):
    fake_pool(rows=[TTR, PRL, HERC2], building=["prl"])
    first = _built(client, limit=2)
    assert [s["slug"] for s in first["proteins"]] == ["ttr"]
    assert first["next"].endswith("|prl")
    assert [s["slug"] for s in _built(client, limit=2, before=first["next"])["proteins"]] == \
        ["herc2"]


def test_a_cursor_keeps_the_microsecond_in_utc(client, fake_pool):
    later = T0.astimezone(timezone(timedelta(hours=2)))
    fake_pool(rows=[_row(resolved_at=later), PRL])
    assert _built(client, limit=1)["next"] == "2026-10-09T08:45:20.800935Z|ttr"


def test_a_cursor_this_service_did_not_give_is_refused(client, fake_pool):
    pool = fake_pool(rows=[TTR])
    for before in ("nonsense", "2026-10-09T08:45:20.800935Z", "2026-10-09T08:45:20.800935Z|",
                   "yesterday|ttr", "2026-10-09 08:45:20|ttr", "|ttr"):
        response = client.get("/proteins/built", params={"before": before})
        assert response.status_code == 422, before
    # Refused before anything was asked.
    assert pool.statements == []


def test_the_limit_is_bounded(client, fake_pool):
    pool = fake_pool(rows=[])
    assert client.get("/proteins/built?limit=0").status_code == 422
    assert client.get("/proteins/built?limit=201").status_code == 422
    assert client.get("/proteins/built?limit=200").status_code == 200
    client.get("/proteins/built")
    assert _pages(pool)[-1]["limit"] == 51
