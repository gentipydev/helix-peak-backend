"""Tests for `POST /proteins/resolve` and `GET /proteins/resolve/{gene}`.

The trade `test_catalog.py` and `test_suggest.py` make: a fake pool answers
each statement by its shape, so what the routes promise -- every state, when a
request is written, when the resolver is woken, the daily cap, and an
unreadable database reported as unavailable -- is testable without a socket.
The statements themselves run against Postgres in
`pipeline/resolver/test_worker_pg.py`.
"""

import pytest

from app import db, resolves
from app.config import settings

BRCA1 = ("P38398", "BRCA1", True, None)
UNBUILDABLE = ("Q13625", "TP53BP2", False,
               "MANE Select encodes isoform Q13625-3, and UniProt numbers its features "
               "on the canonical sequence.")


class FakeCursor:
    def __init__(self, rows):
        self._rows = rows

    def fetchone(self):
        return self._rows[0] if self._rows else None

    def fetchall(self):
        return list(self._rows)


class FakePool:
    def __init__(self, protein=None, index=None, latest=None, today=0, queued=True,
                 failing=False):
        self.protein = protein
        self.index = index
        self.latest = latest
        self.today = today
        self.queued = queued
        self.failing = failing
        self.inserted = []

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
        if flat.startswith("select slug from protein"):
            return FakeCursor([(self.protein,)] if self.protein else [])
        if flat.startswith("select uniprot, gene, buildable"):
            return FakeCursor([self.index] if self.index else [])
        if flat.startswith("select state, reason from resolve_request"):
            return FakeCursor([self.latest] if self.latest else [])
        if flat.startswith("select count(*) from resolve_request"):
            return FakeCursor([(self.today,)])
        if flat.startswith("insert into resolve_request"):
            self.inserted.append(params)
            return FakeCursor([(1,)] if self.queued else [])
        raise AssertionError(f"unexpected statement: {flat}")


@pytest.fixture
def pool(monkeypatch):
    def _install(**kwargs):
        fake = FakePool(**kwargs)
        monkeypatch.setattr(db, "pool", fake)
        return fake
    return _install


@pytest.fixture
def woken(monkeypatch):
    calls = []
    monkeypatch.setattr(resolves, "wake", lambda: calls.append(True))
    return calls


def _ask(client, gene="BRCA1", **extra):
    return client.post("/proteins/resolve", json={"gene": gene, **extra})


# ------------------------------------------------------------------ POST


def test_a_gene_the_catalog_lists_is_ready_where_it_shipped(client, pool, woken):
    fake = pool(protein="insulin", index=("P01308", "INS", True, None))
    response = _ask(client, "INS")
    assert response.status_code == 200
    assert response.json() == {"slug": "insulin", "state": "ready", "reason": None}
    assert fake.inserted == [] and woken == []


def test_a_new_ask_is_queued_and_wakes_the_resolver(client, pool, woken):
    fake = pool(index=BRCA1)
    response = _ask(client, "brca1")
    assert response.status_code == 202
    assert response.json() == {"slug": "brca1", "state": "pending", "reason": None}
    assert fake.inserted == [("BRCA1", "P38398", "brca1")]
    assert woken == [True]


def test_an_ask_already_in_flight_is_pending_and_writes_nothing(client, pool, woken):
    fake = pool(index=BRCA1, latest=("running", None))
    response = _ask(client)
    assert response.status_code == 202
    assert response.json()["state"] == "pending"
    assert fake.inserted == [] and woken == []


def test_a_race_lost_to_another_reader_is_pending_without_a_second_wake(client, pool, woken):
    fake = pool(index=BRCA1, queued=False)
    response = _ask(client)
    assert response.status_code == 202
    assert fake.inserted == [("BRCA1", "P38398", "brca1")]
    assert woken == []


def test_a_refusal_stands(client, pool, woken):
    fake = pool(index=BRCA1, latest=("refused", "Exons alone are 81,000 bp, over the budget."))
    response = _ask(client)
    assert response.status_code == 200
    assert response.json() == {"slug": None, "state": "refused",
                               "reason": "Exons alone are 81,000 bp, over the budget."}
    assert fake.inserted == [] and woken == []


def test_a_failed_request_is_asked_again(client, pool, woken):
    fake = pool(index=BRCA1, latest=("failed", "rest.uniprot.org did not answer"))
    response = _ask(client)
    assert response.status_code == 202
    assert fake.inserted == [("BRCA1", "P38398", "brca1")]
    assert woken == [True]


def test_an_unbuildable_protein_says_why(client, pool, woken):
    fake = pool(index=UNBUILDABLE)
    response = _ask(client, "TP53BP2")
    assert response.status_code == 200
    assert response.json() == {"slug": None, "state": "unavailable", "reason": UNBUILDABLE[3]}
    assert fake.inserted == []


def test_a_gene_the_index_does_not_hold_is_not_found(client, pool):
    pool()
    response = _ask(client, "NOTAGENE")
    assert response.status_code == 404
    assert response.json()["detail"] == \
        "No reviewed human protein made by 'NOTAGENE' in the index."


def test_another_organism_is_not_found(client, pool):
    fake = pool(index=BRCA1)
    assert _ask(client, taxon=10090).status_code == 404
    assert fake.inserted == []


def test_an_empty_gene_is_not_a_request(client, pool):
    pool(index=BRCA1)
    assert _ask(client, "").status_code == 422


def test_a_days_requests_are_capped(client, pool, woken, monkeypatch):
    monkeypatch.setattr(settings, "resolves_per_day", 50)
    fake = pool(index=BRCA1, today=50)
    response = _ask(client)
    assert response.status_code == 429
    assert response.json()["detail"] == \
        "The service builds 50 proteins a day, and today's are taken. Ask again tomorrow."
    assert fake.inserted == [] and woken == []


def test_the_cap_never_stops_an_answer_that_costs_nothing(client, pool, monkeypatch):
    monkeypatch.setattr(settings, "resolves_per_day", 0)
    pool(protein="insulin", index=("P01308", "INS", True, None))
    assert _ask(client, "INS").status_code == 200


def test_an_unreadable_database_is_unavailable_not_absent(client, pool, woken):
    pool(failing=True)
    response = _ask(client)
    assert response.status_code == 503
    assert woken == []


def test_no_database_is_unavailable(client, monkeypatch):
    monkeypatch.setattr(db, "pool", None)
    assert _ask(client).status_code == 503


# ------------------------------------------------------------------- GET


def test_a_read_never_writes(client, pool, woken):
    fake = pool(index=BRCA1)
    response = client.get("/proteins/resolve/BRCA1")
    assert response.status_code == 200
    assert response.json() == {"slug": "brca1", "state": "buildable", "reason": None}
    assert fake.inserted == [] and woken == []


@pytest.mark.parametrize("latest, said", [
    (("queued", None), {"slug": "brca1", "state": "pending", "reason": None}),
    (("running", None), {"slug": "brca1", "state": "pending", "reason": None}),
    (("failed", "UniProt did not answer"),
     {"slug": "brca1", "state": "failed", "reason": "UniProt did not answer"}),
    (("refused", "Too long."), {"slug": None, "state": "refused", "reason": "Too long."}),
])
def test_a_read_says_where_a_request_is(client, pool, latest, said):
    pool(index=BRCA1, latest=latest)
    assert client.get("/proteins/resolve/BRCA1").json() == said


def test_a_read_of_a_resolved_protein_is_ready(client, pool):
    pool(protein="brca1", index=BRCA1, latest=("done", None))
    assert client.get("/proteins/resolve/BRCA1").json() == \
        {"slug": "brca1", "state": "ready", "reason": None}


def test_a_read_of_an_unknown_gene_is_not_found(client, pool):
    pool()
    assert client.get("/proteins/resolve/NOTAGENE").status_code == 404


def test_a_read_of_an_unreadable_database_is_unavailable(client, pool):
    pool(failing=True)
    assert client.get("/proteins/resolve/BRCA1").status_code == 503


# ------------------------------------------------------------------ wake


class FakeResponse:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_waking_posts_to_the_endpoint_with_its_proxy_token(monkeypatch):
    sent = []
    monkeypatch.setattr(settings, "modal_wake_url", "https://w--helix-peak-resolver-wake.modal.run")
    monkeypatch.setattr(settings, "modal_key", "wk-key")
    monkeypatch.setattr(settings, "modal_secret", "ws-secret")
    monkeypatch.setattr(resolves.urllib.request, "urlopen",
                        lambda request, timeout: sent.append(request) or FakeResponse())
    resolves.wake()
    (request,) = sent
    assert request.full_url == "https://w--helix-peak-resolver-wake.modal.run"
    assert request.get_method() == "POST"
    assert request.get_header("Modal-key") == "wk-key"
    assert request.get_header("Modal-secret") == "ws-secret"


def test_waking_unconfigured_does_nothing(monkeypatch):
    monkeypatch.setattr(settings, "modal_wake_url", None)

    def never(*args, **kwargs):
        raise AssertionError("woke with no URL")

    monkeypatch.setattr(resolves.urllib.request, "urlopen", never)
    resolves.wake()


def test_a_wake_that_fails_is_logged_and_never_raised(monkeypatch, caplog):
    monkeypatch.setattr(settings, "modal_wake_url", "https://w--helix-peak-resolver-wake.modal.run")

    def unreachable(*args, **kwargs):
        raise OSError("modal.run did not answer")

    monkeypatch.setattr(resolves.urllib.request, "urlopen", unreachable)
    resolves.wake()
    assert "Could not wake the resolver" in caplog.text
