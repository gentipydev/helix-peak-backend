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
                 failing=False, request_build=None, bake_build=None, ask=None, bake=None,
                 resumable=False, stoppable=True):
        self.protein = protein
        self.index = index
        self.latest = latest
        self.today = today
        self.queued = queued
        self.failing = failing
        # `_REQUEST_BUILD`'s row and `_BAKE_BUILD`'s, where there is one.
        self.request_build = request_build
        self.bake_build = bake_build
        # `_LATEST_ASK`'s row and `_LATEST_BAKE`'s: (id, state, asker).
        self.ask = ask
        self.bake = bake
        # Whether `_RESUME` finds a stopped protein to score again, and whether
        # a stop's update still finds what it stops.
        self.resumable = resumable
        self.stoppable = stoppable
        self.inserted = []
        self.built = []
        self.written = []

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
        if flat.startswith("select r.state, extract(epoch"):
            self.built.append(("request", params))
            return FakeCursor([self.request_build] if self.request_build else [])
        if flat.startswith("select j.state, case when j.progress_at"):
            self.built.append(("bake", params))
            return FakeCursor([self.bake_build] if self.bake_build else [])
        if flat.startswith("select id, state, asker from resolve_request"):
            return FakeCursor([self.ask] if self.ask else [])
        if flat.startswith("select id, state, asker from bake_job"):
            return FakeCursor([self.bake] if self.bake else [])
        if flat.startswith("update resolve_request set state = 'stopped'"):
            self.written.append(("request stopped", params))
            return FakeCursor([(12.0,)] if self.stoppable else [])
        if flat.startswith("update bake_job set state = 'stopped'"):
            self.written.append(("bake stopped", params))
            return FakeCursor([(7,)] if self.stoppable else [])
        if flat.startswith("update protein_track set state = 'absent'"):
            self.written.append(("track absent", params))
            return FakeCursor([])
        if flat.startswith("insert into bake_job"):
            self.written.append(("resumed", params))
            return FakeCursor([(8,)] if self.resumable else [])
        if flat.startswith("update protein_track set state = 'pending'"):
            self.written.append(("track pending", params))
            return FakeCursor([])
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
    assert response.json() == {"slug": "insulin", "state": "ready", "reason": None, "build": None}
    assert fake.inserted == [] and woken == []


def test_a_new_ask_is_queued_and_wakes_the_resolver(client, pool, woken):
    fake = pool(index=BRCA1)
    response = _ask(client, "brca1")
    assert response.status_code == 202
    assert response.json() == {"slug": "brca1", "state": "pending", "reason": None, "build": None}
    assert fake.inserted == [("BRCA1", "P38398", "brca1", None)]
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
    assert fake.inserted == [("BRCA1", "P38398", "brca1", None)]
    assert woken == []


def test_a_refusal_stands(client, pool, woken):
    fake = pool(index=BRCA1, latest=("refused", "Exons alone are 81,000 bp, over the budget."))
    response = _ask(client)
    assert response.status_code == 200
    assert response.json() == {"slug": None, "state": "refused",
                               "reason": "Exons alone are 81,000 bp, over the budget.",
                               "build": None}
    assert fake.inserted == [] and woken == []


def test_a_failed_request_is_asked_again(client, pool, woken):
    fake = pool(index=BRCA1, latest=("failed", "rest.uniprot.org did not answer"))
    response = _ask(client)
    assert response.status_code == 202
    assert fake.inserted == [("BRCA1", "P38398", "brca1", None)]
    assert woken == [True]


def test_an_unbuildable_protein_says_why(client, pool, woken):
    fake = pool(index=UNBUILDABLE)
    response = _ask(client, "TP53BP2")
    assert response.status_code == 200
    assert response.json() == {"slug": None, "state": "unavailable", "reason": UNBUILDABLE[3],
                               "build": None}
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
    assert response.json() == {"slug": "brca1", "state": "buildable", "reason": None, "build": None}
    assert fake.inserted == [] and woken == []


@pytest.mark.parametrize("latest, said", [
    (("queued", None), {"slug": "brca1", "state": "pending", "reason": None, "build": None}),
    (("running", None), {"slug": "brca1", "state": "pending", "reason": None, "build": None}),
    (("failed", "UniProt did not answer"),
     {"slug": "brca1", "state": "failed", "reason": "UniProt did not answer", "build": None}),
    (("refused", "Too long."),
     {"slug": None, "state": "refused", "reason": "Too long.", "build": None}),
])
def test_a_read_says_where_a_request_is(client, pool, latest, said):
    pool(index=BRCA1, latest=latest)
    assert client.get("/proteins/resolve/BRCA1").json() == said


def test_a_read_of_a_resolved_protein_is_ready(client, pool):
    pool(protein="brca1", index=BRCA1, latest=("done", None))
    assert client.get("/proteins/resolve/BRCA1").json() == \
        {"slug": "brca1", "state": "ready", "reason": None, "build": None}


# ------------------------------------------------------------------ build


def _report(step, elapsed, **said):
    return {"step": step, "ahead": None, "scored": None, "residues": None, "left": None,
            "elapsed": elapsed, "transcript": None, "reason": None, "stopped": False,
            "ran": None, **said}


def test_a_new_ask_says_how_many_proteins_are_ahead_of_it(client, pool, woken):
    fake = pool(index=BRCA1, request_build=("queued", 0.2, 2))
    response = _ask(client)
    assert response.status_code == 202
    assert response.json()["build"] == _report("queued", 0.2, ahead=2)
    # Read after the row is written, so the ask counts itself in.
    assert fake.inserted and fake.built == [("request", ("BRCA1",))]


def test_a_request_a_worker_holds_is_making_its_record(client, pool):
    pool(index=BRCA1, latest=("running", None), request_build=("running", 3.5, 0))
    assert client.get("/proteins/resolve/BRCA1").json()["build"] == _report("record", 3.5)


def test_a_refused_request_ends_at_its_record_and_says_why_once(client, pool):
    pool(index=BRCA1, latest=("refused", "Too long."), request_build=("refused", 4.0, 0))
    said = client.get("/proteins/resolve/BRCA1").json()
    assert said["reason"] == "Too long."
    assert said["build"] == _report("record", 4.0)


def _bake(job="running", done=None, total=None, left=None, ahead=None, elapsed=60.0,
          track="pending", reason=None, ran=None):
    return (job, done, total, left, ahead, elapsed, "NM_007294.4", track, reason, ran)


@pytest.mark.parametrize("bake, report", [
    (_bake(job="queued", ahead=1),
     _report("scoring", 60.0, ahead=1, transcript="NM_007294.4")),
    # Claimed, and the model still loading: no line from the scorer yet.
    (_bake(),
     _report("scoring", 60.0, transcript="NM_007294.4")),
    (_bake(done=450, total=1863, left=612.5, elapsed=331.0, ran=290.0),
     _report("scoring", 331.0, scored=450, residues=1863, left=612.5,
             transcript="NM_007294.4", ran=290.0)),
    (_bake(done=1863, total=1863, left=0.0, elapsed=900.0, ran=860.0),
     _report("check", 900.0, scored=1863, residues=1863, left=0.0,
             transcript="NM_007294.4", ran=860.0)),
    (_bake(job="done", done=1863, total=1863, elapsed=905.0, track="ready"),
     _report("done", 905.0, scored=1863, residues=1863, transcript="NM_007294.4")),
])
def test_a_resolved_protein_says_how_far_esm2_has_got(client, pool, bake, report):
    fake = pool(protein="brca1", index=BRCA1, latest=("done", None), bake_build=bake)
    said = client.get("/proteins/resolve/BRCA1").json()
    assert (said["state"], said["slug"]) == ("ready", "brca1")
    assert said["build"] == report
    assert fake.built == [("bake", ("brca1",))]


def test_a_track_the_gate_refused_ends_at_the_check_with_its_reason(client, pool):
    reason = "ESM-2 prefers the residue that is there ... The scores are not drawn."
    pool(protein="sln", index=BRCA1, bake_build=_bake(
        job="failed", done=31, total=31, elapsed=20.0, track="refused", reason=reason))
    assert client.get("/proteins/resolve/SLN").json()["build"] == _report(
        "check", 20.0, scored=31, residues=31, transcript="NM_007294.4", reason=reason)


def test_a_scorer_that_never_got_going_ends_at_the_scoring(client, pool):
    reason = "Scoring failed 3 times: The scorer exited with status 1."
    pool(protein="brca1", index=BRCA1, bake_build=_bake(
        job="failed", elapsed=30.0, track="refused", reason=reason))
    assert client.get("/proteins/resolve/BRCA1").json()["build"] == _report(
        "scoring", 30.0, transcript="NM_007294.4", reason=reason)


def test_one_of_the_twenty_has_no_build(client, pool):
    # `_BAKE_BUILD` reads only proteins with no reading order, so it finds no row.
    fake = pool(protein="insulin", index=("P01308", "INS", True, None))
    assert client.get("/proteins/resolve/INS").json()["build"] is None
    assert fake.built == [("bake", ("insulin",))]


def test_a_protein_nobody_asked_for_has_no_build_and_reads_nothing_more(client, pool):
    fake = pool(index=BRCA1)
    assert client.get("/proteins/resolve/BRCA1").json()["build"] is None
    assert fake.built == []


# ------------------------------------------------------------------- stop


def _stop(client, gene="BRCA1", asker="install-a"):
    return client.post(f"/proteins/resolve/{gene}/stop", json={"asker": asker})


def test_an_ask_carries_who_asked(client, pool, woken):
    fake = pool(index=BRCA1)
    _ask(client, asker="install-a")
    assert fake.inserted == [("BRCA1", "P38398", "brca1", "install-a")]


def test_a_queued_request_is_withdrawn_by_the_phone_that_asked(client, pool):
    fake = pool(index=BRCA1, latest=("queued", None), ask=(41, "queued", "install-a"))
    response = _stop(client)
    assert response.status_code == 200
    assert response.json() == {
        "slug": None, "state": "stopped", "reason": None,
        "build": _report("queued", 12.0, stopped=True, reason="Stopped before it was built."),
    }
    assert fake.written == [("request stopped", (41,))]


@pytest.mark.parametrize("asker", ["install-b", None])
def test_only_the_phone_that_asked_may_stop_it(client, pool, asker):
    # Another install, or a request from before installs were named.
    fake = pool(index=BRCA1, latest=("queued", None), ask=(41, "queued", asker))
    response = _stop(client)
    assert response.status_code == 403
    assert response.json()["detail"] == "Only the phone that asked for it can stop it."
    assert fake.written == []


def test_a_request_a_worker_holds_is_left_to_it(client, pool):
    fake = pool(index=BRCA1, latest=("running", None), ask=(41, "running", "install-a"))
    response = _stop(client)
    assert response.status_code == 409
    assert response.json()["detail"] == "Its gene record is being written. Stop it in a moment."
    assert fake.written == []


def test_scoring_is_stopped_and_the_protein_stays_built_without_it(client, pool):
    fake = pool(protein="brca1", index=BRCA1, ask=(41, "done", "install-a"),
                bake=(7, "running", "install-a"),
                bake_build=_bake(job="stopped", done=450, total=1863, elapsed=400.0,
                                 track="absent"))
    response = _stop(client)
    assert response.status_code == 200
    said = response.json()
    assert (said["state"], said["slug"]) == ("stopped", "brca1")
    assert said["build"] == _report(
        "scoring", 400.0, scored=450, residues=1863, transcript="NM_007294.4",
        stopped=True, reason="Stopped at 450 of 1,863 residues.")
    assert fake.written == [("bake stopped", (7,)), ("track absent", ("brca1",))]


def test_a_stop_that_comes_as_the_track_lands_says_so(client, pool):
    fake = pool(protein="brca1", index=BRCA1, bake=(7, "running", "install-a"),
                stoppable=False)
    response = _stop(client)
    assert response.status_code == 409
    assert response.json()["detail"] == "Its ESM-2 track has just been made."
    assert fake.written == [("bake stopped", (7,))]


@pytest.mark.parametrize("protein, bake", [
    (None, None),                                  # nobody asked
    ("brca1", (7, "done", "install-a")),           # scored already
    ("insulin", None),                             # one of the twenty
])
def test_nothing_under_way_is_nothing_to_stop(client, pool, protein, bake):
    pool(protein=protein, index=BRCA1, bake=bake)
    response = _stop(client)
    assert response.status_code == 409
    assert response.json()["detail"] == "Nothing is building it."


def test_a_stop_for_a_gene_the_index_does_not_hold_is_not_found(client, pool):
    pool()
    assert _stop(client, "NOTAGENE").status_code == 404


def test_a_stop_says_who_it_is_from(client, pool):
    pool(index=BRCA1)
    assert client.post("/proteins/resolve/BRCA1/stop", json={}).status_code == 422


def test_asked_again_a_stopped_protein_is_scored_afresh(client, pool, woken):
    fake = pool(protein="brca1", index=BRCA1, resumable=True,
                bake_build=_bake(job="queued", ahead=0, elapsed=0.1))
    response = _ask(client, asker="install-a")
    assert response.status_code == 200
    said = response.json()
    assert (said["state"], said["build"]["step"]) == ("ready", "scoring")
    assert fake.written == [("resumed", {"slug": "brca1", "asker": "install-a"}),
                            ("track pending", ("brca1",))]
    assert fake.inserted == [] and woken == [True]


def test_asked_again_a_protein_with_nothing_stopped_writes_nothing(client, pool, woken):
    fake = pool(protein="brca1", index=BRCA1, resumable=False)
    assert _ask(client).status_code == 200
    assert fake.written == [("resumed", {"slug": "brca1", "asker": None})]
    assert woken == []


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
