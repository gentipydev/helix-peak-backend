"""The worker's SQL against a real Postgres, from request to ready track.

Opt-in: skipped unless `RESOLVER_TEST_DATABASE_URL` names a server this may
create a scratch database on, such as a local one:

    RESOLVER_TEST_DATABASE_URL=postgresql://postgres@127.0.0.1:54329/postgres \\
        .venv/bin/pytest pipeline/resolver/test_worker_pg.py

The scratch database gets every migration in order, so this also proves 0009
applies on top of 0001-0008. Nothing reaches NCBI, UniProt or storage: the
record and the entry are the insulin fixtures `test_resolve.py` uses, storage
is a dict, and ESM-2 is a function that returns a track. What is real is every
statement `store.py` runs, and the service reading the rows they write.
"""

from __future__ import annotations

import dataclasses
import json
import os
import re
import sys
import threading
import time
import uuid
from pathlib import Path
from urllib.error import URLError

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

ADMIN_URL = os.environ.get("RESOLVER_TEST_DATABASE_URL")
pytestmark = pytest.mark.skipif(
    not ADMIN_URL, reason="set RESOLVER_TEST_DATABASE_URL to run the worker against Postgres")

os.environ.setdefault("NCBI_EMAIL", "tests@example.com")

from pipeline.resolver import local_worker, store, worker  # noqa: E402
from pipeline.resolver.test_local_worker import (  # noqa: E402
    SCORES_A_TRACK, fake_python, make_settings,
)
from pipeline.resolver.test_resolve import INS, _body, _mutated  # noqa: E402
from pipeline.resolver.test_worker import _another_genes_record  # noqa: E402

MIGRATIONS = sorted((BACKEND / "migrations").glob("0*.sql"))


@pytest.fixture(scope="module")
def database_url():
    import psycopg
    from psycopg.conninfo import conninfo_to_dict, make_conninfo

    name = f"helixpeak_resolver_{uuid.uuid4().hex[:8]}"
    with psycopg.connect(ADMIN_URL, autocommit=True) as admin:
        admin.execute(f"create database {name}")
    url = make_conninfo(**{**conninfo_to_dict(ADMIN_URL), "dbname": name})
    with psycopg.connect(url, autocommit=True) as conn:
        for migration in MIGRATIONS:
            conn.execute(migration.read_text())
    yield url
    with psycopg.connect(ADMIN_URL, autocommit=True) as admin:
        admin.execute(f"drop database {name} with (force)")


@pytest.fixture
def conn(database_url):
    with store.connect(database_url) as connection:
        connection.execute(
            "truncate resolve_request, bake_job, protein_track, protein_alias, protein, "
            "protein_index_term, protein_index, protein_index_release restart identity cascade")
        connection.execute(
            """
            insert into protein_index (uniprot, gene, name, length, annotation_score, existence,
                refseq_nuc, refseq_prot, chrom_acc, chrom_start, chrom_end, chrom_strand,
                refseqgene, buildable)
            values ('P01308', 'INS', 'Insulin', 110, 5, 1, 'NM_000207.3', 'NP_000198.1',
                    'NC_000011.10', 2159779, 2161209, -1, 'NG_007114.1', true)
            """)
        connection.execute(
            "insert into protein_index_release (uniprot, mane, entries, buildable) "
            "values ('2026_03', 'v1.5', 1, 1)")
        yield connection


class Storage:
    def __init__(self):
        self.objects: dict[tuple[str, str], bytes] = {}

    def put(self, bucket, path, payload, content_type):
        self.objects[(bucket, path)] = payload

    def get(self, bucket, path):
        return self.objects[(bucket, path)]


def _ask(conn, gene="INS", uniprot="P01308"):
    conn.execute("insert into resolve_request (gene, uniprot, slug) values (%s, %s, %s)",
                 (gene, uniprot, gene.lower()))


def _one(conn, sql, *params):
    return conn.execute(sql, params).fetchone()


def _constraint(target, record: bytes) -> bytes:
    """What the scorer writes, as far as the upload gate reads it."""
    return (json.dumps({
        "gene": target.gene, "uniprot": target.uniprot,
        "model": "facebook/esm2_t33_650M_UR50D", "revision": "test",
        "method": "masked_marginals", "normalization": "minmax",
        "score_units": "natural_log_ratio_to_wildtype", "entropy_vocabulary": "full",
        "vocabulary_size": 33, "context": {"mode": "full", "residues": target.aa},
    }) + "\n").encode()


def test_a_request_becomes_a_protein_with_its_record_ready(conn, offline):
    storage = Storage()
    _ask(conn)
    outcome = worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))

    assert outcome == {"id": 1, "gene": "INS", "state": "done", "slug": "ins"}
    assert _one(conn, "select state, slug, attempts, resolver_version, reason "
                      "from resolve_request") == ("done", "ins", 1, 1, None)
    assert _one(conn, "select catalog_order, resolver_version, display, residues "
                      "from protein where slug = 'ins'") == (None, 1, "Insulin", 110)
    assert _one(conn, "select count(*) from protein_alias where slug = 'ins'") == (5,)

    bucket, path, sha, state = _one(
        conn, "select bucket, object_path, sha256, state from protein_track "
              "where slug = 'ins' and kind = 'record'")
    assert (bucket, state) == ("tracks", "ready")
    assert path == f"record/ins.{sha[:12]}.json"
    assert json.loads(storage.get(bucket, path))["gene"] == "INS"

    assert _one(conn, "select state from protein_track "
                      "where slug = 'ins' and kind = 'constraint'") == ("pending",)
    assert _one(conn, "select kind, state from bake_job where slug = 'ins'") == \
        ("constraint", "queued")


def test_the_service_serves_what_the_resolver_wrote(conn, offline, database_url, monkeypatch):
    from fastapi.testclient import TestClient
    from psycopg_pool import ConnectionPool

    from app import db
    from app.config import settings
    from app.main import app

    # A ready track is only served with somewhere to fetch it from.
    monkeypatch.setattr(settings, "supabase_url", "https://project.supabase.co")
    _ask(conn)
    assert worker.resolve_next(
        conn, Storage(), fetch_entry=lambda accession: offline(_body()))["state"] == "done"
    with ConnectionPool(database_url, kwargs={"autocommit": True}, open=True) as pool:
        monkeypatch.setattr(db, "pool", pool)
        client = TestClient(app)

        detail = client.get("/protein/ins")
        assert detail.status_code == 200
        body = detail.json()
        assert (body["slug"], body["gene"], body["catalog_order"]) == ("ins", "INS", None)
        assert body["provenance"]["prose"] == "templated"
        assert body["tracks"]["record"] == "ready" and body["tracks"]["constraint"] == "pending"

        tracks = client.get("/protein/ins/tracks").json()
        assert tracks["record"]["state"] == "ready"
        assert tracks["record"]["url"].startswith(
            "https://project.supabase.co/storage/v1/object/public/tracks/record/ins.")
        assert tracks["constraint"]["state"] == "pending"

        # The list is the curated rows only; a resolved protein never joins it.
        assert client.get("/catalog").json()["proteins"] == []


def test_a_queued_constraint_is_scored_and_lands_ready(conn, offline):
    storage = Storage()
    _ask(conn)
    worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))

    assert worker.score_next(conn, storage, score=_constraint) == \
        {"id": 1, "slug": "ins", "state": "ready"}
    state, provenance = _one(conn, "select state, provenance from protein_track "
                                   "where slug = 'ins' and kind = 'constraint'")
    assert state == "ready" and provenance["model"] == "facebook/esm2_t33_650M_UR50D"
    assert _one(conn, "select state, attempts from bake_job") == ("done", 1)
    assert worker.score_next(conn, storage, score=_constraint) is None


def test_a_scorer_refusal_is_said_on_the_track(conn, offline):
    storage = Storage()
    _ask(conn)
    worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))

    def misaligned(target, record):
        raise ValueError("Alignment gate FAILED: 40.0% before, 45.0% after. No JSON written.")

    outcome = worker.score_next(conn, storage, score=misaligned)
    assert outcome["state"] == "refused"
    state, reason = _one(conn, "select state, reason from protein_track "
                               "where slug = 'ins' and kind = 'constraint'")
    assert state == "refused"
    assert reason.startswith("ESM-2 prefers the residue that is there to the one before it 40.0%")
    assert "to the one after it 45.0%" in reason and reason.endswith("The scores are not drawn.")
    assert _one(conn, "select state from bake_job") == ("failed",)


def test_a_scorer_that_breaks_is_tried_three_times_then_said(conn, offline):
    storage = Storage()
    _ask(conn)
    worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))

    def broken(target, record):
        raise RuntimeError("CUDA out of memory")

    states = [worker.score_next(conn, storage, score=broken)["state"] for _ in range(3)]
    assert states == ["queued", "queued", "refused"]
    assert _one(conn, "select state, reason from protein_track where kind = 'constraint'") == \
        ("refused", "Scoring failed 3 times: CUDA out of memory")


def test_a_refused_protein_writes_nothing_but_the_reason(conn, offline):
    _ask(conn)
    outcome = worker.resolve_next(
        conn, Storage(), fetch_entry=lambda accession: offline(_mutated(_body(), [30, 60, 70, 80])))
    assert outcome["state"] == "refused"
    state, reason, version = _one(conn, "select state, reason, resolver_version from resolve_request")
    assert state == "refused" and "at 4 residues" in reason and version == 1
    assert _one(conn, "select count(*) from protein") == (0,)


def test_a_record_the_upload_gate_declines_is_refused_once(conn, offline, monkeypatch):
    monkeypatch.setattr(worker, "resolve", _another_genes_record)
    storage = Storage()
    _ask(conn)
    outcome = worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))

    assert outcome["state"] == "refused"
    state, reason, attempts, version = _one(
        conn, "select state, reason, attempts, resolver_version from resolve_request")
    assert (state, attempts, version) == ("refused", 1, 1) and "expected 'INS'" in reason
    assert _one(conn, "select count(*) from protein") == (0,)
    assert storage.objects == {}
    # Nothing is left on the queue to try again.
    assert worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body())) is None


def test_an_unreachable_uniprot_is_retried_then_failed(conn, offline):
    _ask(conn)

    def unreachable(accession):
        raise URLError("rest.uniprot.org did not answer")

    states = [worker.resolve_next(conn, Storage(), fetch_entry=unreachable)["state"]
              for _ in range(3)]
    assert states == ["queued", "queued", "failed"]
    state, reason, attempts = _one(conn, "select state, reason, attempts from resolve_request")
    assert (state, attempts) == ("failed", 3) and "did not answer" in reason


def test_a_gene_the_catalog_already_lists_opens_where_it_is(conn, offline):
    conn.execute(
        """
        insert into protein (slug, gene, uniprot, accession, display, summary, residues, exons,
            chains, bridges, regions, disulfides, provenance, resolver_version, catalog_order)
        values ('insulin', 'INS', 'P01308', 'NG_007114', 'Insulin', 'Curated.', 110, 3, 3, 3,
                '[]', '[]', '{}', 0, 0)
        """)
    _ask(conn)
    outcome = worker.resolve_next(conn, Storage(), fetch_entry=lambda accession: offline(_body()))
    assert outcome == {"id": 1, "gene": "INS", "state": "done", "slug": "insulin"}
    assert _one(conn, "select slug from resolve_request") == ("insulin",)
    assert _one(conn, "select count(*) from protein") == (1,)


def test_the_resolver_never_writes_over_a_curated_row(conn, offline):
    from pipeline.resolver.resolve import resolve

    conn.execute(
        """
        insert into protein (slug, gene, uniprot, accession, display, summary, residues, exons,
            chains, bridges, regions, disulfides, provenance, resolver_version, catalog_order)
        values ('ins', 'OTHER', 'P00000', 'NG_000000', 'Curated', 'Curated.', 1, 1, 1, 0,
                '[]', '[]', '{}', 0, 7)
        """)
    _ask(conn)
    request = store.claim_request(conn)
    resolution = resolve(INS, offline(_body()))
    record = worker.stage(Storage(), "record", resolution.target, resolution.record)
    with pytest.raises(RuntimeError, match="curated protein"):
        store.write_resolution(conn, request["id"], resolution, record, 1)
    assert _one(conn, "select display, catalog_order from protein where slug = 'ins'") == \
        ("Curated", 7)
    assert _one(conn, "select count(*) from protein_track") == (0,)


def test_a_bake_for_a_curated_protein_is_left_for_the_hand_bakes(conn):
    conn.execute(
        """
        insert into protein (slug, gene, uniprot, accession, display, summary, residues, exons,
            chains, bridges, regions, disulfides, provenance, resolver_version, catalog_order)
        values ('insulin', 'INS', 'P01308', 'NG_007114', 'Insulin', 'Curated.', 110, 3, 3, 3,
                '[]', '[]', '{}', 0, 0)
        """)
    conn.execute("insert into bake_job (slug, kind) values ('insulin', 'constraint')")
    assert store.claim_bake(conn, "constraint") is None
    assert _one(conn, "select state from bake_job") == ("queued",)


def test_two_asks_for_one_gene_are_one_request(conn):
    import psycopg

    _ask(conn)
    with pytest.raises(psycopg.errors.UniqueViolation):
        _ask(conn, gene="ins")


def test_a_dead_workers_claims_go_back_on_the_queue(conn, offline):
    storage = Storage()
    _ask(conn)
    worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))
    _ask(conn, gene="TP53", uniprot="P04637")
    store.claim_request(conn)
    store.claim_bake(conn, "constraint")
    conn.execute("update resolve_request set started_at = now() - interval '2 hours' "
                 "where state = 'running'")
    conn.execute("update bake_job set started_at = now() - interval '2 hours'")

    store.reap(conn)
    assert _one(conn, "select state, started_at from resolve_request where gene = 'TP53'") == \
        ("queued", None)
    assert _one(conn, "select state from bake_job") == ("queued",)


def test_a_sweep_resolves_and_says_what_waits_for_a_gpu(conn, offline):
    _ask(conn)
    summary = worker.sweep(conn, Storage(), fetch_entry=lambda accession: offline(_body()))
    assert [r["state"] for r in summary["resolved"]] == ["done"]
    assert summary["constraint_queued"] == 1


# ------------------------------------------------- the service's endpoints


@pytest.fixture
def api(database_url, monkeypatch):
    """The service on this database, with the resolver's wake recorded."""
    from fastapi.testclient import TestClient
    from psycopg_pool import ConnectionPool

    from app import db, resolves
    from app.main import app

    woken = []
    monkeypatch.setattr(resolves, "wake", lambda: woken.append(True))
    with ConnectionPool(database_url, kwargs={"autocommit": True}, open=True) as pool:
        monkeypatch.setattr(db, "pool", pool)
        client = TestClient(app)
        client.woken = woken
        yield client


def test_two_readers_asking_make_one_request_and_one_wake(conn, api):
    first = api.post("/proteins/resolve", json={"gene": "INS"})
    second = api.post("/proteins/resolve", json={"gene": "ins"})
    assert (first.status_code, second.status_code) == (202, 202)
    assert first.json() == second.json() == {"slug": "ins", "state": "pending", "reason": None}
    assert _one(conn, "select count(*), min(gene), min(uniprot) from resolve_request") == \
        (1, "INS", "P01308")
    assert api.woken == [True]


def test_ask_resolve_and_read_ready(conn, api, offline):
    assert api.get("/proteins/resolve/INS").json()["state"] == "buildable"
    api.post("/proteins/resolve", json={"gene": "INS"})
    assert api.get("/proteins/resolve/INS").json()["state"] == "pending"

    worker.resolve_next(conn, Storage(), fetch_entry=lambda accession: offline(_body()))
    assert api.get("/proteins/resolve/INS").json() == \
        {"slug": "ins", "state": "ready", "reason": None}
    assert api.post("/proteins/resolve", json={"gene": "INS"}).status_code == 200


def test_a_refusal_stands_and_a_failure_is_asked_again(conn, api, offline):
    api.post("/proteins/resolve", json={"gene": "INS"})
    worker.resolve_next(conn, Storage(),
                        fetch_entry=lambda accession: offline(_mutated(_body(), [30, 60, 70, 80])))
    refused = api.post("/proteins/resolve", json={"gene": "INS"})
    assert refused.status_code == 200 and refused.json()["state"] == "refused"
    assert _one(conn, "select count(*) from resolve_request") == (1,)

    conn.execute("update resolve_request set state = 'failed', reason = 'NCBI did not answer'")
    assert api.get("/proteins/resolve/INS").json()["state"] == "failed"
    assert api.post("/proteins/resolve", json={"gene": "INS"}).status_code == 202
    assert _one(conn, "select count(*) from resolve_request where state = 'queued'") == (1,)


def test_the_daily_cap_counts_the_days_requests(conn, api, monkeypatch):
    from app.config import settings

    monkeypatch.setattr(settings, "resolves_per_day", 1)
    conn.execute("insert into resolve_request (gene, uniprot, slug, state, reason) "
                 "values ('TP53', 'P04637', 'tp53', 'refused', 'Not today.')")
    assert api.post("/proteins/resolve", json={"gene": "INS"}).status_code == 429
    conn.execute("update resolve_request set requested_at = now() - interval '2 days'")
    assert api.post("/proteins/resolve", json={"gene": "INS"}).status_code == 202


# ------------------------------------------------- the worker on a Mac


@pytest.fixture
def mac(database_url, offline, tmp_path, monkeypatch):
    """`local_worker.cycle` on this database: storage a dict, UniProt and NCBI the
    insulin fixtures, `caffeinate` a process that only waits. Call it with a
    name and a `fake_python` body, and it runs one cycle scored by that script."""
    storage = Storage()
    sweep = worker.sweep
    monkeypatch.setattr(store, "TrackStorage", lambda url, key: storage)
    monkeypatch.setattr(
        worker, "sweep", lambda conn, storage, **said: sweep(
            conn, storage, fetch_entry=lambda accession: offline(_body()), **said))
    monkeypatch.setattr(local_worker, "CAFFEINATE",
                        (sys.executable, "-c", "import time; time.sleep(120)"))
    monkeypatch.setattr(local_worker, "_CHILD_POLL", 0.02)

    def run(name, body, stop=None, **check):
        folder = tmp_path / name
        folder.mkdir()
        settings = dataclasses.replace(
            make_settings(tmp_path, fake_python(folder, body, **check)),
            database_url=database_url)
        stop = stop or threading.Event()
        return local_worker.cycle(settings, local_worker.Scorer(settings, stop), stop)

    run.storage = storage
    return run


def _bake_and_track(conn):
    return (_one(conn, "select state, attempts, error, started_at from bake_job"),
            _one(conn, "select state from protein_track where kind = 'constraint'"))


def test_the_macs_worker_takes_a_request_to_a_scored_protein(conn, mac):
    _ask(conn)
    assert store.queued_requests(conn) == 1

    done = mac("scores", SCORES_A_TRACK)
    assert [(o["gene"], o["state"]) for o in done["resolved"]] == [("INS", "done")]
    assert [(o["slug"], o["state"]) for o in done["scored"]] == [("ins", "ready")]
    assert done["idle"] is False

    assert _one(conn, "select state, slug, attempts from resolve_request") == ("done", "ins", 1)
    assert _one(conn, "select state, attempts, error from bake_job") == ("done", 1, None)
    state, provenance = _one(conn, "select state, provenance from protein_track "
                                   "where kind = 'constraint'")
    assert state == "ready" and provenance["model"] == "facebook/esm2_t33_650M_UR50D"
    assert sorted(path.split("/")[0] for _, path in mac.storage.objects) == \
        ["constraint", "record"]
    assert store.queued_requests(conn) == 0

    # Nothing is left: the next cycle only reaps.
    assert mac("again", SCORES_A_TRACK) == {"resolved": [], "scored": [], "idle": True}


def test_a_stop_while_scoring_puts_the_bake_back_once(conn, mac, tmp_path, monkeypatch):
    # The scorer says when it has started, and the stop comes then: mid-score,
    # however long resolving took.
    started = tmp_path / "scoring-has-started"
    monkeypatch.setenv("SCORING_STARTED", str(started))
    stop = threading.Event()

    def stop_once_scoring():
        while not started.exists():
            time.sleep(0.01)
        stop.set()

    threading.Thread(target=stop_once_scoring, daemon=True).start()
    _ask(conn)
    done = mac("slow", """
        import time
        Path(os.environ["SCORING_STARTED"]).write_text("")
        time.sleep(120)
    """, stop=stop)

    assert [o["state"] for o in done["resolved"]] == ["done"]
    assert [(o["state"], o["reason"]) for o in done["scored"]] == \
        [("queued", local_worker.STOPPED)]
    # Back on the queue at once, not left running for the reaper's 90 minutes.
    # Claimed once: a worker that is leaving does not take it again, which
    # would be its second and third tries, and then its refusal.
    assert _bake_and_track(conn) == (("queued", 1, local_worker.STOPPED, None), ("pending",))

    done = mac("scores", SCORES_A_TRACK)
    assert [o["state"] for o in done["scored"]] == ["ready"]
    assert _one(conn, "select state, attempts, error from bake_job") == ("done", 2, None)
    assert _one(conn, "select state from protein_track where kind = 'constraint'") == ("ready",)


def test_a_scorer_that_cannot_start_claims_no_bake(conn, mac):
    _ask(conn)
    done = mac("no-torch", "sys.exit(1)", check=1, check_says="No module named 'torch'\n")
    assert [o["state"] for o in done["resolved"]] == ["done"] and done["scored"] == []
    assert _bake_and_track(conn) == (("queued", 0, None, None), ("pending",))


def test_a_scorer_that_breaks_is_tried_once_a_cycle(conn, mac):
    _ask(conn)
    states = []
    for attempt in ("first", "second", "third"):
        done = mac(attempt, "sys.stderr.write('RuntimeError: MPS backend out of memory')\n"
                            "sys.exit(1)")
        states.append([o["state"] for o in done["scored"]])
        states.append(_one(conn, "select attempts from bake_job")[0])
    # One try a cycle, a poll apart, not three in the second it takes to fail.
    assert states == [["queued"], 1, ["queued"], 2, ["refused"], 3]
    assert _one(conn, "select state, reason from protein_track where kind = 'constraint'") == (
        "refused", "Scoring failed 3 times: The scorer exited with status 1: "
                   "RuntimeError: MPS backend out of memory")


def test_the_queue_is_read_as_lines(conn, mac):
    assert local_worker.queue(conn).splitlines() == [
        "requests  0 queued, 0 running (0 done, 0 refused, 0 failed)",
        "bakes     0 queued, 0 running (0 done, 0 failed)",
        "no request yet"]

    _ask(conn)
    assert local_worker.queue(conn).splitlines()[0] == \
        "requests  1 queued, 0 running (0 done, 0 refused, 0 failed)"
    mac("no-torch", "sys.exit(1)", check=1)
    lines = local_worker.queue(conn).splitlines()
    assert lines[:4] == [
        "requests  0 queued, 0 running (1 done, 0 refused, 0 failed)",
        "bakes     1 queued, 0 running (0 done, 0 failed)",
        "  ins constraint: queued, claimed 0 time(s)",
        "last requests"]
    assert re.fullmatch(r"  \d{4}-\d\d-\d\d \d\d:\d\d  INS +done", lines[4]) and len(lines) == 5
    # It only read: the queue is as it was.
    assert _bake_and_track(conn) == (("queued", 0, None, None), ("pending",))


def test_the_queue_command_prints_the_queue_and_minds_no_reader_leaving(
        conn, database_url, tmp_path):
    import subprocess

    _ask(conn)
    home = tmp_path / "home"
    home.mkdir()
    # Every key named, so the repository's `.env` is never what is read.
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith(("RESOLVER_", "HELIXPEEK_"))}
    environment.update(
        DATABASE_URL=database_url, SUPABASE_URL="https://project.invalid",
        SUPABASE_SERVICE_KEY="sb_secret_not_a_real_key", NCBI_EMAIL="tests@example.com",
        HOME=str(home))
    command = [sys.executable, "-m", "pipeline.resolver.local_worker", "--queue"]

    done = subprocess.run(command, cwd=BACKEND, env=environment, capture_output=True,
                          text=True, timeout=120)
    assert (done.returncode, done.stderr) == (0, "")
    lines = done.stdout.splitlines()
    assert lines[:3] == ["requests  1 queued, 0 running (0 done, 0 refused, 0 failed)",
                         "bakes     0 queued, 0 running (0 done, 0 failed)", "last requests"]
    assert re.fullmatch(r"  \d{4}-\d\d-\d\d \d\d:\d\d  INS +queued", lines[3])
    # It read, and that is all: no log, no directory, and the request untouched.
    assert not (home / "Library" / "Logs").exists()
    assert not (home / "Library" / "Application Support").exists()
    assert _one(conn, "select state, attempts from resolve_request") == ("queued", 0)

    # `status | head`: whoever was reading has gone before it prints.
    piped = subprocess.Popen(command, cwd=BACKEND, env=environment, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, text=True)
    piped.stdout.close()
    assert piped.wait(timeout=120) == 0
    assert piped.stderr.read() == ""
    piped.stderr.close()
