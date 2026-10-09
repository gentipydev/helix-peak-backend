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

from pipeline import upload_tracks  # noqa: E402
from pipeline.resolver import local_worker, store, worker  # noqa: E402
from pipeline.resolver.test_local_worker import (  # noqa: E402
    MAKES_A_MODEL, SCORES_A_TRACK, fake_model_python, fake_python, make_settings,
)
from pipeline.resolver.test_resolve import INS, _body, _mutated  # noqa: E402
from pipeline.resolver.test_worker import _another_genes_record, a_model  # noqa: E402

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


def _constraint(target, record: bytes, report=None) -> bytes:
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

    def misaligned(target, record, report):
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

    def broken(target, record, report):
        raise RuntimeError("CUDA out of memory")

    states = [worker.score_next(conn, storage, score=broken)["state"] for _ in range(3)]
    assert states == ["queued", "queued", "refused"]
    assert _one(conn, "select state, reason from protein_track where kind = 'constraint'") == \
        ("refused", "Scoring failed 3 times: CUDA out of memory")


def _service(database_url, monkeypatch):
    from psycopg_pool import ConnectionPool

    from app import db

    pool = ConnectionPool(database_url, kwargs={"autocommit": True}, open=True)
    monkeypatch.setattr(db, "pool", pool)
    return pool


def test_a_reader_sees_each_step_of_a_build_as_the_worker_takes_it(
        conn, offline, database_url, monkeypatch):
    from app import resolves

    storage = Storage()
    with _service(database_url, monkeypatch):
        said, queued = resolves.request("INS")
        assert queued and said["state"] == "pending"
        build = said["build"]
        assert (build["step"], build["ahead"]) == ("queued", 0)
        assert 0 <= build["elapsed"] < 60

        worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))
        build = resolves.current("INS")["build"]
        # Its bake is queued, and nobody else's is ahead of it.
        assert (build["step"], build["ahead"], build["transcript"]) == \
            ("scoring", 0, "NM_000207.3")

        seen = []

        def scoring(target, record, report):
            seen.append(resolves.current("INS")["build"])  # claimed, nothing scored yet
            report(50, 110, 30.0)
            seen.append(resolves.current("INS")["build"])
            report(110, 110, 0.0)
            seen.append(resolves.current("INS")["build"])
            return _constraint(target, record)

        assert worker.score_next(conn, storage, score=scoring)["state"] == "ready"
        loading, halfway, checking = seen
        assert (loading["step"], loading["scored"], loading["left"]) == ("scoring", None, None)
        assert (halfway["step"], halfway["scored"], halfway["residues"]) == ("scoring", 50, 110)
        # The scorer's estimate, aged by the moment since it was written.
        assert 29.0 < halfway["left"] <= 30.0
        assert (checking["step"], checking["scored"]) == ("check", 110)

        said = resolves.current("INS")
        assert said["state"] == "ready"
        assert (said["build"]["step"], said["build"]["left"]) == ("done", None)
        assert said["build"]["elapsed"] >= halfway["elapsed"]


def test_progress_from_a_run_that_was_put_back_does_not_count(
        conn, offline, database_url, monkeypatch):
    from app import resolves

    storage = Storage()
    _ask(conn)
    worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))

    def breaks_halfway(target, record, report):
        report(50, 110, 30.0)
        raise RuntimeError("MPS backend out of memory")

    assert worker.score_next(conn, storage, score=breaks_halfway)["state"] == "queued"
    with _service(database_url, monkeypatch):
        build = resolves.current("INS")["build"]
        assert (build["step"], build["scored"], build["left"], build["ahead"]) == \
            ("scoring", None, None, 0)

        def claimed(target, record, report):
            # Claimed again: the first run's 50 of 110 is not this run's.
            assert resolves.current("INS")["build"]["scored"] is None
            return _constraint(target, record)

        assert worker.score_next(conn, storage, score=claimed)["state"] == "ready"


def test_a_build_refused_at_the_gate_says_so_at_the_check(conn, offline, database_url, monkeypatch):
    from app import resolves

    storage = Storage()
    _ask(conn)
    worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))

    def misaligned(target, record, report):
        report(110, 110, 0.0)
        raise ValueError("Alignment gate FAILED: 40.0% before, 45.0% after. No JSON written.")

    assert worker.score_next(conn, storage, score=misaligned)["state"] == "refused"
    with _service(database_url, monkeypatch):
        build = resolves.current("INS")["build"]
    assert (build["step"], build["scored"], build["residues"]) == ("check", 110, 110)
    assert build["reason"].startswith("ESM-2 prefers the residue that is there")


def test_proteins_ahead_are_counted_in_the_order_the_worker_takes_them(
        conn, offline, database_url, monkeypatch):
    from app import resolves

    storage = Storage()
    _ask(conn)
    worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))
    # INS waits for ESM-2; two more are asked for after it.
    conn.execute("insert into protein_index (uniprot, gene, name, length, annotation_score, "
                 "existence, buildable) values ('P61769', 'B2M', 'Beta-2-microglobulin', 119, "
                 "5, 1, true), ('P06213', 'INSR', 'Insulin receptor', 1382, 5, 1, true)")
    with _service(database_url, monkeypatch):
        assert resolves.request("B2M")[0]["build"]["ahead"] == 1
        assert resolves.request("INSR")[0]["build"]["ahead"] == 2
        assert resolves.current("INS")["build"]["ahead"] == 0


def _builds(conn):
    from app import suggest

    return set(conn.execute(suggest._BUILDS, {"genes": ["ins", "b2m"]}).fetchall())


def test_suggestions_know_which_proteins_are_being_built(conn, offline):
    assert _builds(conn) == set()
    _ask(conn)
    assert _builds(conn) == {("ins", "building")}  # asked for, not resolved
    storage = Storage()
    worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))
    assert _builds(conn) == {("ins", "building")}  # resolved, its ESM-2 bake to come
    worker.score_next(conn, storage, score=_constraint)
    assert _builds(conn) == set()


def test_the_built_list_is_every_resolved_protein_newest_first(
        conn, offline, database_url, monkeypatch):
    from fastapi.testclient import TestClient

    from app.main import app

    # A curated protein, never listed, and an older build whose gene has a
    # second index row, listed once.
    conn.execute(
        """
        insert into protein_index (uniprot, gene, name, length, annotation_score, existence,
            buildable, unavailable_reason)
        values ('P68871', 'HBB', 'Hemoglobin subunit beta', 147, 5, 1, true, null),
               ('O00241', 'SIRPB1', 'Signal-regulatory protein beta-1', 398, 5, 1, true, null),
               ('Q5TFQ8', 'SIRPB1', 'Signal-regulatory protein beta-1 isoform 3', 400, 3, 1,
                false, 'MANE Select encodes another protein.')
        """)
    conn.execute(
        """
        insert into protein (slug, gene, uniprot, accession, display, summary, residues, exons,
            chains, bridges, regions, disulfides, provenance, resolver_version, catalog_order,
            resolved_at)
        values ('hbb', 'HBB', 'P68871', 'NG_000007', 'Hemoglobin (beta chain)', 'Curated.',
                147, 3, 1, 0, '[]', '[]', '{}', 0, 3, now()),
               ('sirpb1', 'SIRPB1', 'O00241', 'NG_000000', 'Signal-regulatory protein beta-1',
                'Built.', 398, 6, 1, 0, '[]', '[]', '{}', 1, null, now() - interval '1 hour')
        """)
    storage = Storage()
    _ask(conn)
    worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))
    with _service(database_url, monkeypatch):
        client = TestClient(app)
        # Insulin's row is written and its ESM-2 bake is to come: not yet.
        body = client.get("/proteins/built").json()
        assert [(s["slug"], s["uniprot"]) for s in body["proteins"]] == [("sirpb1", "O00241")]

        worker.score_next(conn, storage, score=_constraint)
        first = client.get("/proteins/built", params={"limit": 1}).json()
        assert [(s["slug"], s["status"], s["name"], s["length"]) for s in first["proteins"]] == \
            [("ins", "ready", "Insulin", 110)]
        rest = client.get("/proteins/built", params={"limit": 1, "before": first["next"]}).json()
        assert [s["slug"] for s in rest["proteins"]] == ["sirpb1"]
        assert rest["next"] is None


def test_a_build_is_stopped_by_the_phone_that_asked_and_scored_again_when_asked(
        conn, offline, database_url, monkeypatch):
    from app import resolves

    storage = Storage()
    with _service(database_url, monkeypatch):
        # Stopped while queued: nothing was built, and asking again starts over.
        resolves.request("INS", asker="install-a")
        with pytest.raises(resolves.NotYours):
            resolves.stop("INS", "install-b")
        said = resolves.stop("INS", "install-a")
        assert (said["state"], said["build"]["step"], said["build"]["stopped"]) == \
            ("stopped", "queued", True)
        assert resolves.current("INS")["state"] == "buildable"
        assert worker.resolve_next(conn, storage) is None
        with pytest.raises(resolves.NothingToStop):
            resolves.stop("INS", "install-a")

        resolves.request("INS", asker="install-a")
        worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))
        # The bake the resolution queued is the asker's too.
        assert _one(conn, "select asker from bake_job where state = 'queued'") == ("install-a",)

        # Stopped mid-scoring: the reader's stop lands, and the scorer hears it
        # at its next line.
        def scoring(target, record, report):
            report(25, 110, 30.0)
            said = resolves.stop("INS", "install-a")
            assert said["build"]["reason"] == "Stopped at 25 of 110 residues."
            report(50, 110, 20.0)
            raise AssertionError("the scorer went on after it was stopped")

        assert worker.score_next(conn, storage, score=scoring) == \
            {"id": 1, "slug": "ins", "state": "stopped"}
        assert _one(conn, "select state, error, progress_done from bake_job") == \
            ("stopped", "Stopped by the reader who asked.", 25)
        assert _one(conn, "select state from protein_track "
                          "where slug = 'ins' and kind = 'constraint'") == ("absent",)
        assert worker.score_next(conn, storage, score=_constraint) is None
        assert _builds(conn) == {("ins", "stopped")}
        said = resolves.current("INS")
        assert (said["state"], said["build"]["stopped"]) == ("ready", True)

        # Asked again, by anyone: its scoring is queued afresh, theirs to stop.
        said, queued = resolves.request("INS", asker="install-b")
        assert queued and said["build"]["step"] == "scoring"
        assert _builds(conn) == {("ins", "building")}
        assert resolves.request("INS", asker="install-c")[1] is False  # once
        assert worker.score_next(conn, storage, score=_constraint)["state"] == "ready"
        assert resolves.current("INS")["build"]["step"] == "done"
        assert _builds(conn) == set()


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
    # The build's seconds go on between the two; everything else is one answer.
    said = [{**answer.json(), "build": {**answer.json()["build"], "elapsed": None}}
            for answer in (first, second)]
    assert said[0] == said[1]
    assert (said[0]["slug"], said[0]["state"], said[0]["reason"]) == ("ins", "pending", None)
    assert (said[0]["build"]["step"], said[0]["build"]["ahead"]) == ("queued", 0)
    assert _one(conn, "select count(*), min(gene), min(uniprot) from resolve_request") == \
        (1, "INS", "P01308")
    assert api.woken == [True]


def test_ask_resolve_and_read_ready(conn, api, offline):
    assert api.get("/proteins/resolve/INS").json()["state"] == "buildable"
    api.post("/proteins/resolve", json={"gene": "INS"})
    assert api.get("/proteins/resolve/INS").json()["state"] == "pending"

    worker.resolve_next(conn, Storage(), fetch_entry=lambda accession: offline(_body()))
    said = api.get("/proteins/resolve/INS").json()
    assert (said["slug"], said["state"], said["reason"]) == ("ins", "ready", None)
    # Ready to open its row, and its ESM-2 track still to come.
    assert said["build"]["step"] == "scoring"
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


def test_a_protein_whose_scoring_broke_is_scored_once_it_is_queued_again_by_hand(conn, mac):
    # The two statements the README gives an operator ("Operating it").
    _ask(conn)
    for attempt in ("first", "second", "third"):
        mac(attempt, "sys.exit(1)")
    assert _one(conn, "select state from protein_track where kind = 'constraint'") == ("refused",)
    assert _one(conn, "select state from bake_job") == ("failed",)

    conn.execute("update protein_track set state = 'pending', reason = null, updated_at = now() "
                 "where slug = 'ins' and kind = 'constraint' and state = 'refused'")
    conn.execute("insert into bake_job (slug, kind) values ('ins', 'constraint')")

    done = mac("scores", SCORES_A_TRACK)
    assert [(o["slug"], o["state"]) for o in done["scored"]] == [("ins", "ready")]
    assert _one(conn, "select state, reason from protein_track where kind = 'constraint'") == \
        ("ready", None)
    assert conn.execute("select state from bake_job order by id").fetchall() == \
        [("failed",), ("done",)]


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


# ------------------------------------------------- a model made on demand


def _modelled(target, record: bytes):
    """What a baker hands over, as far as the worker's gate reads it."""
    return a_model()


def _resolved_with_a_model_queued(conn, offline, storage, asker=None):
    conn.execute("insert into resolve_request (gene, uniprot, slug, asker) "
                 "values ('INS', 'P01308', 'ins', %s)", (asker,))
    return worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()),
                               structures=True)


def _structure(conn):
    return (_one(conn, "select state, attempts, error from bake_job where kind = 'structure'"),
            _one(conn, "select state, reason from protein_track where kind = 'structure'"))


def test_by_default_no_model_is_queued(conn, offline):
    # Modal's worker, and every worker before this one: a track is never left
    # pending where nothing will bake it.
    _ask(conn)
    summary = worker.sweep(conn, Storage(), fetch_entry=lambda accession: offline(_body()))
    assert (summary["constraint_queued"], summary["structure_queued"]) == (1, 0)
    assert _one(conn, "select count(*) from protein_track where kind = 'structure'") == (0,)
    assert _one(conn, "select count(*) from bake_job where kind = 'structure'") == (0,)
    assert worker.structure_next(conn, Storage(), model=_modelled) is None


def test_a_worker_that_makes_models_queues_one_beside_the_scoring(conn, offline):
    outcome = _resolved_with_a_model_queued(conn, offline, Storage(), asker="install-a")
    assert outcome["state"] == "done"
    assert conn.execute("select kind, state, asker from bake_job order by kind").fetchall() == \
        [("constraint", "queued", "install-a"), ("structure", "queued", "install-a")]
    assert conn.execute("select kind, state, format from protein_track "
                        "where kind in ('constraint', 'structure') order by kind").fetchall() == \
        [("constraint", "pending", "json"), ("structure", "pending", "fsceneb")]
    assert _one(conn, "select structure from protein where slug = 'ins'") == (None,)


def test_a_queued_model_is_made_and_lands_ready(conn, offline):
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    job = _one(conn, "select id from bake_job where kind = 'structure'")[0]

    assert worker.structure_next(conn, storage, model=_modelled) == \
        {"id": job, "slug": "ins", "state": "ready"}
    built = a_model()
    bucket, path, sha, fmt, provenance, state = _one(
        conn, "select bucket, object_path, sha256, format, provenance, state "
              "from protein_track where slug = 'ins' and kind = 'structure'")
    assert (bucket, fmt, state) == ("models", "fsceneb", "ready")
    assert path == f"structure/ins.{sha[:12]}.fsceneb"
    assert storage.get(bucket, path) == built.scene
    assert storage.get(bucket, provenance["glb"]["path"]) == built.glb
    assert (provenance["entry"], provenance["licence"]) == ("AF-P01308-F1", "CC BY 4.0")
    # The fold page's words and chains, on the protein's own row.
    assert _one(conn, "select structure from protein where slug = 'ins'") == \
        ({"chrome": built.chrome, "chains": built.chains},)
    assert _structure(conn) == (("done", 1, None), ("ready", None))
    assert worker.structure_next(conn, storage, model=_modelled) is None
    # The scoring is its own bake, and has not moved.
    assert _one(conn, "select state from bake_job where kind = 'constraint'") == ("queued",)
    assert worker.score_next(conn, storage, score=_constraint)["state"] == "ready"


def test_the_service_serves_a_model_made_on_demand(conn, offline, database_url, monkeypatch):
    from fastapi.testclient import TestClient
    from psycopg_pool import ConnectionPool

    from app import db, resolves
    from app.config import settings
    from app.main import app

    monkeypatch.setattr(settings, "supabase_url", "https://project.supabase.co")
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    with ConnectionPool(database_url, kwargs={"autocommit": True}, open=True) as pool:
        monkeypatch.setattr(db, "pool", pool)
        client = TestClient(app)

        # Queued, and then being made: the track says a model is on its way,
        # and the build a reader watches is still its ESM-2 track's alone.
        assert client.get("/protein/ins").json()["tracks"]["structure"] == "pending"
        assert resolves.current("INS")["build"]["step"] == "scoring"

        worker.structure_next(conn, storage, model=_modelled)
        built = a_model()
        body = client.get("/protein/ins").json()
        assert body["chains"] == built.chains and body["structure"] == built.chrome
        assert body["tracks"]["structure"] == "ready"
        track = client.get("/protein/ins/tracks").json()["structure"]
        assert (track["state"], track["format"], track["bytes"]) == \
            ("ready", "fsceneb", len(built.scene))
        assert track["url"] == ("https://project.supabase.co/storage/v1/object/public/models/"
                                f"structure/ins.{track['sha256'][:12]}.fsceneb")
        assert track["provenance"]["glb"]["bytes"] == len(built.glb)
        assert resolves.current("INS")["build"]["step"] == "scoring"
        assert client.get("/catalog").json()["proteins"] == []


def test_a_protein_with_no_model_says_why_on_its_track(conn, offline):
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    why = "AlphaFold DB has no model of proteins over 2,700 residues; this one has 4,834."

    def none(target, record):
        raise worker.Unmodelled(why)

    outcome = worker.structure_next(conn, storage, model=none)
    assert (outcome["state"], outcome["reason"]) == ("refused", why)
    assert _structure(conn) == (("failed", 1, why), ("refused", why))
    assert _one(conn, "select structure from protein where slug = 'ins'") == (None,)
    assert [path for _, path in storage.objects if path.startswith("structure/")] == []
    # A refusal is final: nothing is left to claim.
    assert worker.structure_next(conn, storage, model=_modelled) is None


def test_a_bake_that_breaks_is_tried_three_times_and_a_reader_is_not_told_how(conn, offline):
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)

    def broken(target, record):
        raise RuntimeError("the scene importer failed (255):\nUnhandled exception")

    states = [worker.structure_next(conn, storage, model=broken) for _ in range(3)]
    assert [outcome["state"] for outcome in states] == ["queued", "queued", "refused"]
    said = ("AlphaFold's model could not be made: the bake broke each of the 3 times "
            "it was tried.")
    assert states[2]["reason"] == said
    # The track carries the reader's sentence, and the job what broke.
    assert _structure(conn) == (
        ("failed", 3, "the scene importer failed (255): Unhandled exception"),
        ("refused", said))


def test_a_model_the_gate_declines_is_not_stored_and_is_tried_again(conn, offline):
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    wrong = lambda target, record: a_model(chains=[{"node": "chainA", "tint": "mature1"}])  # noqa: E731
    outcome = worker.structure_next(conn, storage, model=wrong)
    assert outcome["state"] == "queued" and "names no tint" in outcome["reason"]
    assert [path for _, path in storage.objects if path.startswith("structure/")] == []
    assert _structure(conn)[1] == ("pending", None)


def test_a_model_is_never_written_onto_a_curated_row(conn):
    conn.execute(
        """
        insert into protein (slug, gene, uniprot, accession, display, summary, residues, exons,
            chains, bridges, regions, disulfides, structure, provenance, resolver_version,
            catalog_order)
        values ('insulin', 'INS', 'P01308', 'NG_007114', 'Insulin', 'Curated.', 110, 3, 3, 3,
                '[]', '[]', '{"chrome": {"pdb": "3I40"}, "chains": []}', '{}', 0, 0)
        """)
    conn.execute("insert into bake_job (slug, kind) values ('insulin', 'structure')")
    # Its job is left for the hand bakes, as its scoring would be.
    assert worker.structure_next(conn, Storage(), model=_modelled) is None
    assert _one(conn, "select state from bake_job") == ("queued",)

    built = a_model()
    row, _ = upload_tracks.structure_row("insulin", built.glb, built.scene, built.provenance)
    with pytest.raises(RuntimeError, match="curated protein"):
        store.finish_structure(conn, 1, row, {"chrome": built.chrome, "chains": built.chains})
    assert _one(conn, "select structure from protein where slug = 'insulin'") == \
        ({"chrome": {"pdb": "3I40"}, "chains": []},)
    assert _one(conn, "select count(*) from protein_track") == (0,)
    assert _one(conn, "select state from bake_job") == ("queued",)


def test_a_dead_modellers_claim_goes_back_long_before_a_scorers(conn, offline):
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    store.claim_bake(conn, "structure")
    store.claim_bake(conn, "constraint")
    # Half an hour: past a model's patience, well inside a scoring's.
    conn.execute("update bake_job set started_at = now() - interval '30 minutes'")

    store.reap(conn)
    assert conn.execute("select kind, state from bake_job order by kind").fetchall() == \
        [("constraint", "running"), ("structure", "queued")]

    # Dead on its last try, its track says so in the modeller's words.
    conn.execute("update bake_job set state = 'running', attempts = 3, "
                 "started_at = now() - interval '30 minutes' where kind = 'structure'")
    store.reap(conn)
    assert _structure(conn) == (
        ("failed", 3, "The worker stopped before it finished."),
        ("refused", "The model's bake stopped before it finished, every time it was tried."))
    assert _one(conn, "select state from protein_track where kind = 'constraint'") == ("pending",)


def test_a_stop_ends_the_scoring_and_leaves_the_model_to_be_made(
        conn, offline, database_url, monkeypatch):
    from app import resolves

    storage = Storage()
    with _service(database_url, monkeypatch):
        _resolved_with_a_model_queued(conn, offline, storage, asker="install-a")
        said = resolves.stop("INS", "install-a")
        assert (said["state"], said["build"]["stopped"]) == ("stopped", True)
        assert conn.execute("select kind, state from bake_job order by kind").fetchall() == \
            [("constraint", "stopped"), ("structure", "queued")]
        # A model takes seconds, as the record did: it is made all the same,
        # and the protein opens with its fold and without its scores.
        assert worker.structure_next(conn, storage, model=_modelled)["state"] == "ready"
        assert worker.score_next(conn, storage, score=_constraint) is None
        assert _builds(conn) == {("ins", "stopped")}


@pytest.fixture
def mac_with_models(database_url, offline, tmp_path, monkeypatch):
    """`mac`, for a worker that makes models: one cycle scored by one script
    and modelled by another."""
    storage = Storage()
    sweep = worker.sweep
    monkeypatch.setattr(store, "TrackStorage", lambda url, key: storage)
    monkeypatch.setattr(
        worker, "sweep", lambda conn, storage, **said: sweep(
            conn, storage, fetch_entry=lambda accession: offline(_body()), **said))
    monkeypatch.setattr(local_worker, "CAFFEINATE",
                        (sys.executable, "-c", "import time; time.sleep(120)"))
    monkeypatch.setattr(local_worker, "_CHILD_POLL", 0.02)

    def run(name, scores, models, **check):
        folder = tmp_path / name
        (folder / "scorer").mkdir(parents=True)
        (folder / "modeller").mkdir()
        settings = dataclasses.replace(
            make_settings(tmp_path, fake_python(folder / "scorer", scores)),
            database_url=database_url, structures=True,
            structure_python=fake_model_python(folder / "modeller", models, **check))
        stop = threading.Event()
        return local_worker.cycle(settings, local_worker.Scorer(settings, stop), stop)

    run.storage = storage
    return run


def test_the_macs_worker_makes_the_model_before_it_scores(conn, mac_with_models, tmp_path):
    # The scorer notes whether the model was there when it started.
    seen = tmp_path / "model-first"
    scores = SCORES_A_TRACK + f"""
    Path({str(seen)!r}).write_text(os.environ["MODEL_STORED"])
"""
    _ask(conn)
    storage = mac_with_models.storage
    os.environ["MODEL_STORED"] = "unset"
    try:
        put = storage.put

        def noting(bucket, path, payload, content_type):
            if path.endswith(".fsceneb"):
                os.environ["MODEL_STORED"] = "stored"
            put(bucket, path, payload, content_type)

        storage.put = noting
        done = mac_with_models("first", scores, MAKES_A_MODEL)
    finally:
        del os.environ["MODEL_STORED"]

    assert [(o["gene"], o["state"]) for o in done["resolved"]] == [("INS", "done")]
    assert [(o["slug"], o["state"]) for o in done["modelled"]] == [("ins", "ready")]
    assert [(o["slug"], o["state"]) for o in done["scored"]] == [("ins", "ready")]
    assert seen.read_text() == "stored"
    assert conn.execute("select kind, state, attempts from bake_job order by kind").fetchall() == \
        [("constraint", "done", 1), ("structure", "done", 1)]
    assert sorted(path.split("/")[0] for _, path in storage.objects) == \
        ["constraint", "record", "structure", "structure"]
    chrome = _one(conn, "select structure from protein where slug = 'ins'")[0]["chrome"]
    assert chrome["pdb"] == "AF-P01308-F1"
    assert mac_with_models("again", scores, MAKES_A_MODEL) == \
        {"resolved": [], "scored": [], "idle": True, "modelled": []}


def test_a_modeller_that_cannot_start_claims_no_model_and_the_scoring_goes_on(
        conn, mac_with_models):
    _ask(conn)
    done = mac_with_models("no-pymol", SCORES_A_TRACK, "sys.exit(1)", check=1,
                           check_says="PyMOL (/opt/homebrew/bin/pymol) could not be run\n")
    assert done["modelled"] == [] and [o["state"] for o in done["scored"]] == ["ready"]
    # Left as it was queued, with none of its three tries used.
    assert _structure(conn) == (("queued", 0, None), ("pending", None))


def test_a_protein_alphafold_has_no_model_of_is_refused_by_the_macs_worker(
        conn, mac_with_models):
    _ask(conn)
    done = mac_with_models("no-model", SCORES_A_TRACK, """
        (out / "refusal.txt").parent.mkdir(parents=True, exist_ok=True)
        (out / "refusal.txt").write_text("AlphaFold DB holds no model of UniProt P01308.")
        sys.exit(3)
    """)
    assert [(o["state"], o["reason"]) for o in done["modelled"]] == \
        [("refused", "AlphaFold DB holds no model of UniProt P01308.")]
    assert _structure(conn)[1] == ("refused", "AlphaFold DB holds no model of UniProt P01308.")
    assert [o["state"] for o in done["scored"]] == ["ready"]
