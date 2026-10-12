"""A fold's publication, its generations and its searches, against a real Postgres.

Opt-in, as `test_worker_pg.py` is, whose scratch database (every migration in
order, 0012 included), fixtures and stand-ins these share: storage is a dict,
the modeller a function, and what is real is every statement `store.py` runs
and the service reading the rows they write.

    RESOLVER_TEST_DATABASE_URL=postgresql://postgres@127.0.0.1:54329/postgres \\
        .venv/bin/pytest pipeline/resolver/test_fold_pg.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.resolver import local_worker, store, worker  # noqa: E402
from pipeline.resolver.resolve import target_of  # noqa: E402
from pipeline.resolver.test_worker import a_fold, a_model, a_selection  # noqa: E402
from pipeline.resolver.test_worker_pg import (  # noqa: E402,F401 -- fixtures, used by name
    Storage, _builds, _constraint, _one, _resolved_with_a_model_queued, _service, conn,
    database_url, pytestmark,
)

POLICY = worker.pdb_policy.POLICY_VERSION


def _modeller(made):
    """A modeller that hands over `made` (or raises it), and notes what it was
    asked, its target's slug and the prior it was given."""
    asked = []

    def model(target, record, ask):
        asked.append((target.slug, ask.prior))
        if isinstance(made, BaseException):
            raise made
        return made

    model.asked = asked
    return model


def _fold(conn) -> dict:
    """What a reader is served of INS's fold, and what the tables record."""
    tracks = {kind: (state, sha256, reason) for kind, state, sha256, reason in conn.execute(
        "select kind, state, sha256, reason from protein_track where slug = 'ins' "
        "and kind in ('structure', 'folding')").fetchall()}
    return {
        "structure": _one(conn, "select structure from protein where slug = 'ins'")[0],
        "tracks": tracks,
        "generations": conn.execute(
            "select id, source, state, active from fold_generation where slug = 'ins' "
            "order by id").fetchall(),
        "selections": conn.execute(
            "select outcome, candidate from fold_selection where slug = 'ins' "
            "order by id").fetchall(),
        "jobs": conn.execute(
            "select state, attempts, refresh from bake_job where slug = 'ins' "
            "and kind = 'structure' order by id").fetchall(),
    }


def _published_alphafold_before_0012(conn, offline, storage) -> worker.Modelled:
    """INS built as it was before 0012: its AlphaFold model written by the old
    publication, with no generation and no folding track."""
    _resolved_with_a_model_queued(conn, offline, storage)
    job = store.claim_bake(conn, "structure")
    target = target_of(store.protein(conn, "ins"))
    built = a_model()
    row, structure = worker.stage_structure(storage, target, built)
    store.finish_structure(conn, job["id"], row, structure)
    return built


# ------------------------------------------------- publication


def test_an_experimental_fold_is_published_with_its_animation_whole(
        conn, offline, database_url, monkeypatch):
    from fastapi.testclient import TestClient
    from psycopg_pool import ConnectionPool

    from app import db
    from app.config import settings
    from app.main import app

    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    folded = a_fold()
    outcome = worker.structure_next(conn, storage, model=_modeller(folded))
    assert outcome == {"id": 2, "slug": "ins", "state": "ready", "source": "pdb",
                       "generation": 1}

    said = _fold(conn)
    assert said["structure"] == {"chrome": folded.model.chrome, "chains": folded.model.chains}
    scene = hashlib.sha256(folded.model.scene).hexdigest()
    folding = hashlib.sha256(folded.folding).hexdigest()
    assert said["tracks"] == {"structure": ("ready", scene, None),
                              "folding": ("ready", folding, None)}
    assert said["generations"] == [(1, "pdb", "ready", True)]
    assert said["selections"] == [("chosen", "1ABC:A@25,B@90")]
    assert said["jobs"] == [("done", 1, False)]
    # What the generation keeps is what was written, and what made it.
    generation = _one(conn, "select structure_track, folding_track, audit, record_sha256, "
                            "regions, policy, selection_id from fold_generation")
    structure_track, folding_track, audit, record_sha256, regions, policy, selection = generation
    assert structure_track["sha256"] == scene and folding_track["sha256"] == folding
    assert audit["path"].startswith("fold_audit/ins.")
    assert storage.get(audit["bucket"], audit["path"]) == folded.audit
    assert record_sha256 == store.ready_track(conn, "ins", "record")["sha256"]
    assert (regions, policy, selection) == (store.regions_of(conn, "ins"), POLICY, 1)

    monkeypatch.setattr(settings, "supabase_url", "https://project.supabase.co")
    with ConnectionPool(database_url, kwargs={"autocommit": True}, open=True) as pool:
        monkeypatch.setattr(db, "pool", pool)
        client = TestClient(app)
        tracks = client.get("/protein/ins/tracks").json()
        assert tracks["folding"]["state"] == "ready"
        assert tracks["folding"]["url"] == (
            "https://project.supabase.co/storage/v1/object/public/tracks/"
            f"folding/ins.{folding[:12]}.json")
        # The pair the app holds together: the track names the model it was made in.
        assert json.loads(folded.folding)["frame"]["glb_sha256"] == \
            tracks["structure"]["provenance"]["glb"]["sha256"] == \
            tracks["folding"]["provenance"]["frame"]["glb_sha256"]
        assert tracks["structure"]["provenance"]["audit"]["sha256"] == \
            hashlib.sha256(folded.audit).hexdigest()
        body = client.get("/protein/ins").json()
        assert body["structure"]["pdb"] == "1ABC" and body["tracks"]["folding"] == "ready"


def test_a_bake_taken_from_its_worker_writes_nothing(conn, offline):
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)

    def taken(target, record, ask):
        # The reaper put the bake back, and another worker took it.
        conn.execute("update bake_job set attempts = attempts + 1 where kind = 'structure'")
        return a_fold()

    before = _fold(conn)
    assert worker.structure_next(conn, storage, model=taken)["state"] == "superseded"
    assert _fold(conn) == {**before, "jobs": [("running", 2, False)]}
    # The same for a bake put back on the queue and not yet taken.
    conn.execute("update bake_job set state = 'queued', attempts = 1 where kind = 'structure'")
    job = store.claim_bake(conn, "structure")
    conn.execute("update bake_job set state = 'queued', started_at = null where id = %s",
                 (job["id"],))
    record = store.ready_track(conn, "ins", "record")["sha256"]
    with pytest.raises(store.Superseded):
        store.publish_fold(conn, job, {"source": "none", "reason": "No model.",
                                       "record_sha256": record,
                                       "regions": store.regions_of(conn, "ins")})
    assert _fold(conn)["tracks"] == before["tracks"] and _fold(conn)["generations"] == []
    # Nor may it put itself back, or say it failed.
    assert store.end_fold_bake(conn, job, error="broke", requeue=False, reason="No.") is False
    assert _one(conn, "select state, error from bake_job where kind = 'structure'") == \
        ("queued", None)


def test_a_protein_changed_while_its_fold_was_made_writes_nothing(conn, offline):
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    record_sha = store.ready_track(conn, "ins", "record")["sha256"]

    def changed(sql):
        def model(target, record, ask):
            conn.execute(sql)
            return a_fold()
        return model

    before = _fold(conn)
    # Its record written again with other bytes.
    outcome = worker.structure_next(conn, storage, model=changed(
        "update protein_track set sha256 = 'another' where slug = 'ins' and kind = 'record'"))
    assert outcome["state"] == "superseded" and "changed" in outcome["reason"]
    assert _fold(conn)["tracks"] == before["tracks"] and _fold(conn)["generations"] == []
    # Or its regions.
    conn.execute("update bake_job set state = 'queued' where kind = 'structure'")
    conn.execute("update protein_track set sha256 = %s where kind = 'record'", (record_sha,))
    outcome = worker.structure_next(conn, storage, model=changed(
        "update protein set provenance = jsonb_set(provenance, '{regions_from,0,end}', '23') "
        "where slug = 'ins'"))
    assert outcome["state"] == "superseded"
    assert _fold(conn)["generations"] == [] and _fold(conn)["selections"] == []


def test_an_interrupted_publication_leaves_the_fold_drawn_before_it_whole(
        conn, offline, monkeypatch):
    storage = Storage()
    built = _published_alphafold_before_0012(conn, offline, storage)
    before = _fold(conn)
    assert store.queue_structure(conn, "ins", refresh=True, policy=POLICY) == "queued"

    put = store._put_track

    def failing(conn, row):
        if row["kind"] == "folding":
            raise RuntimeError("the connection went away mid-publication")
        put(conn, row)

    monkeypatch.setattr(store, "_put_track", failing)
    outcome = worker.structure_next(conn, storage, model=_modeller(a_fold()))
    assert outcome["state"] == "queued" and "went away" in outcome["reason"]
    after = _fold(conn)
    # Nothing of it: the words, the model and its generation are as they were.
    assert (after["structure"], after["tracks"], after["generations"], after["selections"]) == \
        (before["structure"], before["tracks"], [], [])
    assert after["structure"]["chrome"]["pdb"] == built.chrome["pdb"]
    assert after["jobs"][-1] == ("queued", 1, True)


def test_a_curated_row_is_never_written_put_back_or_queued(conn):
    conn.execute(
        """
        insert into protein (slug, gene, uniprot, accession, display, summary, residues, exons,
            chains, bridges, regions, disulfides, structure, provenance, resolver_version,
            catalog_order)
        values ('insulin', 'INS', 'P01308', 'NG_007114', 'Insulin', 'Curated.', 110, 3, 3, 3,
                '[]', '[]', '{"chrome": {"pdb": "3I40"}, "chains": []}', '{}', 0, 0)
        """)
    conn.execute("insert into protein_track (slug, kind, state, bucket, object_path, sha256, "
                 "format, provenance) values ('insulin', 'record', 'ready', 'tracks', 'r', "
                 "'abc', 'json', '{}')")
    conn.execute("insert into bake_job (slug, kind, state, attempts) "
                 "values ('insulin', 'structure', 'running', 1)")
    job = {"id": 1, "slug": "insulin", "kind": "structure", "attempts": 1}
    with pytest.raises(RuntimeError, match="curated protein"):
        store.publish_fold(conn, job, {"source": "none", "reason": "No.", "record_sha256": "abc",
                                       "regions": "x"})
    with pytest.raises(RuntimeError, match="baked by hand"):
        store.queue_structure(conn, "insulin")
    with pytest.raises(RuntimeError, match="not a protein resolved on demand"):
        store.roll_back_fold(conn, "insulin")
    assert _one(conn, "select structure from protein where slug = 'insulin'") == \
        ({"chrome": {"pdb": "3I40"}, "chains": []},)
    assert _one(conn, "select count(*) from fold_generation") == (0,)
    assert _one(conn, "select count(*) from protein_track where kind <> 'record'") == (0,)


# ------------------------------------------------- generations


def test_a_fold_is_put_back_and_forward_whole(conn, offline):
    storage = Storage()
    built = _published_alphafold_before_0012(conn, offline, storage)
    alphafold = _fold(conn)
    assert store.queue_structure(conn, "ins", refresh=True, policy=POLICY) == "queued"
    worker.structure_next(conn, storage, model=_modeller(a_fold()))
    experimental = _fold(conn)
    # The AlphaFold model it drew before is kept as a generation of its own.
    assert experimental["generations"] == [(1, "alphafold", "ready", False),
                                           (2, "pdb", "ready", True)]
    assert experimental["tracks"]["folding"][0] == "ready"

    # Back to AlphaFold's: its words, its model, and no animation.
    with pytest.raises(RuntimeError, match="queued or running"):
        conn.execute("insert into bake_job (slug, kind) values ('ins', 'structure')")
        store.roll_back_fold(conn, "ins")
    conn.execute("update bake_job set state = 'done' where state = 'queued'")
    assert store.roll_back_fold(conn, "ins") == {"from": 2, "to": 1, "source": "alphafold",
                                                 "state": "ready"}
    back = _fold(conn)
    assert back["structure"] == alphafold["structure"] == \
        {"chrome": built.chrome, "chains": built.chains}
    assert back["tracks"] == {"structure": alphafold["tracks"]["structure"],
                              "folding": ("absent", None, None)}
    assert [(g[0], g[3]) for g in back["generations"]] == [(1, True), (2, False)]

    # And forward again, to the experimental fold, by its generation.
    assert store.roll_back_fold(conn, "ins", to=2)["source"] == "pdb"
    forward = _fold(conn)
    assert (forward["structure"], forward["tracks"]) == \
        (experimental["structure"], experimental["tracks"])
    with pytest.raises(RuntimeError, match="draws now"):
        store.roll_back_fold(conn, "ins", to=2)
    with pytest.raises(RuntimeError, match="no generation 9"):
        store.roll_back_fold(conn, "ins", to=9)


def test_a_refusal_stood_on_before_is_kept_when_its_fold_is_queued_again(conn, offline):
    """A protein with no model before 0012 (AlphaFold DB had none of it), queued
    for a search: its refusal is a generation before its track says a model is
    on its way, so the sentence it showed can be put back."""
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    why = "AlphaFold DB has no model of proteins over 2,700 residues; this one has 4,834."
    job = store.claim_bake(conn, "structure")
    store.refuse_bake(conn, job["id"], "ins", "structure", why)
    refused = _fold(conn)
    assert refused["tracks"] == {"structure": ("refused", None, why)}
    assert refused["generations"] == []

    assert store.queue_structure(conn, "ins", policy=POLICY) == "queued"
    queued = _fold(conn)
    assert queued["tracks"] == {"structure": ("pending", None, None)}
    assert queued["generations"] == [(1, "none", "refused", False)]

    worker.structure_next(conn, storage, model=_modeller(a_fold()))
    assert _fold(conn)["generations"] == [(1, "none", "refused", False),
                                          (2, "pdb", "ready", True)]
    # Put back: no model, no words for one, and the sentence it showed.
    assert store.roll_back_fold(conn, "ins") == {"from": 2, "to": 1, "source": "none",
                                                 "state": "refused"}
    back = _fold(conn)
    assert back["structure"] is None
    assert back["tracks"] == {"structure": ("refused", None, why),
                              "folding": ("absent", None, None)}


def test_a_protein_with_one_fold_has_nothing_to_put_back(conn, offline):
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    worker.structure_next(conn, storage, model=_modeller(a_fold()))
    with pytest.raises(RuntimeError, match="no earlier fold"):
        store.roll_back_fold(conn, "ins")


# ------------------------------------------------- the search between attempts


def test_a_search_cut_short_keeps_how_far_it_got_and_spends_no_try(conn, offline):
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    said = a_selection("incomplete", refused=["1ABC:A@25", "2DEF:B@3"])
    model = _modeller(worker.Unfinished(said))
    outcome = worker.structure_next(conn, storage, model=model)
    assert outcome["state"] == "queued" and "goes on next cycle" in outcome["reason"]
    assert _fold(conn)["jobs"] == [("queued", 0, False)]
    assert _fold(conn)["selections"] == [("incomplete", None)]
    assert _fold(conn)["tracks"]["structure"] == ("pending", None, None)

    # The next attempt is told where it got to.
    model = _modeller(a_fold())
    worker.structure_next(conn, storage, model=model)
    (_, prior), = model.asked
    assert prior["progress"] == {"id": 1, "policy": 1, "inputs": "the-inputs-digest",
                                 "search": "the-ids-digest", "refused": ["1ABC:A@25", "2DEF:B@3"]}
    assert prior["refresh"] is False and prior["none"] is None
    assert _fold(conn)["jobs"] == [("done", 1, False)]
    # Once decided, no progress is left to go on from.
    assert store.fold_prior(conn, "ins")["progress"] is None


def test_a_search_that_found_nothing_is_remembered_until_it_expires(conn, offline):
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    folded = worker.Folded("alphafold", a_model(), selection=a_selection("none"))
    worker.structure_next(conn, storage, model=_modeller(folded))
    assert _fold(conn)["generations"] == [(1, "alphafold", "ready", True)]
    assert "folding" not in _fold(conn)["tracks"]
    none = store.fold_prior(conn, "ins")["none"]
    assert (none["id"], none["policy"], none["inputs"], none["entities"]) == \
        (1, 1, "the-inputs-digest", 3)
    days = _one(conn, "select extract(day from expires_at - decided_at) from fold_selection")[0]
    assert int(days) == 28
    # A search relied on is not recorded again: the fold names the one it used.
    store.queue_structure(conn, "ins", policy=POLICY)
    cached = worker.Folded("alphafold", a_model(), selection={"cached": none["id"]})
    worker.structure_next(conn, storage, model=_modeller(cached))
    assert _fold(conn)["selections"] == [("none", None)]
    assert _one(conn, "select selection_id from fold_generation where active") == (1,)
    # Expired, it is searched again.
    conn.execute("update fold_selection set expires_at = now() - interval '1 second'")
    assert store.fold_prior(conn, "ins")["none"] is None


def test_an_experimental_fold_stands_until_a_search_is_asked_afresh(conn, offline):
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    worker.structure_next(conn, storage, model=_modeller(a_fold()))
    published = _fold(conn)
    assert store.queue_structure(conn, "ins", policy=POLICY) == "kept"
    # A bake asked for all the same, by hand, keeps it and asks no modeller.
    conn.execute("insert into bake_job (slug, kind) values ('ins', 'structure')")
    untouched = _modeller(AssertionError("the modeller was asked"))
    assert worker.structure_next(conn, storage, model=untouched)["state"] == "kept"
    assert _fold(conn)["generations"] == published["generations"]
    assert _fold(conn)["jobs"][-1] == ("done", 1, False)
    # Asked afresh, it searches, and publishes what it finds.
    assert store.queue_structure(conn, "ins", refresh=True, policy=POLICY) == "queued"
    assert store.queue_structure(conn, "ins", refresh=True, policy=POLICY) == "waiting"
    model = _modeller(a_fold(folding={"secondary_structure": "PyMOL 3.1.0 dss, 1ABC"}))
    assert worker.structure_next(conn, storage, model=model)["state"] == "ready"
    assert model.asked[0][1]["refresh"] is True
    assert [g[1:] for g in _fold(conn)["generations"]] == [
        ("pdb", "ready", False), ("pdb", "ready", True)]


# ------------------------------------------------- a backfill


def test_a_backfill_keeps_the_model_drawn_until_another_is_published(conn, offline):
    storage = Storage()
    _published_alphafold_before_0012(conn, offline, storage)
    before = _fold(conn)
    assert store.queue_structure(conn, "ins", policy=POLICY) == "queued"
    # Its model stays ready while the new one is made.
    assert _fold(conn)["tracks"] == before["tracks"]

    # Neither source has one now: it keeps the one it had, and nothing is published.
    none = worker.Unmodelled("AlphaFold DB holds no model of UniProt P01308.",
                             a_selection("none"))
    outcome = worker.structure_next(conn, storage, model=_modeller(none))
    assert outcome["state"] == "refused"
    after = _fold(conn)
    assert (after["structure"], after["tracks"]) == (before["structure"], before["tracks"])
    assert after["generations"] == [(1, "alphafold", "ready", False)]
    assert after["selections"] == [("none", None)]
    assert after["jobs"][-1] == ("failed", 1, False)

    # Nor where its bake breaks, every time it is tried.
    store.queue_structure(conn, "ins", refresh=True, policy=POLICY)
    for _ in range(3):
        outcome = worker.structure_next(conn, storage, model=_modeller(RuntimeError("PyMOL")))
    assert outcome["reason"] == ("The fold could not be made this time: its bake broke each "
                                 "of the 3 times it was tried.")
    assert _fold(conn)["tracks"] == before["tracks"]


def test_a_new_build_drawn_from_neither_source_says_why_and_keeps_its_search(conn, offline):
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    why = ("The PDB holds no experimental structure of it. AlphaFold DB holds no model of "
           "UniProt P01308.")
    worker.structure_next(conn, storage, model=_modeller(
        worker.Unmodelled(why, a_selection("none", search={"entities": 0}))))
    said = _fold(conn)
    assert said["tracks"] == {"structure": ("refused", None, why)}
    assert said["structure"] is None
    assert said["generations"] == [(1, "none", "refused", True)]
    assert said["selections"] == [("none", None)]


# ------------------------------------------------- what a reader is told


def test_a_build_is_held_at_its_structure_while_its_first_fold_is_made(
        conn, offline, database_url, monkeypatch):
    from app import resolves

    storage = Storage()
    with _service(database_url, monkeypatch):
        _resolved_with_a_model_queued(conn, offline, storage)
        assert worker.score_next(conn, storage, score=_constraint)["state"] == "ready"
        # Scored, and its fold still queued: held, with nothing to end it.
        build = resolves.current("INS")["build"]
        assert (build["step"], build["reason"], build["fold"]) == ("structure", None, "pending")
        assert _builds(conn) == {("ins", "building")}

        # Made: done, with its fold.
        worker.structure_next(conn, storage, model=_modeller(a_fold()))
        build = resolves.current("INS")["build"]
        assert (build["step"], build["fold"]) == ("done", "ready")
        assert _builds(conn) == set()


def test_a_build_is_held_for_its_fold_no_longer_than_the_hold(
        conn, offline, database_url, monkeypatch):
    from app import resolves

    storage = Storage()
    with _service(database_url, monkeypatch):
        _resolved_with_a_model_queued(conn, offline, storage)
        worker.score_next(conn, storage, score=_constraint)
        conn.execute("update bake_job set requested_at = now() - interval '21 minutes' "
                     "where kind = 'structure'")
        # Queued all this while: no worker is making it, and the build opens.
        build = resolves.current("INS")["build"]
        assert (build["step"], build["fold"]) == ("done", "pending")
        assert _builds(conn) == set()
        # One being made is held, however long ago it was asked for.
        store.claim_bake(conn, "structure")
        assert resolves.current("INS")["build"]["step"] == "structure"


def test_a_build_made_without_models_is_never_held(conn, offline, database_url, monkeypatch):
    from app import resolves

    from pipeline.resolver.test_worker_pg import _ask, _body

    storage = Storage()
    with _service(database_url, monkeypatch):
        _ask(conn)
        worker.resolve_next(conn, storage, fetch_entry=lambda accession: offline(_body()))
        worker.score_next(conn, storage, score=_constraint)
        build = resolves.current("INS")["build"]
        assert (build["step"], build["fold"]) == ("done", None)


def test_a_stopped_build_is_not_held_for_its_fold(conn, offline, database_url, monkeypatch):
    from app import resolves

    storage = Storage()
    with _service(database_url, monkeypatch):
        _resolved_with_a_model_queued(conn, offline, storage, asker="install-a")
        resolves.stop("INS", "install-a")
        build = resolves.current("INS")["build"]
        assert (build["step"], build["stopped"], build["fold"]) == ("scoring", True, "pending")
        assert _builds(conn) == {("ins", "stopped")}


def test_the_queue_says_what_each_protein_draws_and_how_its_search_ended(conn, offline):
    storage = Storage()
    _resolved_with_a_model_queued(conn, offline, storage)
    worker.structure_next(conn, storage, model=_modeller(a_fold()))
    lines = local_worker.queue(conn).splitlines()
    assert lines[2] == "folds     1 pdb ready; searches 1 chosen"


# ------------------------------------------------- the operator's two commands


def _script(name: str):
    import importlib

    scripts = str(BACKEND / "scripts")
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    return importlib.import_module(name)


def test_the_backfill_and_the_rollback_do_what_they_say(
        conn, offline, database_url, monkeypatch, capsys):
    from app.config import settings

    monkeypatch.setattr(settings, "database_url", database_url)
    queue_structures, fold_rollback = _script("queue_structures"), _script("fold_rollback")
    storage = Storage()
    _published_alphafold_before_0012(conn, offline, storage)

    assert queue_structures.main(["--dry-run", "--all"]) == 0
    assert capsys.readouterr().out == \
        "ins: draws AlphaFold DB AF-P01308-F1; would queue a structure bake\n"
    assert _one(conn, "select count(*) from bake_job where kind = 'structure' "
                      "and state = 'queued'") == (0,)
    assert queue_structures.main(["INS"]) == 0
    assert capsys.readouterr().out == "ins: draws AlphaFold DB AF-P01308-F1; queued\n"
    assert queue_structures.main(["INS"]) == 0
    assert capsys.readouterr().out == "ins: draws AlphaFold DB AF-P01308-F1; waiting\n"
    assert queue_structures.main(["TTR"]) == 1
    assert "TTR: no protein built on demand" in capsys.readouterr().err

    worker.structure_next(conn, storage, model=_modeller(a_fold()))
    assert queue_structures.main(["--dry-run", "INS"]) == 0
    assert capsys.readouterr().out == \
        "ins: draws PDB 1ABC, animated; would queue nothing: its experimental fold stands\n"

    assert fold_rollback.main(["INS", "--list"]) == 0
    listed = capsys.readouterr().out.splitlines()
    assert [line[:6] for line in listed] == ["*    2", "     1"]
    assert listed[0].endswith("pdb       PDB 1ABC, animated")
    assert listed[1].endswith("alphafold AlphaFold DB AF-P01308-F1")
    assert fold_rollback.main(["INS", "--dry-run"]) == 0
    assert capsys.readouterr().out.splitlines()[-1] == \
        "ins: would put back generation 1 (alphafold)"
    assert _fold(conn)["tracks"]["folding"][0] == "ready"

    assert fold_rollback.main(["INS"]) == 0
    assert capsys.readouterr().out.splitlines()[0] == \
        "ins: generation 1 (alphafold, ready) put back over 2"
    assert _fold(conn)["tracks"]["folding"] == ("absent", None, None)
    assert fold_rollback.main(["INS", "--to", "2"]) == 0
    capsys.readouterr()
    assert _fold(conn)["tracks"]["folding"][0] == "ready"
    assert fold_rollback.main(["TTR"]) == 1
