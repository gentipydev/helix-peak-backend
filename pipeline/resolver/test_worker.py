"""The worker's own decisions, offline. Its SQL is `test_worker_pg.py`'s."""

from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import replace
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import upload_tracks  # noqa: E402
from pipeline.resolver import scoring, store, worker  # noqa: E402
from pipeline.resolver.resolve import resolve  # noqa: E402
from pipeline.resolver.test_resolve import INS, _body  # noqa: E402


class Storage:
    def __init__(self):
        self.objects = {}

    def put(self, bucket, path, payload, content_type):
        self.objects[(bucket, path, content_type)] = payload


def test_a_track_is_stored_where_the_uploader_would_store_it(offline):
    resolution = resolve(INS, offline(_body()))
    storage = Storage()
    row = worker.stage(storage, "record", resolution.target, resolution.record)

    digest = hashlib.sha256(resolution.record).hexdigest()
    assert row == {
        "slug": "ins", "kind": "record", "bucket": "tracks",
        "object_path": f"record/ins.{digest[:12]}.json",
        "bytes": len(resolution.record), "sha256": digest, "format": "json",
        "provenance": {"source": "NCBI Entrez", "accession": "NG_007114.1",
                       "protein_id": "NP_000198.1", "transcript_id": "NM_000207.3",
                       "built_by": "pipeline/mock/build_gene_record.py"},
    }
    assert storage.objects == {
        ("tracks", row["object_path"], "application/json"): resolution.record,
    }


def test_a_payload_the_upload_gate_refuses_is_never_stored(offline):
    resolution = resolve(INS, offline(_body()))
    storage = Storage()
    wrong = json.dumps({**json.loads(resolution.record), "gene": "TH"}).encode()
    try:
        worker.stage(storage, "record", resolution.target, wrong)
    except ValueError as exc:
        assert "expected 'INS'" in str(exc)
    else:
        raise AssertionError("a record for another gene was staged")
    assert storage.objects == {}


def _another_genes_record(row, entry, **said):
    """What `resolve` returns, carrying a record the upload gate declines."""
    resolution = resolve(row, entry, **said)
    wrong = json.dumps({**json.loads(resolution.record), "gene": "TH"}).encode()
    return replace(resolution, record=wrong)


def _one_request(monkeypatch):
    """`store`'s side of `resolve_next` with no database: one request for INS,
    and what the worker did with it, as (finished, requeued)."""
    finished, requeued = [], []
    monkeypatch.setattr(store, "claim_request", lambda conn: {
        "id": 7, "gene": "INS", "uniprot": "P01308", "slug": "ins", "attempts": 1})
    monkeypatch.setattr(store, "protein_for_gene", lambda conn, gene: None)
    monkeypatch.setattr(store, "index_row", lambda conn, uniprot, gene: INS)
    monkeypatch.setattr(store, "mane_release", lambda conn: "v1.5")
    monkeypatch.setattr(store, "finish_request", lambda conn, request_id, state, **said:
                        finished.append((request_id, state, said)))
    monkeypatch.setattr(store, "requeue_request", lambda conn, request_id:
                        requeued.append(request_id))
    return finished, requeued


def test_a_record_the_upload_gate_declines_is_refused_and_not_tried_again(offline, monkeypatch):
    monkeypatch.setattr(worker, "resolve", _another_genes_record)
    finished, requeued = _one_request(monkeypatch)
    storage = Storage()
    outcome = worker.resolve_next(None, storage, fetch_entry=lambda accession: offline(_body()))

    assert outcome["state"] == "refused" and "expected 'INS'" in outcome["reason"]
    assert finished == [(7, "refused", {"reason": outcome["reason"], "resolver_version": 1})]
    assert requeued == []
    assert storage.objects == {}


def test_a_record_that_cannot_be_stored_is_tried_again(offline, monkeypatch):
    # Only the gate's own verdict is a refusal. Storage answering wrongly, even
    # with a ValueError (a storage URL that is not one), is worth another try.
    class NoStorage:
        def put(self, bucket, path, payload, content_type):
            raise ValueError("unknown url type: 'project.supabase.co/storage/v1'")

    finished, requeued = _one_request(monkeypatch)
    outcome = worker.resolve_next(None, NoStorage(),
                                  fetch_entry=lambda accession: offline(_body()))

    assert outcome["state"] == "queued" and "unknown url type" in outcome["reason"]
    assert (finished, requeued) == ([], [7])


def test_a_failed_alignment_gate_is_told_as_a_sentence():
    said = worker._refusal(ValueError(
        "Alignment gate FAILED: 48.2% before, 51.0% after. No JSON written. "
        "Inspect token offsets and how scores are filed before any UI work."))
    assert said == (
        "ESM-2 prefers the residue that is there to the one before it 48.2% of the time, "
        "and to the one after it 51.0%. More than half is needed to check that its scores "
        "sit on their own residues, and a protein the model knows little about falls short. "
        "The scores are not drawn.")
    # What was measured, and no claim that anything was misfiled. As long as
    # any other reason is allowed to be, at its longest.
    assert "filed" not in said
    assert len(said.replace("48.2%", "100.0%").replace("51.0%", "100.0%")) <= 300


def test_any_other_refusal_is_its_own_message_on_one_line():
    assert worker._refusal(ValueError("ins asset is 1,000 bytes,\n over its budget")) == \
        "ins asset is 1,000 bytes, over its budget"
    assert len(worker._said(RuntimeError("x" * 1000))) == 300


def test_the_worker_scores_with_the_function_the_scorers_environment_imports(monkeypatch):
    # One function, not a copy of it: what Modal calls in the worker's own
    # process is what an environment holding only the scorer calls, with its
    # printed progress passed on as it prints it.
    assert worker.score_with_esm is scoring.score_with_esm
    assert worker.score_next.__kwdefaults__["score"] is worker.score_in_process
    assert worker.score_all.__kwdefaults__["score"] is worker.score_in_process

    called = []
    monkeypatch.setattr(worker, "score_with_esm",
                        lambda target, record: called.append((target, record)) or b"track")
    assert worker.score_in_process("target", b"record", lambda *said: None) == b"track"
    assert called == [("target", b"record")]


def test_the_scorers_progress_lines_are_read_as_it_prints_them():
    assert worker.progress_of("") is None
    assert worker.progress_of("ins: 110 residues, one full-length pass; vocab=33.") is None
    # The last line wins, and the seconds left are its rate over what is left:
    # 283.3 s for 450 residues, with 388 still to score.
    said = "Scored 25/838 (15.9s, 8.6 min left)\nScored 450/838 (283.3s, 4.1 min left)\n"
    done, total, left = worker.progress_of(said)
    assert (done, total) == (450, 838)
    assert left == pytest.approx(283.3 / 450 * 388)
    assert worker.progress_of("Scored 838/838 (533.0s, 0.0 min left)") == (838, 838, 0.0)
    # Nonsense is not progress.
    assert worker.progress_of("Scored 0/838 (0.0s, 0.0 min left)") is None
    assert worker.progress_of("Scored 900/838 (1.0s, 0.0 min left)") is None


def test_progress_printed_in_this_process_is_told_and_still_printed(monkeypatch, capsys):
    def scorer(target, record):
        print("oca2: 838 residues, one full-length pass; vocab=33.", flush=True)
        print("Scored 25/838 (15.9s, 8.6 min left)", flush=True)
        # A line written in pieces is read once it is whole.
        sys.stdout.write("Scored 50/838 ")
        sys.stdout.write("(31.8s, 8.4 min left)\n")
        return b"track"

    monkeypatch.setattr(worker, "score_with_esm", scorer)
    told = []
    assert worker.score_in_process("target", b"record", lambda *said: told.append(said)) == \
        b"track"
    assert [(done, total) for done, total, _ in told] == [(25, 838), (50, 838)]
    assert "Scored 50/838 (31.8s, 8.4 min left)" in capsys.readouterr().out


def test_progress_that_cannot_be_written_costs_the_bake_nothing():
    class Unreachable:
        def execute(self, *args):
            raise OSError("server closed the connection unexpectedly")

    report = worker._Reporter(Unreachable(), 7)
    report(450, 838, 245.0)
    report.check()


class _Rows:
    """A connection whose bake is in `state`: what the reporter's two
    statements read back."""

    def __init__(self, state):
        self.state = state

    def execute(self, sql, params=None):
        flat = " ".join(sql.split())
        running = self.state == "running"

        class Cursor:
            def fetchone(cursor):
                if flat.startswith("update bake_job set progress_done"):
                    return (7,) if running else None
                return (self.state,)

        return Cursor()


def test_a_bake_the_reader_stopped_is_stopped_at_its_next_line_or_check():
    running = worker._Reporter(_Rows("running"), 7)
    running(450, 838, 245.0)
    running.check()

    stopped = worker._Reporter(_Rows("stopped"), 7)
    with pytest.raises(worker.Stopped):
        stopped(475, 838, 230.0)
    with pytest.raises(worker.Stopped):
        stopped.check()


def test_a_stop_printed_through_is_raised_out_of_the_scorer(monkeypatch):
    def scorer(target, record):
        print("Scored 25/838 (15.9s, 8.6 min left)", flush=True)
        raise AssertionError("the scorer went on after it was stopped")

    def stopped(*said):
        raise worker.Stopped()

    monkeypatch.setattr(worker, "score_with_esm", scorer)
    with pytest.raises(worker.Stopped):
        worker.score_in_process("target", b"record", stopped)


# ------------------------------------------------------------ a model


def a_glb(*nodes: str) -> bytes:
    """The least a `.glb` can be and still name its nodes: its JSON chunk."""
    import struct

    document = json.dumps({
        "asset": {"version": "2.0"}, "scenes": [{"nodes": list(range(len(nodes)))}],
        "nodes": [{"name": name, "mesh": index} for index, name in enumerate(nodes)],
        "meshes": [{"primitives": []} for _ in nodes],
    }).encode()
    document += b" " * (-len(document) % 4)
    return (struct.pack("<4sII", b"glTF", 2, 20 + len(document))
            + struct.pack("<I4s", len(document), b"JSON") + document)


def a_model(*nodes: str, **changes) -> worker.Modelled:
    """What a baker hands over for a model of these nodes (`alphafold.Built`)."""
    nodes = nodes or ("plddtVeryHigh", "plddtLow", "bonds")
    said = dict(
        glb=a_glb(*nodes), scene=b"a compiled scene",
        chrome={"pdb": "AF-P01308-F1", "modelled": [25, 110], "label": "the predicted fold",
                "count": 86, "unit": "residues",
                "sentence": "AlphaFold prediction, mean pLDDT 71.0. 60% of residues at 70 or over.",
                "semantics": "AlphaFold's predicted fold of Insulin. Drag to turn it."},
        chains=[{"node": node, "tint": "cysteine" if node == "bonds" else node}
                for node in nodes],
        provenance={"pdb": "AF-P01308-F1", "nodes": list(nodes), "source": "AlphaFold DB",
                    "entry": "AF-P01308-F1", "model_version": 6, "mean_plddt": 71.0,
                    "licence": "CC BY 4.0", "span": [25, 110], "sampling": 8},
    )
    said.update(changes)
    return worker.Modelled(**said)


def test_a_model_is_stored_where_the_uploader_would_store_one(offline):
    target = resolve(INS, offline(_body())).target
    storage = Storage()
    built = a_model()
    row, structure = worker.stage_structure(storage, target, built)

    glb, scene = (hashlib.sha256(payload).hexdigest() for payload in (built.glb, built.scene))
    assert row == {
        "slug": "ins", "kind": "structure", "bucket": "models",
        # The row is the scene the phone loads; the `.glb` it came from is beside it.
        "object_path": f"structure/ins.{scene[:12]}.fsceneb",
        "bytes": len(built.scene), "sha256": scene, "format": "fsceneb",
        "provenance": {**built.provenance, "glb": {
            "path": f"structure/ins.{glb[:12]}.glb", "sha256": glb, "bytes": len(built.glb)}},
    }
    assert storage.objects == {
        ("models", f"structure/ins.{glb[:12]}.glb", "model/gltf-binary"): built.glb,
        ("models", row["object_path"], "application/octet-stream"): built.scene,
    }
    # The protein row's `structure` column, as the twenty carry it.
    assert structure == {"chrome": built.chrome, "chains": built.chains}


def test_the_uploader_names_a_models_objects_as_it_named_the_twentys():
    # Insulin's stored row, from its two digests: `structure_row` is the code
    # that branch of the uploader ran inline.
    row, uploads = upload_tracks.structure_row(
        "insulin", b"glTF-not-the-real-bytes", b"not-the-real-scene",
        {"pdb": "3I40", "nodes": ["chainA", "chainB"]})
    glb, scene = (hashlib.sha256(payload).hexdigest()
                  for payload in (b"glTF-not-the-real-bytes", b"not-the-real-scene"))
    assert [(path, mime) for path, _, mime in uploads] == [
        (f"structure/insulin.{glb[:12]}.glb", "model/gltf-binary"),
        (f"structure/insulin.{scene[:12]}.fsceneb", "application/octet-stream")]
    assert (row["bucket"], row["format"], row["object_path"], row["sha256"]) == (
        "models", "fsceneb", uploads[1][0], scene)
    assert list(row["provenance"]) == ["pdb", "nodes", "glb"]


@pytest.mark.parametrize("changes, why", [
    ({"glb": b"not a model"}, "not a .glb"),
    ({"scene": b""}, "no compiled scene"),
    ({"chrome": {"pdb": "AF-P01308-F1", "label": "the predicted fold"}}, "not the seven"),
    ({"chains": []}, "names no tint"),
    ({"chains": [{"node": "plddtVeryHigh", "tint": "mauve"}]}, "names no tint"),
    ({"chains": [{"node": "plddtVeryHigh", "tint": "plddtVeryHigh"}]}, "its row names"),
    ({"provenance": {"pdb": "AF-P99999-F1"}}, "another entry"),
])
def test_a_model_the_app_would_trip_on_is_never_stored(offline, changes, why):
    target = resolve(INS, offline(_body())).target
    storage = Storage()
    with pytest.raises(worker.Declined, match=why):
        worker.stage_structure(storage, target, a_model(**changes))
    assert storage.objects == {}


def test_a_model_is_made_for_the_protein_as_its_record_and_row_have_it(offline):
    resolution = resolve(INS, offline(_body()))
    job = worker.model_job(resolution.target, resolution.record)
    protein = json.loads(resolution.record)["protein"]["translation"]
    assert job == {
        "slug": "ins", "accession": "P01308", "display": "Insulin", "protein": protein,
        # The kept regions alone: no signal peptide, no C-peptide, no cut site.
        "kept": [[25, 54], [90, 110]],
        "disulfides": [[31, 96], [43, 109], [95, 100]],
        # As far as the resolver let the record differ from UniProt.
        "allowed": 3,
    }
    assert len(protein) == 110
    long = replace(resolution.target, regions=())
    record = json.dumps({"protein": {"translation": "M" * 2500}}).encode()
    assert worker.model_job(long, record)["allowed"] == 25
    assert worker.model_job(long, record)["kept"] == []


def test_the_worker_and_the_bake_agree_on_the_words_and_the_tints():
    alphafold = pytest.importorskip("pipeline.structure.alphafold")   # needs numpy
    assert worker._TINTS == {name for name, _ in alphafold.BANDS} | {"cysteine"}
    job = alphafold.Job(slug="x", accession="P00000", display="X", protein="M" * 9,
                        kept=((1, 9),))
    entry = alphafold.Entry("P00000", "AF-P00000-F1", 6, "https://example.invalid/m.pdb")
    chrome, chains = alphafold.describe(
        job, entry, (1, 9), alphafold.Confidence(80.0, (0.5, 0.5, 0.0, 0.0)), [(2, 8)],
        ["plddtVeryHigh", "plddtConfident", alphafold.BONDS])
    assert tuple(chrome) == worker._CHROME
    assert {chain["tint"] for chain in chains} <= worker._TINTS


def _one_bake(monkeypatch, attempts=1):
    """`store`'s side of `structure_next` with no database: one bake for INS,
    and what the worker did with it."""
    did = []
    monkeypatch.setattr(store, "claim_bake", lambda conn, kind: {
        "id": 4, "slug": "ins", "kind": kind, "attempts": attempts})
    monkeypatch.setattr(store, "protein", lambda conn, slug: None)
    monkeypatch.setattr(store, "ready_track", lambda conn, slug, kind: None)
    for name in ("finish_structure", "refuse_bake", "requeue_bake"):
        monkeypatch.setattr(store, name, lambda conn, *said, _name=name, **more:
                            did.append((_name, said, more)))
    return did


def test_a_protein_with_no_record_has_no_model_to_be_held_to(monkeypatch):
    did = _one_bake(monkeypatch)
    outcome = worker.structure_next(None, Storage(), model=lambda target, record: a_model())
    assert outcome == {"id": 4, "slug": "ins", "state": "refused",
                       "reason": "The protein has no ready record to hold a model to."}
    assert did == [("refuse_bake", (4, "ins", "structure", outcome["reason"]), {})]
