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


# ------------------------------------------------------------ the evidence


class Stored(Storage):
    """Storage that serves what was put, as the public read every client makes."""

    def get(self, bucket, path):
        return next(payload for (stored, named, _), payload in self.objects.items()
                    if (stored, named) == (bucket, path))


def _row(storage, kind, payload):
    """A ready track's row, its bytes stored where the row names them."""
    digest = hashlib.sha256(payload).hexdigest()
    path = f"{kind}/ins.{digest[:12]}.json"
    storage.put("tracks", path, payload, "application/json")
    return {"bucket": "tracks", "object_path": path, "sha256": digest}


def an_avi(target, record: bytes, **changes) -> bytes:
    """An AVI track the upload gate takes for this record
    (`check_assets.check_impact`): its letters, one run over every drawn base,
    three scores a base, and exons above their introns' edges, which are above
    the introns' interiors."""
    from pipeline.check_assets import bases

    mock = json.loads(record)
    start, end = mock["location"]["start"], mock["location"]["end"]
    exons = sorted((exon["start"], exon["end"]) for exon in mock["exons"])
    exonic = {local for low, high in exons for local in range(low, high + 1)}
    edges = set()
    for (_, before), (after, _) in zip(exons, exons[1:]):
        low, high = before + 1, after - 1
        edge = min(8, (high - low + 1) // 2)
        edges |= set(range(low, low + edge)) | set(range(high - edge + 1, high + 1))
    asset = {
        "gene": target.gene, "uniprot": target.uniprot, "accession": target.source.accession,
        "assembly": "GRCh38", "annotation": "GENCODE v46", "chromosome": "chr11",
        "transcript": "ENST00000381330.5", "orientation": -1, "complemented": True,
        "scorer": "AVI_SCORE", "score_units": "phred", "alt_order": "ACGT minus wildtype",
        "high_phred": 20.0, "middle_phred": 10.0, "start": start,
        "sequence": bases(mock, list(range(start, end + 1))),
        "generation": {"drawn_bases": end - start + 1},
        "runs": [{"local": start, "genomic": 2_161_209, "step": -1, "length": end - start + 1}],
        "positions": {str(local): [20.0 if local in exonic else 12.0 if local in edges else 2.0] * 3
                      for local in range(start, end + 1)},
    }
    asset.update(changes)
    return json.dumps(asset, separators=(",", ":")).encode()


def a_clinvar(target, record: bytes, avi: bytes, **changes) -> bytes:
    """A ClinVar track the upload gate takes, placed by this AVI track's map:
    two records searched, neither a single-base substitution."""
    from pipeline.clinvar.bake_clinvar import SCOPE

    mock, impact = json.loads(record), json.loads(avi)
    asset = {
        "schema_version": 1, "source": "NCBI ClinVar", "gene": target.gene,
        "accession": target.source.accession, "assembly": "GRCh38",
        "chromosome": impact["chromosome"], "scope": SCOPE, "start": impact["start"],
        "sequence": impact["sequence"], "protein_sequence": mock["protein"]["translation"],
        "complemented": impact["complemented"], "runs": impact["runs"],
        "retrieved_at": "2026-10-09T12:00:00+00:00", "query": f"{target.gene}[gene]",
        "searched_records": 2, "excluded": {"not_a_single_base_substitution": 2},
        "source_batches": [{"file": "batch-0.xml", "sha256": "0" * 64}],
        "traits": {}, "variants": [],
    }
    asset.update(changes)
    return (json.dumps(asset, indent=2) + "\n").encode()


@pytest.fixture
def data(monkeypatch, tmp_path):
    """`paths.DATA` as the worker's own directory, wherever a module took it."""
    from pipeline import check_assets, paths

    folder = tmp_path / "data"
    for module in (paths, upload_tracks, check_assets):
        monkeypatch.setattr(module, "DATA", folder)
    return folder


def _evidence(monkeypatch, resolution, kind, storage, impact=None, attempts=1,
              impact_track=("absent", None)):
    """`store`'s side of `impact_next` and `clinvar_next` with no database: one
    bake of `kind` for INS, its record stored ready (and its AVI track, where
    one is given), and what the worker did with it."""
    did = []
    rows = {"record": _row(storage, "record", resolution.record)}
    if impact is not None:
        rows["impact"] = _row(storage, "impact", impact)
    monkeypatch.setattr(store, "claim_bake", lambda conn, wanted: {
        "id": 5, "slug": "ins", "kind": wanted, "attempts": attempts} if wanted == kind else None)
    monkeypatch.setattr(store, "protein", lambda conn, slug: {"slug": slug})
    monkeypatch.setattr(worker, "target_of", lambda protein: resolution.target)
    monkeypatch.setattr(store, "ready_track", lambda conn, slug, wanted: rows.get(wanted))
    monkeypatch.setattr(store, "track_state", lambda conn, slug, wanted: impact_track)
    for name in ("finish_bake", "refuse_bake", "requeue_bake"):
        monkeypatch.setattr(store, name, lambda conn, *said, _name=name, **more:
                            did.append((_name, said, more)))
    return did


def _raising(error):
    def bake(*said):
        raise error
    return bake


def test_an_avi_track_is_stored_ready_as_the_uploader_would_store_it(offline, monkeypatch, data):
    resolution = resolve(INS, offline(_body()))
    storage = Stored()
    did = _evidence(monkeypatch, resolution, "impact", storage)
    payload = an_avi(resolution.target, resolution.record)
    handed = []
    outcome = worker.impact_next(
        None, storage, avi=lambda target, record: handed.append((target, record)) or payload)

    assert outcome == {"id": 5, "slug": "ins", "state": "ready"}
    assert handed == [(resolution.target, resolution.record)]
    [(name, (job, track), more)] = did
    digest = hashlib.sha256(payload).hexdigest()
    assert (name, job, more) == ("finish_bake", 5, {})
    assert (track["kind"], track["bucket"], track["object_path"], track["sha256"]) == \
        ("impact", "tracks", f"impact/ins.{digest[:12]}.json", digest)
    assert (track["provenance"]["scorer"], track["provenance"]["chromosome"]) == \
        ("AVI_SCORE", "chr11")
    assert storage.get("tracks", track["object_path"]) == payload
    # The gate held it to the record, laid out in the worker's own directory.
    assert (data / resolution.target.mock_asset).read_bytes() == resolution.record


def test_an_avi_gate_is_a_readers_sentence_on_the_track_and_its_words_on_the_job(
        offline, monkeypatch, data):
    resolution = resolve(INS, offline(_body()))
    storage = Stored()
    did = _evidence(monkeypatch, resolution, "impact", storage)
    gate = ("INS: the coordinate map is wrong — 1,071 of 1,431 bases (74.8%) disagree with the "
            "reference:\n      local 4986 -> 2161209: draws G, reference A")
    outcome = worker.impact_next(None, storage, avi=_raising(worker.Unplaced(gate)))

    said = ("AlphaGenome's scores were not placed on INS's bases: the coordinate map is wrong — "
            "1,071 of 1,431 bases (74.8%) disagree with the reference.")
    assert outcome == {"id": 5, "slug": "ins", "state": "refused", "reason": said}
    assert did == [("refuse_bake", (5, "ins", "impact", said), {"error": worker._said(gate)})]
    assert [path for (_, path, _) in storage.objects if path.startswith("impact/")] == []


@pytest.mark.parametrize("changes, why", [
    ({"sequence": "A" * 1431}, "the impact track's sequence is not the record's drawn letters"),
    ({"positions": {}}, "only 0 of 1,431 drawn bases are scored"),
    ({"uniprot": "P99999"}, "impact track says uniprot='P99999', expected 'P01308'"),
])
def test_an_avi_track_the_gate_declines_is_refused_and_never_stored(
        offline, monkeypatch, data, changes, why):
    resolution = resolve(INS, offline(_body()))
    storage = Stored()
    did = _evidence(monkeypatch, resolution, "impact", storage)
    payload = an_avi(resolution.target, resolution.record, **changes)
    outcome = worker.impact_next(None, storage, avi=lambda target, record: payload)

    assert outcome["state"] == "refused"
    assert outcome["reason"] == f"AlphaGenome's scores were not placed on INS's bases: {why}."
    [(name, said, more)] = did
    assert name == "refuse_bake" and why in more["error"]
    assert [path for (_, path, _) in storage.objects if path.startswith("impact/")] == []


def test_an_atlas_that_does_not_answer_is_tried_again_and_a_reader_is_not_told_how(
        offline, monkeypatch, data):
    resolution = resolve(INS, offline(_body()))
    broke = RuntimeError("The AVI bake exited with status 1: chr11:2159779-2161209 failed "
                         "after 7 attempts: RESOURCE_EXHAUSTED")
    storage = Stored()
    did = _evidence(monkeypatch, resolution, "impact", storage, attempts=1)
    outcome = worker.impact_next(None, storage, avi=_raising(broke))
    assert (outcome["state"], did) == ("queued", [("requeue_bake", (5, worker._said(broke)), {})])

    storage = Stored()
    did = _evidence(monkeypatch, resolution, "impact", storage, attempts=store.MAX_ATTEMPTS)
    outcome = worker.impact_next(None, storage, avi=_raising(broke))
    said = ("AlphaGenome's scores could not be fetched: the bake broke each of the 3 times it "
            "was tried.")
    assert outcome == {"id": 5, "slug": "ins", "state": "refused", "reason": said}
    assert did == [("refuse_bake", (5, "ins", "impact", said), {"error": worker._said(broke)})]


def test_a_protein_with_no_record_has_no_bases_for_either(monkeypatch):
    did = _one_bake(monkeypatch)
    assert worker.impact_next(None, Stored(), avi=_raising(AssertionError()))["reason"] == \
        "The protein has no ready record to place AlphaGenome's scores on."
    assert worker.clinvar_next(None, Stored(), clinvar=_raising(AssertionError()))["reason"] == \
        "The protein has no ready record to place ClinVar's records on."
    assert [(name, said[2]) for name, said, _ in did] == \
        [("refuse_bake", "impact"), ("refuse_bake", "clinvar")]


def test_clinvar_is_placed_by_the_avi_tracks_map(offline, monkeypatch, data):
    resolution = resolve(INS, offline(_body()))
    avi = an_avi(resolution.target, resolution.record)
    storage = Stored()
    did = _evidence(monkeypatch, resolution, "clinvar", storage, impact=avi)

    def clinvar(target):
        # The two it reads, where the baker reads them.
        assert (data / target.mock_asset).read_bytes() == resolution.record
        assert (data / target.impact_asset).read_bytes() == avi
        return a_clinvar(target, resolution.record, avi)

    assert worker.clinvar_next(None, storage, clinvar=clinvar) == \
        {"id": 5, "slug": "ins", "state": "ready"}
    [(name, (job, track), _)] = did
    assert (name, job, track["kind"]) == ("finish_bake", 5, "clinvar")
    assert track["object_path"].startswith("clinvar/ins.")
    assert (track["provenance"]["source"], track["provenance"]["searched_records"]) == \
        ("NCBI ClinVar", 2)


def test_clinvar_is_refused_where_avi_was_and_says_why_in_avis_words(offline, monkeypatch, data):
    resolution = resolve(INS, offline(_body()))
    unplaced = "AlphaGenome's scores were not placed on INS's bases: record has 3 exons, MANE " \
               "Select has 4."
    did = _evidence(monkeypatch, resolution, "clinvar", Stored(),
                    impact_track=("refused", unplaced))
    outcome = worker.clinvar_next(None, Stored(), clinvar=_raising(AssertionError(
        "ClinVar was fetched with no map to place it by")))
    said = ("ClinVar's records are placed by AlphaGenome's coordinate map, which INS has none "
            "of. " + unplaced)
    assert outcome == {"id": 5, "slug": "ins", "state": "refused", "reason": said}
    assert did == [("refuse_bake", (5, "ins", "clinvar", said), {})]

    # With no AVI track at all, there is nothing more to say.
    did = _evidence(monkeypatch, resolution, "clinvar", Stored())
    assert worker.clinvar_next(None, Stored(), clinvar=_raising(AssertionError()))["reason"] == \
        "ClinVar's records are placed by AlphaGenome's coordinate map, which INS has none of."


def test_clinvars_verdict_is_a_refusal_and_ncbi_not_answering_is_tried_again(
        offline, monkeypatch, data):
    resolution = resolve(INS, offline(_body()))
    avi = an_avi(resolution.target, resolution.record)

    storage = Stored()
    did = _evidence(monkeypatch, resolution, "clinvar", storage, impact=avi)
    outcome = worker.clinvar_next(None, storage, clinvar=_raising(ValueError(
        "Unaccounted records")))
    said = "ClinVar's records were not placed on INS's bases: Unaccounted records."
    assert (outcome["state"], outcome["reason"]) == ("refused", said)
    assert did == [("refuse_bake", (5, "ins", "clinvar", said), {"error": "Unaccounted records"})]

    storage = Stored()
    did = _evidence(monkeypatch, resolution, "clinvar", storage, impact=avi)
    wrong = a_clinvar(resolution.target, resolution.record, avi, runs=[])
    outcome = worker.clinvar_next(None, storage, clinvar=lambda target: wrong)
    assert outcome["reason"] == ("ClinVar's records were not placed on INS's bases: mapping "
                                 "differs from the gene/AVI track.")
    assert [path for (_, path, _) in storage.objects if path.startswith("clinvar/")] == []

    missed = worker.Unfetched("ClinVar did not answer in full: Batch returned 99 of the 100 "
                              "records asked for")
    storage = Stored()
    did = _evidence(monkeypatch, resolution, "clinvar", storage, impact=avi)
    outcome = worker.clinvar_next(None, storage, clinvar=_raising(missed))
    assert (outcome["state"], did) == ("queued", [("requeue_bake", (5, worker._said(missed)), {})])


def test_clinvar_is_fetched_beside_the_worker_and_its_raw_responses_go_with_the_bake(
        offline, monkeypatch, data, tmp_path, capsys):
    from pipeline.clinvar import bake_clinvar

    target = resolve(INS, offline(_body())).target
    seen = {}

    def fetch(gene, cache, replay):
        # The baker's own, as it runs: its cache, and its batches.
        seen.update(cache=cache, replay=replay)
        cache.mkdir(parents=True)
        (cache / "batch-0.xml").write_text("<ClinVarResult-Set/>")
        seen["batch"] = bake_clinvar.fetch_batch("efetch", ["1"])
        return "root", {"batches": []}

    def bake(target, replay):
        bake_clinvar.fetch(target.gene, Path(bake_clinvar.__file__).parent / "cache" / target.gene,
                           replay)
        print("INS: 0 mapped SNVs")
        path = data / f"assets/clinvar/{target.slug}_clinvar.json"
        path.parent.mkdir(parents=True)
        path.write_bytes(b"the track\n")

    monkeypatch.setattr(bake_clinvar, "fetch", fetch)
    monkeypatch.setattr(bake_clinvar, "fetch_batch", lambda efetch, ids: b"<batch/>")
    monkeypatch.setattr(bake_clinvar, "bake", bake)
    cache = tmp_path / "fetching" / "clinvar"
    assert worker.clinvar_in_process(target, cache) == b"the track\n"

    # Never beside the twenty's raw responses, and gone once the bake is over.
    assert (seen["cache"], seen["replay"], seen["batch"]) == (cache / "INS", False, b"<batch/>")
    assert not (cache / "INS").exists()
    assert bake_clinvar.fetch is fetch
    assert capsys.readouterr() == ("", "")

    # A stop is heard before the next batch, and NCBI's own failure is no verdict.
    with pytest.raises(worker.Unfetched, match="stopped while ClinVar"):
        worker.clinvar_in_process(target, cache, stopped=lambda: True)
    monkeypatch.setattr(bake_clinvar, "fetch", _raising(ValueError(
        "Search changed during pagination; retry")))
    with pytest.raises(worker.Unfetched, match="did not answer in full: Search changed"):
        worker.clinvar_in_process(target, cache)
    assert not isinstance(worker.Unfetched("x"), ValueError)
