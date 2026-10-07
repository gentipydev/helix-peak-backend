"""The worker's own decisions, offline. Its SQL is `test_worker_pg.py`'s."""

from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import replace
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.resolver import store, worker  # noqa: E402
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
