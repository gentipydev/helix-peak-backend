"""The worker's own decisions, offline. Its SQL is `test_worker_pg.py`'s."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.resolver import worker  # noqa: E402
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


def test_a_failed_alignment_gate_is_told_as_a_sentence():
    said = worker._refusal(ValueError(
        "Alignment gate FAILED: 48.2% before, 51.0% after. No JSON written. "
        "Inspect token offsets and how scores are filed before any UI work."))
    assert said == ("ESM-2's scores failed the check that each is filed under its own residue "
                    "(48.2% and 51.0% where more than half is needed), so they are not drawn.")


def test_any_other_refusal_is_its_own_message_on_one_line():
    assert worker._refusal(ValueError("ins asset is 1,000 bytes,\n over its budget")) == \
        "ins asset is 1,000 bytes, over its budget"
    assert len(worker._said(RuntimeError("x" * 1000))) == 300
