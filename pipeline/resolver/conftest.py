"""What the resolver's tests share: a builder that cannot reach the network,
and no locus inputs snapshot unless a test makes one.

The record is the real NG_007114 under `tests/fixtures/`, put where the
builder's own flat-file cache looks; UniProt is whichever entry body a test
serves; Entrez raises if anything calls it. The pinned snapshot the locus is
made from (`pipeline/locus/inputs.py`) is looked for where there is none, so a
protein resolved here has no locus track, whether or not the machine has
built the real one; a test of the locus step builds a small one
(`pipeline/locus/testing.py`).
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.mock import build_gene_record as builder  # noqa: E402
from pipeline.uniprot import Entry  # noqa: E402

RECORD = BACKEND / "tests" / "fixtures" / "ng_007114.gb"


@pytest.fixture(autouse=True)
def no_locus_inputs(monkeypatch, tmp_path):
    """The locus inputs snapshot looked for where there is none."""
    from pipeline.locus import inputs

    monkeypatch.setenv(inputs.ENV, str(tmp_path / "no-locus-inputs.jsonl.gz"))
    monkeypatch.setattr(inputs, "_LOADED", {})


@pytest.fixture
def offline(monkeypatch, tmp_path):
    """Serve an entry body as UniProt P01308, with NG_007114 in the cache."""
    cache = tmp_path / "genbank"
    cache.mkdir()
    shutil.copy(RECORD, cache / "NG_007114.1.gb")
    monkeypatch.setattr(builder, "CACHE", cache)

    def no_entrez(**kwargs):
        raise AssertionError(f"Entrez was called: {kwargs}")

    monkeypatch.setattr(builder.Entrez, "efetch", no_entrez)

    def serve(body: dict) -> Entry:
        monkeypatch.setattr(builder, "canonical_sequence",
                            lambda accession: body["sequence"]["value"])
        return Entry(accession="P01308", body=body, release="2026_03",
                     release_date="2026-09-02", retrieved="2026-10-02")

    return serve
