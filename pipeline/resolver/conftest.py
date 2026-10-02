"""What the resolver's tests share: a builder that cannot reach the network.

The record is the real NG_007114 under `tests/fixtures/`, put where the
builder's own flat-file cache looks; UniProt is whichever entry body a test
serves; Entrez raises if anything calls it.
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
