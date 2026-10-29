from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[3]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

pytest.importorskip("torch", reason="pipeline/denovo/venv has torch; .venv does not")

from pipeline.denovo.fetch_peptides import (
    MATURE, PROPEPTIDE, parse_fasta, products_from_entry, products_from_json,
)

ENTRY = {
    "primaryAccession": "P00000",
    "sequence": {"value": "MAAAKWKLFKKIGGGRRPQVTRR"},
    "features": [
        {"type": "Peptide", "location": {"start": {"value": 4, "modifier": "EXACT"},
                                         "end": {"value": 15, "modifier": "EXACT"}}},
        {"type": "Propeptide", "location": {"start": {"value": 16, "modifier": "EXACT"},
                                            "end": {"value": 23, "modifier": "EXACT"}}},
        {"type": "Chain", "location": {"start": {"value": 1, "modifier": "EXACT"},
                                       "end": {"value": 23, "modifier": "EXACT"}}},
    ],
}


def test_a_product_is_cut_on_inclusive_one_based_bounds():
    assert products_from_entry(ENTRY) == [
        (MATURE, "AKWKLFKKIGGG"),
        (PROPEPTIDE, "RRPQVTRR"),
    ]


def test_only_the_kinds_asked_for_are_cut():
    assert products_from_entry(ENTRY, kinds=(MATURE,)) == [(MATURE, "AKWKLFKKIGGG")]


def test_a_chain_is_not_a_peptide():
    assert all(kind != "Chain" for kind, _ in products_from_entry(ENTRY))


def test_an_unpinned_boundary_is_skipped():
    entry = dict(ENTRY, features=[
        {"type": "Peptide", "location": {"start": {"value": None, "modifier": "UNKNOWN"},
                                         "end": {"value": 9, "modifier": "EXACT"}}},
        {"type": "Peptide", "location": {"start": {"value": 1, "modifier": "EXACT"},
                                         "end": {"value": 5, "modifier": "UNSURE"}}},
    ])
    assert products_from_entry(entry) == []


def test_a_boundary_past_the_sequence_is_skipped():
    entry = dict(ENTRY, features=[
        {"type": "Peptide", "location": {"start": {"value": 4, "modifier": "EXACT"},
                                         "end": {"value": 400, "modifier": "EXACT"}}},
    ])
    assert products_from_entry(entry) == []


def test_products_from_json_reads_the_results_list():
    text = json.dumps({"results": [ENTRY, ENTRY]})
    assert len(products_from_json(text, kinds=(MATURE,))) == 2


def test_parse_fasta_joins_wrapped_lines():
    assert parse_fasta(">sp|P1|A\nMKA\nLL\n\n>sp|P2|B\nQQ\n") == ["MKALL", "QQ"]


def test_parsing_never_reaches_the_network(monkeypatch):
    import urllib.request

    def forbidden(*args, **kwargs):
        raise AssertionError("parsing must not call out")

    monkeypatch.setattr(urllib.request, "urlopen", forbidden)
    assert products_from_entry(ENTRY)
    assert parse_fasta(">a\nMKALL\n")
