from __future__ import annotations

import sys
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[3]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.denovo.alphabet import MAX_LENGTH, MIN_LENGTH
from pipeline.denovo.curated_peptides import (
    in_range, length, offered, pieces, report, sources, too_short, whole, withheld,
)
from pipeline.targets import TARGETS


def _target(slug):
    return next(t for t in TARGETS if t.slug == slug)


def test_length_is_inclusive():
    insulin = _target("insulin")
    a_chain = next(r for r in insulin.regions if r.label == "A chain")
    assert length(a_chain) == 21


def test_in_range_uses_the_model_s_own_bounds():
    assert not in_range(MIN_LENGTH - 1)
    assert in_range(MIN_LENGTH)
    assert in_range(MAX_LENGTH)
    assert not in_range(MAX_LENGTH + 1)


def test_pieces_drops_what_the_precursor_loses():
    kept = [label for label, _ in pieces(_target("insulin"))]
    assert kept == ["B chain", "A chain"]
    everything = [label for label, _ in pieces(_target("insulin"), kept_only=False)]
    assert "Signal peptide" in everything and "C-peptide" in everything


def test_offered_is_empty_without_mature_peptides():
    glucagon = _target("glucagon")
    assert not glucagon.mature_peptides
    assert pieces(glucagon)
    assert offered(glucagon) == []


def test_insulin_and_relaxin_are_the_cleaved_pairs():
    found = {t.slug: dict(p) for t, p in sources()}
    assert found["insulin"] == {"B chain": 30, "A chain": 21}
    assert found["relaxin"] == {"B chain": 29, "A chain": 24}


def test_sources_are_ordered_by_how_much_they_offer():
    counts = [len(found) for _, found in sources()]
    assert counts == sorted(counts, reverse=True)


def test_no_curated_protein_is_itself_in_range():
    assert whole() == []


def test_the_nonapeptides_are_under_the_floor():
    small = {t.slug: dict(p) for t, p in too_short()}
    assert small["oxytocin"]["Oxytocin"] == 9
    assert small["vasopressin"]["Vasopressin"] == 9


def test_glucagon_is_the_one_held_back():
    assert [t.slug for t, _ in withheld()] == ["glucagon"]


def test_report_renders_both_ways():
    assert "GeneRecord.peptides" in report()
    assert report(markdown=True).startswith("| protein |")
