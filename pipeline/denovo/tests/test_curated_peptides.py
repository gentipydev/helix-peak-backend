from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[3]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.denovo.alphabet import MAX_LENGTH, MIN_LENGTH
from pipeline.denovo.curated_peptides import (
    LEADER, MATURE, by_protein, collect, entries, leaders, missing, modelable, of,
    read_record, record_path, report, under_floor, withheld, write,
)
from pipeline.targets import TARGETS

stored = pytest.mark.skipif(
    bool(missing()),
    reason="run pipeline/fetch_tracks.py --kind record to fill pipeline/data",
)


def _target(slug):
    return next(t for t in TARGETS if t.slug == slug)


def test_modelable_uses_the_model_s_own_bounds():
    assert not modelable("A" * (MIN_LENGTH - 1))
    assert modelable("A" * MIN_LENGTH)
    assert modelable("A" * MAX_LENGTH)
    assert not modelable("A" * (MAX_LENGTH + 1))


def test_modelable_refuses_a_residue_the_tokenizer_cannot_read():
    assert not modelable("ACDEFGHIKX")


def test_record_path_is_under_the_data_directory():
    assert record_path(_target("insulin")).name == "gene_ins.json"


def test_entries_reads_the_peptides_and_the_pieces_beside_them():
    payload = {
        "peptides": [{"product": "B chain", "translation": "AAAAAAAAAA"}],
        "signal_peptide": {"product": None, "translation": "CCCCCCCCCC"},
        "proprotein": {"product": "pro", "translation": "DDDDDDDDDD"},
    }
    found = entries(payload)
    assert [e[0] for e in found] == [MATURE, "signal_peptide", "proprotein"]
    assert found[1][1] == "signal_peptide"


def test_entries_of_a_record_with_nothing_is_empty():
    assert entries({"peptides": [], "protein": None}) == []


@stored
def test_insulin_carries_its_three_chains():
    found = {label: s for _, label, s in of(_target("insulin"), with_leaders=True)}
    assert len(found["insulin A chain"]) == 21
    assert len(found["insulin B chain"]) == 30
    assert len(found["C-peptide"]) == 31
    assert len(found["signal_peptide"]) == 24
    assert found["insulin A chain"] == "GIVEQCCTSICSLYQLENYCN"


@stored
def test_a_leader_is_not_a_product():
    assert {e[0] for e in of(_target("insulin"))} == {MATURE}
    assert of(_target("lysozyme")) == []
    assert [e[0] for e in of(_target("lysozyme"), with_leaders=True)] == [LEADER]


@stored
def test_the_leaders_are_set_aside_not_lost():
    held = leaders()
    assert len(held) == 31
    assert all(e[0] == LEADER for _, e in held)
    assert "lysozyme" in {t.slug for t, _ in held}


@stored
def test_only_three_proteins_carry_a_peptide_the_model_can_train_on():
    assert [t.slug for t, _ in by_protein()] == ["insulin", "relaxin", "vasopressin"]
    assert len(collect()) == 6


@stored
def test_collect_is_unique_and_every_sequence_is_modelable():
    unique = collect()
    assert all(modelable(s) for s in unique)
    assert len(unique) == sum(len(found) for _, found in by_protein())


@stored
def test_a_record_with_no_peptide_is_named_against_the_table():
    held = {t.slug for t, _ in withheld()}
    assert {"p53", "dystrophin", "glucagon"} <= held


@stored
def test_the_nonapeptides_are_seen_and_refused():
    small = {t.slug: [len(s) for _, _, s in found] for t, found in under_floor()}
    assert small["oxytocin"] == [9]
    assert small["vasopressin"] == [9]


@stored
def test_write_round_trips(tmp_path):
    out = tmp_path / "peptides.txt"
    count = write(out)
    lines = out.read_text().split()
    assert count == len(lines) == len(set(lines))
    assert all(modelable(line) for line in lines)


@stored
def test_the_report_says_why_the_leaders_are_out():
    assert "31 signal peptides set aside" in report()
    assert "set aside" not in report(with_leaders=True)


@stored
def test_report_renders_both_ways():
    assert "unique peptides" in report()
    assert report(markdown=True).startswith("| protein |")


def _unfetched():
    return dataclasses.replace(_target("insulin"), gene="NOWHERE")


def test_an_unfetched_record_is_named_not_guessed_at():
    absent = _unfetched()
    assert absent.mock_asset == "assets/mock/gene_nowhere.json"
    assert missing([absent]) == [absent]
    with pytest.raises(FileNotFoundError, match="fetch_tracks.py"):
        read_record(absent)


def test_report_says_what_to_run_when_a_record_is_not_fetched():
    said = report(targets=[_unfetched()])
    assert "fetch_tracks.py --kind record" in said
    assert "1 of 1 records" in said
