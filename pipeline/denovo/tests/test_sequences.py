from __future__ import annotations

import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[3]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

torch = pytest.importorskip("torch", reason="pipeline/denovo/venv has torch; .venv does not")

from pipeline.denovo.alphabet import BOS_ID, EOS_ID, MAX_LENGTH, MIN_LENGTH, PAD_ID
from pipeline.denovo.sequences import (
    FRAME_WIDTH, batches, build_dataset, clean, is_modelable, read_fasta, read_peptides,
)

TOY = ["A" * 12, "C" * 20, "D" * 31]


def test_only_the_twenty_within_range_are_modelable():
    assert is_modelable("A" * MIN_LENGTH)
    assert is_modelable("A" * MAX_LENGTH)
    assert not is_modelable("A" * (MIN_LENGTH - 1))
    assert not is_modelable("A" * (MAX_LENGTH + 1))
    for odd in "XUBZO":
        assert not is_modelable("A" * 11 + odd)


def test_clean_drops_duplicates_and_keeps_order():
    assert clean(["A" * 12, "a" * 12, "C" * 20]) == ["A" * 12, "C" * 20]


def test_clean_drops_what_cannot_be_tokenised():
    assert clean(["A" * 12, "X" * 12, "short"]) == ["A" * 12]


def test_dataset_frames_every_sequence_once():
    frames, lengths = build_dataset(TOY)
    assert frames.shape == (3, FRAME_WIDTH)
    assert lengths.tolist() == [14, 22, 33]
    for row, sequence in zip(frames, TOY):
        assert row[0].item() == BOS_ID
        assert row[len(sequence) + 1].item() == EOS_ID
        assert (row[len(sequence) + 2 :] == PAD_ID).all()


def test_batches_cover_every_sequence_exactly_once():
    frames, lengths = build_dataset(TOY * 4)
    seen = 0
    for x, y in batches(frames, lengths, batch_size=5):
        assert x.shape == y.shape
        seen += x.shape[0]
    assert seen == 12


def test_a_batch_is_trimmed_to_its_own_longest_frame():
    frames, lengths = build_dataset(TOY)
    widths = sorted(x.shape[1] for x, _ in batches(frames, lengths, batch_size=1))
    assert widths == [13, 21, 32]


def test_targets_are_the_inputs_shifted_by_one():
    frames, lengths = build_dataset(["A" * 12])
    (x, y), = list(batches(frames, lengths, batch_size=1))
    assert x[0, 1:].tolist() == y[0, :-1].tolist()


def test_bucketing_puts_like_lengths_together():
    frames, lengths = build_dataset(["A" * 10, "C" * 50, "D" * 11, "E" * 49])
    spans = [x.shape[1] for x, _ in batches(frames, lengths, batch_size=2)]
    assert sorted(spans) == [12, 51]


def test_read_fasta_and_read_peptides(tmp_path):
    fasta = tmp_path / "p.fasta"
    fasta.write_text(">one\nMKA\nLL\n>two\nQQ\n")
    assert read_fasta(fasta) == ["MKALL", "QQ"]
    plain = tmp_path / "p.txt"
    plain.write_text("# a comment\nMKALL\n\nQQ\n")
    assert read_peptides(plain) == ["MKALL", "QQ"]
