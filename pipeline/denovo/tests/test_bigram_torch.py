from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[3]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

torch = pytest.importorskip("torch", reason="pipeline/denovo/venv has torch; .venv does not")

from pipeline.denovo.alphabet import BOS_ID, EOS_ID, PAD_ID, TOKEN_TO_ID
from pipeline.denovo.bigram import counts, cross_entropy, probabilities
from pipeline.denovo.bigram_torch import (
    loss_of, masked, observed_rows, pairs, train,
)

TOY = ["MAM", "MAA"]
M, A, W = TOKEN_TO_ID["M"], TOKEN_TO_ID["A"], TOKEN_TO_ID["W"]


def test_pairs_are_the_adjacent_tokens():
    first, second = pairs(["MAM"])
    assert first.tolist() == [BOS_ID, M, A, M]
    assert second.tolist() == [M, A, M, EOS_ID]
    assert len(first) == len(second)


def test_pairs_run_every_sequence_together():
    first, second = pairs(TOY)
    assert len(first) == sum(len(s) + 1 for s in TOY)


def test_masked_blocks_what_cannot_come_next():
    blocked = masked(torch.zeros(23, 23))
    assert torch.isinf(blocked[:, PAD_ID]).all()
    assert torch.isinf(blocked[:, BOS_ID]).all()
    assert torch.isfinite(blocked[:, EOS_ID]).all()
    assert torch.isfinite(blocked[:, M]).all()


def test_rows_are_distributions_after_masking():
    probs = masked(torch.randn(23, 23)).softmax(dim=1)
    assert torch.allclose(probs.sum(dim=1), torch.ones(23))
    assert probs[:, PAD_ID].sum() == 0.0
    assert probs[:, BOS_ID].sum() == 0.0


def test_training_beats_the_uniform_baseline():
    first, second = pairs(TOY)
    logits = train(first, second, steps=200, learning_rate=10.0)
    assert loss_of(logits, first, second) < math.log(20)


def test_training_learns_the_order_of_the_counts():
    first, second = pairs(TOY)
    probs = masked(train(first, second, steps=200, learning_rate=10.0)).softmax(dim=1)
    assert probs[M][A] > probs[M][W]
    assert probs[BOS_ID][M] > probs[BOS_ID][W]


def test_it_converges_on_the_counting_table():
    first, second = pairs(TOY)
    logits = train(first, second, steps=4000, learning_rate=10.0)
    counted = probabilities(counts(TOY), smoothing=0.0)
    assert loss_of(logits, first, second) == pytest.approx(
        cross_entropy(counted, TOY), abs=0.02
    )
    learned = masked(logits).softmax(dim=1)
    gap = (learned - torch.tensor(counted))[observed_rows(first)].abs().max().item()
    assert gap < 0.05
