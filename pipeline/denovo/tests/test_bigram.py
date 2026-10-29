from __future__ import annotations

import math
import sys
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[3]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.denovo.alphabet import BOS_ID, EOS_ID, PAD_ID, TOKEN_TO_ID
from pipeline.denovo.bigram import (
    counts, cross_entropy, probabilities, split,
)

TOY = ["MAM", "MAA"]
M, A, W = TOKEN_TO_ID["M"], TOKEN_TO_ID["A"], TOKEN_TO_ID["W"]


def test_counts_are_the_framed_pairs():
    table = counts(TOY)
    assert table[BOS_ID][M] == 2
    assert table[M][A] == 2
    assert table[A][M] == 1
    assert table[A][A] == 1
    assert table[M][EOS_ID] == 1
    assert table[A][EOS_ID] == 1
    assert sum(map(sum, table)) == sum(len(s) + 1 for s in TOY)


def test_rows_are_distributions_over_what_can_come_next():
    probs = probabilities(counts(TOY))
    assert math.isclose(sum(probs[BOS_ID]), 1.0)
    assert math.isclose(sum(probs[M]), 1.0)
    assert probs[M][PAD_ID] == 0.0
    assert probs[M][BOS_ID] == 0.0


def test_smoothing_keeps_unseen_pairs_possible():
    assert probabilities(counts(TOY), smoothing=1.0)[M][W] > 0.0
    assert probabilities(counts(TOY), smoothing=0.0)[M][W] == 0.0


def test_a_counted_pair_outscores_an_unseen_one():
    probs = probabilities(counts(TOY))
    assert probs[M][A] > probs[M][W]


def test_cross_entropy_is_lower_than_uniform():
    probs = probabilities(counts(TOY))
    assert cross_entropy(probs, TOY) < math.log(20)


def test_split_is_disjoint_and_repeatable():
    sequences = ["A" * n for n in range(1, 11)]
    train, validation = split(sequences, held_back=3)
    assert len(validation) == 3
    assert sorted(train + validation) == sorted(sequences)
    assert split(sequences, held_back=3) == (train, validation)
