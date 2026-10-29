from __future__ import annotations

import math
import random
from typing import Iterable, List, Sequence, Tuple

from pipeline.denovo.alphabet import (
    BOS_ID, EOS_ID, PAD_ID, TOKEN_TO_ID, VOCAB_SIZE,
)

BLOCKED = (PAD_ID, BOS_ID)
ALLOWED = tuple(i for i in range(VOCAB_SIZE) if i not in BLOCKED)


def framed(sequence: str) -> List[int]:
    return [BOS_ID] + [TOKEN_TO_ID[letter] for letter in sequence] + [EOS_ID]


def adjacent(sequences: Iterable[str]) -> Iterable[Tuple[int, int]]:
    for sequence in sequences:
        ids = framed(sequence)
        for first, second in zip(ids, ids[1:]):
            yield first, second


def counts(sequences: Sequence[str]) -> List[List[int]]:
    table = [[0] * VOCAB_SIZE for _ in range(VOCAB_SIZE)]
    for first, second in adjacent(sequences):
        table[first][second] += 1
    return table


def probabilities(table: Sequence[Sequence[int]], smoothing: float = 1.0) -> List[List[float]]:
    rows = [[0.0] * VOCAB_SIZE for _ in range(VOCAB_SIZE)]
    for i, counted in enumerate(table):
        total = sum(counted[j] for j in ALLOWED) + smoothing * len(ALLOWED)
        if total <= 0:
            continue
        for j in ALLOWED:
            rows[i][j] = (counted[j] + smoothing) / total
    return rows


def cross_entropy(probs: Sequence[Sequence[float]], sequences: Sequence[str]) -> float:
    total = 0.0
    seen = 0
    for first, second in adjacent(sequences):
        p = probs[first][second]
        total += -math.log(p) if p > 0.0 else math.inf
        seen += 1
    if seen == 0:
        raise ValueError("no pairs to score")
    return total / seen


def split(sequences: Sequence[str], held_back: int = 2, seed: int = 0) -> Tuple[List[str], List[str]]:
    shuffled = list(sequences)
    random.Random(seed).shuffle(shuffled)
    return shuffled[held_back:], shuffled[:held_back]
