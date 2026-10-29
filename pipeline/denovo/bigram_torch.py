from __future__ import annotations

from typing import Sequence, Tuple

import torch
import torch.nn.functional as F

from pipeline.denovo.alphabet import VOCAB_SIZE
from pipeline.denovo.bigram import BLOCKED, framed


def pairs(sequences: Sequence[str]) -> Tuple[torch.Tensor, torch.Tensor]:
    first, second = [], []
    for sequence in sequences:
        ids = framed(sequence)
        first.extend(ids[:-1])
        second.extend(ids[1:])
    return (
        torch.tensor(first, dtype=torch.long),
        torch.tensor(second, dtype=torch.long),
    )


def masked(logits: torch.Tensor) -> torch.Tensor:
    blocked = logits.clone()
    for token in BLOCKED:
        blocked[:, token] = float("-inf")
    return blocked


def observed_rows(first: torch.Tensor) -> torch.Tensor:
    return torch.unique(first)


def loss_of(logits: torch.Tensor, first: torch.Tensor, second: torch.Tensor) -> torch.Tensor:
    return F.cross_entropy(masked(logits)[first], second)


def train(
    first: torch.Tensor,
    second: torch.Tensor,
    steps: int = 1000,
    learning_rate: float = 10.0,
) -> torch.Tensor:
    logits = torch.zeros(VOCAB_SIZE, VOCAB_SIZE, requires_grad=True)
    optimizer = torch.optim.SGD([logits], lr=learning_rate)
    for _ in range(steps):
        loss = loss_of(logits, first, second)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    return logits.detach()
