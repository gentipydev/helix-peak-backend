from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[3]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

torch = pytest.importorskip("torch", reason="pipeline/denovo/venv has torch; .venv does not")

from pipeline.denovo.alphabet import BOS_ID, EOS_ID, PAD_ID
from pipeline.denovo.helix_peek_LM import HelixPeekLM
from pipeline.denovo.sequences import build_dataset
from pipeline.denovo.train import (
    Config, _token_losses, evaluate, lr_schedule, param_groups, run,
)

CPU = torch.device("cpu")
SMALL = Config(d_model=32, nhead=2, num_layers=2, dim_feedforward=64,
               batch_size=4, max_epochs=2, patience=2, warmup_steps=2)
TOY = ["A" * 12, "C" * 14, "D" * 16, "E" * 18, "F" * 20, "G" * 22, "H" * 24, "I" * 26]


def test_only_matrices_are_decayed():
    model = HelixPeekLM(d_model=32, nhead=2, num_layers=2, dim_feedforward=64)
    decay, plain = param_groups(model, 0.1)
    assert decay["weight_decay"] == 0.1
    assert plain["weight_decay"] == 0.0
    decayed = {id(p) for p in decay["params"]}
    assert id(model.token_embedding.weight) not in decayed
    assert id(model.positional_embedding) not in decayed
    assert id(model.lm_head.bias) not in decayed
    assert id(model.lm_head.weight) in decayed


def test_the_schedule_warms_up_then_anneals_to_eta_min():
    config = Config(learning_rate=1e-3, eta_min=1e-5, warmup_steps=10)
    schedule = lr_schedule(config, total_steps=100)
    assert schedule(0) == pytest.approx(0.1)
    assert schedule(9) == pytest.approx(1.0)
    assert schedule(9) > schedule(50) > schedule(99)
    assert schedule(100) == pytest.approx(config.eta_min / config.learning_rate)


def test_a_padded_position_costs_nothing_and_stays_finite():
    logits = torch.randn(1, 4, 23, requires_grad=True)
    targets = torch.tensor([[5, EOS_ID, PAD_ID, PAD_ID]])
    losses = _token_losses(logits, targets)
    assert losses[2].item() == 0.0 and losses[3].item() == 0.0
    assert losses[0].item() > 0.0
    losses.sum().backward()
    assert torch.isfinite(logits.grad).all()


def test_no_probability_is_spent_on_what_cannot_come_next():
    logits = torch.zeros(1, 1, 23)
    targets = torch.tensor([[5]])
    assert _token_losses(logits, targets)[0].item() == pytest.approx(math.log(21))


def test_the_blocked_columns_never_learn():
    logits = torch.randn(1, 3, 23, requires_grad=True)
    _token_losses(logits, torch.tensor([[5, 6, EOS_ID]])).sum().backward()
    assert logits.grad[0, :, [PAD_ID, BOS_ID]].abs().max().item() == 0.0


def test_evaluate_scores_residues_and_the_ending_apart():
    model = HelixPeekLM(d_model=32, nhead=2, num_layers=2, dim_feedforward=64)
    frames, lengths = build_dataset(TOY)
    scores = evaluate(model, frames, lengths, batch_size=4, device=CPU)
    assert set(scores) == {"loss", "residue_loss", "eos_loss"}
    assert all(math.isfinite(v) for v in scores.values())
    assert min(scores["residue_loss"], scores["eos_loss"]) <= scores["loss"]
    assert scores["loss"] <= max(scores["residue_loss"], scores["eos_loss"])


def test_a_short_run_learns_something_and_keeps_its_best_weights():
    result = run(SMALL, TOY, TOY[:4], device=CPU, log=False)
    assert result["best"]["loss"] < math.log(21)
    assert result["state_dict"] is not None
    assert 1 <= result["best_epoch"] <= SMALL.max_epochs
    assert result["history"][0]["epoch"] == 1


def test_a_run_is_repeatable_from_its_seed():
    first = run(SMALL, TOY, TOY[:4], device=CPU, log=False)["best"]["loss"]
    second = run(SMALL, TOY, TOY[:4], device=CPU, log=False)["best"]["loss"]
    assert first == pytest.approx(second, abs=1e-9)
