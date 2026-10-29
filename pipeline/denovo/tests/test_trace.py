from __future__ import annotations

import base64
import math
import sys
from dataclasses import replace
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[3]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

torch = pytest.importorskip("torch", reason="pipeline/denovo/venv has torch; .venv does not")
np = pytest.importorskip("numpy")

from pipeline.denovo import check_trace
from pipeline.denovo.helix_peek_LM import HelixPeekLM
from pipeline.denovo.trace import (
    TOY, TOY_PEPTIDES, Peptide, Source, attention, deterministic, floor, group_names,
    group_of, milestones, point_names, record, serialize,
)
from pipeline.denovo.train import run

CPU = torch.device("cpu")
TOY_SEQUENCES = sorted((p.sequence for p in TOY_PEPTIDES), key=lambda s: (len(s), s))


@pytest.fixture(scope="module")
def toy():
    return record(TOY, TOY_PEPTIDES, name="toy")


def stored(spec):
    raw = base64.b64decode(spec["data"])
    return np.frombuffer(raw, dtype="<" + spec["dtype"]).reshape(spec["shape"]).astype(np.float64)


def test_the_recorder_trains_as_the_engine_does(toy):
    config = replace(TOY, patience=TOY.max_epochs)
    with deterministic():
        history = run(config, TOY_SEQUENCES, TOY_SEQUENCES, device=CPU, log=False)["history"]
    recorded = stored(toy["steps"]["loss"]).tolist()
    assert len(history) == len(recorded) == TOY.max_epochs
    assert [h["train_loss"] for h in history] == pytest.approx(recorded, rel=1e-6)


def test_two_records_are_the_same_bytes(toy):
    assert serialize(record(TOY, TOY_PEPTIDES, name="toy")) == serialize(toy)


def test_what_it_records_is_a_run_the_app_can_play(toy):
    summary = check_trace.validate(serialize(toy), name="toy")
    assert summary["steps"] == TOY.max_epochs
    assert summary["checkpoints"][0] == 1 and summary["checkpoints"][-1] == TOY.max_epochs
    assert summary["targets"] == sum(len(s) + 1 for s in TOY_SEQUENCES)


def test_the_committed_fixture_is_a_toy_run():
    fixture = Path(__file__).resolve().parent / "fixtures" / "toy_trace.json"
    assert check_trace.validate(fixture.read_bytes(), name="toy")["steps"] == TOY.max_epochs


def test_the_peptides_are_kept_shortest_first_with_their_sources(toy):
    assert [p["sequence"] for p in toy["peptides"]] == TOY_SEQUENCES
    assert toy["peptides"][0]["sources"] == [{"slug": "toy", "gene": "TOY", "label": "toy A"}]


def test_a_checkpoint_falls_wherever_the_loss_first_passes_a_level():
    losses = [3.5, 3.1, 2.9, 2.5, 0.9, 0.05]
    reference = {"uniform": 3.0445, "unigram": 2.95, "bigram": 2.6}
    assert milestones(losses, reference) == [
        (1, ["first"]),
        (3, ["uniform", "unigram"]),
        (4, ["bigram"]),
        (5, ["one_nat"]),
        (6, ["tenth_nat", "last"]),
    ]


def test_a_level_never_passed_marks_nothing():
    assert milestones([3.5, 3.4], {"uniform": 3.0, "unigram": 2.9, "bigram": 2.0}) == [
        (1, ["first"]), (2, ["last"])]


def test_the_floor_is_what_no_model_can_know_from_the_prefix():
    assert floor(["A" * 10, "C" * 10]) == pytest.approx(2 * math.log(2) / 22)
    assert floor(["AC" * 5, "AD" * 5]) == pytest.approx(2 * math.log(2) / 22)
    assert floor(["A" * 10]) == 0.0


def test_no_step_is_below_the_floor(toy):
    assert stored(toy["steps"]["loss"]).min() >= toy["floor"]


def test_the_attention_worked_out_again_is_the_layers_own():
    torch.manual_seed(0)
    model = HelixPeekLM(d_model=16, nhead=2, num_layers=1, dim_feedforward=32, dropout=0.0)
    seen = {}

    def hook(module, args, kwargs, output):
        seen["weights"], seen["error"] = attention(module, args[0], kwargs.get("attn_mask"), output[0])

    layer = model.transformer_encoder.layers[0]
    layer.self_attn.register_forward_hook(hook, with_kwargs=True)
    model(torch.tensor([[1, 5, 9, 13, 2]]))
    weights = seen["weights"][0]
    assert seen["error"] < 1e-5
    assert torch.allclose(weights.sum(dim=-1), torch.ones(2, 5))
    assert weights.triu(diagonal=1).abs().max().item() == 0.0


def test_each_attention_row_reads_only_what_came_before_and_sums_to_one(toy):
    spans = [len(s) + 1 for s in TOY_SEQUENCES]
    for checkpoint in toy["checkpoints"]:
        maps = stored(checkpoint["attention"]) * checkpoint["attention"]["scale"]
        assert maps.shape[-1] == sum(n * (n + 1) // 2 for n in spans)
        at = 0
        for span in spans:
            for i in range(span):
                row = maps[..., at:at + i + 1]
                assert np.allclose(row.sum(axis=-1), 1.0, atol=(i + 1) / 255 + 1e-3)
                at += i + 1


def test_the_update_is_how_far_the_head_moved(toy):
    weights = stored(toy["steps"]["head_weight"])
    scales = toy["steps"]["head_weight"]["scales"]
    for checkpoint in toy["checkpoints"]:
        step = checkpoint["step"]
        if step == toy["steps"]["count"]:
            continue
        moved = weights[step] * scales[step] - weights[step - 1] * scales[step - 1]
        update = stored(checkpoint["head_update"]) * checkpoint["head_update"]["scale"]
        tolerance = scales[step] + scales[step - 1] + checkpoint["head_update"]["scale"]
        assert np.abs(moved - update).max() <= tolerance


def test_the_probabilities_are_over_what_can_come_next(toy):
    for checkpoint in toy["checkpoints"]:
        probs = stored(checkpoint["probs"]) * checkpoint["probs"]["scale"]
        assert probs.shape[1] == len(toy["allowed"]) == 21
        assert np.allclose(probs.sum(axis=1), 1.0, atol=21 / 255)


def test_points_and_groups_follow_the_layers():
    assert point_names(2) == ["x0", "m1", "h1", "m2", "h2", "hf"]
    assert len(group_names(4)) == 16
    model = HelixPeekLM(d_model=16, nhead=2, num_layers=3, dim_feedforward=32)
    assert {group_of(name) for name, _ in model.named_parameters()} == set(group_names(3))


def test_a_run_that_needs_two_batches_a_step_is_refused():
    with pytest.raises(ValueError, match="one batch"):
        record(replace(TOY, batch_size=1), TOY_PEPTIDES, name="toy")


def test_a_peptide_the_model_cannot_take_is_refused():
    short = Peptide("MKA", (Source("toy", "TOY", "short"),))
    with pytest.raises(ValueError, match="not a peptide"):
        record(TOY, TOY_PEPTIDES + (short,), name="toy")


def test_the_settings_a_record_needs_are_put_back():
    threads = torch.get_num_threads()
    with deterministic():
        assert torch.are_deterministic_algorithms_enabled()
        assert torch.get_num_threads() == 1
    assert torch.get_num_threads() == threads
