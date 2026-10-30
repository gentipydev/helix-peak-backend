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
    SHOWCASE, TOY, TOY_PEPTIDES, Peptide, Source, attention, deterministic, floor, group_names,
    group_of, inside, milestones, point_names, readout_names, record, serialize,
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


def scaled(spec):
    return stored(spec) * spec.get("scale", 1.0)


def gelu(x):
    return 0.5 * x * (1 + np.vectorize(math.erf)(x / math.sqrt(2)))


def rows_of(spans):
    """Each peptide's places' rows of its lower triangle, as (start, length)."""
    at = 0
    for span in spans:
        for i in range(span):
            yield at, i + 1
            at += i + 1


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


def test_the_embedding_tables_are_the_ones_the_step_read(toy):
    torch.manual_seed(TOY.seed)
    model = HelixPeekLM(d_model=TOY.d_model, nhead=TOY.nhead, num_layers=TOY.num_layers,
                        dim_feedforward=TOY.dim_feedforward, dropout=TOY.dropout)
    first = toy["checkpoints"][0]
    for key, weights in (("token_embedding", model.token_embedding.weight),
                         ("position_embedding", model.positional_embedding[0])):
        spec = first[key]
        assert spec["shape"] == list(weights.shape)
        # Kept whole: the weights as the step found them, float for float.
        assert spec["dtype"] == "f4"
        assert np.array_equal(stored(spec), weights.detach().numpy().astype(np.float32))


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


def test_every_weight_is_in_one_group(toy):
    model = HelixPeekLM(d_model=TOY.d_model, nhead=TOY.nhead, num_layers=TOY.num_layers,
                        dim_feedforward=TOY.dim_feedforward)
    assert sum(toy["group_sizes"]) == toy["model"]["parameters"] == sum(
        p.numel() for p in model.parameters())
    assert len(toy["group_sizes"]) == len(toy["groups"])


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


def test_the_showcase_is_two_blocks_of_two_heads_sixteen_wide():
    assert (SHOWCASE.num_layers, SHOWCASE.nhead, SHOWCASE.d_model, SHOWCASE.dim_feedforward) == (2, 2, 16, 64)
    assert SHOWCASE.dropout == 0.0


def test_the_stream_is_kept_whole_and_x0_is_the_two_tables_added(toy):
    spans = [len(s) + 1 for s in TOY_SEQUENCES]
    for checkpoint in toy["checkpoints"]:
        assert checkpoint["points"]["dtype"] == "f4"
        x0 = stored(checkpoint["points"])[0].astype(np.float32)
        token = stored(checkpoint["token_embedding"]).astype(np.float32)
        place = stored(checkpoint["position_embedding"]).astype(np.float32)
        ids = [toy["batch"]["x"][c][i] for c, span in enumerate(spans) for i in range(span)]
        places = [i for span in spans for i in range(span)]
        assert np.array_equal(x0, token[ids] + place[places])


def test_each_block_keeps_what_it_worked_out(toy):
    targets, width, ff = toy["targets"], TOY.d_model, TOY.dim_feedforward
    triangle = sum(n * (n + 1) // 2 for n in (len(s) + 1 for s in TOY_SEQUENCES))
    shapes = {"query": ("i1", [targets, width]), "key": ("i1", [targets, width]),
              "value": ("i1", [targets, width]), "mixed": ("i1", [targets, width]),
              "scores": ("f4", [TOY.nhead, triangle]), "ff_input": ("f4", [targets, width]),
              "hidden": ("i1", [targets, ff]), "ff_in": ("f4", [ff, width]),
              "ff_in_bias": ("f4", [ff]), "ff_out": ("f4", [width, ff]), "ff_out_bias": ("f4", [width])}
    for checkpoint in toy["checkpoints"]:
        assert len(checkpoint["blocks"]) == TOY.num_layers
        for block in checkpoint["blocks"]:
            assert {key: (spec["dtype"], spec["shape"]) for key, spec in block.items()} == shapes


def test_the_scores_softmax_to_the_attention_maps_and_mix_the_values(toy):
    spans = [len(s) + 1 for s in TOY_SEQUENCES]
    size = TOY.d_model // TOY.nhead
    for checkpoint in toy["checkpoints"]:
        maps = scaled(checkpoint["attention"])
        for b, block in enumerate(checkpoint["blocks"]):
            scores, value, mixed = scaled(block["scores"]), scaled(block["value"]), scaled(block["mixed"])
            starts = [sum(spans[:c]) for c in range(len(spans))]
            for h in range(TOY.nhead):
                rows = list(rows_of(spans))
                for (at, n), target in zip(rows, [starts[c] + i for c, span in enumerate(spans)
                                                  for i in range(span)]):
                    row = scores[h, at:at + n]
                    shares = np.exp(row - row.max())
                    shares /= shares.sum()
                    assert np.abs(shares - maps[b, h, at:at + n]).max() <= 0.5 / 255 + 1e-6
                    first = target - (n - 1)
                    mix = shares @ value[first:target + 1, h * size:(h + 1) * size]
                    tolerance = (block["value"]["scale"] + block["mixed"]["scale"]) / 2 + 1e-5
                    assert np.abs(mix - mixed[target, h * size:(h + 1) * size]).max() <= tolerance


def test_the_neurons_are_gelu_of_their_input_and_add_what_the_stream_gained(toy):
    for checkpoint in toy["checkpoints"]:
        points = stored(checkpoint["points"])
        for b, block in enumerate(checkpoint["blocks"]):
            neurons = gelu(scaled(block["ff_input"]) @ scaled(block["ff_in"]).T + scaled(block["ff_in_bias"]))
            assert np.abs(neurons - scaled(block["hidden"])).max() <= block["hidden"]["scale"] / 2 + 1e-5
            added = neurons @ scaled(block["ff_out"]).T + scaled(block["ff_out_bias"])
            assert np.abs(added - (points[2 * b + 2] - points[2 * b + 1])).max() < 1e-5


def test_the_feed_forward_weights_are_the_ones_the_step_read(toy):
    torch.manual_seed(TOY.seed)
    model = HelixPeekLM(d_model=TOY.d_model, nhead=TOY.nhead, num_layers=TOY.num_layers,
                        dim_feedforward=TOY.dim_feedforward, dropout=TOY.dropout)
    first = toy["checkpoints"][0]
    for block, layer in zip(first["blocks"], model.transformer_encoder.layers):
        assert np.array_equal(stored(block["ff_in"]), layer.linear1.weight.detach().numpy())
        assert np.array_equal(stored(block["ff_out_bias"]), layer.linear2.bias.detach().numpy())


def test_the_readout_is_the_output_layer_read_early(toy):
    torch.manual_seed(TOY.seed)
    model = HelixPeekLM(d_model=TOY.d_model, nhead=TOY.nhead, num_layers=TOY.num_layers,
                        dim_feedforward=TOY.dim_feedforward, dropout=TOY.dropout)
    assert readout_names(2) == ["x0", "m1", "h1", "m2", "h2"]
    for checkpoint in toy["checkpoints"]:
        read = scaled(checkpoint["readout"])
        assert read.shape == (2 * TOY.num_layers + 1, toy["targets"], 21)
        assert np.allclose(read.sum(axis=-1), 1.0, atol=21 / 255)
        # Its last point is the model's own output.
        assert np.array_equal(stored(checkpoint["readout"])[-1], stored(checkpoint["probs"]))
    # At step 1 the weights are the initial ones: x0 through the final norm
    # and the output layer, as the recorder reads it.
    first = toy["checkpoints"][0]
    x0 = torch.from_numpy(stored(first["points"])[0].astype(np.float32))
    from pipeline.denovo.trace import readout
    again = readout(model, x0.unsqueeze(0))[0].numpy()
    assert np.abs(again - scaled(first["readout"])[0]).max() <= 0.5 / 255 + 1e-6


def test_inside_is_the_layers_own_attention_taken_apart():
    torch.manual_seed(0)
    model = HelixPeekLM(d_model=16, nhead=2, num_layers=1, dim_feedforward=32, dropout=0.0)
    seen = {}

    def hook(module, args, kwargs, output):
        seen["parts"], seen["error"] = inside(module, args[0], kwargs.get("attn_mask"), output[0])

    layer = model.transformer_encoder.layers[0]
    layer.self_attn.register_forward_hook(hook, with_kwargs=True)
    model(torch.tensor([[1, 5, 9, 13, 2]]))
    parts = seen["parts"]
    assert seen["error"] < 1e-5
    for key in ("query", "key", "value", "mixed"):
        assert parts[key].shape == (1, 5, 16)
    assert parts["scores"].shape == parts["weights"].shape == (1, 2, 5, 5)
    # Softmax of the scores, under the mask, is the weights.
    lower = torch.ones(5, 5, dtype=torch.bool).tril()
    again = parts["scores"].masked_fill(~lower, float("-inf")).softmax(dim=-1)
    assert torch.allclose(again, parts["weights"])
