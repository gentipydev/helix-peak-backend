from __future__ import annotations

import base64
import json
import struct
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[3]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.denovo import check_trace
from pipeline.denovo.check_trace import Refused, validate

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "toy_trace.json"


def load() -> dict:
    return json.loads(FIXTURE.read_bytes())


def dump(data: dict) -> bytes:
    return json.dumps(data).encode()


def rewrite(spec: dict, change) -> None:
    raw = bytearray(base64.b64decode(spec["data"]))
    change(raw)
    spec["data"] = base64.b64encode(bytes(raw)).decode("ascii")


def test_the_fixture_is_a_run_the_app_can_play():
    summary = validate(FIXTURE.read_bytes(), name="toy")
    assert summary["name"] == "toy"
    assert summary["peptides"] == 2
    assert summary["targets"] == 25
    assert summary["checkpoints"][0] == 1 and summary["checkpoints"][-1] == summary["steps"]
    assert summary["loss"][0] > summary["loss"][1] >= summary["floor"]


def test_another_run_is_not_this_one():
    with pytest.raises(Refused, match="'toy', not 'curated'"):
        validate(FIXTURE.read_bytes(), name="curated")


def test_a_schema_this_does_not_read_is_refused():
    for schema in (1, 2, 3, 4, 5):
        data = load()
        data["schema"] = schema
        with pytest.raises(Refused, match="schema"):
            validate(dump(data))


def test_a_vocabulary_not_the_models_is_refused():
    data = load()
    data["vocab"] = data["vocab"][::-1]
    with pytest.raises(Refused, match="vocabulary"):
        validate(dump(data))


def test_a_cut_short_array_is_refused():
    data = load()
    rewrite(data["steps"]["loss"], lambda raw: raw.__delitem__(slice(-4, None)))
    with pytest.raises(Refused, match="loss holds"):
        validate(dump(data))


def test_a_loss_under_the_floor_is_refused():
    data = load()
    data["floor"] = 10.0
    with pytest.raises(Refused, match="floor"):
        validate(dump(data))


def test_targets_must_average_to_their_steps_loss():
    data = load()
    rewrite(data["steps"]["token_loss"], lambda raw: raw.__setitem__(slice(0, 8), bytes(8)))
    with pytest.raises(Refused, match="targets average"):
        validate(dump(data))


def test_an_attention_row_must_sum_to_one():
    data = load()
    # The first row of the first map is one cell wide: <bos> reads only itself.
    rewrite(data["checkpoints"][0]["attention"], lambda raw: raw.__setitem__(0, 0))
    with pytest.raises(Refused, match="attention has a row"):
        validate(dump(data))


def test_the_groups_must_hold_every_weight():
    data = load()
    data["group_sizes"][0] += 1
    with pytest.raises(Refused, match="sizes"):
        validate(dump(data))


def test_the_batch_must_frame_the_peptides():
    data = load()
    data["batch"]["x"][0][1] = 0
    with pytest.raises(Refused, match="does not frame"):
        validate(dump(data))


def test_a_checkpoint_without_its_embedding_tables_is_refused():
    data = load()
    del data["checkpoints"][0]["position_embedding"]
    with pytest.raises(Refused, match="position_embedding"):
        validate(dump(data))


def test_the_padding_row_of_the_token_embedding_is_zero():
    data = load()
    rewrite(data["checkpoints"][0]["token_embedding"], lambda raw: raw.__setitem__(0, 5))
    with pytest.raises(Refused, match="padding row"):
        validate(dump(data))


def test_the_tables_must_be_kept_whole():
    data = load()
    spec = data["checkpoints"][0]["position_embedding"]
    raw = base64.b64decode(spec["data"])
    values = struct.unpack("<%df" % (len(raw) // 4), raw)
    scale = max(abs(v) for v in values) / 127
    quantized = struct.pack("<%db" % len(values), *(round(v / scale) for v in values))
    spec.update(dtype="i1", scale=scale, data=base64.b64encode(quantized).decode("ascii"))
    with pytest.raises(Refused, match="not kept whole"):
        validate(dump(data))


def test_x0_must_be_each_ids_token_row_plus_its_places():
    data = load()
    width = data["model"]["d_model"]
    # <bos>'s row, which every peptide's first place reads: its first number's
    # top byte, so the float moves far.
    rewrite(data["checkpoints"][0]["token_embedding"],
            lambda raw: raw.__setitem__(4 * width + 3, raw[4 * width + 3] ^ 0x40))
    with pytest.raises(Refused, match="token row plus"):
        validate(dump(data))


def test_the_checkpoints_must_run_from_the_first_step_to_the_last():
    data = load()
    data["checkpoints"] = data["checkpoints"][1:]
    with pytest.raises(Refused, match="checkpoints"):
        validate(dump(data))


def test_a_payload_over_the_limit_is_refused(monkeypatch):
    monkeypatch.setattr(check_trace, "LIMIT", 1000)
    with pytest.raises(Refused, match="over the 1000"):
        validate(FIXTURE.read_bytes())


def test_main_says_why_it_refused(tmp_path, capsys):
    bad = tmp_path / "trace.json"
    bad.write_bytes(b"{}")
    assert check_trace.main([str(bad)]) == 1
    assert "REFUSED" in capsys.readouterr().err
    assert check_trace.main([str(FIXTURE)]) == 0


def test_the_stream_must_be_kept_whole():
    data = load()
    spec = data["checkpoints"][0]["points"]
    raw = base64.b64decode(spec["data"])
    values = struct.unpack("<%df" % (len(raw) // 4), raw)
    scale = max(abs(v) for v in values) / 127
    quantized = struct.pack("<%db" % len(values), *(round(v / scale) for v in values))
    spec.update(dtype="i1", scales=[scale] * spec["shape"][0], data=base64.b64encode(quantized).decode("ascii"))
    with pytest.raises(Refused, match="points is not kept whole"):
        validate(dump(data))


def test_a_checkpoint_without_its_blocks_is_refused():
    data = load()
    del data["checkpoints"][0]["blocks"]
    with pytest.raises(Refused, match="blocks"):
        validate(dump(data))


def test_a_row_that_is_not_the_stream_evened_out_is_refused():
    data = load()
    rewrite(data["checkpoints"][-1]["blocks"][0]["attn_input"], lambda raw: raw.__setitem__(2, raw[2] ^ 0x40))
    with pytest.raises(Refused, match="not the stream evened out"):
        validate(dump(data))


def test_a_layernorm_scale_that_is_not_the_one_used_is_refused():
    data = load()
    rewrite(data["checkpoints"][-1]["blocks"][1]["norm1_weight"], lambda raw: raw.__setitem__(2, raw[2] ^ 0x40))
    with pytest.raises(Refused, match="not the stream evened out"):
        validate(dump(data))


def test_queries_keys_and_values_must_be_the_weights_sums_of_the_evened_row():
    for key, named in (("attn_in", "query"), ("key", "key"), ("value", "value")):
        data = load()
        # A mantissa bit of the first number of a query's weights, or of the
        # first key or value: a change of up to a half, never a NaN.
        rewrite(data["checkpoints"][-1]["blocks"][0][key], lambda raw: raw.__setitem__(2, raw[2] ^ 0x40))
        with pytest.raises(Refused, match="%s is not its weights' sum" % named):
            validate(dump(data))


def test_scores_that_are_not_query_times_key_are_refused():
    data = load()
    # The second place's two scores (after <bos>'s one): swap their bytes.
    spec = data["checkpoints"][-1]["blocks"][0]["scores"]
    rewrite(spec, lambda raw: raw.__setitem__(slice(4, 12), raw[8:12] + raw[4:8]))
    with pytest.raises(Refused, match="query times key"):
        validate(dump(data))


def test_a_map_that_is_not_softmax_of_the_scores_is_refused():
    data = load()
    # Block 1, head 1, the second place's row of two shares: swapped, it still
    # sums to one.
    spec = data["checkpoints"][-1]["attention"]
    raw = base64.b64decode(spec["data"])
    assert raw[1] != raw[2]
    rewrite(spec, lambda raw: raw.__setitem__(slice(1, 3), bytes([raw[2], raw[1]])))
    with pytest.raises(Refused, match="do not softmax"):
        validate(dump(data))


def test_a_mix_that_is_not_the_shares_of_the_values_is_refused():
    data = load()
    # A mantissa bit of the first number of the mix: a change of up to a half.
    rewrite(data["checkpoints"][-1]["blocks"][1]["mixed"], lambda raw: raw.__setitem__(2, raw[2] ^ 0x40))
    with pytest.raises(Refused, match="shares of the values"):
        validate(dump(data))


def quantized(spec: dict) -> None:
    """`spec` stored as int8 with one scale, as schema 4 kept the mix."""
    raw = base64.b64decode(spec["data"])
    values = struct.unpack("<%df" % (len(raw) // 4), raw)
    scale = max(abs(v) for v in values) / 127
    packed = struct.pack("<%db" % len(values), *(round(v / scale) for v in values))
    spec.update(dtype="i1", scale=scale, data=base64.b64encode(packed).decode("ascii"))


def test_the_mix_must_be_kept_whole():
    data = load()
    quantized(data["checkpoints"][0]["blocks"][0]["mixed"])
    with pytest.raises(Refused, match="mixed is not kept whole"):
        validate(dump(data))


def test_a_block_without_its_output_weights_is_refused():
    for key in ("attn_out", "attn_out_bias"):
        data = load()
        del data["checkpoints"][0]["blocks"][1][key]
        with pytest.raises(Refused, match="%s is not kept whole" % key):
            validate(dump(data))


def test_the_output_weights_must_be_kept_whole():
    data = load()
    quantized(data["checkpoints"][-1]["blocks"][0]["attn_out"])
    with pytest.raises(Refused, match="attn_out is not kept whole"):
        validate(dump(data))


def test_output_weights_that_do_not_add_what_the_stream_gained_are_refused():
    for b, key in ((0, "attn_out"), (1, "attn_out_bias")):
        data = load()
        # A mantissa bit of the first weight, or of the first bias (trained by
        # the last step, so not zero): a change of up to a half.
        rewrite(data["checkpoints"][-1]["blocks"][b][key], lambda raw: raw.__setitem__(2, raw[2] ^ 0x40))
        with pytest.raises(Refused, match="attention does not add what the stream gained"):
            validate(dump(data))


def test_neurons_that_are_not_gelu_of_their_input_are_refused():
    data = load()
    rewrite(data["checkpoints"][-1]["blocks"][0]["hidden"],
            lambda raw: raw.__setitem__(3, (raw[3] + 40) % 128))
    with pytest.raises(Refused, match="GELU"):
        validate(dump(data))


def test_a_feed_forward_that_does_not_add_what_the_stream_gained_is_refused():
    data = load()
    rewrite(data["checkpoints"][-1]["blocks"][1]["ff_out_bias"],
            lambda raw: raw.__setitem__(3, raw[3] ^ 0x40))
    with pytest.raises(Refused, match="what the stream gained"):
        validate(dump(data))


def test_a_readout_row_must_sum_to_one():
    data = load()
    rewrite(data["checkpoints"][0]["readout"], lambda raw: raw.__setitem__(slice(0, 21), bytes(21)))
    with pytest.raises(Refused, match="readout has a row"):
        validate(dump(data))


def test_the_last_readout_must_be_the_models_probabilities():
    data = load()
    spec = data["checkpoints"][-1]["readout"]
    targets, allowed = spec["shape"][1], spec["shape"][2]
    last = (spec["shape"][0] - 1) * targets * allowed
    # The last point's first row, its two largest shares swapped: still sums to one.

    def swap(raw):
        row = list(raw[last:last + allowed])
        first = max(range(allowed), key=row.__getitem__)
        second = max((i for i in range(allowed) if i != first), key=row.__getitem__)
        row[first], row[second] = row[second], row[first]
        raw[last:last + allowed] = bytes(row)

    rewrite(spec, swap)
    with pytest.raises(Refused, match="not the model's probabilities"):
        validate(dump(data))


def test_a_checkpoint_without_its_output_layer_is_refused():
    data = load()
    del data["checkpoints"][0]["head_weight"]
    with pytest.raises(Refused, match="head_weight is not kept whole"):
        validate(dump(data))


def test_the_output_layer_and_its_logits_must_be_kept_whole():
    for key in ("head_weight", "logits"):
        data = load()
        quantized(data["checkpoints"][-1][key])
        with pytest.raises(Refused, match="%s is not kept whole" % key):
            validate(dump(data))


def test_an_output_layer_not_the_one_the_step_found_is_refused():
    data = load()
    # A mantissa bit of the first weight: a change of up to a half, far past
    # the half step the step's own int8 copy of it allows.
    rewrite(data["checkpoints"][-1]["head_weight"], lambda raw: raw.__setitem__(2, raw[2] ^ 0x40))
    with pytest.raises(Refused, match="not the output layer the step found"):
        validate(dump(data))


def test_logits_that_are_not_the_last_row_times_the_output_layer_are_refused():
    data = load()
    rewrite(data["checkpoints"][-1]["logits"], lambda raw: raw.__setitem__(2, raw[2] ^ 0x40))
    with pytest.raises(Refused, match="logits are not the last row times the output layer"):
        validate(dump(data))


def test_probabilities_that_are_not_softmax_of_the_logits_are_refused():
    data = load()
    checkpoint = data["checkpoints"][-1]
    allowed = checkpoint["probs"]["shape"][1]
    targets = checkpoint["probs"]["shape"][0]
    # The first row's two largest shares swapped, in the probabilities and in
    # the last readout alike: each still sums to one and the two still agree,
    # but neither is softmax of the logits any more.

    def swap(at):
        def change(raw):
            row = list(raw[at:at + allowed])
            top = max(range(allowed), key=row.__getitem__)
            second = max((i for i in range(allowed) if i != top), key=row.__getitem__)
            row[top], row[second] = row[second], row[top]
            raw[at:at + allowed] = bytes(row)
        return change

    rewrite(checkpoint["probs"], swap(0))
    readout = checkpoint["readout"]
    rewrite(readout, swap((readout["shape"][0] - 1) * targets * allowed))
    with pytest.raises(Refused, match="probabilities are not softmax of the logits"):
        validate(dump(data))
