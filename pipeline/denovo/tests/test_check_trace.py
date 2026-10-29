from __future__ import annotations

import base64
import json
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
    data = load()
    data["schema"] = 2
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
