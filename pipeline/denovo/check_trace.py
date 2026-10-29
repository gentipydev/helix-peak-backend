"""Hold a recorded training run (`trace.py`) to its own shape before it is stored.

    python3 pipeline/denovo/check_trace.py [PATH]

Standard library only, and Python 3.9, so the uploader's environment, which
has no torch, can read what the recorder wrote. A payload that fails here is
never published: the app parses what this accepts.
"""

from __future__ import annotations

import base64
import binascii
import json
import math
import struct
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

HERE = Path(__file__).resolve().parent
if str(HERE.parents[1]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[1]))

from pipeline.denovo.alphabet import AMINO_ACIDS, MAX_LENGTH, MIN_LENGTH, PAD_ID, TOKENS
from pipeline.denovo.bigram import ALLOWED, framed

SCHEMA = 1
LIMIT = 4_000_000
DEFAULT = HERE / "runs" / "curated" / "trace.json"

_FORMATS = {"i1": ("b", 1), "u1": ("B", 1), "u2": ("H", 2), "f4": ("f", 4)}

# A row of quantised probabilities sums to 255 give or take half a step a cell.
_ROUNDING = 0.5


class Refused(ValueError):
    pass


def numbers(spec: Dict, shape: Sequence[int], what: str) -> List[float]:
    """The values `spec` holds, as stored (before any scale), checked against `shape`."""
    if not isinstance(spec, dict) or spec.get("dtype") not in _FORMATS:
        raise Refused("%s is not an array this reads" % what)
    if list(spec.get("shape") or []) != list(shape):
        raise Refused("%s is shaped %s, not %s" % (what, spec.get("shape"), list(shape)))
    letter, size = _FORMATS[spec["dtype"]]
    try:
        raw = base64.b64decode(spec.get("data") or "", validate=True)
    except (binascii.Error, ValueError):
        raise Refused("%s is not base64" % what) from None
    count = math.prod(shape)
    if len(raw) != count * size:
        raise Refused("%s holds %d bytes, not %d" % (what, len(raw), count * size))
    return list(struct.unpack("<%d%s" % (count, letter), raw))


def _finite(values: Sequence[float], what: str, low: Optional[float] = None) -> None:
    for value in values:
        if not math.isfinite(value) or (low is not None and value < low):
            raise Refused("%s holds %r" % (what, value))


def _scale(spec: Dict, what: str) -> float:
    scale = spec.get("scale")
    if not isinstance(scale, (int, float)) or not math.isfinite(scale) or scale <= 0:
        raise Refused("%s has no usable scale" % what)
    return float(scale)


def _scales(spec: Dict, count: int, what: str) -> List[float]:
    scales = spec.get("scales")
    if not isinstance(scales, list) or len(scales) != count:
        raise Refused("%s needs %d scales" % (what, count))
    for scale in scales:
        if not isinstance(scale, (int, float)) or not math.isfinite(scale) or scale <= 0:
            raise Refused("%s has a scale of %r" % (what, scale))
    return [float(s) for s in scales]


def _rows_sum_to_one(values: Sequence[float], widths: Sequence[int], what: str) -> None:
    at = 0
    for width in widths:
        total = sum(values[at:at + width])
        if abs(total - 255) > _ROUNDING * width + 1:
            raise Refused("%s has a row summing to %d of 255" % (what, total))
        at += width


def _peptides(data: Dict) -> List[str]:
    peptides = data.get("peptides")
    if not isinstance(peptides, list) or not peptides:
        raise Refused("no peptides")
    sequences = []
    for peptide in peptides:
        sequence = peptide.get("sequence") if isinstance(peptide, dict) else None
        if (not isinstance(sequence, str)
                or not MIN_LENGTH <= len(sequence) <= MAX_LENGTH
                or not set(AMINO_ACIDS).issuperset(sequence)):
            raise Refused("%r is not a peptide the model takes" % (sequence,))
        sources = peptide.get("sources")
        if not sources or not all(
                isinstance(s, dict) and all(isinstance(s.get(k), str) and s.get(k)
                                            for k in ("slug", "gene", "label"))
                for s in sources):
            raise Refused("%s names no source" % sequence)
        sequences.append(sequence)
    if sequences != sorted(sequences, key=lambda s: (len(s), s)) or len(set(sequences)) != len(sequences):
        raise Refused("the peptides are not each once, shortest first")
    return sequences


def _batch(data: Dict, sequences: Sequence[str]) -> None:
    # The padded frames, cut as `sequences.batches` cuts them: inputs drop the
    # frame's last column and targets its first, so a shorter peptide's
    # inputs carry its <eos> too, as a position whose target is padding.
    width = max(len(s) for s in sequences) + 2
    frames = [framed(s) + [PAD_ID] * (width - len(s) - 2) for s in sequences]
    batch = data.get("batch") or {}
    for side, cut in (("x", slice(None, -1)), ("y", slice(1, None))):
        rows = batch.get(side)
        if not isinstance(rows, list) or len(rows) != len(sequences):
            raise Refused("the batch's %s does not hold one row a peptide" % side)
        for row, frame, sequence in zip(rows, frames, sequences):
            if row != frame[cut]:
                raise Refused("the batch's %s does not frame %s" % (side, sequence))


def validate(payload: bytes, name: Optional[str] = None) -> Dict:
    """What the run is, if it is one the app can play; `Refused` saying why if not."""
    if len(payload) > LIMIT:
        raise Refused("%d bytes, over the %d the app is asked to fetch" % (len(payload), LIMIT))
    try:
        data = json.loads(payload)
    except ValueError:
        raise Refused("not JSON") from None
    if not isinstance(data, dict) or data.get("schema") != SCHEMA:
        raise Refused("not a schema %d run" % SCHEMA)
    if name is not None and data.get("name") != name:
        raise Refused("the run is %r, not %r" % (data.get("name"), name))
    if data.get("vocab") != list(TOKENS):
        raise Refused("the vocabulary is not the model's")
    if data.get("allowed") != list(ALLOWED):
        raise Refused("the tokens that can come next are not the model's")

    model = data.get("model") or {}
    for key in ("d_model", "nhead", "num_layers", "dim_feedforward", "parameters"):
        if not isinstance(model.get(key), int) or model[key] <= 0:
            raise Refused("the model has no %s" % key)
    width, heads, layers = model["d_model"], model["nhead"], model["num_layers"]
    if width % heads:
        raise Refused("%d does not split into %d heads" % (width, heads))

    sequences = _peptides(data)
    _batch(data, sequences)
    spans = [len(s) + 1 for s in sequences]
    targets = sum(spans)
    if data.get("targets") != targets:
        raise Refused("the run counts %r targets, not %d" % (data.get("targets"), targets))
    vocab, allowed = len(TOKENS), len(ALLOWED)
    points = data.get("points") or []
    groups = data.get("groups") or []
    if len(points) != 2 * layers + 2 or len(groups) != 3 * layers + 4:
        raise Refused("the points or groups do not fit %d layers" % layers)

    reference = data.get("baselines") or {}
    for key in ("uniform", "unigram", "bigram"):
        if not isinstance(reference.get(key), (int, float)) or not reference[key] > 0:
            raise Refused("no %s baseline" % key)
    if abs(reference["uniform"] - math.log(allowed)) > 1e-6:
        raise Refused("the uniform baseline is not ln %d" % allowed)
    floor = data.get("floor")
    if not isinstance(floor, (int, float)) or not math.isfinite(floor) or floor < 0:
        raise Refused("no floor")

    steps = data.get("steps") or {}
    count = steps.get("count")
    if not isinstance(count, int) or count < 1:
        raise Refused("no steps")
    if (data.get("config") or {}).get("max_epochs") != count:
        raise Refused("the run is not one step an epoch")
    loss = numbers(steps.get("loss"), [count], "loss")
    _finite(loss, "loss", low=0)
    for key in ("residue_loss", "eos_loss", "grad_norm"):
        _finite(numbers(steps.get(key), [count], key), key, low=0)
    rates = numbers(steps.get("lr"), [count], "lr")
    _finite(rates, "lr")
    clip = numbers(steps.get("clip"), [count], "clip")
    if not all(0 < rate for rate in rates) or not all(0 < c <= 1 for c in clip):
        raise Refused("a learning rate or clip is out of range")
    for key in ("group_grad", "group_update"):
        _finite(numbers(steps.get(key), [count, len(groups)], key), key, low=0)
    if any(value < floor - 1e-4 for value in loss):
        raise Refused("a step's loss is under the floor no model can pass")
    token_spec = steps.get("token_loss") or {}
    token_scale = _scale(token_spec, "token_loss")
    token_loss = numbers(token_spec, [count, targets], "token_loss")
    for step in range(count):
        mean = sum(token_loss[step * targets:(step + 1) * targets]) * token_scale / targets
        if abs(mean - loss[step]) > 2e-3:
            raise Refused("step %d's targets average %.4f, its loss is %.4f"
                          % (step + 1, mean, loss[step]))
    head = steps.get("head_weight") or {}
    numbers(head, [count, vocab, width], "head_weight")
    _scales(head, count, "head_weight")
    _finite(numbers(steps.get("head_bias"), [count, vocab], "head_bias"), "head_bias")

    checkpoints = data.get("checkpoints")
    if not isinstance(checkpoints, list) or not checkpoints:
        raise Refused("no checkpoints")
    at = [c.get("step") for c in checkpoints]
    if at != sorted(set(at)) or at[0] != 1 or at[-1] != count:
        raise Refused("the checkpoints are not steps 1 to %d in order" % count)
    triangle = sum(span * (span + 1) // 2 for span in spans)
    rows = [i + 1 for span in spans for i in range(span)]
    for checkpoint in checkpoints:
        what = "step %d" % checkpoint["step"]
        if not checkpoint.get("why") or not all(isinstance(w, str) for w in checkpoint["why"]):
            raise Refused("%s says no why" % what)
        stream = checkpoint.get("points") or {}
        numbers(stream, [len(points), targets, width], what + " points")
        _scales(stream, len(points), what + " points")
        _finite(numbers(checkpoint.get("point_grad"), [len(points), targets], what + " point_grad"),
                what + " point_grad", low=0)
        maps = checkpoint.get("attention") or {}
        _scale(maps, what + " attention")
        values = numbers(maps, [layers, heads, triangle], what + " attention")
        for start in range(0, len(values), triangle):
            _rows_sum_to_one(values[start:start + triangle], rows, what + " attention")
        probs = checkpoint.get("probs") or {}
        _scale(probs, what + " probs")
        _rows_sum_to_one(numbers(probs, [targets, allowed], what + " probs"),
                         [allowed] * targets, what + " probs")
        for key, shape in (("logits", [targets, vocab]), ("head_grad", [vocab, width]),
                           ("head_update", [vocab, width])):
            spec = checkpoint.get(key) or {}
            numbers(spec, shape, "%s %s" % (what, key))
            _scale(spec, "%s %s" % (what, key))

    return {
        "schema": SCHEMA,
        "name": data["name"],
        "peptides": len(sequences),
        "targets": targets,
        "parameters": model["parameters"],
        "steps": count,
        "checkpoints": at,
        "loss": [loss[0], loss[-1]],
        "floor": floor,
        "provenance": data.get("provenance") or {},
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = list(argv if argv is not None else sys.argv[1:])
    path = Path(args[0]) if args else DEFAULT
    try:
        summary = validate(path.read_bytes())
    except Refused as exc:
        print("%s: REFUSED -- %s" % (path, exc), file=sys.stderr)
        return 1
    print("%s: %s, %d steps, checkpoints %s, loss %.4f -> %.4f"
          % (path, summary["name"], summary["steps"],
             ", ".join(str(s) for s in summary["checkpoints"]),
             summary["loss"][0], summary["loss"][1]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
