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

SCHEMA = 4
LIMIT = 4_000_000
DEFAULT = HERE / "runs" / "curated" / "trace.json"

_FORMATS = {"i1": ("b", 1), "u1": ("B", 1), "u2": ("H", 2), "f4": ("f", 4)}

# A row of quantised probabilities sums to 255 give or take half a step a cell.
_ROUNDING = 0.5

# GELU's least value, at x = -0.7518: no neuron of the feed-forward goes below it.
_GELU_LEAST = -0.16997

# Float32 against the same sum worked out in Python's doubles.
_WHOLE = 1e-4

# The engine's LayerNorm's epsilon (`nn.TransformerEncoderLayer`'s default).
_NORM_EPS = 1e-5


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


def _whole(spec: Dict, shape: Sequence[int], what: str) -> List[float]:
    """An array kept whole, as float32, and finite."""
    if not isinstance(spec, dict) or spec.get("dtype") != "f4":
        raise Refused("%s is not kept whole" % what)
    values = numbers(spec, shape, what)
    _finite(values, what)
    return values


def _scaled(spec: Dict, shape: Sequence[int], what: str) -> List[float]:
    """An int8 array with its one scale, as the values it stands for."""
    if not isinstance(spec, dict) or spec.get("dtype") != "i1":
        raise Refused("%s is not int8" % what)
    values = numbers(spec, shape, what)
    scale = _scale(spec, what)
    return [v * scale for v in values]


def _gelu(x: float) -> float:
    return 0.5 * x * (1 + math.erf(x / math.sqrt(2)))


def _blocks(checkpoint: Dict, stream: Sequence[float], attention: Sequence[float], layers: int, heads: int,
            width: int, ff: int, targets: int, spans: Sequence[int], triangle: int, what: str) -> None:
    """What each block worked out, as the stages draw it: its row evened out
    by its first LayerNorm, that norm's scale and shift, the weights that make
    the queries, keys and values, those, its scores, each head's mix of values,
    and its feed-forward's input, neurons and weights. Held to each other and
    to the stream: the evened row is the stream before the block, less its
    mean, over its spread, scaled and shifted; the queries, keys and values are
    the weights' sums of it; each score is a query times a key, number by
    number, summed, over the root of a head's size; softmax of the scores is the
    block's attention map, each head's mix is its shares of the values, every
    neuron is GELU of its weighted input, and what the feed-forward adds is
    what the stream gained across it."""
    blocks = checkpoint.get("blocks")
    if not isinstance(blocks, list) or len(blocks) != layers:
        raise Refused("%s does not keep its %d blocks" % (what, layers))
    size = width // heads
    for b, block in enumerate(blocks):
        named = "%s block %d" % (what, b + 1)
        if not isinstance(block, dict):
            raise Refused("%s is not an object" % named)
        gain = _whole(block.get("norm1_weight"), [width], named + " norm1_weight")
        shift = _whole(block.get("norm1_bias"), [width], named + " norm1_bias")
        evened = _whole(block.get("attn_input"), [targets, width], named + " attn_input")
        before = stream[2 * b * targets * width:]
        for t in range(targets):
            row = before[t * width:(t + 1) * width]
            mean = sum(row) / width
            spread = math.sqrt(sum((v - mean) ** 2 for v in row) / width + _NORM_EPS)
            for d in range(width):
                if abs((row[d] - mean) / spread * gain[d] + shift[d] - evened[t * width + d]) > _WHOLE:
                    raise Refused("%s attn_input is not the stream evened out" % named)
        made = _whole(block.get("attn_in"), [3 * width, width], named + " attn_in")
        made_bias = _whole(block.get("attn_in_bias"), [3 * width], named + " attn_in_bias")
        query = _whole(block.get("query"), [targets, width], named + " query")
        key = _whole(block.get("key"), [targets, width], named + " key")
        value = _whole(block.get("value"), [targets, width], named + " value")
        for part, kept in enumerate((query, key, value)):
            for t in range(targets):
                x = evened[t * width:(t + 1) * width]
                for d in range(width):
                    out = part * width + d
                    total = sum(made[out * width + e] * x[e] for e in range(width)) + made_bias[out]
                    if abs(total - kept[t * width + d]) > _WHOLE:
                        raise Refused("%s %s is not its weights' sum of the evened row"
                                      % (named, ("query", "key", "value")[part]))
        mixed = _scaled(block.get("mixed"), [targets, width], named + " mixed")
        mixed_step = _scale(block["mixed"], named + " mixed")
        scores = _whole(block.get("scores"), [heads, triangle], named + " scores")
        root = math.sqrt(size)
        for h in range(heads):
            at, first = h * triangle, 0
            for span in spans:
                for i in range(span):
                    q = query[(first + i) * width + h * size:(first + i) * width + (h + 1) * size]
                    for j in range(i + 1):
                        k = key[(first + j) * width + h * size:(first + j) * width + (h + 1) * size]
                        if abs(sum(a * c for a, c in zip(q, k)) / root - scores[at + j]) > _WHOLE:
                            raise Refused("%s scores are not query times key over the root of a head" % named)
                    at += i + 1
                first += span
        maps = attention[b * heads * triangle:(b + 1) * heads * triangle]
        for h in range(heads):
            at, first = h * triangle, 0
            for span in spans:
                for i in range(span):
                    row = scores[at:at + i + 1]
                    top = max(row)
                    exps = [math.exp(v - top) for v in row]
                    total = sum(exps)
                    shares = [e / total for e in exps]
                    for j, share in enumerate(shares):
                        if abs(share * 255 - maps[at + j]) > _ROUNDING + 1e-3:
                            raise Refused("%s scores do not softmax to its attention" % named)
                    target = first + i
                    for d in range(h * size, (h + 1) * size):
                        mix = sum(share * value[(first + j) * width + d] for j, share in enumerate(shares))
                        if abs(mix - mixed[target * width + d]) > mixed_step / 2 + _WHOLE:
                            raise Refused("%s mix is not its shares of the values" % named)
                    at += i + 1
                first += span
        ff_input = _whole(block.get("ff_input"), [targets, width], named + " ff_input")
        hidden = _scaled(block.get("hidden"), [targets, ff], named + " hidden")
        hidden_step = _scale(block["hidden"], named + " hidden")
        if min(hidden) < _GELU_LEAST - hidden_step / 2:
            raise Refused("%s has a neuron under GELU's least" % named)
        w_in = _whole(block.get("ff_in"), [ff, width], named + " ff_in")
        b_in = _whole(block.get("ff_in_bias"), [ff], named + " ff_in_bias")
        w_out = _whole(block.get("ff_out"), [width, ff], named + " ff_out")
        b_out = _whole(block.get("ff_out_bias"), [width], named + " ff_out_bias")
        middle, after = stream[(2 * b + 1) * targets * width:], stream[(2 * b + 2) * targets * width:]
        for t in range(targets):
            x = ff_input[t * width:(t + 1) * width]
            neurons = []
            for u in range(ff):
                neuron = _gelu(sum(w_in[u * width + d] * x[d] for d in range(width)) + b_in[u])
                if abs(neuron - hidden[t * ff + u]) > hidden_step / 2 + _WHOLE:
                    raise Refused("%s neurons are not GELU of their input" % named)
                neurons.append(neuron)
            for d in range(width):
                added = sum(w_out[d * ff + u] * neurons[u] for u in range(ff)) + b_out[d]
                gained = after[t * width + d] - middle[t * width + d]
                if abs(added - gained) > _WHOLE:
                    raise Refused("%s feed-forward does not add what the stream gained" % named)


def _readout(checkpoint: Dict, points: int, targets: int, allowed: int, what: str) -> None:
    """The stream read out early at every point but the last (the final norm
    is the readout's own): each row a distribution over what can come next,
    and the last point's the model's own probabilities."""
    spec = checkpoint.get("readout") or {}
    _scale(spec, what + " readout")
    values = numbers(spec, [points - 1, targets, allowed], what + " readout")
    _rows_sum_to_one(values, [allowed] * ((points - 1) * targets), what + " readout")
    probs = numbers(checkpoint.get("probs") or {}, [targets, allowed], what + " probs")
    last = values[(points - 2) * targets * allowed:]
    if any(abs(a - b) > 1 for a, b in zip(last, probs)):
        raise Refused("%s readout at the last point is not the model's probabilities" % what)


def _embedded(checkpoint: Dict, x0: Sequence[float], inputs: Sequence[Sequence[int]],
              spans: Sequence[int], vocab: int, context: int, width: int, what: str) -> None:
    """The tables the step read its ids through, whole (float32): the token
    embedding, its padding row zero, and the position embedding; and the
    stream's first point, whole too, each target's token row plus its
    place's."""
    tables = []
    for key, rows in (("token_embedding", vocab), ("position_embedding", context)):
        spec = checkpoint.get(key) or {}
        if spec.get("dtype") != "f4":
            raise Refused("%s %s is not kept whole" % (what, key))
        values = numbers(spec, [rows, width], "%s %s" % (what, key))
        _finite(values, "%s %s" % (what, key))
        tables.append(values)
    token, place = tables
    if any(token[PAD_ID * width:(PAD_ID + 1) * width]):
        raise Refused("%s's padding row of the token embedding is not zero" % what)
    target = 0
    for row, span in zip(inputs, spans):
        for i in range(span):
            for d in range(width):
                made = token[row[i] * width + d] + place[i * width + d]
                if abs(made - x0[target * width + d]) > 1e-5:
                    raise Refused("%s x0 is not each id's token row plus its place's" % what)
            target += 1


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
    if not isinstance(model.get("max_len"), int) or model["max_len"] <= 0:
        raise Refused("the model has no max_len")
    width, heads, layers = model["d_model"], model["nhead"], model["num_layers"]
    if width % heads:
        raise Refused("%d does not split into %d heads" % (width, heads))

    sequences = _peptides(data)
    _batch(data, sequences)
    spans = [len(s) + 1 for s in sequences]
    if max(spans) > model["max_len"]:
        raise Refused("a peptide is longer than the %d places the model holds" % model["max_len"])
    targets = sum(spans)
    if data.get("targets") != targets:
        raise Refused("the run counts %r targets, not %d" % (data.get("targets"), targets))
    vocab, allowed = len(TOKENS), len(ALLOWED)
    points = data.get("points") or []
    groups = data.get("groups") or []
    if len(points) != 2 * layers + 2 or len(groups) != 3 * layers + 4:
        raise Refused("the points or groups do not fit %d layers" % layers)
    sizes = data.get("group_sizes")
    if (not isinstance(sizes, list) or len(sizes) != len(groups)
            or not all(isinstance(n, int) and n > 0 for n in sizes)
            or sum(sizes) != model["parameters"]):
        raise Refused("the groups' sizes do not add up to the model's parameters")

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
        flowing = _whole(checkpoint.get("points"), [len(points), targets, width], what + " points")
        _embedded(checkpoint, flowing[:targets * width], data["batch"]["x"], spans,
                  vocab, model["max_len"], width, what)
        _finite(numbers(checkpoint.get("point_grad"), [len(points), targets], what + " point_grad"),
                what + " point_grad", low=0)
        maps = checkpoint.get("attention") or {}
        _scale(maps, what + " attention")
        values = numbers(maps, [layers, heads, triangle], what + " attention")
        for start in range(0, len(values), triangle):
            _rows_sum_to_one(values[start:start + triangle], rows, what + " attention")
        _blocks(checkpoint, flowing, values, layers, heads, width, model["dim_feedforward"], targets,
                spans, triangle, what)
        _readout(checkpoint, len(points), targets, allowed, what)
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
