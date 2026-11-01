"""Record one training run of the model, number by number, for the app to play.

    python pipeline/denovo/trace.py                      # runs/curated/trace.json
    python pipeline/denovo/trace.py --check              # record again, compare, write nothing
    python pipeline/denovo/trace.py --check --against https://helix-peak-backend.onrender.com
    python pipeline/denovo/trace.py --toy --out PATH     # the tiny run the tests read

The run is the engine's own: `train.build` and `train.train_step`, driven as
`train.run` drives them, on the curated proteins' peptides, with nothing held
out. Read-only hooks watch it. Every step's loss, learning rate, gradient norms,
per-target losses and output-layer weights are kept, and at each checkpoint so
is everything else: the embedding tables the step read its ids through, the
residual stream at every layer, its gradient, every attention map, what each
block worked out on the way (its row evened out by its first LayerNorm, that
norm's scale and shift, the weights that make the queries, keys and values from
it, the queries, keys, values and scores, each head's mix of values, the output
weights that mix the heads back into one row, its feed-forward's input, neurons
and weights), the stream read
out early through the final norm and the output layer, the output layer's
weights as the step found them, the logits, and the update the step made.

It is recorded twice. The first pass picks the checkpoints from the losses;
the second keeps the activations at them and must train to the same losses,
float for float, or nothing is written.
"""

from __future__ import annotations

import argparse
import base64
import contextlib
import hashlib
import json
import math
import sys
import urllib.error
import urllib.request
from collections import Counter
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Callable, Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import torch
import torch.nn.functional as F

HERE = Path(__file__).resolve().parent
if str(HERE.parents[1]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[1]))

from pipeline.denovo.alphabet import EOS_ID, MAX_CONTEXT, TOKENS
from pipeline.denovo.baselines import baselines
from pipeline.denovo.bigram import ALLOWED, framed
from pipeline.denovo.bigram_torch import masked
from pipeline.denovo.helix_peek_LM import parameter_count
from pipeline.denovo.sequences import batches, build_dataset, is_modelable
from pipeline.denovo.tokenizer import decode
from pipeline.denovo.train import Config, _token_losses, build, train_step

SCHEMA = 6
NAME = "curated"
OUT = HERE / "runs" / NAME / "trace.json"
CPU = torch.device("cpu")

# The engine's model and optimiser at a size of its own, with dropout off:
# the forward pass drawn is then the forward pass trained, number for number.
# The page writes out every number of an embedding row, and 16 of them fit
# across a phone; it walks one block's attention and feed-forward, so two
# blocks of two heads; at the swept 1e-3 so narrow a model never passes one
# nat in the run's 200 steps, at 1e-2 it passes every milestone.
SHOWCASE = replace(Config(), dropout=0.0, d_model=16, dim_feedforward=64, learning_rate=1e-2,
                   num_layers=2, nhead=2)

# Small enough to commit as a fixture, and fast enough to pass every milestone
# up to one nat in its twenty steps.
TOY = replace(SHOWCASE, d_model=8, nhead=2, num_layers=2, dim_feedforward=16,
              max_epochs=20, warmup_steps=2, learning_rate=5e-2)

# A target's loss as uint16, in 1/4096 of a nat: sixteen nats of range.
LOSS_SCALE = 1.0 / 4096

ATTENTION_TOLERANCE = 1e-4


@dataclass(frozen=True)
class Source:
    slug: str
    gene: str
    label: str


@dataclass(frozen=True)
class Peptide:
    sequence: str
    sources: Tuple[Source, ...]


TOY_PEPTIDES = (
    Peptide("GIVEQCCTSIC", (Source("toy", "TOY", "toy A"),)),
    Peptide("QLYSALANKCCH", (Source("toy", "TOY", "toy B"),)),
)


def point_names(layers: int) -> List[str]:
    names = ["x0"]
    for i in range(1, layers + 1):
        names += ["m%d" % i, "h%d" % i]
    return names + ["hf"]


def group_names(layers: int) -> List[str]:
    names = ["token_embedding", "position_embedding"]
    for i in range(1, layers + 1):
        names += ["block%d.attention" % i, "block%d.feed_forward" % i, "block%d.norms" % i]
    return names + ["final_norm", "head"]


def group_of(parameter: str) -> str:
    if parameter.startswith("token_embedding."):
        return "token_embedding"
    if parameter == "positional_embedding":
        return "position_embedding"
    if parameter.startswith("transformer_encoder.layers."):
        _, _, index, part = parameter.split(".")[:4]
        kind = {"self_attn": "attention", "linear1": "feed_forward",
                "linear2": "feed_forward", "norm1": "norms", "norm2": "norms"}[part]
        return "block%d.%s" % (int(index) + 1, kind)
    if parameter.startswith("transformer_encoder.norm."):
        return "final_norm"
    if parameter.startswith("lm_head."):
        return "head"
    raise ValueError("no group for parameter %r" % parameter)


def floor(sequences: Sequence[str]) -> float:
    """The least loss any model reading whole prefixes can reach on these
    sequences: the entropy of what follows each prefix, per target."""
    following: Dict[Tuple[int, ...], Counter] = {}
    targets = 0
    for sequence in sequences:
        ids = framed(sequence)
        for i in range(1, len(ids)):
            following.setdefault(tuple(ids[:i]), Counter())[ids[i]] += 1
            targets += 1
    cost = 0.0
    for counts in following.values():
        seen = sum(counts.values())
        cost += sum(-count * math.log(count / seen) for count in counts.values())
    return cost / targets


def milestones(losses: Sequence[float], reference: Dict[str, float]) -> List[Tuple[int, List[str]]]:
    """Step 1, the first step below each baseline, below one nat and below a
    tenth of one, and the last."""
    marks: Dict[int, List[str]] = {1: ["first"]}
    levels = (
        ("uniform", reference["uniform"]),
        ("unigram", reference["unigram"]),
        ("bigram", reference["bigram"]),
        ("one_nat", 1.0),
        ("tenth_nat", 0.1),
    )
    for name, level in levels:
        step = next((i + 1 for i, loss in enumerate(losses) if loss < level), None)
        if step is not None:
            marks.setdefault(step, []).append(name)
    marks.setdefault(len(losses), []).append("last")
    return sorted(marks.items())


def inside(
    module: torch.nn.MultiheadAttention,
    query: torch.Tensor,
    mask: Optional[torch.Tensor],
    output: torch.Tensor,
) -> Tuple[Dict[str, torch.Tensor], float]:
    """What a layer's attention worked out, again from its input and its
    weights as they stood: its input (the row its block's first LayerNorm
    evened out), the queries, keys and values (heads side by side),
    every score before softmax (heads apart), the weights, and each head's mix
    of values before the output projection; and how far the output they give
    is from the layer's own (a check, not a correction)."""
    with torch.no_grad():
        batch, length, width = query.shape
        heads = module.num_heads
        size = width // heads
        made = F.linear(query, module.in_proj_weight, module.in_proj_bias).chunk(3, dim=-1)
        q, k, v = (t.reshape(batch, length, heads, size).transpose(1, 2) for t in made)
        scores = q @ k.transpose(-2, -1) / math.sqrt(size)
        if mask is None:
            mask = torch.triu(torch.ones(length, length, dtype=torch.bool), diagonal=1)
        if mask.dtype == torch.bool:
            masked_scores = scores.masked_fill(mask, float("-inf"))
        else:
            masked_scores = scores + mask
        weights = masked_scores.softmax(dim=-1)
        mixed = (weights @ v).transpose(1, 2).reshape(batch, length, width)
        out = F.linear(mixed, module.out_proj.weight, module.out_proj.bias)
        parts = {"attn_input": query.detach(), "query": made[0], "key": made[1], "value": made[2],
                 "scores": scores, "weights": weights, "mixed": mixed}
        return parts, (out - output).abs().max().item()


def attention(
    module: torch.nn.MultiheadAttention,
    query: torch.Tensor,
    mask: Optional[torch.Tensor],
    output: torch.Tensor,
) -> Tuple[torch.Tensor, float]:
    """The attention weights a layer used, worked out again (`inside`), and
    how far the output they give is from the layer's own."""
    parts, error = inside(module, query, mask, output)
    return parts["weights"], error


class _Capture:
    """Hooks that only look: the logits at every step and, while `on`, every
    point of the residual stream with its gradient kept, every attention map,
    and what each block worked out on the way (`inside`, the feed-forward's
    input and its neurons)."""

    def __init__(self, model: torch.nn.Module):
        self.on = False
        self.logits: Optional[torch.Tensor] = None
        self.points: Dict[str, torch.Tensor] = {}
        self.attention: Dict[int, torch.Tensor] = {}
        self.blocks: Dict[int, Dict[str, torch.Tensor]] = {}
        self.error = 0.0
        layers = model.transformer_encoder.layers
        self.handles = [layers[0].register_forward_pre_hook(self._input("x0"))]
        for i, layer in enumerate(layers, start=1):
            self.handles.append(layer.norm2.register_forward_pre_hook(self._input("m%d" % i)))
            self.handles.append(layer.norm2.register_forward_hook(self._kept(i, "ff_input")))
            # The neurons as linear2 reads them: GELU's output, dropout being off.
            self.handles.append(layer.linear2.register_forward_pre_hook(self._kept_input(i, "hidden")))
            self.handles.append(layer.register_forward_hook(self._output("h%d" % i)))
            self.handles.append(layer.self_attn.register_forward_hook(
                self._attention(i), with_kwargs=True))
        self.handles.append(model.transformer_encoder.norm.register_forward_hook(self._output("hf")))
        self.handles.append(model.lm_head.register_forward_hook(self._logits))

    def clear(self) -> None:
        self.points, self.attention, self.blocks = {}, {}, {}

    def remove(self) -> None:
        for handle in self.handles:
            handle.remove()

    def _input(self, name: str) -> Callable:
        def hook(module, args):
            if self.on:
                args[0].retain_grad()
                self.points[name] = args[0]
        return hook

    def _output(self, name: str) -> Callable:
        def hook(module, args, output):
            if self.on:
                output.retain_grad()
                self.points[name] = output
        return hook

    def _kept(self, layer: int, name: str) -> Callable:
        def hook(module, args, output):
            if self.on:
                self.blocks.setdefault(layer, {})[name] = output.detach()
        return hook

    def _kept_input(self, layer: int, name: str) -> Callable:
        def hook(module, args):
            if self.on:
                self.blocks.setdefault(layer, {})[name] = args[0].detach()
        return hook

    def _attention(self, layer: int) -> Callable:
        def hook(module, args, kwargs, output):
            if self.on:
                parts, error = inside(module, args[0], kwargs.get("attn_mask"), output[0])
                self.attention[layer] = parts["weights"]
                self.blocks.setdefault(layer, {}).update(parts)
                self.error = max(self.error, error)
        return hook

    def _logits(self, module, args, output) -> None:
        self.logits = output.detach()


def _norm(tensors: Sequence[torch.Tensor]) -> float:
    return math.sqrt(sum(float(t.double().pow(2).sum()) for t in tensors))


def readout_names(layers: int) -> List[str]:
    """The points read out early: every point of the stream before the final
    norm, which the readout itself applies."""
    return point_names(layers)[:-1]


def readout(model: torch.nn.Module, point: torch.Tensor) -> torch.Tensor:
    """What the output layer would say if the model stopped at `point`: the
    final norm and the output layer applied to it, with the weights as they
    stand, masked and normalised as `probs` is, over the tokens that can come
    next. Worked out with the modules' weights, not by calling them, so no
    hook sees it."""
    with torch.no_grad():
        norm, head = model.transformer_encoder.norm, model.lm_head
        normed = F.layer_norm(point, norm.normalized_shape, norm.weight, norm.bias, norm.eps)
        logits = F.linear(normed, head.weight, head.bias)
        probabilities = masked(logits.reshape(-1, logits.shape[-1])).softmax(dim=-1)
        return probabilities.reshape(logits.shape)[..., list(ALLOWED)]


@contextlib.contextmanager
def deterministic() -> Iterator[None]:
    """One thread and deterministic kernels, so a run records the same bytes
    every time on the same machine; what was set before is put back."""
    was_deterministic = torch.are_deterministic_algorithms_enabled()
    threads = torch.get_num_threads()
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(1)
    try:
        yield
    finally:
        torch.use_deterministic_algorithms(was_deterministic)
        torch.set_num_threads(threads)


def _run(config: Config, sequences: Sequence[str], capture_at: frozenset) -> Dict:
    torch.manual_seed(config.seed)
    frames, lengths = build_dataset(sequences)
    if math.ceil(len(sequences) / config.batch_size) != 1:
        raise ValueError("a recorded run takes all its peptides in one batch: %d peptides, batch %d"
                         % (len(sequences), config.batch_size))
    model, optimizer, scheduler = build(config, config.max_epochs, CPU)
    capture = _Capture(model)
    groups = group_names(config.num_layers)
    members: Dict[str, List[torch.nn.Parameter]] = {name: [] for name in groups}
    for name, parameter in model.named_parameters():
        members[group_of(name)].append(parameter)
    names = point_names(config.num_layers)
    head = model.lm_head
    canonical = {sequence: i for i, sequence in enumerate(sequences)}
    spans = [len(sequence) + 1 for sequence in sequences]

    kept: Dict[str, list] = {key: [] for key in (
        "loss", "residue_loss", "eos_loss", "lr", "grad_norm", "clip",
        "group_grad", "group_update", "token_loss", "head_weight", "head_bias")}
    checkpoints: Dict[int, Dict] = {}
    batch = None
    try:
        for epoch in range(1, config.max_epochs + 1):
            generator = torch.Generator().manual_seed(config.seed * 10_000 + epoch)
            for x, y in batches(frames, lengths, config.batch_size, generator):
                step = epoch
                order = [canonical[decode(row.tolist())] for row in x]
                rows = [order.index(c) for c in range(len(sequences))]

                def real(t: torch.Tensor) -> torch.Tensor:
                    return torch.cat([t[rows[c], :spans[c]] for c in range(len(sequences))])

                if batch is None:
                    batch = (x[rows].tolist(), y[rows].tolist())
                capture.on = step in capture_at
                capture.clear()
                lr = optimizer.param_groups[0]["lr"]
                before = {name: [p.detach().clone() for p in members[name]] for name in groups}
                weight, bias = head.weight.detach().clone(), head.bias.detach().clone()
                # Each block's weights as the step found them: its first
                # LayerNorm's scale and shift, the weights that make the
                # queries, keys and values and those that mix the heads' answers
                # back into one row, which the page writes out; and its
                # feed-forward, the edges the page draws between its neurons.
                fed = [{
                    "norm1_weight": layer.norm1.weight.detach().clone(),
                    "norm1_bias": layer.norm1.bias.detach().clone(),
                    "attn_in": layer.self_attn.in_proj_weight.detach().clone(),
                    "attn_in_bias": layer.self_attn.in_proj_bias.detach().clone(),
                    "attn_out": layer.self_attn.out_proj.weight.detach().clone(),
                    "attn_out_bias": layer.self_attn.out_proj.bias.detach().clone(),
                    "ff_in": layer.linear1.weight.detach().clone(),
                    "ff_in_bias": layer.linear1.bias.detach().clone(),
                    "ff_out": layer.linear2.weight.detach().clone(),
                    "ff_out_bias": layer.linear2.bias.detach().clone(),
                } for layer in model.transformer_encoder.layers] if capture.on else []
                seen: Dict[str, object] = {}

                def after_backward() -> None:
                    seen["group_grad"] = [
                        _norm([p.grad for p in members[name] if p.grad is not None])
                        for name in groups]
                    if capture.on:
                        seen["point_grad"] = {
                            name: tensor.grad.detach().norm(dim=-1)
                            for name, tensor in capture.points.items()}
                        # Before the weights move: the step's own norm and head.
                        seen["readout"] = [readout(model, capture.points[name].detach())
                                           for name in readout_names(config.num_layers)]

                loss, _, norm = train_step(
                    model, optimizer, scheduler, x, y, config.grad_clip, after_backward)

                per_target = real(_token_losses(capture.logits, y).reshape(y.shape))
                targets = real(y)
                ends = targets == EOS_ID
                kept["loss"].append(float(loss))
                kept["residue_loss"].append(float(per_target[~ends].mean()))
                kept["eos_loss"].append(float(per_target[ends].mean()))
                kept["lr"].append(float(lr))
                kept["grad_norm"].append(float(norm))
                kept["clip"].append(min(1.0, config.grad_clip / (float(norm) + 1e-6)))
                kept["group_grad"].append(seen["group_grad"])
                kept["group_update"].append([
                    _norm([p.detach() - b for p, b in zip(members[name], before[name])])
                    for name in groups])
                kept["token_loss"].append(per_target.numpy().copy())
                kept["head_weight"].append(weight.numpy())
                kept["head_bias"].append(bias.numpy())

                if capture.on:
                    logits = capture.logits
                    probabilities = masked(logits.reshape(-1, logits.shape[-1])).softmax(dim=-1)
                    probabilities = probabilities.reshape(logits.shape)[..., list(ALLOWED)]

                    def triangles(squares: torch.Tensor) -> np.ndarray:
                        """Each head's batch × length × length, as each
                        peptide's lower triangle in a row: what reads what."""
                        heads = []
                        for h in range(config.nhead):
                            pieces = []
                            for c in range(len(sequences)):
                                square = squares[rows[c], h, :spans[c], :spans[c]]
                                pieces.append(square.numpy()[np.tril_indices(spans[c])])
                            heads.append(np.concatenate(pieces))
                        return np.array(heads)

                    blocks = []
                    for layer in range(1, config.num_layers + 1):
                        kept_block = capture.blocks[layer]
                        blocks.append({
                            **{key: real(kept_block[key]).numpy()
                               for key in ("attn_input", "query", "key", "value", "mixed",
                                           "ff_input", "hidden")},
                            "scores": triangles(kept_block["scores"]),
                            **{key: tensor.numpy() for key, tensor in fed[layer - 1].items()},
                        })
                    checkpoints[step] = {
                        # As the step found them: the tables its x0 was made from.
                        "token_embedding": before["token_embedding"][0].numpy(),
                        "position_embedding": before["position_embedding"][0][0].numpy(),
                        "points": [real(capture.points[name].detach()).numpy() for name in names],
                        "point_grad": [real(seen["point_grad"][name]).numpy() for name in names],
                        "attention": np.array([triangles(capture.attention[layer])
                                               for layer in range(1, config.num_layers + 1)]),
                        "blocks": blocks,
                        "readout": np.array([real(r).numpy() for r in seen["readout"]]),
                        # The output layer as the step found it: its logits
                        # are the last point times it, plus its bias.
                        "head_weight": weight.numpy(),
                        "logits": real(logits).numpy(),
                        "probs": real(probabilities).numpy(),
                        "head_grad": head.weight.grad.detach().numpy().copy(),
                        "head_update": (head.weight.detach() - weight).numpy(),
                    }
    finally:
        capture.remove()
    if capture.error > ATTENTION_TOLERANCE:
        raise RuntimeError("the attention worked out again is %.2e from the layer's own" % capture.error)
    return {"steps": kept, "checkpoints": checkpoints, "batch": batch,
            "parameters": parameter_count(model),
            "group_sizes": [sum(p.numel() for p in members[name]) for name in groups]}


def _array(values, dtype: str, **extra) -> Dict:
    data = np.ascontiguousarray(values, dtype=np.dtype("<" + dtype))
    return {"dtype": dtype, "shape": list(data.shape),
            "data": base64.b64encode(data.tobytes()).decode("ascii"), **extra}


def _signed(values) -> Tuple[np.ndarray, float]:
    values = np.asarray(values, dtype=np.float64)
    top = float(np.abs(values).max())
    scale = top / 127 if top > 0 else 1.0
    return np.clip(np.rint(values / scale), -127, 127).astype(np.int8), scale


def _unit(values) -> np.ndarray:
    return np.clip(np.rint(np.asarray(values, dtype=np.float64) * 255), 0, 255).astype(np.uint8)


def _int8(values) -> Dict:
    quantized, scale = _signed(values)
    return _array(quantized, "i1", scale=scale)


def _int8_each(slices) -> Dict:
    pairs = [_signed(s) for s in slices]
    return _array(np.stack([q for q, _ in pairs]), "i1", scales=[s for _, s in pairs])


def record(
    config: Config,
    peptides: Sequence[Peptide],
    name: str = NAME,
    records: Optional[Dict[str, str]] = None,
) -> Dict:
    """The run, recorded: a payload `serialize` writes and `check_trace` reads."""
    peptides = sorted(peptides, key=lambda p: (len(p.sequence), p.sequence))
    sequences = [p.sequence for p in peptides]
    if len(set(sequences)) != len(sequences):
        raise ValueError("a peptide is given twice")
    for sequence in sequences:
        if not is_modelable(sequence):
            raise ValueError("%r is not a peptide the model takes" % sequence)

    with deterministic():
        plain = _run(config, sequences, frozenset())
        reference = baselines(sequences, sequences)
        marks = milestones(plain["steps"]["loss"], reference)
        played = _run(config, sequences, frozenset(step for step, _ in marks))
    if played["steps"]["loss"] != plain["steps"]["loss"]:
        raise RuntimeError("the run with its activations recorded trained differently")

    steps = played["steps"]
    targets = sum(len(sequence) + 1 for sequence in sequences)
    checkpoints = []
    for step, why in marks:
        kept = played["checkpoints"][step]
        checkpoints.append({
            "step": step,
            "why": why,
            # Whole, as float32: the page writes every number of them out
            # to two places, and their sum is x0 itself.
            "token_embedding": _array(kept["token_embedding"], "f4"),
            "position_embedding": _array(kept["position_embedding"], "f4"),
            # Whole too: the page writes each point's numbers and what each
            # part added, and they add up.
            "points": _array(np.stack(kept["points"]), "f4"),
            "point_grad": _array(np.stack(kept["point_grad"]), "f4"),
            "attention": _array(_unit(kept["attention"]), "u1", scale=1.0 / 255),
            # What each block worked out: drawn as colours where it is int8,
            # written out where it is whole. The evened row, the weights that
            # make the queries, keys and values from it, and those, are whole:
            # the page writes them out, and a score is their sum. So are the
            # heads' mix and the output weights that make what attention
            # adds from it, which the page writes out too.
            "blocks": [{
                "norm1_weight": _array(block["norm1_weight"], "f4"),
                "norm1_bias": _array(block["norm1_bias"], "f4"),
                "attn_input": _array(block["attn_input"], "f4"),
                "attn_in": _array(block["attn_in"], "f4"),
                "attn_in_bias": _array(block["attn_in_bias"], "f4"),
                "query": _array(block["query"], "f4"),
                "key": _array(block["key"], "f4"),
                "value": _array(block["value"], "f4"),
                "scores": _array(block["scores"], "f4"),
                "mixed": _array(block["mixed"], "f4"),
                "attn_out": _array(block["attn_out"], "f4"),
                "attn_out_bias": _array(block["attn_out_bias"], "f4"),
                "ff_input": _array(block["ff_input"], "f4"),
                "hidden": _int8(block["hidden"]),
                "ff_in": _array(block["ff_in"], "f4"),
                "ff_in_bias": _array(block["ff_in_bias"], "f4"),
                "ff_out": _array(block["ff_out"], "f4"),
                "ff_out_bias": _array(block["ff_out_bias"], "f4"),
            } for block in kept["blocks"]],
            # The stream read out early (`readout`), at every point before
            # the final norm: the last is `probs`.
            "readout": _array(_unit(kept["readout"]), "u1", scale=1.0 / 255),
            # The output layer's weights as the step found them, and the
            # logits, whole: the page writes a logit out as the last row's
            # numbers times a token's weights, summed, plus its bias, and
            # softmax of the logits makes its shares.
            "head_weight": _array(kept["head_weight"], "f4"),
            "logits": _array(kept["logits"], "f4"),
            "probs": _array(_unit(kept["probs"]), "u1", scale=1.0 / 255),
            "head_grad": _int8(kept["head_grad"]),
            "head_update": _int8(kept["head_update"]),
        })

    token_loss = np.clip(np.rint(np.stack(steps["token_loss"]) / LOSS_SCALE), 0, 65535)
    return {
        "schema": SCHEMA,
        "name": name,
        "vocab": list(TOKENS),
        "allowed": list(ALLOWED),
        "model": {
            "name": "HelixPeekLM",
            "d_model": config.d_model,
            "nhead": config.nhead,
            "num_layers": config.num_layers,
            "dim_feedforward": config.dim_feedforward,
            "max_len": MAX_CONTEXT,
            "parameters": played["parameters"],
        },
        "config": asdict(config),
        "peptides": [
            {"sequence": p.sequence, "sources": [asdict(s) for s in p.sources]}
            for p in peptides
        ],
        "batch": {"x": played["batch"][0], "y": played["batch"][1]},
        "targets": targets,
        "baselines": reference,
        "floor": floor(sequences),
        "points": point_names(config.num_layers),
        "groups": group_names(config.num_layers),
        # How many weights each group holds: a group's norm over the root of
        # its size is the typical weight's, which is what Adam evens out.
        "group_sizes": played["group_sizes"],
        "steps": {
            "count": len(steps["loss"]),
            "loss": _array(steps["loss"], "f4"),
            "residue_loss": _array(steps["residue_loss"], "f4"),
            "eos_loss": _array(steps["eos_loss"], "f4"),
            "lr": _array(steps["lr"], "f4"),
            "grad_norm": _array(steps["grad_norm"], "f4"),
            "clip": _array(steps["clip"], "f4"),
            "group_grad": _array(steps["group_grad"], "f4"),
            "group_update": _array(steps["group_update"], "f4"),
            "token_loss": _array(token_loss, "u2", scale=LOSS_SCALE),
            "head_weight": _int8_each(steps["head_weight"]),
            "head_bias": _array(np.stack(steps["head_bias"]), "f4"),
        },
        "checkpoints": checkpoints,
        "provenance": {
            "engine": "pipeline/denovo/train.py train_step",
            "torch": torch.__version__,
            "numpy": np.__version__,
            "device": "cpu",
            "threads": 1,
            "deterministic": True,
            "peptides_sha256": hashlib.sha256(
                ("\n".join(sorted(sequences)) + "\n").encode()).hexdigest(),
            "records": dict(sorted((records or {}).items())),
        },
    }


def serialize(payload: Dict) -> bytes:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def curated() -> Tuple[List[Peptide], Dict[str, str]]:
    """The curated proteins' peptides, as `curated_peptides.collect` finds them
    on the stored records, and the digest of each record read."""
    from pipeline.denovo.curated_peptides import collect, record_path
    from pipeline.targets import TARGETS

    by_slug = {target.slug: target for target in TARGETS}
    peptides = [
        Peptide(sequence, tuple(Source(slug, by_slug[slug].gene, label)
                                for slug, _, label in found))
        for sequence, found in collect().items()
    ]
    slugs = sorted({source.slug for p in peptides for source in p.sources})
    records = {slug: hashlib.sha256(record_path(by_slug[slug]).read_bytes()).hexdigest()
               for slug in slugs}
    return peptides, records


def served_sha256(service: str, name: str) -> Optional[str]:
    url = "%s/peptide-lm/runs/%s" % (service.rstrip("/"), name)
    request = urllib.request.Request(url, headers={"User-Agent": "helix-peek-denovo"})
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            return json.load(response)["sha256"]
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        raise


def _summary(payload: Dict, body: bytes) -> str:
    def numbers(spec: Dict) -> np.ndarray:
        return np.frombuffer(base64.b64decode(spec["data"]), dtype="<" + spec["dtype"])

    loss = numbers(payload["steps"]["loss"])
    reference = payload["baselines"]
    lines = [
        "%d peptides, %d targets | %s parameters | %d steps on cpu"
        % (len(payload["peptides"]), payload["targets"],
           format(payload["model"]["parameters"], ","), payload["steps"]["count"]),
        "loss %.4f -> %.4f (uniform %.4f, unigram %.4f, bigram %.4f, floor %.4f)"
        % (loss[0], loss[-1], reference["uniform"], reference["unigram"],
           reference["bigram"], payload["floor"]),
        "checkpoints: %s" % ", ".join(
            "%d (%s)" % (c["step"], " ".join(c["why"])) for c in payload["checkpoints"]),
        "%s bytes, sha256 %s" % (format(len(body), ","), hashlib.sha256(body).hexdigest()),
    ]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--check", action="store_true",
                        help="record again and compare with --out, or with --against; write nothing")
    parser.add_argument("--against", default=None, metavar="SERVICE",
                        help="with --check: compare with the run the service names")
    parser.add_argument("--toy", action="store_true",
                        help="a tiny model on two toy peptides: the tests' fixture")
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])
    if args.against and not args.check:
        parser.error("--against compares a fresh record with what is served; add --check")

    if args.toy:
        payload = record(TOY, TOY_PEPTIDES, name="toy")
    else:
        peptides, records = curated()
        payload = record(SHOWCASE, peptides, name=NAME, records=records)
    body = serialize(payload)
    digest = hashlib.sha256(body).hexdigest()
    print(_summary(payload, body))

    if args.check:
        if args.against:
            other, where = served_sha256(args.against, payload["name"]), args.against
        else:
            other = hashlib.sha256(args.out.read_bytes()).hexdigest() if args.out.is_file() else None
            where = str(args.out)
        if other is None:
            print("\nnothing recorded at %s to compare with" % where)
            return 1
        same = other == digest
        print("\n%s %s" % ("identical to" if same else "DIFFERS from", where))
        return 0 if same else 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_bytes(body)
    print("wrote %s" % args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
