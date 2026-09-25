"""Bake masked-marginal scores offline; publish only after the gate passes.

The sequence is never typed in: it is read out of the protein the walk draws,
the stored record `assets/mock/gene_<gene>.json` under `pipeline/data/` (baked
there, or fetched by `fetch_tracks.py`), so a constraint track cannot come to
describe a different molecule from the one under it. Region boundaries and
disulfide pairs come from `pipeline/targets.py` and are written into the asset,
which is what lets the Dart side name a residue's domain and its bonding
partner without knowing which protein it is looking at.

    pipeline/.esm-venv/bin/python -u pipeline/constraint/score_protein.py --target lysozyme
    pipeline/.esm-venv/bin/python -u pipeline/constraint/score_protein.py --all
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import time

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import BY_SLUG, ESM_CONTEXT_RESIDUES, TARGETS, Target, partition  # noqa: E402

MODEL = "facebook/esm2_t33_650M_UR50D"
REVISION = "08e4846e537177426273712802403f7ba8261b6c"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# What the gate is actually for: a mask or context offset bug that files each
# score under the wrong residue. The token checks on every pass rule out most of
# those; this catches what gets past them, from the scores themselves.
#
# It used to be read off the disulfide cysteines, on the argument that they are
# the most constrained residues a protein has, and the rule was rewritten three
# times, each time because the data was right and the rule was not:
#
#   * "all within the top decile" cannot be satisfied at all by oxytocin,
#     whose precursor has sixteen bonded cysteines in 125 residues.
#   * "all within two places each" then rejected growth hormone over C215,
#     which ranks 28/217 while its other three rank 1, 2 and 4.
#   * "the median in the top decile, none past halfway" then rejected three of
#     the second ten. SOD1's bridge ranks 46 and 51 of 154, under the
#     histidines that hold its copper and zinc; leptin's ranks 24 and 33 of 167,
#     under its core leucines; amylase's C85 and C130 rank 449 and 265 of 511,
#     while its other eight bonded cysteines rank 1-41.
#
# A proxy that is biology can be beaten by biology, and it was no check at all
# on a protein with no disulfides. So alignment is now measured directly, on
# every protein: wherever a residue differs from the one beside it, the model
# has to score the residue that is there above its neighbour more often than
# not, on both sides. The shipped tracks measure 65-100%. A track shifted by one
# position scores each residue's neighbour as its own, which inverts that.
# Bonded cysteine ranks are still printed, as findings.
GATE_ALIGNMENT_FLOOR = 0.5

# Bytes. Generous for a small protein and binding for dystrophin, whose 3,685
# positions are most of what the app ships.
BUDGET_BYTES_PER_RESIDUE = 400
BUDGET_FLOOR = 100_000


def normalize_conservation(entropies: list[float]) -> list[float]:
    """Inverse min-max entropy, full vocabulary, with no domain adjustments."""
    if not entropies or any(not math.isfinite(h) or h < 0 for h in entropies):
        raise ValueError("Entropies must be finite and nonnegative")
    low, high = min(entropies), max(entropies)
    if high == low:
        raise ValueError("Flat entropy cannot produce a min-max constraint track")
    return [(high - h) / (high - low) for h in entropies]


def sequence_of(target: Target) -> str:
    """The protein the app draws, read back out of its own gene record."""
    record = json.loads((DATA / target.mock_asset).read_text())
    sequence = record["protein"]["translation"]
    if len(sequence) != target.aa:
        raise ValueError(
            f"{target.slug}: {target.mock_asset} holds {len(sequence)} residues, "
            f"the table says {target.aa}. Rebuild the gene record first."
        )
    return sequence


def alignment(sequence: str, positions: list[dict]) -> tuple[float, float]:
    """How often the residue that is there outscores the one before it, and the
    one after it.

    Counted only where the two differ, so a run of one amino acid — a polyQ
    tract, a proline-rich stretch — is evidence neither way.
    """
    fractions = []
    for step in (-1, 1):
        wins = total = 0
        for i, p in enumerate(positions):
            j = i + step
            if 0 <= j < len(sequence) and sequence[j] != sequence[i]:
                total += 1
                wins += p["substitutions"][sequence[j]] < 0
        if not total:
            raise ValueError("No residue differs from its neighbour; alignment cannot be measured")
        fractions.append(wins / total)
    return fractions[0], fractions[1]


def gate(target: Target, sequence: str, positions: list[dict]) -> None:
    """Report the shape of the track, and refuse to publish a misaligned one.

    Runs before any file is created. The universal half checks what the Dart
    loader will check anyway, early and with a better message. The alignment
    half is the one only this side can check: whether each score is filed under
    the residue it was measured for. Bonded cysteines are checked to be
    cysteines, and their ranks reported rather than gated; see
    `GATE_ALIGNMENT_FLOOR` for why.
    """
    if len(positions) != len(sequence):
        raise ValueError(f"{len(positions)} positions for {len(sequence)} residues")
    if "".join(p["wildtype"] for p in positions) != sequence:
        raise ValueError("Scored wildtypes are not the sequence")
    for p in positions:
        if not 0.0 <= p["conservation"] <= 1.0:
            raise ValueError(f"Conservation out of range at residue {p['index'] + 1}")
        if p["substitutions"][p["wildtype"]] != 0:
            raise ValueError(f"Wildtype does not score 0 at residue {p['index'] + 1}")

    print("\nRegional means (unmodified inverse min-max entropy):", flush=True)
    for region in partition(target):
        values = [p["conservation"] for p in positions[region["start"] - 1 : region["end"]]]
        print(
            f"  {region['label'][:26]:26} {region['start']:>4}-{region['end']:<4}: "
            f"{sum(values) / len(values):.6f}  ({len(values)} residues)"
        )

    ordered = sorted(positions, key=lambda p: p["entropy"])
    ranks = {p["index"]: rank for rank, p in enumerate(ordered, 1)}
    print("\nMost constrained positions:", flush=True)
    for p in ordered[:15]:
        print(f"  {p['wildtype']}{p['index'] + 1}: {p['conservation']:.6f}")

    bonded = sorted({n for pair in target.disulfides for n in pair})
    for number in bonded:
        if sequence[number - 1] != "C":
            raise ValueError(f"Residue {number} is {sequence[number - 1]}, not a cysteine")

    before, after = alignment(sequence, positions)
    print(
        f"\nAlignment gate (the residue that is there outscores its neighbour's, "
        f"where they differ; more than {GATE_ALIGNMENT_FLOOR:.0%} on each side):\n"
        f"  over the residue before: {before:.1%}\n"
        f"  over the residue after:  {after:.1%}",
        flush=True,
    )
    if before <= GATE_ALIGNMENT_FLOOR or after <= GATE_ALIGNMENT_FLOOR:
        raise ValueError(
            f"Alignment gate FAILED: {before:.1%} before, {after:.1%} after. No JSON "
            "written. Inspect token offsets and how scores are filed before any UI work."
        )
    print("Alignment gate PASSED.", flush=True)

    if not bonded:
        print("\nNo disulfides to report.", flush=True)
        return

    total = len(sequence)
    print(f"\nBonded cysteines ({len(bonded)}, reported rather than gated):", flush=True)
    for number in bonded:
        p = positions[number - 1]
        print(
            f"  C{number:<5d} conservation={p['conservation']:.6f} "
            f"entropy={p['entropy']:.6f} rank={ranks[number - 1]}/{total}",
            flush=True,
        )
    bonded_ranks = sorted(ranks[n - 1] for n in bonded)
    middle = len(bonded_ranks) // 2
    median = (
        bonded_ranks[middle]
        if len(bonded_ranks) % 2
        else (bonded_ranks[middle - 1] + bonded_ranks[middle]) / 2
    )
    print(f"  median rank {median}/{total}", flush=True)


def write_asset(payload: dict, target: Target, destination: Path) -> int:
    encoded = (json.dumps(payload, separators=(",", ":"), allow_nan=False) + "\n").encode()
    budget = max(BUDGET_FLOOR, len(payload["sequence"]) * BUDGET_BYTES_PER_RESIDUE)
    if len(encoded) >= budget:
        raise ValueError(f"{target.slug} asset is {len(encoded):,} bytes, over its {budget:,} budget")
    destination.parent.mkdir(parents=True, exist_ok=True)
    # An interrupted bake must not leave a partial asset behind.
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as temp:
        temp.write(encoded)
        temporary_path = temp.name
    try:
        os.chmod(temporary_path, 0o644)
        os.replace(temporary_path, destination)
    finally:
        Path(temporary_path).unlink(missing_ok=True)
    return len(encoded)


def refresh_metadata(target: Target, destination: Path) -> None:
    """Put the table's regions and disulfides back into an already-scored asset."""
    payload = json.loads(destination.read_text())
    if payload["sequence"] != sequence_of(target):
        raise ValueError(f"{target.slug}: the asset is for a different sequence; re-score it")
    payload["regions"] = partition(target)
    payload["disulfides"] = [list(pair) for pair in target.disulfides]
    size = write_asset(payload, target, destination)
    print(f"{target.slug:14} metadata refreshed: {len(payload['regions'])} regions, {size:,} B")


def score_protein(target: Target, output_path: Path, context: int) -> None:
    """Score every residue in its full context and atomically write a compact asset."""
    sequence = sequence_of(target)
    if not sequence or set(sequence) - set(AMINO_ACIDS):
        raise ValueError("Provide a nonempty sequence of uppercase canonical amino acids")

    import scipy
    from scipy.special import entr
    import torch
    import transformers
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    torch.manual_seed(0)
    device = "mps" if torch.backends.mps.is_available() else "cpu"
    print(f"Loading {MODEL}@{REVISION} on {device} (float32).", flush=True)
    tokenizer = AutoTokenizer.from_pretrained(MODEL, revision=REVISION)
    model = AutoModelForMaskedLM.from_pretrained(
        MODEL, revision=REVISION, dtype=torch.float32,
        attn_implementation="eager", use_safetensors=True,
    ).to(device).eval()
    aa_ids = {aa: tokenizer.convert_tokens_to_ids(aa) for aa in AMINO_ACIDS}
    if tokenizer(sequence[:8], return_tensors="pt")["input_ids"][0].tolist() != [
        tokenizer.cls_token_id, *[aa_ids[aa] for aa in sequence[:8]], tokenizer.eos_token_id
    ]:
        raise ValueError("Tokenizer does not map residues to one token each")

    # ESM-2 positions are rotary, so a longer sequence runs rather than raising
    # — it is just being asked about distances it never saw in training. A
    # protein past the trained context is therefore scored in windows centred
    # on the mask instead: every residue still gets the largest context the
    # model was actually trained to use, and the ones in the middle get it on
    # both sides. Only dystrophin and CFTR need this; every other protein takes
    # one window that is the whole of it, down the identical code path.
    width = min(len(sequence), context)
    windowed = width < len(sequence)
    if context + 2 > model.config.max_position_embeddings:
        raise ValueError(f"Window of {context} exceeds the model's {model.config.max_position_embeddings}")
    print(
        f"{target.slug}: {len(sequence)} residues, "
        f"{'windows of ' + str(width) if windowed else 'one full-length pass'}; "
        f"vocab={model.config.vocab_size}.",
        flush=True,
    )

    positions = []
    started = time.monotonic()
    with torch.inference_mode():
        for i, wildtype in enumerate(sequence):
            # Centred, then slid back inside the sequence at either end.
            start = min(max(i - width // 2, 0), len(sequence) - width)
            window = sequence[start : start + width]
            offset = i - start
            masked = window[:offset] + tokenizer.mask_token + window[offset + 1 :]
            inputs = tokenizer(masked, return_tensors="pt")
            ids = inputs["input_ids"][0]
            mask_indices = (ids == tokenizer.mask_token_id).nonzero(as_tuple=True)[0]
            if mask_indices.numel() != 1:
                raise ValueError(f"Expected exactly one mask at residue {i + 1}")
            mask_pos = mask_indices.item()
            expected = [
                tokenizer.cls_token_id, *[aa_ids[aa] for aa in window], tokenizer.eos_token_id
            ]
            expected[offset + 1] = tokenizer.mask_token_id
            if mask_pos != offset + 1 or ids.tolist() != expected:
                raise ValueError(f"Mask/context token offset mismatch at residue {i + 1}")
            logits = model(**inputs.to(device)).logits[0, mask_pos]
            # Float64 CPU reduction makes tiny entropy differences auditable.
            log_probs = torch.log_softmax(logits.cpu().double(), dim=-1)
            if not torch.isfinite(log_probs).all():
                raise ValueError(f"Nonfinite model output at residue {i + 1}")
            probabilities = log_probs.exp().numpy()
            entropy = float(entr(probabilities).sum())  # All 33 vocabulary tokens.
            wt_lp = log_probs[aa_ids[wildtype]].item()
            positions.append({
                "index": i,
                "wildtype": wildtype,
                "entropy": entropy,
                "substitutions": {
                    aa: round(log_probs[token_id].item() - wt_lp, 3)
                    for aa, token_id in aa_ids.items()
                },
            })
            if (i + 1) % 25 == 0 or i + 1 == len(sequence):
                elapsed = time.monotonic() - started
                rate = elapsed / (i + 1)
                print(
                    f"Scored {i + 1}/{len(sequence)} ({elapsed:.1f}s, "
                    f"{rate * (len(sequence) - i - 1) / 60:.1f} min left)",
                    flush=True,
                )

    for p, conservation in zip(positions, normalize_conservation([p["entropy"] for p in positions])):
        p["conservation"] = conservation
    gate(target, sequence, positions)  # Deliberately before any output file is created.

    payload = {
        "gene": target.gene,
        "sequence": sequence,
        "uniprot": target.uniprot,
        "sequence_source": target.mock_asset,
        "model": MODEL,
        "revision": REVISION,
        "method": "masked_marginals",
        "context": {"mode": "window" if windowed else "full", "residues": width},
        "normalization": "minmax",
        "conservation_definition": "(max_entropy - entropy) / (max_entropy - min_entropy)",
        "entropy_vocabulary": "full",
        "vocabulary_size": model.config.vocab_size,
        "score_units": "natural_log_ratio_to_wildtype",
        "regions": partition(target),
        "disulfides": [list(pair) for pair in target.disulfides],
        "generation": {
            "device": device, "dtype": "float32", "torch": torch.__version__,
            "transformers": transformers.__version__, "scipy": scipy.__version__,
        },
        "positions": [
            {**p, "entropy": round(p["entropy"], 9), "conservation": round(p["conservation"], 6)}
            for p in positions
        ],
    }
    size = write_asset(payload, target, Path(output_path))
    print(f"Wrote {output_path}: {len(positions)} positions, {size:,} bytes.\n", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", action="append", help="slug; repeatable")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--context", type=int, default=ESM_CONTEXT_RESIDUES)
    parser.add_argument("--output", default=None, help="Override, for a one-off bake")
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Rewrite the region and disulfide tables of an existing asset without "
             "re-scoring. They come from the table, not the model, so a change to a "
             "region boundary does not need four hours of GPU to land.",
    )
    args = parser.parse_args()

    if args.target:
        chosen = [BY_SLUG[slug] for slug in args.target]
        # Named on purpose, so refused out loud: a track for a row the table
        # calls unscored is one `check_assets.py` would then reject.
        unscored = [t.slug for t in chosen if not t.scored]
        if unscored:
            raise SystemExit(f"{unscored} are not scored in targets.py; set scored=True first")
    elif args.all:
        # Every protein the table says is scored. The rest are not missing a
        # track; they do not have one yet, and the app draws them that way.
        chosen = [t for t in TARGETS if t.scored]
    else:
        raise SystemExit("Pass --target <slug> or --all")

    for target in chosen:
        destination = Path(args.output or DATA / target.constraint_asset)
        if args.metadata_only:
            refresh_metadata(target, destination)
        else:
            score_protein(target, destination, args.context)
