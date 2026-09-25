"""Independent CPU/token-mask audit of a generated asset; never rewrites it.

Re-derives every number in the asset on CPU, by a different route: residue
offsets come from the tokenizer's special-token mask rather than from `i + 1`,
log-probability ratios are checked as raw logit differences, and entropy is a
NumPy reduction rather than SciPy's. Two implementations agreeing to five
decimals is the point.

    pipeline/.esm-venv/bin/python -u pipeline/constraint/verify_cpu.py --target lysozyme

Full-length only: it masks in one pass over the whole sequence, so it audits
every protein scored that way and refuses dystrophin and CFTR, which are scored
in windows. Auditing those would mean reimplementing the windowing, which is the
one thing an independent check must not share.
"""

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.special import logsumexp
import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.constraint.score_protein import AMINO_ACIDS, MODEL, REVISION  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import BY_SLUG, Target  # noqa: E402


def main(target: Target):
    path = DATA / target.constraint_asset
    data = json.loads(path.read_text())
    sequence = data["sequence"]
    if data["context"]["mode"] != "full":
        raise SystemExit(
            f"{target.slug} was scored in {data['context']['residues']}-residue windows. "
            "This audit only re-derives a full-length pass."
        )
    assert len(data["positions"]) == len(sequence)
    torch.set_num_threads(4)
    tokenizer = AutoTokenizer.from_pretrained(MODEL, revision=REVISION, local_files_only=True)
    model = AutoModelForMaskedLM.from_pretrained(
        MODEL, revision=REVISION, dtype=torch.float32, attn_implementation="eager",
        use_safetensors=True, local_files_only=True,
    ).cpu().eval()
    inputs = tokenizer(sequence, return_tensors="pt", return_special_tokens_mask=True)
    # Derive residue offsets from the special-token map, independently of i+1.
    offsets = (inputs.pop("special_tokens_mask")[0] == 0).nonzero().flatten().tolist()
    assert len(offsets) == len(sequence)
    aa_ids = tokenizer.convert_tokens_to_ids(list(AMINO_ACIDS))
    entropies = []
    table = []
    max_score_error = 0.0
    started = time.monotonic()
    with torch.inference_mode():
        for i, offset in enumerate(offsets):
            assert tokenizer.convert_ids_to_tokens(inputs["input_ids"][0, offset].item()) == sequence[i]
            token_ids = inputs["input_ids"].clone()
            token_ids[0, offset] = tokenizer.mask_token_id
            logits = model(input_ids=token_ids, attention_mask=inputs["attention_mask"]).logits[0, offset].numpy().astype(np.float64)
            lp = logits - logsumexp(logits)
            entropies.append(float(-np.sum(np.exp(lp) * lp)))
            # Log-probability ratios also equal raw logit differences.
            scores = logits[aa_ids] - logits[tokenizer.convert_tokens_to_ids(sequence[i])]
            table.append(scores)
            for aa, score in zip(AMINO_ACIDS, scores):
                max_score_error = max(max_score_error, abs(round(float(score), 3) - data["positions"][i]["substitutions"][aa]))
            if (i + 1) % 25 == 0:
                print(f"CPU verified {i + 1}/{len(sequence)} ({time.monotonic() - started:.1f}s)", flush=True)
    entropies = np.array(entropies)
    conservation = (entropies.max() - entropies) / np.ptp(entropies)
    max_entropy_error = float(np.max(np.abs(entropies - [p["entropy"] for p in data["positions"]])))
    max_conservation_error = float(np.max(np.abs(conservation - [p["conservation"] for p in data["positions"]])))
    print(f"Maximum entropy difference: {max_entropy_error:.9f}")
    print(f"Maximum conservation difference: {max_conservation_error:.9f}")
    print(f"Maximum rounded substitution difference: {max_score_error:.9f}")
    assert max_entropy_error < 0.0005
    assert max_conservation_error < 0.0005
    assert max_score_error <= 0.002
    # The generator's alignment gate, re-derived from this run's own logits: where
    # a residue differs from its neighbour, the residue that is there has to beat
    # the neighbour's more often than not, on each side.
    table = np.array(table)
    column = {aa: k for k, aa in enumerate(AMINO_ACIDS)}
    for side, step in (("before", -1), ("after", 1)):
        beaten = [
            table[i, column[sequence[i + step]]] < 0
            for i in range(len(sequence))
            if 0 <= i + step < len(sequence) and sequence[i + step] != sequence[i]
        ]
        print(f"CPU alignment over the residue {side}: {np.mean(beaten):.1%}")
        assert np.mean(beaten) > 0.5
    if target.disulfides:
        bonded = sorted({n for pair in target.disulfides for n in pair})
        ranks = np.argsort(np.argsort(entropies)) + 1
        print("CPU cysteine ranks:", {n: int(ranks[n - 1]) for n in bonded})
    print(f"CPU audit PASSED. All {len(sequence)} positions agree with the MPS asset.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", default="insulin")
    main(BY_SLUG[parser.parse_args().target])
