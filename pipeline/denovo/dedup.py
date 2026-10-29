from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import torch

HERE = Path(__file__).resolve().parent
if str(HERE.parents[1]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[1]))

from pipeline.denovo.alphabet import AMINO_ACIDS
from pipeline.denovo.sequences import read_peptides

RESIDUE = {letter: i for i, letter in enumerate(AMINO_ACIDS)}
DEFAULT_K = 3
DEFAULT_THRESHOLD = 0.5
DEFAULT_HELD_OUT = 0.1


def best_device() -> torch.device:
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def kmer_codes(sequence: str, k: int = DEFAULT_K) -> List[int]:
    codes, base = [], len(AMINO_ACIDS)
    for i in range(len(sequence) - k + 1):
        code = 0
        for letter in sequence[i : i + k]:
            code = code * base + RESIDUE[letter]
        codes.append(code)
    return codes


def presence(sequences: Sequence[str], k: int = DEFAULT_K, device=None) -> torch.Tensor:
    device = device or best_device()
    width = len(AMINO_ACIDS) ** k
    matrix = torch.zeros((len(sequences), width), dtype=torch.float16, device=device)
    for i, sequence in enumerate(sequences):
        codes = kmer_codes(sequence, k)
        if codes:
            matrix[i, torch.tensor(sorted(set(codes)), device=device)] = 1.0
    return matrix


def cluster(
    sequences: Sequence[str],
    threshold: float = DEFAULT_THRESHOLD,
    k: int = DEFAULT_K,
    device=None,
    block: int = 2048,
) -> List[int]:
    device = device or best_device()
    count = len(sequences)
    order = sorted(range(count), key=lambda i: (-len(sequences[i]), sequences[i]))
    ordered = [sequences[i] for i in order]

    matrix = presence(ordered, k=k, device=device)
    sizes = matrix.sum(dim=1)

    assigned = [-1] * count
    next_id = 0
    for start in range(0, count, block):
        stop = min(start + block, count)
        intersection = (matrix[start:stop] @ matrix.T).float()
        union = sizes[start:stop, None].float() + sizes[None, :].float() - intersection
        similarity = intersection / union.clamp(min=1.0)
        close = (similarity >= threshold).cpu()
        for row in range(stop - start):
            here = start + row
            if assigned[here] != -1:
                continue
            assigned[here] = next_id
            for other in torch.nonzero(close[row]).flatten().tolist():
                if other > here and assigned[other] == -1:
                    assigned[other] = next_id
            next_id += 1

    out = [-1] * count
    for position, original in enumerate(order):
        out[original] = assigned[position]
    return out


def split_by_cluster(
    sequences: Sequence[str],
    clusters: Sequence[int],
    held_out: float = DEFAULT_HELD_OUT,
    seed: int = 0,
) -> Tuple[List[str], List[str]]:
    members = {}
    for sequence, cluster_id in zip(sequences, clusters):
        members.setdefault(cluster_id, []).append(sequence)

    ids = sorted(members)
    shuffled = torch.randperm(len(ids), generator=torch.Generator().manual_seed(seed)).tolist()

    wanted = int(round(len(sequences) * held_out))
    validation, taken = [], 0
    held = set()
    for position in shuffled:
        if taken >= wanted:
            break
        cluster_id = ids[position]
        held.add(cluster_id)
        validation.extend(members[cluster_id])
        taken += len(members[cluster_id])

    train = [s for s, c in zip(sequences, clusters) if c not in held]
    return train, validation


def report(sequences: Sequence[str], clusters: Sequence[int]) -> dict:
    sizes = Counter(Counter(clusters).values())
    biggest = Counter(clusters).most_common(1)
    return {
        "sequences": len(sequences),
        "clusters": len(set(clusters)),
        "singletons": sizes.get(1, 0),
        "largest_cluster": biggest[0][1] if biggest else 0,
        "size_histogram": {str(size): n for size, n in sorted(sizes.items())},
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=HERE / "data" / "peptides.txt")
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("-k", type=int, default=DEFAULT_K)
    parser.add_argument("--held-out", type=float, default=DEFAULT_HELD_OUT)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--report", action="store_true")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    sequences = read_peptides(args.corpus)
    device = best_device()
    print("%d peptides, clustering on %s at Jaccard >= %.2f over %d-mers"
          % (len(sequences), device.type, args.threshold, args.k))

    clusters = cluster(sequences, threshold=args.threshold, k=args.k, device=device)
    summary = report(sequences, clusters)
    print("  %(clusters)d clusters, %(singletons)d singletons, largest holds %(largest_cluster)d"
          % summary)

    train, validation = split_by_cluster(
        sequences, clusters, held_out=args.held_out, seed=args.seed
    )
    print("  split: %d train, %d validation (%.1f%% held out)"
          % (len(train), len(validation), 100 * len(validation) / len(sequences)))

    by_cluster = dict(zip(sequences, clusters))
    shared = {by_cluster[s] for s in train} & {by_cluster[s] for s in validation}
    print("  clusters spanning the split: %d" % len(shared))
    assert not shared, "a cluster landed on both sides of the split"

    if args.report:
        print(json.dumps(summary["size_histogram"], indent=2))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps({
            "threshold": args.threshold, "k": args.k, "seed": args.seed,
            "summary": summary,
            "train": train, "validation": validation,
        }, indent=2) + "\n")
        print("wrote %s" % args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
