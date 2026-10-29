from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Dict, Optional, Sequence

HERE = Path(__file__).resolve().parent
if str(HERE.parents[1]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[1]))

from pipeline.denovo.alphabet import VOCAB_SIZE
from pipeline.denovo.bigram import (
    ALLOWED, adjacent, counts, cross_entropy, probabilities,
)
from pipeline.denovo.sequences import read_peptides

CORPUS = HERE / "data" / "peptides.txt"
SMOOTHING = 0.1


def uniform_table() -> list:
    even = 1.0 / len(ALLOWED)
    row = [even if i in ALLOWED else 0.0 for i in range(VOCAB_SIZE)]
    return [list(row) for _ in range(VOCAB_SIZE)]


def unigram_table(sequences: Sequence[str], smoothing: float = SMOOTHING) -> list:
    tally = [0] * VOCAB_SIZE
    for _, second in adjacent(sequences):
        tally[second] += 1
    total = sum(tally[j] for j in ALLOWED) + smoothing * len(ALLOWED)
    row = [
        (tally[j] + smoothing) / total if j in ALLOWED else 0.0
        for j in range(VOCAB_SIZE)
    ]
    return [list(row) for _ in range(VOCAB_SIZE)]


def baselines(
    train: Sequence[str],
    validation: Sequence[str],
    smoothing: float = SMOOTHING,
) -> Dict[str, float]:
    return {
        "uniform": cross_entropy(uniform_table(), validation),
        "unigram": cross_entropy(unigram_table(train, smoothing), validation),
        "bigram": cross_entropy(
            probabilities(counts(train), smoothing=smoothing), validation
        ),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--held-out", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--random-split", action="store_true",
                        help="split by sequence instead of by cluster, to show the difference")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    sequences = read_peptides(args.corpus)

    if args.random_split:
        from pipeline.denovo.bigram import split
        held = int(round(len(sequences) * args.held_out))
        train, validation = split(sequences, held_back=held, seed=args.seed)
        how = "random, by sequence"
    else:
        from pipeline.denovo.dedup import cluster, split_by_cluster
        clusters = cluster(sequences, threshold=args.threshold)
        train, validation = split_by_cluster(
            sequences, clusters, held_out=args.held_out, seed=args.seed
        )
        how = "by cluster at Jaccard >= %.2f" % args.threshold

    scores = baselines(train, validation)
    print("%d peptides, split %s: %d train, %d validation"
          % (len(sequences), how, len(train), len(validation)))
    print()
    print("%-10s %12s %12s" % ("model", "nats/token", "perplexity"))
    for name, value in scores.items():
        print("%-10s %12.4f %12.2f" % (name, value, math.exp(value)))

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps({
            "split": how, "seed": args.seed,
            "train": len(train), "validation": len(validation),
            "baselines": scores,
        }, indent=2) + "\n")
        print("\nwrote %s" % args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
