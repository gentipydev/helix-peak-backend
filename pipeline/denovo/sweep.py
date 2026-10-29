from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import torch

HERE = Path(__file__).resolve().parent
if str(HERE.parents[1]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[1]))

from pipeline.denovo.baselines import baselines
from pipeline.denovo.train import (
    Config, RUNS, best_device, corpus_split, run,
)

CORPUS = HERE / "data" / "peptides.txt"

SHAPES: Sequence[Tuple[int, int, int]] = (
    (64, 4, 256),
    (96, 4, 384),
    (128, 2, 512),
    (128, 3, 512),
    (128, 4, 512),
)
DROPOUTS = (0.1, 0.2, 0.3)
RATES = (3e-4, 1e-3, 3e-3)
SEEDS = (0, 1, 2)


def _key(config: Config) -> str:
    return json.dumps(asdict(config), sort_keys=True)


class Sweep:
    def __init__(self, train_sequences, validation_sequences, device, epochs: int):
        self.train_sequences = train_sequences
        self.validation_sequences = validation_sequences
        self.device = device
        self.epochs = epochs
        self.done: Dict[str, Dict] = {}
        self.rows: List[Dict] = []

    def measure(self, config: Config) -> Dict:
        config = replace(config, max_epochs=self.epochs)
        key = _key(config)
        if key in self.done:
            return self.done[key]
        result = run(
            config, self.train_sequences, self.validation_sequences,
            device=self.device, log=False,
        )
        result.pop("state_dict", None)
        result.pop("history", None)
        self.done[key] = result
        row = {
            "d_model": config.d_model,
            "layers": config.num_layers,
            "ff": config.dim_feedforward,
            "dropout": config.dropout,
            "lr": config.learning_rate,
            "seed": config.seed,
            "parameters": result["parameters"],
            "best_epoch": result["best_epoch"],
            "val_loss": result["best"]["loss"],
            "residue_loss": result["best"]["residue_loss"],
            "eos_loss": result["best"]["eos_loss"],
            "seconds": result["seconds"],
        }
        self.rows.append(row)
        print(
            "  d%-4d L%-2d ff%-4d drop%.1f lr%.0e seed%d | %7d params | "
            "val %.4f @ep%-3d | %5.1fs"
            % (row["d_model"], row["layers"], row["ff"], row["dropout"], row["lr"],
               row["seed"], row["parameters"], row["val_loss"], row["best_epoch"],
               row["seconds"]),
            flush=True,
        )
        return result

    def over_seeds(self, config: Config) -> Dict:
        losses, params = [], 0
        for seed in SEEDS:
            result = self.measure(replace(config, seed=seed))
            losses.append(result["best"]["loss"])
            params = result["parameters"]
        return {
            "config": config,
            "parameters": params,
            "mean": statistics.mean(losses),
            "stdev": statistics.stdev(losses) if len(losses) > 1 else 0.0,
            "losses": losses,
        }


def choose(candidates: List[Dict]) -> Dict:
    spreads = [c["stdev"] for c in candidates if c["stdev"] > 0]
    noise = statistics.mean(spreads) if spreads else 0.0
    best = min(candidates, key=lambda c: c["mean"])
    within = [c for c in candidates if c["mean"] <= best["mean"] + noise]
    return min(within, key=lambda c: c["parameters"])


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--device", default=None)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--held-out", type=float, default=0.1)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--out", type=Path, default=RUNS / "sweep")
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    device = best_device(args.device)
    train_sequences, validation_sequences, clusters = corpus_split(
        args.corpus, args.threshold, args.held_out, seed=0, device=device
    )
    reference = baselines(train_sequences, validation_sequences)
    print("corpus %d clusters | %d train, %d validation | device %s"
          % (clusters, len(train_sequences), len(validation_sequences), device.type))
    print("baselines: uniform %.4f, unigram %.4f, bigram %.4f\n"
          % (reference["uniform"], reference["unigram"], reference["bigram"]))

    sweep = Sweep(train_sequences, validation_sequences, device, args.epochs)
    started = time.perf_counter()

    print("stage A -- capacity (dropout 0.2, lr 1e-3)")
    stage_a = [
        sweep.over_seeds(Config(
            d_model=d, num_layers=layers, dim_feedforward=ff,
            dropout=0.2, learning_rate=1e-3,
        ))
        for d, layers, ff in SHAPES
    ]
    shape = choose(stage_a)["config"]
    print("  -> %d/%d/%d\n" % (shape.d_model, shape.num_layers, shape.dim_feedforward))

    print("stage B -- dropout")
    stage_b = [sweep.over_seeds(replace(shape, dropout=d)) for d in DROPOUTS]
    with_dropout = choose(stage_b)["config"]
    print("  -> dropout %.1f\n" % with_dropout.dropout)

    print("stage C -- learning rate")
    stage_c = [sweep.over_seeds(replace(with_dropout, learning_rate=r)) for r in RATES]
    winner = choose(stage_c)
    final = winner["config"]
    print("  -> lr %.0e\n" % final.learning_rate)

    elapsed = time.perf_counter() - started
    print("%d runs in %.1f minutes" % (len(sweep.rows), elapsed / 60))
    print("\nchosen: d_model %d, %d layers, ff %d, dropout %.1f, lr %.0e"
          % (final.d_model, final.num_layers, final.dim_feedforward,
             final.dropout, final.learning_rate))
    print("  %d parameters | validation %.4f +/- %.4f over %d seeds"
          % (winner["parameters"], winner["mean"], winner["stdev"], len(SEEDS)))
    print("  bigram is %.4f, so the transformer buys %.4f nats"
          % (reference["bigram"], reference["bigram"] - winner["mean"]))

    args.out.mkdir(parents=True, exist_ok=True)
    columns = ["d_model", "layers", "ff", "dropout", "lr", "seed", "parameters",
               "best_epoch", "val_loss", "residue_loss", "eos_loss", "seconds"]
    lines = ["\t".join(columns)]
    for row in sweep.rows:
        lines.append("\t".join(str(row[c]) for c in columns))
    (args.out / "runs.tsv").write_text("\n".join(lines) + "\n")

    def summarise(stage):
        return [{
            "d_model": c["config"].d_model, "layers": c["config"].num_layers,
            "ff": c["config"].dim_feedforward, "dropout": c["config"].dropout,
            "lr": c["config"].learning_rate, "parameters": c["parameters"],
            "mean": c["mean"], "stdev": c["stdev"], "losses": c["losses"],
        } for c in stage]

    (args.out / "sweep.json").write_text(json.dumps({
        "baselines": reference,
        "clusters": clusters,
        "train": len(train_sequences),
        "validation": len(validation_sequences),
        "epochs": args.epochs,
        "seeds": list(SEEDS),
        "stage_a_capacity": summarise(stage_a),
        "stage_b_dropout": summarise(stage_b),
        "stage_c_learning_rate": summarise(stage_c),
        "chosen": asdict(final),
        "chosen_mean": winner["mean"],
        "chosen_stdev": winner["stdev"],
        "minutes": elapsed / 60,
    }, indent=2) + "\n")
    print("\nwrote %s and %s" % (args.out / "runs.tsv", args.out / "sweep.json"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
