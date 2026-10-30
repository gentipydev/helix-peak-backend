from __future__ import annotations

import argparse
import copy
import json
import math
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import torch
import torch.nn.functional as F

HERE = Path(__file__).resolve().parent
if str(HERE.parents[1]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[1]))

from pipeline.denovo.alphabet import EOS_ID, PAD_ID, VOCAB_SIZE
from pipeline.denovo.bigram_torch import masked
from pipeline.denovo.dedup import cluster, split_by_cluster
from pipeline.denovo.helix_peek_LM import HelixPeekLM, parameter_count
from pipeline.denovo.sequences import batches, build_dataset, read_peptides

CORPUS = HERE / "data" / "peptides.txt"
RUNS = HERE / "runs"


@dataclass
class Config:
    d_model: int = 96
    nhead: int = 4
    num_layers: int = 4
    dim_feedforward: int = 384
    dropout: float = 0.2
    learning_rate: float = 1e-3
    weight_decay: float = 0.1
    eta_min: float = 1e-5
    warmup_steps: int = 100
    batch_size: int = 64
    max_epochs: int = 200
    patience: int = 10
    grad_clip: float = 1.0
    seed: int = 0


def best_device(name: Optional[str] = None) -> torch.device:
    if name:
        return torch.device(name)
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def param_groups(model: torch.nn.Module, weight_decay: float) -> List[Dict]:
    decay, plain = [], []
    for name, parameter in model.named_parameters():
        if not parameter.requires_grad:
            continue
        if parameter.ndim >= 2 and "embedding" not in name:
            decay.append(parameter)
        else:
            plain.append(parameter)
    return [
        {"params": decay, "weight_decay": weight_decay},
        {"params": plain, "weight_decay": 0.0},
    ]


def lr_schedule(config: Config, total_steps: int):
    warmup = max(1, config.warmup_steps)
    floor = config.eta_min / config.learning_rate

    def multiplier(step: int) -> float:
        if step < warmup:
            return (step + 1) / warmup
        progress = (step - warmup) / max(1, total_steps - warmup)
        progress = min(1.0, progress)
        return floor + (1.0 - floor) * 0.5 * (1.0 + math.cos(math.pi * progress))

    return multiplier


def _token_losses(logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
    return F.cross_entropy(
        masked(logits.reshape(-1, VOCAB_SIZE)),
        targets.reshape(-1),
        reduction="none",
        ignore_index=PAD_ID,
    )


def build(
    config: Config, total_steps: int, device: torch.device
) -> Tuple[HelixPeekLM, torch.optim.Optimizer, torch.optim.lr_scheduler.LambdaLR]:
    model = HelixPeekLM(
        d_model=config.d_model,
        nhead=config.nhead,
        num_layers=config.num_layers,
        dim_feedforward=config.dim_feedforward,
        dropout=config.dropout,
    ).to(device)
    optimizer = torch.optim.AdamW(
        param_groups(model, config.weight_decay), lr=config.learning_rate
    )
    scheduler = torch.optim.lr_scheduler.LambdaLR(
        optimizer, lr_schedule(config, total_steps)
    )
    return model, optimizer, scheduler


def train_step(
    model: torch.nn.Module,
    optimizer: torch.optim.Optimizer,
    scheduler: torch.optim.lr_scheduler.LambdaLR,
    x: torch.Tensor,
    y: torch.Tensor,
    grad_clip: float,
    after_backward: Optional[Callable[[], None]] = None,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """`after_backward` sees the gradients before they are clipped and the
    weights before they move; `trace.py` records them there."""
    per_token = _token_losses(model(x), y)
    real = (y.reshape(-1) != PAD_ID)
    tokens = real.sum()
    loss = (per_token * real).sum() / tokens.clamp(min=1)

    optimizer.zero_grad(set_to_none=True)
    loss.backward()
    if after_backward is not None:
        after_backward()
    norm = torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
    optimizer.step()
    scheduler.step()
    return loss.detach(), tokens, norm


@torch.no_grad()
def evaluate(
    model: torch.nn.Module,
    frames: torch.Tensor,
    lengths: torch.Tensor,
    batch_size: int,
    device: torch.device,
) -> Dict[str, float]:
    model.eval()
    totals = torch.zeros(3, device=device)
    counts = torch.zeros(3, device=device)
    for x, y in batches(frames, lengths, batch_size):
        x, y = x.to(device), y.to(device)
        per_token = _token_losses(model(x), y)
        flat = y.reshape(-1)
        real = flat != PAD_ID
        ends = flat == EOS_ID
        residues = real & ~ends
        for i, mask in enumerate((real, residues, ends)):
            totals[i] += (per_token * mask).sum()
            counts[i] += mask.sum()
    model.train()
    mean = (totals / counts.clamp(min=1)).tolist()
    return {"loss": mean[0], "residue_loss": mean[1], "eos_loss": mean[2]}


def run(
    config: Config,
    train_sequences: Sequence[str],
    validation_sequences: Sequence[str],
    device: Optional[torch.device] = None,
    log: bool = True,
) -> Dict:
    device = device or best_device()
    torch.manual_seed(config.seed)

    train_frames, train_lengths = build_dataset(train_sequences)
    val_frames, val_lengths = build_dataset(validation_sequences)

    steps_per_epoch = math.ceil(len(train_sequences) / config.batch_size)
    model, optimizer, scheduler = build(
        config, steps_per_epoch * config.max_epochs, device
    )

    if log:
        print(
            "%s | %d params | %d train, %d validation | %d steps/epoch"
            % (device.type, parameter_count(model), len(train_sequences),
               len(validation_sequences), steps_per_epoch)
        )

    history: List[Dict] = []
    best = {"loss": float("inf")}
    best_state, best_epoch, waited = None, 0, 0
    started = time.perf_counter()

    for epoch in range(1, config.max_epochs + 1):
        generator = torch.Generator().manual_seed(config.seed * 10_000 + epoch)
        loss_sum = torch.zeros((), device=device)
        token_sum = torch.zeros((), device=device)

        for x, y in batches(train_frames, train_lengths, config.batch_size, generator):
            x, y = x.to(device), y.to(device)
            loss, tokens, _ = train_step(
                model, optimizer, scheduler, x, y, config.grad_clip
            )
            loss_sum += loss * tokens
            token_sum += tokens

        train_loss = (loss_sum / token_sum.clamp(min=1)).item()
        scores = evaluate(model, val_frames, val_lengths, config.batch_size, device)
        history.append({"epoch": epoch, "train_loss": train_loss,
                        "lr": scheduler.get_last_lr()[0], **scores})

        if log:
            print(
                "epoch %3d | train %.4f | val %.4f (residues %.4f, eos %.4f) | lr %.2e"
                % (epoch, train_loss, scores["loss"], scores["residue_loss"],
                   scores["eos_loss"], scheduler.get_last_lr()[0])
            )

        if scores["loss"] < best["loss"] - 1e-4:
            best, best_epoch, waited = scores, epoch, 0
            best_state = copy.deepcopy({k: v.cpu() for k, v in model.state_dict().items()})
        else:
            waited += 1
            if waited >= config.patience:
                if log:
                    print("stopped at epoch %d; best was epoch %d" % (epoch, best_epoch))
                break

    if best_state is None:
        raise RuntimeError(
            "training never improved on an infinite loss -- the first epoch scored %r"
            % (history[0]["loss"] if history else None)
        )

    return {
        "config": asdict(config),
        "device": device.type,
        "parameters": parameter_count(model),
        "best_epoch": best_epoch,
        "best": best,
        "epochs_run": len(history),
        "seconds": time.perf_counter() - started,
        "history": history,
        "state_dict": best_state,
    }


def corpus_split(
    path: Path,
    threshold: float,
    held_out: float,
    seed: int,
    device,
    random_split: bool = False,
):
    sequences = read_peptides(path)
    clusters = cluster(sequences, threshold=threshold, device=device)
    if random_split:
        from pipeline.denovo.bigram import split

        held = int(round(len(sequences) * held_out))
        train, validation = split(sequences, held_back=held, seed=seed)
    else:
        train, validation = split_by_cluster(
            sequences, clusters, held_out=held_out, seed=seed
        )
    return train, validation, len(set(clusters))


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--tag", default="default")
    parser.add_argument("--device", default=None)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--held-out", type=float, default=0.1)
    parser.add_argument("--random-split", action="store_true",
                        help="split by sequence rather than by cluster, to show what leaks")
    for field, value in asdict(Config()).items():
        parser.add_argument("--" + field.replace("_", "-"), type=type(value), default=value)
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    config = Config(**{f: getattr(args, f) for f in asdict(Config())})
    device = best_device(args.device)

    train_sequences, validation_sequences, clusters = corpus_split(
        args.corpus, args.threshold, args.held_out, config.seed, device,
        random_split=args.random_split,
    )
    print("corpus %s: %d clusters, split %s"
          % (args.corpus, clusters,
             "at random by sequence" if args.random_split else "by cluster"))

    result = run(config, train_sequences, validation_sequences, device=device)

    out = RUNS / args.tag
    out.mkdir(parents=True, exist_ok=True)
    torch.save(result.pop("state_dict"), out / "best.pt")
    (out / "history.json").write_text(json.dumps(result, indent=2) + "\n")
    print("\nbest validation %.4f at epoch %d, in %.1fs"
          % (result["best"]["loss"], result["best_epoch"], result["seconds"]))
    print("wrote %s" % (out / "best.pt"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
