from __future__ import annotations

from typing import Iterable, Iterator, List, Optional, Sequence, Tuple

import torch
from torch.nn.utils.rnn import pad_sequence

from pipeline.denovo.alphabet import (
    AMINO_ACIDS, MAX_LENGTH, MIN_LENGTH, PAD_ID,
)
from pipeline.denovo.bigram import split
from pipeline.denovo.tokenizer import encode

CANONICAL = frozenset(AMINO_ACIDS)

FRAME_WIDTH = MAX_LENGTH + 2


def is_modelable(sequence: str) -> bool:
    return (
        MIN_LENGTH <= len(sequence) <= MAX_LENGTH
        and CANONICAL.issuperset(sequence)
    )


def clean(sequences: Iterable[str]) -> List[str]:
    seen, kept = set(), []
    for sequence in sequences:
        sequence = sequence.strip().upper()
        if not is_modelable(sequence) or sequence in seen:
            continue
        seen.add(sequence)
        kept.append(sequence)
    return kept


def read_fasta(path) -> List[str]:
    sequences, current = [], []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                    current = []
            else:
                current.append(line.strip())
    if current:
        sequences.append("".join(current))
    return sequences


def read_peptides(path) -> List[str]:
    with open(path) as handle:
        return [
            line.strip()
            for line in handle
            if line.strip() and not line.startswith("#")
        ]


def load_and_process_sequences(file_path, min_len=MIN_LENGTH, max_len=MAX_LENGTH) -> List[str]:
    sequences = [
        s for s in read_fasta(file_path)
        if min_len <= len(s.strip()) <= max_len
    ]
    processed = clean(sequences)
    print(
        "Loaded %d clean, natural peptide sequences from %s."
        % (len(processed), file_path)
    )
    return processed


def load_fake_sequences() -> List[str]:
    hardcoded_peptides = [
        "GLFDIVKKVVGAFGSIGKFLHSAKK",
        "KWKLFKKIPKFLHAKKF",
        "KWKLFKKIPKFLHAKKFAKAFV",
        "GLLDIVKKVVGAIGSIFKFLKK",
        "VNWKKILGKIIKVVK",
        "VNWKKILGKIIKVVKALAKFVAL",
        "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES",
    ]
    return hardcoded_peptides


def build_dataset(sequences: Sequence[str]) -> Tuple[torch.Tensor, torch.Tensor]:
    frames = torch.full((len(sequences), FRAME_WIDTH), PAD_ID, dtype=torch.long)
    lengths = torch.empty(len(sequences), dtype=torch.long)
    for i, sequence in enumerate(sequences):
        ids = encode(sequence, bos=True, eos=True)
        frames[i, : len(ids)] = torch.tensor(ids, dtype=torch.long)
        lengths[i] = len(ids)
    return frames, lengths


def batches(
    frames: torch.Tensor,
    lengths: torch.Tensor,
    batch_size: int,
    generator: Optional[torch.Generator] = None,
    bucket: bool = True,
) -> Iterator[Tuple[torch.Tensor, torch.Tensor]]:
    count = frames.shape[0]
    order = (
        torch.randperm(count, generator=generator)
        if generator is not None
        else torch.arange(count)
    )
    if bucket:
        order = order[torch.argsort(lengths[order], stable=True)]

    chunks = [order[i : i + batch_size] for i in range(0, count, batch_size)]
    if generator is not None:
        shuffled = torch.randperm(len(chunks), generator=generator).tolist()
        chunks = [chunks[i] for i in shuffled]

    for chunk in chunks:
        width = int(lengths[chunk].max())
        frame = frames[chunk, :width]
        yield frame[:, :-1], frame[:, 1:]


def prepare_padded_sequence_batch(sequences: Sequence[str]) -> Tuple[torch.Tensor, torch.Tensor]:
    batch_x, batch_y = [], []
    for seq in sequences:
        ids = encode(seq, bos=True, eos=True)
        batch_x.append(torch.tensor(ids[:-1], dtype=torch.long))
        batch_y.append(torch.tensor(ids[1:], dtype=torch.long))
    x_padded = pad_sequence(batch_x, batch_first=True, padding_value=PAD_ID)
    y_padded = pad_sequence(batch_y, batch_first=True, padding_value=PAD_ID)
    return x_padded, y_padded
