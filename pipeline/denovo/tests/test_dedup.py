from __future__ import annotations

import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[3]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

torch = pytest.importorskip("torch", reason="pipeline/denovo/venv has torch; .venv does not")

from pipeline.denovo.alphabet import AMINO_ACIDS
from pipeline.denovo.dedup import (
    cluster, kmer_codes, presence, report, split_by_cluster,
)

TWINS = [
    "KWKLFKKIPKFLHAKKFAKAF",
    "KWKLFKKIPKFLHAKKFAKAG",
    "KWKLFKKIPKFLHAKKFAKAW",
]
STRANGERS = ["GGGGGPPPPPDDDDDEEEEEQ", "GLFDIVKKVVGAFGSIGKFLH"]
CPU = torch.device("cpu")


def jaccard(a: str, b: str, k: int = 3) -> float:
    first, second = set(kmer_codes(a, k)), set(kmer_codes(b, k))
    return len(first & second) / len(first | second)


def test_kmer_codes_are_positional_base_twenty():
    base = len(AMINO_ACIDS)
    assert kmer_codes("AAA") == [0]
    assert kmer_codes("AAC") == [1]
    assert kmer_codes("ACA") == [base]
    assert len(kmer_codes("A" * 10)) == 8


def test_presence_marks_each_kmer_once():
    matrix = presence(["AAAA"], device=CPU)
    assert matrix.sum().item() == 1.0
    assert matrix[0, 0].item() == 1.0


def test_the_matmul_jaccard_matches_the_sets():
    sequences = TWINS + STRANGERS
    matrix = presence(sequences, device=CPU)
    sizes = matrix.sum(dim=1)
    intersection = (matrix @ matrix.T).float()
    union = sizes[:, None].float() + sizes[None, :].float() - intersection
    computed = intersection / union.clamp(min=1.0)
    for i, a in enumerate(sequences):
        for j, b in enumerate(sequences):
            assert computed[i][j].item() == pytest.approx(jaccard(a, b), abs=1e-5)


def test_near_twins_land_in_one_cluster():
    assigned = cluster(TWINS + STRANGERS, threshold=0.5, device=CPU)
    assert len(set(assigned[:3])) == 1
    assert len(set(assigned)) == 3


def test_a_stricter_threshold_never_makes_fewer_clusters():
    sequences = TWINS + STRANGERS
    loose = len(set(cluster(sequences, threshold=0.3, device=CPU)))
    tight = len(set(cluster(sequences, threshold=0.9, device=CPU)))
    assert tight >= loose


def test_blocking_does_not_change_the_clusters():
    sequences = (TWINS + STRANGERS) * 3
    whole = cluster(sequences, threshold=0.5, device=CPU, block=4096)
    blocked = cluster(sequences, threshold=0.5, device=CPU, block=2)
    assert whole == blocked


def test_no_cluster_spans_the_split():
    sequences = TWINS + STRANGERS
    assigned = cluster(sequences, threshold=0.5, device=CPU)
    train, validation = split_by_cluster(sequences, assigned, held_out=0.4, seed=0)
    assert set(train) | set(validation) == set(sequences)
    assert not set(train) & set(validation)
    where = dict(zip(sequences, assigned))
    assert not {where[s] for s in train} & {where[s] for s in validation}


def test_report_counts_what_it_says():
    sequences = TWINS + STRANGERS
    summary = report(sequences, cluster(sequences, threshold=0.5, device=CPU))
    assert summary["sequences"] == 5
    assert summary["clusters"] == 3
    assert summary["singletons"] == 2
    assert summary["largest_cluster"] == 3
