from __future__ import annotations

import sys
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[3]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.constraint.score_protein import AMINO_ACIDS as SCORED
from pipeline.denovo import alphabet


def test_ids_are_pinned():
    assert alphabet.TOKENS[:3] == ("<pad>", "<bos>", "<eos>")
    assert "".join(alphabet.TOKENS[3:]) == "ACDEFGHIKLMNPQRSTVWY"
    assert (alphabet.PAD_ID, alphabet.BOS_ID, alphabet.EOS_ID) == (0, 1, 2)
    assert alphabet.VOCAB_SIZE == 23
    assert all(alphabet.TOKENS[i] == t for t, i in alphabet.TOKEN_TO_ID.items())


def test_the_twenty_are_the_scorers():
    assert alphabet.AMINO_ACIDS == SCORED
