from __future__ import annotations

import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[3]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.denovo.alphabet import AMINO_ACIDS, BOS_ID, EOS_ID
from pipeline.denovo.tokenizer import UnknownResidue, decode, encode

INSULIN = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"


def test_framing_is_on_by_default():
    ids = encode(INSULIN)
    assert ids[0] == BOS_ID
    assert ids[-1] == EOS_ID
    assert len(ids) == len(INSULIN) + 2


def test_framing_can_be_left_off():
    assert encode(INSULIN, bos=False, eos=False) == encode(INSULIN)[1:-1]


def test_round_trip():
    assert decode(encode(INSULIN)) == INSULIN
    assert decode(encode(AMINO_ACIDS)) == AMINO_ACIDS
    assert decode(encode("")) == ""


def test_specials_are_dropped_unless_asked_for():
    assert decode(encode(INSULIN), specials=True) == "<bos>" + INSULIN + "<eos>"


def test_a_letter_outside_the_twenty_raises():
    with pytest.raises(UnknownResidue) as caught:
        encode("MALUWM")
    assert "position 3" in str(caught.value)

    for letter in "XBZUO.-":
        with pytest.raises(UnknownResidue):
            encode(letter)
