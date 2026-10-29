from __future__ import annotations

from typing import Iterable, List

from pipeline.denovo.alphabet import (
    BOS_ID, EOS_ID, SPECIALS, TOKEN_TO_ID, TOKENS,
)


class UnknownResidue(ValueError):
    pass


def encode(sequence: str, bos: bool = True, eos: bool = True) -> List[int]:
    ids = []
    for i, letter in enumerate(sequence):
        try:
            ids.append(TOKEN_TO_ID[letter])
        except KeyError:
            raise UnknownResidue(
                "position %d is %r, which is not one of the twenty" % (i, letter)
            ) from None
    if bos:
        ids.insert(0, BOS_ID)
    if eos:
        ids.append(EOS_ID)
    return ids


def decode(ids: Iterable[int], specials: bool = False) -> str:
    out = []
    for i in ids:
        token = TOKENS[i]
        if token in SPECIALS and not specials:
            continue
        out.append(token)
    return "".join(out)
