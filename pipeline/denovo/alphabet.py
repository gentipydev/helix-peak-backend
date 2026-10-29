from __future__ import annotations

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

PAD = "<pad>"
BOS = "<bos>"
EOS = "<eos>"
SPECIALS = (PAD, BOS, EOS)

TOKENS = SPECIALS + tuple(AMINO_ACIDS)
VOCAB_SIZE = len(TOKENS)

TOKEN_TO_ID = {token: i for i, token in enumerate(TOKENS)}
ID_TO_TOKEN = {i: token for token, i in TOKEN_TO_ID.items()}

PAD_ID = TOKEN_TO_ID[PAD]
BOS_ID = TOKEN_TO_ID[BOS]
EOS_ID = TOKEN_TO_ID[EOS]

MIN_LENGTH = 10
MAX_LENGTH = 50
MAX_CONTEXT = MAX_LENGTH + 1
