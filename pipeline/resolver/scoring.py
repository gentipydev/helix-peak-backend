"""ESM-2 over one protein's record: the one step of the worker that needs torch.

Kept apart from `worker.py` so that an environment holding only the scorer can
import it: `worker.py` imports the store, the resolver and Biopython, and
`pipeline/.esm-venv`, which baked the twenty, has none of them. This imports
`pipeline.targets` and, when called, the scorer; nothing else may be added.
`worker.score_with_esm` is this function, and `score_local.py` is how that
environment is asked for it.
"""

from __future__ import annotations

from pipeline.targets import Target


def score_with_esm(target: Target, record: bytes) -> bytes:
    """The constraint track, from the unchanged scorer, for the record given.

    The scorer reads the protein out of the record where the bakes keep it
    (`paths.DATA`) and writes its track beside it; both are put there and read
    back here. Torch is imported inside the scorer, so this module is cheap to
    import where no GPU is.
    """
    from pipeline.constraint import score_protein
    from pipeline.targets import ESM_CONTEXT_RESIDUES

    source = score_protein.DATA / target.mock_asset
    source.parent.mkdir(parents=True, exist_ok=True)
    source.write_bytes(record)
    destination = score_protein.DATA / target.constraint_asset
    score_protein.score_protein(target, destination, ESM_CONTEXT_RESIDUES)
    return destination.read_bytes()
