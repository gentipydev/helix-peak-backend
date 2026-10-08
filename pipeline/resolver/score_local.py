"""`scoring.score_with_esm` for one protein, as a process of the scorer's environment.

    pipeline/.esm-venv/bin/python -u -m pipeline.resolver.score_local \\
        target.pkl record.json track.json refusal.txt

Run from the repository root by `local_worker.py`, the resolver's worker on a
Mac. That worker needs psycopg and Biopython; the scorer's environment has
neither and is given neither, because it is the environment that baked the
twenty. So the worker runs beside it and hands over the two things the scorer
needs, the `Target` and the record's bytes, and this calls the function Modal
calls in the worker's own process.

How it ends is the whole protocol:

- exit 0: the track's bytes are in `track.json`.
- exit 3: the scorer raised `ValueError` (its alignment gate, its budget, a
  sequence that is not the record's) and the message is in `refusal.txt`. The
  worker raises it again, and `worker.score_next` tells it as a refusal, as it
  does where the scorer runs in its own process.
- anything else: it broke, stderr says how, and the worker tries again.

A bake that breaks three times is refused for good, so the worker first asks
whether this environment can score at all, and claims nothing while it cannot:

    pipeline/.esm-venv/bin/python -u -m pipeline.resolver.score_local --check

exits 0 where the scorer's packages import and its weights are in the Hugging
Face cache, and otherwise says on stderr which is missing.

The `Target` arrives pickled. Both ends are this repository on one machine,
and the file was written a moment ago by the worker into a directory of its
own, so nothing is unpickled that somebody else could have written. JSON would
mean rebuilding `Target`, `Source` and `Region` field for field, and a field
missed there scores a different protein without saying so.
"""

from __future__ import annotations

import os
import pickle
import sys
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.resolver import scoring  # noqa: E402

REFUSED = 3

# What `score_protein` loads, by the names Hugging Face stores them under: the
# files Modal's scorer image downloads (`modal_app.py`).
_WEIGHTS = ("config.json", "model.safetensors", "special_tokens_map.json",
            "tokenizer_config.json", "vocab.txt")


def check() -> int:
    """Whether a protein could be scored here, without scoring one."""
    import scipy  # noqa: F401
    import torch  # noqa: F401
    import transformers  # noqa: F401
    from huggingface_hub import try_to_load_from_cache

    from pipeline.constraint.score_protein import MODEL, REVISION

    missing = [name for name in _WEIGHTS
               if not isinstance(try_to_load_from_cache(MODEL, name, revision=REVISION), str)]
    if missing:
        print(f"{MODEL}@{REVISION} is not whole in the Hugging Face cache: "
              f"{', '.join(missing)} missing.", file=sys.stderr)
        return 1
    return 0


def main(arguments: list[str]) -> int:
    if arguments == ["--check"]:
        return check()
    if len(arguments) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    if not os.environ.get("HELIXPEEK_DATA"):
        # `paths.DATA` would be `pipeline/data/`, which holds the twenty's
        # stored tracks: the bytes every sha256 proof is read against.
        print("HELIXPEEK_DATA is not set, so the scorer would write under pipeline/data/. "
              "The worker names a directory of its own.", file=sys.stderr)
        return 2
    target_file, record_file, track_file, refusal_file = (Path(p) for p in arguments)
    target = pickle.loads(target_file.read_bytes())
    record = record_file.read_bytes()
    # Only the scorer's own `ValueError` is a refusal. One raised reading the
    # two files above is this process breaking, and is left to exit 1.
    try:
        track = scoring.score_with_esm(target, record)
    except ValueError as exc:
        refusal_file.write_text(str(exc), encoding="utf-8")
        return REFUSED
    track_file.write_bytes(track)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
