"""The resolver on Modal: requests resolved on a CPU, ESM-2 scored on a GPU.

From the repository root, with the Modal CLI logged in to the workspace
(`pip install -r pipeline/resolver/requirements.txt`, then `modal token new`):

    modal deploy pipeline/resolver/modal_app.py

It reads one secret from that workspace, `helix-peak-resolver`: the uploader's
three credentials (DATABASE_URL, SUPABASE_URL, SUPABASE_SERVICE_KEY) and NCBI's
contact address (NCBI_EMAIL). `pipeline/resolver/README.md` has the rest of the
runbook.

Three functions, and what starts each:

- `sweep` resolves queued requests (`worker.sweep`) and, where constraint bakes
  wait, starts `score`. A schedule runs it every five minutes, so a request is
  never left longer than that, even when nobody wakes it.
- `wake` is the URL the service calls once it queues a request (its
  MODAL_WAKE_URL), so a reader does not wait for the schedule. It starts
  `sweep` and does nothing else. Proxy auth guards it: a caller needs one of
  the workspace's proxy auth tokens, which the service holds as MODAL_KEY and
  MODAL_SECRET.
- `score` bakes queued constraint tracks with the scorer the twenty were baked
  with, `facebook/esm2_t33_650M_UR50D` at its pinned revision, on one L4. One
  container at a time, so a busy day costs the hours one GPU is busy.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import modal
from modal.file_pattern_matcher import NON_PYTHON_FILES

# Deploying, this file is run from the repository, which is where `app` and
# `pipeline` are imported from. In a container both are mounted beside it.
_HERE = Path(__file__).resolve()
if len(_HERE.parents) > 2 and (_HERE.parents[2] / "pipeline").is_dir():
    sys.path.insert(0, str(_HERE.parents[2]))

from pipeline.constraint.score_protein import MODEL, REVISION  # noqa: E402

APP_NAME = "helix-peak-resolver"
SECRET_NAME = "helix-peak-resolver"
WAKE_LABEL = "helix-peak-resolver-wake"
SCORER_GPU = "L4"

LOCK = Path(__file__).resolve().parents[1] / "constraint" / "requirements-lock.txt"

# Directories under `pipeline/` that are never code: bake environments, baked
# and fetched tracks, test fixtures. A local `.esm-venv` alone is all of torch.
_NOT_CODE = {".venv", "venv", ".esm-venv", "data", "fixtures", "__pycache__"}


def not_code(path: Path) -> bool:
    """What `add_local_python_source` leaves out: anything not a .py file, and
    anything under a directory that holds no code of ours.

    Modal hands this each file's path relative to the package, so a directory
    above the checkout that happens to be called `data` changes nothing.
    """
    return NON_PYTHON_FILES(path) or any(part in _NOT_CODE for part in Path(path).parts)


_ENV = {
    # Where the bakers read and write. Nothing survives the container, and
    # nothing needs to: the outputs go to storage and the rows.
    "HELIXPEEK_DATA": "/tmp/helixpeek/data",
    "HELIXPEEK_GB_CACHE": "/tmp/helixpeek/genbank",
    "PYTHONUNBUFFERED": "1",
}
_WORKER_PACKAGES = ("biopython==1.85", "psycopg[binary]==3.2.13")

resolver_image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(*_WORKER_PACKAGES, "fastapi[standard]==0.128.8")
    .env(_ENV)
    .add_local_python_source("app", "pipeline", ignore=not_code)
)

scorer_image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install_from_requirements(str(LOCK))
    .pip_install(*_WORKER_PACKAGES)
    # The weights go into the image, at the revision the scorer pins, so a
    # bake never waits on Hugging Face and never reads anything else.
    .run_commands(
        "python -c \"from huggingface_hub import snapshot_download; "
        f"snapshot_download('{MODEL}', revision='{REVISION}', "
        "allow_patterns=['*.json', '*.txt', '*.safetensors'])\""
    )
    .env({**_ENV, "HF_HUB_OFFLINE": "1"})
    .add_local_python_source("app", "pipeline", ignore=not_code)
)

app = modal.App(APP_NAME)
secret = modal.Secret.from_name(
    SECRET_NAME,
    required_keys=["DATABASE_URL", "SUPABASE_URL", "SUPABASE_SERVICE_KEY", "NCBI_EMAIL"],
)


def _connect():
    from pipeline.resolver import store
    return (store.connect(os.environ["DATABASE_URL"]),
            store.TrackStorage(os.environ["SUPABASE_URL"], os.environ["SUPABASE_SERVICE_KEY"]))


@app.function(image=resolver_image, secrets=[secret], schedule=modal.Period(minutes=5),
              timeout=15 * 60, max_containers=1)
def sweep() -> dict:
    """Resolve what is queued, and start the GPU if anything waits for it."""
    from pipeline.resolver import worker

    conn, storage = _connect()
    with conn:
        summary = worker.sweep(conn, storage)
    if summary["constraint_queued"]:
        score.spawn()
    print(json.dumps(summary))
    return summary


@app.function(image=scorer_image, secrets=[secret], gpu=SCORER_GPU,
              timeout=60 * 60, max_containers=1)
def score() -> list:
    """Bake every queued constraint track, one protein at a time."""
    from pipeline.resolver import worker

    conn, storage = _connect()
    with conn:
        done = worker.score_all(conn, storage)
    print(json.dumps(done))
    return done


@app.function(image=resolver_image)
@modal.fastapi_endpoint(method="POST", label=WAKE_LABEL, requires_proxy_auth=True)
def wake() -> dict:
    """Start a sweep now rather than at the next five minutes."""
    sweep.spawn()
    return {"woken": True}
