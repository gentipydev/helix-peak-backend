"""The resolver on Modal: requests resolved on a CPU, ESM-2 scored on a GPU.

From the repository root, with the Modal CLI logged in to the workspace
(`pip install -r pipeline/resolver/requirements.txt`, then `modal token new`):

    modal deploy pipeline/resolver/modal_app.py

It reads one secret from that workspace, `helix-peak-resolver`: the uploader's
three credentials (DATABASE_URL, SUPABASE_URL, SUPABASE_SERVICE_KEY) and NCBI's
contact address (NCBI_EMAIL). `pipeline/resolver/README.md` has the rest of the
runbook.

Three functions, and what starts each:

- `sweep` resolves queued requests (`worker.sweep`) and, where a constraint
  bake waits that `score` can take (`store.claimable_bakes`), starts `score`.
  A schedule runs it every five minutes, so a request is never left longer
  than that, even when nobody wakes it.
- `wake` is the URL the service calls once it queues a request (its
  MODAL_WAKE_URL), so a reader does not wait for the schedule. It starts
  `sweep` and does nothing else. Proxy auth guards it: a caller needs one of
  the workspace's proxy auth tokens, which the service holds as MODAL_KEY and
  MODAL_SECRET.
- `score` bakes one queued constraint track with the scorer the twenty were
  baked with, `facebook/esm2_t33_650M_UR50D` at its pinned revision, on one
  L4, and starts itself again while another waits. One container at a time,
  so a busy day costs the hours one GPU is busy.

Both work one protein at a time, as the Mac's worker does (`local_worker.py`).
A request or a bake that goes back on the queue ends the call's work, and the
next sweep, minutes later, tries it again: taken again at once, a network away
for two seconds would spend all three of its tries in one. A bake that only
waits (for its protein's AVI and ClinVar, which the Mac's worker makes, or one
of the twenty's) starts no GPU, which would find nothing to take.
"""

from __future__ import annotations

import json
import os
import socket
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

# How long an idle container is kept: Modal's least. Requests come minutes or
# days apart, so a warm container would be paid for and seldom used.
IDLE = 2
# One protein's scoring, at most. The longest protein a record can be built
# for (24 kb of exons, about 8,000 residues) is estimated at 12-23 minutes on
# an L4. Each call stays well inside the reaper's patience (`store.STALE_BAKE`,
# `STALE_RESOLVE`): a claim held past it would be taken again.
SCORE_TIMEOUT = 30 * 60
SWEEP_TIMEOUT = 15 * 60
# The record builder's Entrez call sets no timeout of its own; the Mac's worker
# sets the same (`local_worker.NETWORK_TIMEOUT`).
NETWORK_TIMEOUT = 120
# A sweep's most requests.
SWEEP_LIMIT = 10

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


def resolve_queued(conn, storage, sweep=None) -> dict:
    """Resolve queued requests one at a time, as `local_worker._resolve` does:
    up to `SWEEP_LIMIT`, ending where none is left or where one went back on
    the queue. Returns `worker.sweep`'s summary, with every request resolved."""
    if sweep is None:
        from pipeline.resolver import worker
        sweep = worker.sweep
    resolved: list = []
    summary: dict = {}
    for _ in range(SWEEP_LIMIT):
        summary = sweep(conn, storage, limit=1)
        resolved.extend(summary["resolved"])
        if not summary["resolved"] or summary["resolved"][-1]["state"] == "queued":
            break
    return {**summary, "resolved": resolved}


def goes_on(done: list) -> bool:
    """Whether `score` may start again after a call that did `done`: it took
    a bake, and did not put it back on the queue."""
    return bool(done) and done[-1]["state"] != "queued"


@app.function(image=resolver_image, secrets=[secret], schedule=modal.Period(minutes=5),
              timeout=SWEEP_TIMEOUT, max_containers=1, scaledown_window=IDLE)
def sweep() -> dict:
    """Resolve what is queued, and start the GPU if a bake waits that it can take."""
    from pipeline.resolver import store

    socket.setdefaulttimeout(NETWORK_TIMEOUT)
    conn, storage = _connect()
    with conn:
        summary = resolve_queued(conn, storage)
        summary["constraint_claimable"] = store.claimable_bakes(conn, "constraint")
    if summary["constraint_claimable"]:
        score.spawn()
    print(json.dumps(summary))
    return summary


@app.function(image=scorer_image, secrets=[secret], gpu=SCORER_GPU,
              timeout=SCORE_TIMEOUT, max_containers=1, scaledown_window=IDLE)
def score() -> list:
    """Bake the oldest constraint track it can take, and start again while
    another waits."""
    from pipeline.resolver import store, worker

    socket.setdefaulttimeout(NETWORK_TIMEOUT)
    conn, storage = _connect()
    with conn:
        done = worker.score_all(conn, storage, limit=1)
        again = goes_on(done) and store.claimable_bakes(conn, "constraint") > 0
    if again:
        score.spawn()
    print(json.dumps(done))
    return done


@app.function(image=resolver_image, scaledown_window=IDLE)
@modal.fastapi_endpoint(method="POST", label=WAKE_LABEL, requires_proxy_auth=True)
def wake() -> dict:
    """Start a sweep now rather than at the next five minutes."""
    sweep.spawn()
    return {"woken": True}
