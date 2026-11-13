"""The resolver on Modal: everything the dev Mac's worker does, so a protein is
built with the Mac off.

From the repository root, with the Modal CLI logged in to the workspace
(`pip install -r pipeline/resolver/requirements.txt`, then `modal token new`):

    modal deploy pipeline/resolver/modal_app.py

It reads one secret from that workspace, `helix-peak-resolver`: the uploader's
three credentials (DATABASE_URL, SUPABASE_URL, SUPABASE_SERVICE_KEY), NCBI's
contact address (NCBI_EMAIL) and the AlphaGenome Atlas's key
(ALPHAGENOME_API_KEY), which of everything the functions start only the AVI
bake is handed. `pipeline/resolver/README.md` has the rest of the runbook.

Five functions, and what starts each:

- `sweep` resolves queued requests (`worker.sweep`) and starts each of the
  three bakers below where a bake of its kind waits that it can take
  (`store.claimable_bakes`). A schedule runs it every five minutes, so a
  request is never left longer than that, even when nobody wakes it.
- `wake` is the URL the service calls once it queues a request (its
  MODAL_WAKE_URL), so a reader does not wait for the schedule. It starts
  `sweep` and does nothing else. Proxy auth guards it: a caller needs one of
  the workspace's proxy auth tokens, which the service holds as MODAL_KEY and
  MODAL_SECRET.
- `fold` makes a protein's fold as the Mac's worker makes it, with the Mac's
  own driver (`local_worker.Modeller`, `model_local.py`): an experimental
  entry chosen by rule, baked with PyMOL 3.1.0 and held to the app's own
  check, else AlphaFold DB's model. Its image holds the structure bake's
  environment, PyMOL, Flutter and the part of the app's checkout that the
  scene importer and the app's check read.
- `evidence` makes a protein's AVI track (`local_worker.Evidencer`,
  `impact_local.py`, in the AVI bake's environment, with the AlphaGenome
  skill's GENCODE lookup under `uv`) and then a ClinVar track
  (`local_worker.ClinVarBaker`), and starts `score` once a scoring no longer
  waits for them (`store.WAITS_FOR`).
- `score` bakes one queued constraint track with the scorer the twenty were
  baked with, `facebook/esm2_t33_650M_UR50D` at its pinned revision, on one
  GPU: a T4, else an L4, else an A10 (`SCORER_GPU`). One container at a
  time, so a busy day costs the hours one GPU is busy.

Each works one protein a call and starts itself again while another waits,
as the Mac's worker does (`local_worker.py`). A request or a bake that goes
back on the queue ends the call's work, and the next sweep, minutes later,
tries it again: taken again at once, a network away for two seconds would
spend all three of its tries in one. A bake that only waits starts nothing.

`fold` and `evidence` claim nothing while their image could not make what
they are for (`Modeller.unready`, `Evidencer.unready`), exactly as the Mac's
worker: otherwise a broken image would refuse every protein's fold or AVI
after three tries. Whether the sweep asks for folds and evidence at all is
`RESOLVER_STRUCTURES` and `RESOLVER_EVIDENCE`, read where `modal deploy` runs
and read as the Mac's `.env` says them (on, unless 0, false, no or off):

    RESOLVER_STRUCTURES=0 RESOLVER_EVIDENCE=0 modal deploy pipeline/resolver/modal_app.py

The fold's image holds a copy of the app's checkout as it was at the deploy:
deploy again after a change to `lib/shared/folding/`,
`test/fold_bundle_check.dart` or the flutter_scene the app pins.
"""

from __future__ import annotations

import json
import logging
import os
import socket
import sys
import threading
from pathlib import Path
from typing import Mapping

import modal
from modal.file_pattern_matcher import NON_PYTHON_FILES

# Deploying, this file is run from the repository, which is where `app` and
# `pipeline` are imported from. In a container both are mounted beside it.
_HERE = Path(__file__).resolve()
if len(_HERE.parents) > 2 and (_HERE.parents[2] / "pipeline").is_dir():
    sys.path.insert(0, str(_HERE.parents[2]))

from pipeline import paths  # noqa: E402
from pipeline.constraint.score_protein import MODEL, REVISION  # noqa: E402

APP_NAME = "helix-peak-resolver"
SECRET_NAME = "helix-peak-resolver"
WAKE_LABEL = "helix-peak-resolver-wake"
# The scorer's GPUs, in the order Modal tries them: the T4, the cheapest
# ($0.59 an hour), then the L4 ($0.80), then the A10 ($1.10). On 2026-11-13
# Modal gave no L4 for half an hour while a T4 came in seconds, and the first
# protein built whole on Modal (B2M) was scored on a T4. ESM-2 650M in float32
# needs about 3 GB, which each holds; the model, its revision and the alignment
# gate are the same on each. On B2M (nothing written) the A10 took 10.6 s and
# the T4 15.0 s with the model's load, each within 0.001 of the L4 and of the
# Mac's mps at every number, and the track's `generation.device` says `cuda`
# on all.
SCORER_GPU = ["T4", "L4", "A10"]
# The secret's keys: the Mac's worker's four (`local_worker.KEYS`) and the
# Atlas's, which only the AVI bake is handed (`local_worker.Evidencer`).
KEYS = ("DATABASE_URL", "SUPABASE_URL", "SUPABASE_SERVICE_KEY", "NCBI_EMAIL",
        "ALPHAGENOME_API_KEY")

# How long an idle container is kept: Modal's least. Requests come minutes or
# days apart, so a warm container would be paid for and seldom used.
IDLE = 2
# One protein's scoring, at most. The longest protein a record can be built
# for (24 kb of exons, about 8,000 residues) is estimated at 12-23 minutes on
# an L4; on a T4, tried first, a long protein's time is not yet measured. Each
# call stays well inside the reaper's patience (`store.STALE_BAKE`,
# `STALE_RESOLVE`): a claim held past it would be taken again.
SCORE_TIMEOUT = 30 * 60
SWEEP_TIMEOUT = 15 * 60
# One fold: the modeller's check (`local_worker._CHECK_TIMEOUT`, two minutes)
# and one model, which the modeller stops at fifteen
# (`local_worker.MODEL_TIMEOUT`). Under the reaper's patience with a
# structure bake, twenty (`store.STALE_STRUCTURE`).
FOLD_TIMEOUT = 18 * 60
# One protein's evidence: the AVI bake's check, one AVI bake (stopped at fifty
# minutes, `local_worker.EVIDENCE_TIMEOUT`) and one ClinVar bake, a minute or
# two. Under the reaper's patience with either, sixty (`store.STALE_EVIDENCE`).
EVIDENCE_TIMEOUT = 57 * 60
# The record builder's Entrez call sets no timeout of its own; the Mac's worker
# sets the same (`local_worker.NETWORK_TIMEOUT`).
NETWORK_TIMEOUT = 120
# A sweep's most requests.
SWEEP_LIMIT = 10

# A fold runs PyMOL and compiles the app's check; the evidence reads GENCODE's
# 318 MB table, at up to about 4.4 GB, for each gene it looks up.
FOLD_CPU, FOLD_MEMORY = 2.0, 4096
EVIDENCE_CPU, EVIDENCE_MEMORY = 1.0, 6144

PIPELINE = _HERE.parents[1]
LOCK = PIPELINE / "constraint" / "requirements-lock.txt"
# The two bakes' environments, as the dev Mac has them: the structure bake's
# pins, and the AVI bake's environment pinned whole.
STRUCTURE_REQUIREMENTS = PIPELINE / "structure" / "requirements.txt"
IMPACT_LOCK = PIPELINE / "impact" / "requirements-lock.txt"

# Read from outside this repository, as the Mac's worker reads them: the app's
# checkout (`HELIXPEEK_CLIENT`) and the AlphaGenome skill (`ALPHAGENOME_SKILLS`,
# `bake_impact.DEFAULT_SKILL`). Only an image build reads them; importing this
# module never does.
CLIENT = paths.CLIENT
SKILL = paths.SKILLS / "alphagenome-variant-impact-score"
# What a fold needs of the app's checkout, and nothing more. The importer runs
# from its `pubspec.lock` (`alphafold.compile_scene`), and its check of a fold
# (`experimental.app_check`) compiles `test/fold_bundle_check.dart` against
# `lib/`, the assets `pubspec.yaml` declares and the app's build hook. Never
# its `.env` (an empty one stands in for the declared asset), `build/` or
# `.dart_tool/`.
APP_FILES = ("pubspec.yaml", "pubspec.lock", "test/fold_bundle_check.dart")
APP_DIRS = ("lib", "assets", "hook", "test/shared/folding")

# Where each is in the containers.
APP_DIR = "/opt/helix-peek"
FLUTTER_ROOT = "/opt/flutter"
FLUTTER = f"{FLUTTER_ROOT}/bin/flutter"
# The SDK's own binary, as on the Mac (`local_worker.default_dart`).
DART = f"{FLUTTER_ROOT}/bin/cache/dart-sdk/bin/dart"
STRUCTURE_PYTHON = "/opt/structure-venv/bin/python"
PYMOL_ENV = "/opt/pymol"
IMPACT_PYTHON = "/opt/impact-venv/bin/python"
UV = "/usr/local/bin/uv"
SKILLS_DIR = "/opt/skills"
GENCODE = "/opt/alphagenome/gencode.v46.annotation.gtf.gz.feather"
# The worker's own directories: its `data` and `genbank` are the bakers'
# (`_ENV`), and nothing in them outlives a call.
STATE = "/tmp/helixpeek"

# The tools a fold is made with, at the versions the dev Mac makes it with.
FLUTTER_VERSION = "3.47.3"
FLUTTER_ARCHIVE = ("https://storage.googleapis.com/flutter_infra_release/releases/stable/linux/"
                   f"flutter_linux_{FLUTTER_VERSION}-stable.tar.xz")
FLUTTER_SHA256 = "988665565cad9091db1baa54bf6d3868bb40e29719592f3c3a164deefd4208e1"
FLUTTER_REVISION = "e8113bf45620cbeb8aff64947ee4c93e16adb4cf"
ENGINE_REVISION = "06a2e2a110089dff50fe635cffd2a61e1b24fbcd"
DART_VERSION = "3.13.3"
# PyMOL in a conda environment of its own, so its numpy never meets the bake's.
MICROMAMBA = "2.1.1-0"
MICROMAMBA_SHA256 = "0d2dd49fb0171f27a63175e2c21941aafabdd6e46f521c498620b4b14ff604d8"
PYMOL_VERSION = "3.1.0"
PYMOL_BUILD = "py312h01792ee_8"
# The AVI bake's GENCODE lookup: the skill's own environment under `uv`, and
# the table it reads, fetched once here rather than for each protein.
UV_VERSION = "0.12.16"
SKILL_PYTHON = "3.14"
GENCODE_URL = ("https://storage.googleapis.com/alphagenome/reference/gencode/hg38/"
               "gencode.v46.annotation.gtf.gz.feather")
GENCODE_MD5 = "25067f565f7ed6e249b4247e98ecc68f"

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


def switched_on(name: str, environ: Mapping[str, str] = os.environ) -> bool:
    """A switch as the Mac's worker reads one from `.env`
    (`local_worker.load_settings`): on, unless it says 0, false, no or off."""
    return (environ.get(name) or "1").strip().lower() not in ("0", "false", "no", "off")


# Whether a protein this app resolves has its fold, and its AVI and ClinVar,
# queued and made here. Read where `modal deploy` runs, and carried to the
# sweep's containers by its image.
STRUCTURES = switched_on("RESOLVER_STRUCTURES")
EVIDENCE = switched_on("RESOLVER_EVIDENCE")

_ENV = {
    # Where the bakers read and write. Nothing survives the container, and
    # nothing needs to: the outputs go to storage and the rows.
    "HELIXPEEK_DATA": f"{STATE}/data",
    "HELIXPEEK_GB_CACHE": f"{STATE}/genbank",
    "PYTHONUNBUFFERED": "1",
}
_SWITCHES = {"RESOLVER_STRUCTURES": "1" if STRUCTURES else "0",
             "RESOLVER_EVIDENCE": "1" if EVIDENCE else "0"}
_WORKER_PACKAGES = ("biopython==1.85", "psycopg[binary]==3.2.13")
# The worker's numpy, as `.worker-venv` has it through Biopython: a fold's job
# is written with `structure/experimental.py`, which imports it.
_WORKER_NUMPY = "numpy==2.5.3"
_NOT_COPIED = ["**/.DS_Store"]

resolver_image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install(*_WORKER_PACKAGES, "fastapi[standard]==0.128.8")
    .env({**_ENV, **_SWITCHES})
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


def _app_checkout(image: modal.Image) -> modal.Image:
    """The parts of the app's checkout a fold reads (`APP_FILES`, `APP_DIRS`),
    copied into the image so that its build can compile them."""
    for name in APP_FILES:
        image = image.add_local_file(CLIENT / name, f"{APP_DIR}/{name}", copy=True)
    for name in APP_DIRS:
        image = image.add_local_dir(CLIENT / name, f"{APP_DIR}/{name}", copy=True,
                                    ignore=_NOT_COPIED)
    return image


_PYMOL_SAYS = ("printf 'python\\nfrom pymol import cmd\\n"
               "print(\"PYMOL_VERSION \" + cmd.get_version()[0])\\npython end\\n' "
               "> /tmp/pymol-version.pml")

_fold_tools = (
    modal.Image.debian_slim(python_version="3.12")
    # Flutter's own needs on Linux, and curl for what is fetched below.
    .apt_install("git", "curl", "ca-certificates", "xz-utils", "unzip", "libglu1-mesa")
    .pip_install(*_WORKER_PACKAGES, _WORKER_NUMPY)
    # The structure bake's environment (`pipeline/structure/venv` on the Mac).
    .add_local_file(STRUCTURE_REQUIREMENTS, "/opt/structure-requirements.txt", copy=True)
    .run_commands(
        "python -m venv /opt/structure-venv",
        "/opt/structure-venv/bin/pip install --no-cache-dir -r /opt/structure-requirements.txt",
    )
    # PyMOL, the version every model is drawn with (`experimental.PYMOL`), on
    # the path as `pymol`, where `structure/bake.py` looks for it.
    .run_commands(
        "curl -fsSL -o /usr/local/bin/micromamba https://github.com/mamba-org/"
        f"micromamba-releases/releases/download/{MICROMAMBA}/micromamba-linux-64",
        f"echo '{MICROMAMBA_SHA256}  /usr/local/bin/micromamba' | sha256sum -c -",
        "chmod +x /usr/local/bin/micromamba",
        f"MAMBA_ROOT_PREFIX=/opt/mamba micromamba create --yes --prefix {PYMOL_ENV} "
        f"--override-channels -c conda-forge 'pymol-open-source={PYMOL_VERSION}={PYMOL_BUILD}'",
        "MAMBA_ROOT_PREFIX=/opt/mamba micromamba clean --all --yes",
        f"ln -s {PYMOL_ENV}/bin/pymol /usr/local/bin/pymol",
        _PYMOL_SAYS,
        f"pymol -cq /tmp/pymol-version.pml | grep -qx 'PYMOL_VERSION {PYMOL_VERSION}'",
    )
    # Flutter, the release the dev Mac runs, held to its revisions.
    .run_commands(
        f"curl -fsSL -o /tmp/flutter.tar.xz {FLUTTER_ARCHIVE}",
        f"echo '{FLUTTER_SHA256}  /tmp/flutter.tar.xz' | sha256sum -c -",
        "tar --no-same-owner -xJf /tmp/flutter.tar.xz -C /opt && rm /tmp/flutter.tar.xz",
        f"git config --system --add safe.directory {FLUTTER_ROOT}",
        f"{FLUTTER} --disable-analytics || true",
        f"{DART} --disable-analytics || true",
        f"{FLUTTER} --version --machine > /tmp/flutter-version.json",
        f"grep -Eq '\"frameworkRevision\": *\"{FLUTTER_REVISION}\"' /tmp/flutter-version.json "
        f"&& grep -Eq '\"engineRevision\": *\"{ENGINE_REVISION}\"' /tmp/flutter-version.json "
        f"&& grep -Eq '\"dartSdkVersion\": *\"{DART_VERSION}' /tmp/flutter-version.json "
        "|| { cat /tmp/flutter-version.json; exit 1; }",
    )
    .env({"PUB_CACHE": "/opt/pub-cache"})
)

fold_image = (
    _app_checkout(_fold_tools)
    # The app's packages, as its lockfile pins them, and the first compile of
    # its check, the scene importer and the build hooks, made here so that a
    # failure shows at the build and not on a protein.
    .run_commands(
        f": > {APP_DIR}/.env",
        f"cd {APP_DIR} && {FLUTTER} pub get --enforce-lockfile",
        f"cd {APP_DIR} && {FLUTTER} test --no-pub test/fold_bundle_check.dart",
        f"{STRUCTURE_PYTHON} -c \"import trimesh; scene = trimesh.Scene(); "
        "scene.add_geometry(trimesh.creation.box(), geom_name='check', node_name='check'); "
        "scene.export('/tmp/check.glb')\"",
        f"cd {APP_DIR} && {DART} run flutter_scene:import -i /tmp/check.glb "
        "-o /tmp/check.fsceneb && test -s /tmp/check.fsceneb",
        "rm /tmp/check.glb /tmp/check.fsceneb",
    )
    .env(_ENV)
    .add_local_python_source("app", "pipeline", ignore=not_code)
)

evidence_image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("curl", "ca-certificates")
    .pip_install(*_WORKER_PACKAGES, _WORKER_NUMPY, f"uv=={UV_VERSION}")
    # The AVI bake's environment (`pipeline/impact/venv` on the Mac), whole.
    .add_local_file(IMPACT_LOCK, "/opt/impact-requirements-lock.txt", copy=True)
    .run_commands(
        "python -m venv /opt/impact-venv",
        "/opt/impact-venv/bin/pip install --no-cache-dir -r /opt/impact-requirements-lock.txt",
    )
    .env({"ALPHAGENOME_SKILLS": SKILLS_DIR, "ALPHAGENOME_GTF_PATH": GENCODE,
          "UV_PYTHON_INSTALL_DIR": "/opt/uv/python", "UV_CACHE_DIR": "/opt/uv/cache"})
    # The table the skill's GENCODE lookup reads (it reads ALPHAGENOME_GTF_PATH
    # before its cache), the object it would otherwise fetch for each protein.
    .run_commands(
        f"mkdir -p {Path(GENCODE).parent} && curl -fsSL -o {GENCODE} {GENCODE_URL}",
        f"echo '{GENCODE_MD5}  {GENCODE}' | md5sum -c -",
    )
    .add_local_dir(SKILL, f"{SKILLS_DIR}/{SKILL.name}", copy=True,
                   ignore=[".venv", "**/__pycache__", *_NOT_COPIED])
    .run_commands(
        f"cd {SKILLS_DIR}/{SKILL.name} && uv sync --frozen --python {SKILL_PYTHON}",
        f"cd {SKILLS_DIR}/{SKILL.name} && uv run scripts/alphagenome_atlas_avi.py gtf --help "
        "> /dev/null",
    )
    # From here on the skill's environment is as locked, and nothing is fetched.
    .env({"UV_FROZEN": "1", "UV_OFFLINE": "1"})
    .env(_ENV)
    .add_local_python_source("app", "pipeline", ignore=not_code)
)

app = modal.App(APP_NAME)
secret = modal.Secret.from_name(SECRET_NAME, required_keys=list(KEYS))


def _connect():
    """The database, opened as the Mac's worker opens it (`local_worker._conninfo`:
    keepalives, for a claim held minutes with nothing said), and storage."""
    from pipeline.resolver import local_worker, store
    return (store.connect(local_worker._conninfo(os.environ["DATABASE_URL"])),
            store.TrackStorage(os.environ["SUPABASE_URL"], os.environ["SUPABASE_SERVICE_KEY"]))


def _settings():
    """The Mac's worker's settings, for a container (`local_worker.Settings`,
    never `load_settings`, which reads a Mac's `.env` and home): the secret's
    credentials, the worker's directories under `STATE`, and the tools where
    the images put them."""
    from pipeline.resolver import local_worker

    state = Path(STATE)
    return local_worker.Settings(
        database_url=os.environ["DATABASE_URL"],
        supabase_url=os.environ["SUPABASE_URL"],
        service_key=os.environ["SUPABASE_SERVICE_KEY"],
        ncbi_email=os.environ["NCBI_EMAIL"],
        state=state,
        log_file=state / "worker.log",
        structure_python=Path(STRUCTURE_PYTHON),
        dart=Path(DART),
        flutter=Path(FLUTTER),
        app=Path(APP_DIR),
        structures=True,
        evidence=True,
        impact_python=Path(IMPACT_PYTHON),
        uv=Path(UV),
        alphagenome_key=os.environ.get("ALPHAGENOME_API_KEY"),
    )


class _Said(logging.StreamHandler):
    """The handler that puts the Mac's worker's lines in a function's log."""


def _work_here(settings) -> None:
    """This container set up as `local_worker.main` sets up the Mac: the
    bakers' directories named and held to (`take_directories`: never
    `pipeline/data/`), and the drivers' lines in the function's log, with no
    credential in them (`local_worker._Formatter`)."""
    from pipeline.resolver import local_worker

    local_worker.take_directories(settings)
    log = local_worker.log
    if not any(isinstance(handler, _Said) for handler in log.handlers):
        handler = _Said(sys.stdout)
        handler.setFormatter(local_worker._Formatter(settings.secrets))
        log.addHandler(handler)
    log.setLevel(logging.INFO)
    log.propagate = False


def resolve_queued(conn, storage, sweep=None, *, structures: bool = False,
                   evidence: bool = False) -> dict:
    """Resolve queued requests one at a time, as `local_worker._resolve` does:
    up to `SWEEP_LIMIT`, ending where none is left or where one went back on
    the queue. `structures` and `evidence` have each protein resolved queue its
    fold, and its AVI and ClinVar, as well. Returns `worker.sweep`'s summary,
    with every request resolved."""
    if sweep is None:
        from pipeline.resolver import worker
        sweep = worker.sweep
    asked = {}
    if structures:
        asked["structures"] = True
    if evidence:
        asked["evidence"] = True
    resolved: list = []
    summary: dict = {}
    for _ in range(SWEEP_LIMIT):
        summary = sweep(conn, storage, limit=1, **asked)
        resolved.extend(summary["resolved"])
        if not summary["resolved"] or summary["resolved"][-1]["state"] == "queued":
            break
    return {**summary, "resolved": resolved}


def goes_on(done: list) -> bool:
    """Whether a function may start itself again after a call that did `done`:
    it took work, and put none of it back on the queue."""
    return bool(done) and all(outcome["state"] != "queued" for outcome in done)


# The kinds of bake the sweep counts, and starts a function for.
_KINDS = ("constraint", "structure", "impact", "clinvar")


def starts(claimable: Mapping[str, int], structures: bool, evidence: bool) -> list:
    """The functions a sweep starts: each where a bake of its kind waits that
    it can take (`store.claimable_bakes`), and `fold` and `evidence` only
    where this app makes folds and evidence."""
    started = []
    if claimable.get("constraint"):
        started.append("score")
    if structures and claimable.get("structure"):
        started.append("fold")
    if evidence and (claimable.get("impact") or claimable.get("clinvar")):
        started.append("evidence")
    return started


@app.function(image=resolver_image, secrets=[secret], schedule=modal.Period(minutes=5),
              timeout=SWEEP_TIMEOUT, max_containers=1, scaledown_window=IDLE)
def sweep() -> dict:
    """Resolve what is queued, and start each baker for which a bake waits that
    it can take."""
    from pipeline.resolver import store

    socket.setdefaulttimeout(NETWORK_TIMEOUT)
    conn, storage = _connect()
    with conn:
        summary = resolve_queued(conn, storage, structures=STRUCTURES, evidence=EVIDENCE)
        claimable = {kind: store.claimable_bakes(conn, kind) for kind in _KINDS}
    started = starts(claimable, STRUCTURES, EVIDENCE)
    bakers = {"score": score, "fold": fold, "evidence": evidence}
    for name in started:
        bakers[name].spawn()
    summary.update(constraint_claimable=claimable["constraint"], claimable=claimable,
                   structures=STRUCTURES, evidence=EVIDENCE, started=started)
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


# What the fold says while the database has no table to publish a fold in.
NO_FOLDS = "The database cannot record a fold yet: apply migrations/0012_fold_generation.sql."


@app.function(image=fold_image, secrets=[secret], cpu=FOLD_CPU, memory=FOLD_MEMORY,
              timeout=FOLD_TIMEOUT, max_containers=1, scaledown_window=IDLE)
def fold() -> dict:
    """Make the oldest fold it can take, as the Mac's worker makes one, and
    start again while another waits. Claims nothing while no fold could be
    made here."""
    from pipeline.resolver import local_worker, store, worker

    socket.setdefaulttimeout(NETWORK_TIMEOUT)
    settings = _settings()
    _work_here(settings)
    modeller = local_worker.Modeller(settings, threading.Event())
    done: list = []
    again = False
    try:
        # Asked before anything is claimed, as the Mac asks (`local_worker.cycle`).
        unready = modeller.unready()
        if unready is None:
            conn, storage = _connect()
            with conn:
                if not store.folds_ready(conn):
                    unready = NO_FOLDS
                else:
                    done = worker.structure_all(conn, storage, limit=1, model=modeller)
                    again = goes_on(done) and store.claimable_bakes(conn, "structure") > 0
    finally:
        local_worker.clear_workspace(settings)
    if unready:
        local_worker.log.error("no model can be made here now, and none was claimed. %s",
                               unready)
    if again:
        fold.spawn()
    outcome = {"unready": unready, "modelled": done}
    print(json.dumps(outcome))
    return outcome


@app.function(image=evidence_image, secrets=[secret], cpu=EVIDENCE_CPU, memory=EVIDENCE_MEMORY,
              timeout=EVIDENCE_TIMEOUT, max_containers=1, scaledown_window=IDLE)
def evidence() -> dict:
    """Make the oldest AVI track it can take, then the oldest ClinVar track, as
    the Mac's worker makes them; start again while either waits, and start the
    GPU once a scoring no longer waits for them. Claims no AVI bake while the
    Atlas could not be asked."""
    from pipeline.resolver import local_worker, store, worker

    socket.setdefaulttimeout(NETWORK_TIMEOUT)
    settings = _settings()
    _work_here(settings)
    stop = threading.Event()
    evidencer = local_worker.Evidencer(settings, stop)
    clinvar = local_worker.ClinVarBaker(settings, stop)
    placed: list = []
    fetched: list = []
    again = scoring = False
    try:
        # The Atlas is asked (one request) before an AVI bake is claimed.
        # ClinVar needs NCBI alone and is claimed only once its protein's AVI
        # bake is over (`store.WAITS_FOR`), so it goes on either way, as on the
        # Mac (`local_worker.cycle`).
        unready = evidencer.unready()
        conn, storage = _connect()
        with conn:
            if unready is None:
                placed = worker.impact_all(conn, storage, limit=1, avi=evidencer)
            fetched = worker.clinvar_all(conn, storage, limit=1, clinvar=clinvar)
            # What a call of its own could take next: no AVI bake while the
            # Atlas could not be asked.
            waiting = store.claimable_bakes(conn, "clinvar")
            if unready is None:
                waiting += store.claimable_bakes(conn, "impact")
            again = goes_on(placed + fetched) and waiting > 0
            scoring = bool(fetched) and store.claimable_bakes(conn, "constraint") > 0
    finally:
        local_worker.clear_workspace(settings)
    if unready:
        local_worker.log.error("no AVI track can be made here now, and no AVI bake was "
                               "claimed. %s", unready)
    if again:
        evidence.spawn()
    if scoring:
        score.spawn()
    outcome = {"unready": unready, "impact": placed, "clinvar": fetched, "scoring": scoring}
    print(json.dumps(outcome))
    return outcome


@app.function(image=resolver_image, scaledown_window=IDLE)
@modal.fastapi_endpoint(method="POST", label=WAKE_LABEL, requires_proxy_auth=True)
def wake() -> dict:
    """Start a sweep now rather than at the next five minutes."""
    sweep.spawn()
    return {"woken": True}
