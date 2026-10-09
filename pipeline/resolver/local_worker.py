"""The resolver's worker on a Mac, until Modal runs it.

    cd helix-peek-backend
    .worker-venv/bin/python -u -m pipeline.resolver.local_worker           # the loop
    .worker-venv/bin/python -u -m pipeline.resolver.local_worker --once    # one cycle
    .worker-venv/bin/python -u -m pipeline.resolver.local_worker --queue   # what waits; reads only

`scripts/resolver_worker.sh` runs the first as a launchd agent, and the README
beside this file is the runbook.

It is `worker.sweep` and `worker.score_all` as `modal_app.py` runs them, on the
same database and storage, with the same four credentials, read from the
repository's `.env` where Modal reads its secret. A protein built here cannot
be told from one built there, except by the device its track names, and by
its fold: between the two it also makes each protein's model
(`worker.structure_all`), which Modal's image cannot. What differs is how the
work is started and where ESM-2 runs:

- Nothing wakes it. It asks the queue every `RESOLVER_POLL_SECONDS` (60).
- `worker.score_all` scores with whatever it is handed (`worker.Score`). On
  Modal that is `scoring.score_with_esm` in the worker's own process, because
  one image holds torch and psycopg both. Here no environment does: the
  scorer's, `pipeline/.esm-venv`, baked the twenty and is given nothing new.
  So a `Scorer` is handed over instead, which runs `score_local.py` there as a
  subprocess, one protein at a time, and comes back with what
  `score_with_esm` would have: the track's bytes, or the scorer's own
  `ValueError`.
- A model wants numpy, trimesh, PyMOL and the app's scene importer, and the
  worker's environment has none of them. A `Modeller` runs `model_local.py`
  in the structure bake's own environment (`pipeline/structure/venv`), the
  same way. A model takes seconds, so it is made before the scoring: a
  protein's fold page is ready long before its ESM-2 track is.
- The variant evidence comes next, and before the scoring too
  (`worker.impact_all`, then `worker.clinvar_all`): minutes, where ESM-2 can
  take an hour, so an Atlas quota spent shows early. AVI wants the AlphaGenome
  client and its skill's GENCODE lookup, which runs under `uv`, and an
  `Evidencer` runs `impact_local.py` in the AVI bake's own environment
  (`pipeline/impact/venv`), with `ALPHAGENOME_API_KEY` from `.env` and nothing
  else secret. ClinVar wants Biopython only, and is fetched in this process
  (`ClinVarBaker`). `RESOLVER_EVIDENCE=0` in `.env` turns both off.

A laptop is not a container, and four things keep that out of the rows:

- One protein at a time. `sweep` and `score_all` are called with `limit=1`, so
  a stop lands between proteins, and a protein that goes back on the queue
  ends the cycle's work on that queue. The next cycle tries it again, a poll
  later, rather than this one at once with whatever broke it still broken:
  three tries in one second would refuse a protein for good over a network
  that was away for two.
- No bake is claimed while the scorer cannot start (`Scorer.unready`), no
  model while the modeller cannot (`Modeller.unready`), and no AVI bake while
  the Atlas cannot be asked (`Evidencer.unready`).
- The bakers' directories are the worker's own, under Application Support and
  never `pipeline/data/`, and are emptied after each cycle that used them, as
  a container's `/tmp` is. A GenBank reply garbled once is not read twice.
- Nothing waits for ever on a network that went away with the lid: sockets
  time out, and the database connection is probed while a protein is scored.

It keeps the Mac awake (`caffeinate -i`) while it has work, and only then. It
is the only worker on the machine (a lock). SIGTERM or SIGINT stops it between
proteins, ending a scorer that is running so that its bake goes back on the
queue at once. An error that is no protein's (the database out of reach) is
logged and waited out, longer each time, up to ten minutes.

Nothing under `pipeline` is imported when this module is. `paths.DATA` and the
record builder's `CACHE` are fixed the first time they are imported, from the
environment as it is then, and `main` names the worker's directories first.
"""

from __future__ import annotations

import argparse
import json
import logging
import logging.handlers
import math
import os
import pickle
import re
import shutil
import signal
import socket
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Optional
from urllib.parse import unquote

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

log = logging.getLogger("resolver-worker")

# The uploader's three credentials and NCBI's contact address: the four keys of
# Modal's `helix-peak-resolver` secret (`modal_app.py`).
KEYS = ("DATABASE_URL", "SUPABASE_URL", "SUPABASE_SERVICE_KEY", "NCBI_EMAIL")

# What the scorer is started without. It reads a record and writes a track,
# both as files; it has no use for the database, the storage key or the
# Atlas's, which only the AVI bake is given (`Evidencer`).
_CREDENTIALS = ("DATABASE_URL", "SUPABASE_SERVICE_KEY", "ALPHAGENOME_API_KEY")

# `store.STALE_BAKE` is 90 minutes: a sweep takes a bake that has run that long
# for a dead worker's and queues it again. A scoring let run past it could be
# claimed twice, so the timeout stays under it whatever is asked for.
SCORE_TIMEOUT_CEILING = 85 * 60

# A cycle's most, the defaults `worker.sweep` and `worker.score_all` have on Modal.
SWEEP_LIMIT = 10
SCORE_LIMIT = 20

# How long an error that is no protein's is waited out, at most.
MAX_BACKOFF = 600

# Kept from idle sleep while there is work. `-w` ends it with this process, so
# it cannot outlive a worker that is killed.
CAFFEINATE = ("caffeinate", "-i")

# The record builder's Entrez call sets no timeout of its own. Without one, a
# download the Mac slept through waits on a connection nobody is on.
NETWORK_TIMEOUT = 120

# What the worker's connection is opened with, where `DATABASE_URL` does not
# say otherwise. Nothing is said on it from a bake's claim to its track, which
# can be an hour: the probes keep a router from forgetting it, and let a
# connection that died with the network be found dead rather than waited on.
# They are slow on purpose, so that a minute without network costs nothing.
_CONNECTION = {"connect_timeout": "20", "keepalives": "1", "keepalives_idle": "60",
               "keepalives_interval": "30", "keepalives_count": "10"}

_DRIVER = "pipeline.resolver.score_local"
_MODEL_DRIVER = "pipeline.resolver.model_local"
_IMPACT_DRIVER = "pipeline.resolver.impact_local"
# A model is baked in seconds. One still baking after this has hung, and is
# ended well inside `store.STALE_STRUCTURE`, as a scoring is inside the bake's.
MODEL_TIMEOUT = 15 * 60
# An AVI bake takes seconds to minutes, waiting on the Atlas's quota at worst.
# One still going after this is ended inside `store.STALE_EVIDENCE`, and its
# checkpoints carry what it pulled into the next try.
EVIDENCE_TIMEOUT = 50 * 60
_CHECK_TIMEOUT = 120
_CHILD_POLL = 0.5
_PROGRESS_EVERY = 60
# How often a reader watching the protein is told how far ESM-2 has got.
_REPORT_EVERY = 5
_LAST_WORDS = 220
_EX_CONFIG = 78

STOPPED = "The worker was stopped while ESM-2 was scoring."
STOPPED_MODELLING = "The worker was stopped while a model was being made."
STOPPED_PLACING = "The worker was stopped while AlphaGenome's scores were being fetched."

_DEVICE = re.compile(r"^Loading \S+ on (\w+) ", re.MULTILINE)
_BEFORE = re.compile(r"over the residue before:\s+([\d.]+%)")
_AFTER = re.compile(r"over the residue after:\s+([\d.]+%)")
_PROGRESS = re.compile(r"^Scored \d+/\d+ .*$", re.MULTILINE)

# A password, as a database URL carries it and as libpq's keywords do.
_PASSWORDS = (re.compile(r"//[^/@\s:]*:([^@\s]+)@"), re.compile(r"password\s*=\s*'?([^'\s]+)"))


class Misconfigured(Exception):
    """The worker cannot run as it is set up. The message says what to change."""


@dataclass(frozen=True)
class Settings:
    database_url: str
    supabase_url: str
    service_key: str
    ncbi_email: str
    # Everything the worker writes that is neither a row nor a stored object:
    # the bakers' two directories, the scorer's hand-over files, and the lock.
    state: Path
    log_file: Path
    poll: float = 60.0
    score_timeout: float = 3600.0
    esm_python: Path = BACKEND / "pipeline" / ".esm-venv" / "bin" / "python"
    # What a model is made with: the structure bake's interpreter, the Dart
    # that runs the scene importer, and the app's checkout it is run from.
    structure_python: Path = BACKEND / "pipeline" / "structure" / "venv" / "bin" / "python"
    dart: Path = Path("dart")
    app: Path = BACKEND.parent / "helix-peek"
    # Whether this worker makes models at all. Off, no structure bake is
    # queued and none is claimed: the worker is the one Modal runs.
    structures: bool = False
    # Whether it makes AVI and ClinVar tracks, the same way, and what with:
    # the AVI bake's interpreter, the `uv` its GENCODE lookup runs under
    # (launchd's PATH does not name the folder it is in), and the Atlas's key,
    # handed to that bake alone.
    evidence: bool = False
    impact_python: Path = BACKEND / "pipeline" / "impact" / "venv" / "bin" / "python"
    uv: Path = Path("uv")
    alphagenome_key: Optional[str] = None
    # What `load_settings` changed from what was asked, said once the log is open.
    notes: tuple = ()

    @property
    def data(self) -> Path:
        """`HELIXPEEK_DATA`: where the scorer finds the record and leaves its track."""
        return self.state / "data"

    @property
    def genbank(self) -> Path:
        """`HELIXPEEK_GB_CACHE`: the record builder's flat files."""
        return self.state / "genbank"

    @property
    def scoring(self) -> Path:
        return self.state / "scoring"

    @property
    def modelling(self) -> Path:
        """Where a model is baked: the model file, PyMOL's export, the scene."""
        return self.state / "modelling"

    @property
    def fetching(self) -> Path:
        """Where AVI and ClinVar are fetched into: the AVI bake's hand-over,
        and ClinVar's raw responses while its bake runs."""
        return self.state / "fetching"

    @property
    def impact_state(self) -> Path:
        """What an AVI bake keeps between tries, and so is never emptied:
        GENCODE's answers, and the checkpoints of a pull cut short."""
        return self.state / "impact"

    @property
    def lock(self) -> Path:
        return self.state / "worker.lock"

    @property
    def secrets(self) -> tuple:
        """What no line of the log may hold: the storage key, the Atlas's, and
        the database's URL and password. Longest first, so a whole URL goes
        before a part of it."""
        found = [self.service_key, self.database_url]
        if self.alphagenome_key:
            found.append(self.alphagenome_key)
        for pattern in _PASSWORDS:
            match = pattern.search(self.database_url)
            if match:
                found += [match.group(1), unquote(match.group(1))]
        return tuple(sorted({secret for secret in found if len(secret) >= 4},
                            key=len, reverse=True))


def _seconds(said: Mapping[str, str], name: str, default: float) -> float:
    raw = said.get(name)
    if raw is None:
        return default
    try:
        value = float(raw)
    except ValueError:
        value = math.nan
    if not (math.isfinite(value) and value >= 1):
        raise Misconfigured(f"{name} is {raw!r}; it is a number of seconds, 1 or more.")
    return value


def default_dart(home: Path) -> Path:
    """The Dart SDK inside a Flutter checkout in the home directory, where
    there is one, and otherwise whichever `dart` is on the path.

    The SDK's own binary and not `flutter/bin/dart`, which is a script that
    first takes Flutter's start-up lock: a worker should not wait on, or hold
    up, a `flutter run` in the next terminal. Under launchd the path names
    neither, so the first is what the agent finds.
    """
    sdk = home / "flutter" / "bin" / "cache" / "dart-sdk" / "bin" / "dart"
    return sdk if sdk.exists() else Path(shutil.which("dart") or "dart")


def default_uv(home: Path) -> Path:
    """`uv` where its installer puts it in the home directory, and otherwise
    whichever is on the path. Under launchd the path names neither."""
    installed = home / ".local" / "bin" / "uv"
    return installed if installed.exists() else Path(shutil.which("uv") or "uv")


def load_settings(environ: Optional[Mapping[str, str]] = None, env_file: Optional[Path] = None,
                  home: Optional[Path] = None) -> Settings:
    """The worker's settings: the environment's word first, then `.env`'s.

    Under launchd the environment names none of them, so what is read is the
    repository's `.env`, the file the service and the uploader read on this
    machine, and the only place an override can be put for the agent. The
    values go nowhere else: not into a plist, not into a log line.
    """
    from dotenv import dotenv_values

    environ = os.environ if environ is None else environ
    env_file = BACKEND / ".env" if env_file is None else env_file
    home = Path.home() if home is None else home

    said = {key: value for key, value in dotenv_values(env_file).items() if value}
    said.update({key: value for key, value in environ.items() if value})
    missing = [key for key in KEYS if key not in said]
    if missing:
        raise Misconfigured(f"{', '.join(missing)} {'is' if len(missing) == 1 else 'are'} "
                            f"not set in {env_file} or in the environment.")

    notes = []
    timeout = _seconds(said, "RESOLVER_SCORE_TIMEOUT", 3600.0)
    if timeout > SCORE_TIMEOUT_CEILING:
        notes.append(
            f"RESOLVER_SCORE_TIMEOUT is {timeout:g} s. A sweep takes a bake that has run 90 "
            f"minutes for a dead worker's, so scoring stops at {SCORE_TIMEOUT_CEILING} s.")
        timeout = float(SCORE_TIMEOUT_CEILING)
    return Settings(
        database_url=said["DATABASE_URL"],
        supabase_url=said["SUPABASE_URL"],
        service_key=said["SUPABASE_SERVICE_KEY"],
        ncbi_email=said["NCBI_EMAIL"],
        state=home / "Library" / "Application Support" / "HelixPeek" / "resolver-worker",
        log_file=home / "Library" / "Logs" / "HelixPeek" / "resolver-worker.log",
        poll=_seconds(said, "RESOLVER_POLL_SECONDS", 60.0),
        score_timeout=timeout,
        esm_python=Path(said.get("RESOLVER_ESM_PYTHON")
                        or BACKEND / "pipeline" / ".esm-venv" / "bin" / "python"),
        structure_python=Path(said.get("RESOLVER_STRUCTURE_PYTHON")
                              or BACKEND / "pipeline" / "structure" / "venv" / "bin" / "python"),
        dart=Path(said.get("RESOLVER_DART") or default_dart(home)),
        app=Path(said.get("RESOLVER_APP_DIR") or BACKEND.parent / "helix-peek"),
        structures=said.get("RESOLVER_STRUCTURES", "1").strip().lower()
        not in ("0", "false", "no", "off"),
        evidence=said.get("RESOLVER_EVIDENCE", "1").strip().lower()
        not in ("0", "false", "no", "off"),
        impact_python=Path(said.get("RESOLVER_IMPACT_PYTHON")
                           or BACKEND / "pipeline" / "impact" / "venv" / "bin" / "python"),
        uv=Path(said.get("RESOLVER_UV") or default_uv(home)),
        alphagenome_key=said.get("ALPHAGENOME_API_KEY"),
        notes=tuple(notes),
    )


# ------------------------------------------------------------ the scorer


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""


def _last_words(code: int, errors: str, who: str = "The scorer") -> str:
    """Why the scorer broke, inside `worker._said`'s 300 characters, keeping the
    end of what it said: a traceback names the fault on its last line. `who`
    names another child, the modeller."""
    how = f"was ended by signal {-code}" if code < 0 else f"exited with status {code}"
    said = " ".join(errors.split())
    if not said:
        return f"{who} {how} and said nothing."
    if len(said) > _LAST_WORDS:
        said = "..." + said[-_LAST_WORDS:]
    return f"{who} {how}: {said}"


def _facts(target, printed: str, seconds: float) -> str:
    """What the scorer said of one run, on one line: device, size, time, gate."""
    device = _DEVICE.search(printed)
    before, after = _BEFORE.search(printed), _AFTER.search(printed)
    line = (f"ESM-2 on {device.group(1)}" if device else "ESM-2") + \
        f", {target.aa} residues, {seconds:.1f} s"
    if before and after:
        line += f"; alignment {before.group(1)} before, {after.group(1)} after"
    return line


def _end(child: subprocess.Popen) -> None:
    """Make sure a child is gone: asked first, then not asked."""
    if child.poll() is not None:
        return
    child.terminate()
    try:
        child.wait(timeout=5)
    except subprocess.TimeoutExpired:
        child.kill()
        child.wait()


def _end_all(child: subprocess.Popen) -> None:
    """`_end` for a child that leads a session of its own: it and whatever it
    started. The modeller runs PyMOL and the scene importer, and ending the
    modeller alone would leave either running with nobody to read it."""
    for ask in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(child.pid, ask)
        except (ProcessLookupError, PermissionError):
            pass
        try:
            child.wait(timeout=5)
            break
        except subprocess.TimeoutExpired:
            continue
    child.wait()


class Scorer:
    """`worker.Score`, run where torch is: `score_local.py`, one protein a process.

    It raises what `scoring.score_with_esm` would have raised in this process,
    so `worker.score_next` decides as it does on Modal. The scorer's own
    `ValueError` comes back a `ValueError`, which is a refusal. Everything
    else is a `RuntimeError`: the scorer breaking, running past its time, or
    being stopped with the worker. `score_next` puts that bake back on the
    queue, up to `store.MAX_ATTEMPTS` times.
    """

    def __init__(self, settings: Settings, stop: threading.Event):
        self.settings = settings
        self.stop = stop

    def _command(self, *arguments: str) -> list:
        return [str(self.settings.esm_python), "-u", "-m", _DRIVER, *arguments]

    def _environment(self) -> dict:
        environment = {key: value for key, value in os.environ.items()
                       if key not in _CREDENTIALS}
        environment.update(
            HELIXPEEK_DATA=str(self.settings.data),
            HELIXPEEK_GB_CACHE=str(self.settings.genbank),
            # The 650M weights at the scorer's pinned revision are in the
            # Hugging Face cache, and a bake never waits on anything else.
            HF_HUB_OFFLINE="1",
            PYTHONUNBUFFERED="1",
        )
        return environment

    def unready(self) -> Optional[str]:
        """Why no protein could be scored here now. None when one could.

        Asked before a bake is claimed. A scorer that cannot start (a checkout
        without the driver, weights gone from the cache) would break on every
        protein alike, three times each, and leave each one's track refused.
        """
        try:
            done = subprocess.run(
                self._command("--check"), cwd=str(BACKEND), env=self._environment(),
                stdin=subprocess.DEVNULL, capture_output=True, text=True,
                timeout=_CHECK_TIMEOUT)
        except subprocess.TimeoutExpired:
            return f"The scorer's check had not finished after {_CHECK_TIMEOUT} s."
        except OSError as exc:
            return f"The scorer could not be started: {exc}"
        return None if done.returncode == 0 else _last_words(done.returncode, done.stderr)

    def _wait(self, child: subprocess.Popen, slug: str, printed: Path,
              started: float, report=None) -> Optional[str]:
        """Why the scorer was cut short, or None once it has ended by itself.

        The scorer's last progress line goes to the log every
        `_PROGRESS_EVERY` seconds, and to `report` (`worker.Report`), where
        there is one, every `_REPORT_EVERY`: a reader watching the protein
        build is told sooner than the log needs to be. `report` raises
        `worker.Stopped` once the reader who asked has stopped the bake, and
        the scorer is ended on the way out.
        """
        from pipeline.resolver.worker import progress_of

        report_at = started + _PROGRESS_EVERY
        reported = None
        tell_at = started + _REPORT_EVERY
        told = None
        while True:
            try:
                child.wait(timeout=_CHILD_POLL)
                return None
            except subprocess.TimeoutExpired:
                pass
            if self.stop.is_set():
                return STOPPED
            now = time.monotonic()
            if now - started > self.settings.score_timeout:
                return (f"ESM-2 was still scoring after {self.settings.score_timeout / 60:.0f} "
                        f"minutes and was stopped.")
            if report is not None and now >= tell_at:
                tell_at = now + _REPORT_EVERY
                found = progress_of(_read(printed))
                if found is not None and found[:2] != told:
                    told = found[:2]
                    report(*found)
                elif callable(getattr(report, "check", None)):
                    # No new line (the model loading, a long window): a stop
                    # is still heard within the interval.
                    report.check()
            if now >= report_at:
                report_at = now + _PROGRESS_EVERY
                progress = _PROGRESS.findall(_read(printed))
                if progress and progress[-1] != reported:
                    reported = progress[-1]
                    log.info("%s: %s", slug, reported)

    def __call__(self, target, record: bytes, report=None) -> bytes:
        from pipeline.resolver import score_local

        settings = self.settings
        if self.stop.is_set():
            raise RuntimeError(STOPPED)
        settings.scoring.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=settings.scoring, prefix=f"{target.slug}-") as folder:
            work = Path(folder)
            target_file, record_file = work / "target.pkl", work / "record.json"
            track_file, refusal_file = work / "track.json", work / "refusal.txt"
            printed_file, errors_file = work / "stdout.txt", work / "stderr.txt"
            target_file.write_bytes(pickle.dumps(target))
            record_file.write_bytes(record)

            started = time.monotonic()
            with open(printed_file, "wb") as printed, open(errors_file, "wb") as errors:
                child = subprocess.Popen(
                    self._command(str(target_file), str(record_file), str(track_file),
                                  str(refusal_file)),
                    cwd=str(BACKEND), env=self._environment(), stdin=subprocess.DEVNULL,
                    stdout=printed, stderr=errors)
                try:
                    cut = self._wait(child, target.slug, printed_file, started, report)
                finally:
                    _end(child)
            seconds = time.monotonic() - started
            if cut:
                raise RuntimeError(cut)

            code = child.returncode
            if code == score_local.REFUSED and refusal_file.exists():
                log.info("%s: %s", target.slug, _facts(target, _read(printed_file), seconds))
                raise ValueError(refusal_file.read_text(encoding="utf-8"))
            if code == 0 and track_file.exists():
                log.info("%s: %s", target.slug, _facts(target, _read(printed_file), seconds))
                return track_file.read_bytes()

            said = _read(errors_file)
            if said.strip():
                log.error("%s: the scorer's last words:\n%s", target.slug, said[-2000:].rstrip())
            if code == 0:
                raise RuntimeError("The scorer exited with status 0 and wrote no track.")
            if code == score_local.REFUSED:
                raise RuntimeError("The scorer exited with status 3 and gave no reason.")
            raise RuntimeError(_last_words(code, said))


# ------------------------------------------------------------ the modeller


class Modeller:
    """`worker.Model`, run where the structure bake's packages are:
    `model_local.py`, one protein a process.

    It raises what `worker.structure_next` reads: `worker.Unmodelled` where no
    model is drawn for the protein, which is a refusal, and a `RuntimeError`
    for everything else (the bake breaking, running past its time, or being
    stopped with the worker), which puts the bake back on the queue, up to
    `store.MAX_ATTEMPTS` times.
    """

    def __init__(self, settings: Settings, stop: threading.Event):
        self.settings = settings
        self.stop = stop

    def _command(self, *arguments: str) -> list:
        settings = self.settings
        return [str(settings.structure_python), "-u", "-m", _MODEL_DRIVER,
                "--app", str(settings.app), "--dart", str(settings.dart), *arguments]

    def _environment(self) -> dict:
        # A model is made from public data and handed back as files: the
        # baker has no use for the database or the storage key either.
        environment = {key: value for key, value in os.environ.items()
                       if key not in _CREDENTIALS}
        environment.update(PYTHONUNBUFFERED="1")
        return environment

    def unready(self) -> Optional[str]:
        """Why no model could be made here now. None when one could.

        Asked before a bake is claimed, as the scorer is: a baker that cannot
        start (no PyMOL, an app checkout on another flutter_scene) would break
        on every protein alike, and leave each one's fold page refused.
        """
        try:
            done = subprocess.run(
                self._command("--check"), cwd=str(BACKEND), env=self._environment(),
                stdin=subprocess.DEVNULL, capture_output=True, text=True,
                timeout=_CHECK_TIMEOUT)
        except subprocess.TimeoutExpired:
            return f"The modeller's check had not finished after {_CHECK_TIMEOUT} s."
        except OSError as exc:
            return f"The modeller could not be started: {exc}"
        return None if done.returncode == 0 else _last_words(
            done.returncode, done.stderr, "The modeller")

    def _wait(self, child: subprocess.Popen, started: float) -> Optional[str]:
        """Why the bake was cut short, or None once it has ended by itself."""
        while True:
            try:
                child.wait(timeout=_CHILD_POLL)
                return None
            except subprocess.TimeoutExpired:
                pass
            if self.stop.is_set():
                return STOPPED_MODELLING
            if time.monotonic() - started > MODEL_TIMEOUT:
                return (f"The model was still being made after {MODEL_TIMEOUT / 60:.0f} "
                        f"minutes and was stopped.")

    def __call__(self, target, record: bytes):
        from pipeline.resolver import model_local, worker

        settings = self.settings
        if self.stop.is_set():
            raise RuntimeError(STOPPED_MODELLING)
        settings.modelling.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=settings.modelling,
                                         prefix=f"{target.slug}-") as folder:
            work = Path(folder)
            job_file, out = work / "job.json", work / "out"
            printed_file, errors_file = work / "stdout.txt", work / "stderr.txt"
            job_file.write_text(json.dumps(worker.model_job(target, record)), encoding="utf-8")

            started = time.monotonic()
            with open(printed_file, "wb") as printed, open(errors_file, "wb") as errors:
                child = subprocess.Popen(
                    self._command(str(job_file), str(out)), cwd=str(BACKEND),
                    env=self._environment(), stdin=subprocess.DEVNULL,
                    stdout=printed, stderr=errors, start_new_session=True)
                try:
                    cut = self._wait(child, started)
                finally:
                    _end_all(child)
            seconds = time.monotonic() - started
            if cut:
                raise RuntimeError(cut)

            code = child.returncode
            refusal = out / "refusal.txt"
            if code == model_local.REFUSED and refusal.exists():
                raise worker.Unmodelled(refusal.read_text(encoding="utf-8"))
            glb, scene, described = (out / name for name in
                                     ("model.glb", "model.fsceneb", "described.json"))
            if code == 0 and glb.exists() and scene.exists() and described.exists():
                said = json.loads(described.read_text(encoding="utf-8"))
                built = worker.Modelled(glb.read_bytes(), scene.read_bytes(), said["chrome"],
                                        said["chains"], said["provenance"])
                log.info("%s: %s", target.slug, _model_facts(built, seconds))
                return built

            said = _read(errors_file)
            if said.strip():
                log.error("%s: the modeller's last words:\n%s", target.slug,
                          said[-2000:].rstrip())
            if code == 0:
                raise RuntimeError("The modeller exited with status 0 and wrote no model.")
            if code == model_local.REFUSED:
                raise RuntimeError("The modeller exited with status 3 and gave no reason.")
            raise RuntimeError(_last_words(code, said, "The modeller"))


def _model_facts(built, seconds: float) -> str:
    """What the modeller made, on one line: which model, how sure, how large."""
    said = built.provenance
    span = said.get("span") or ["?", "?"]
    return (f"{said.get('entry')} v{said.get('model_version')}, residues {span[0]}-{span[1]}, "
            f"mean pLDDT {said.get('mean_plddt')}, sampling {said.get('sampling')}, "
            f"scene {len(built.scene):,} B, {seconds:.1f} s")


# ------------------------------------------------------------ the evidence

_WINDOWS = re.compile(r"^\s*(\d+) window\(s\), ([\d,]+) bp", re.MULTILINE)
_GATE = re.compile(r"sequence gate: ([\d,]+ of [\d,]+) bases agree")


class Evidencer:
    """`worker.Avi`, run where the AVI bake's packages are: `impact_local.py`,
    one protein a process.

    It raises what `worker.impact_next` reads: `worker.Unplaced` where one of
    the bake's own gates said no, which is a refusal, and a `RuntimeError` for
    everything else (the Atlas or GENCODE's lookup not answering, the bake
    running past its time, or being stopped with the worker), which puts the
    bake back on the queue, up to `store.MAX_ATTEMPTS` times.
    """

    def __init__(self, settings: Settings, stop: threading.Event):
        self.settings = settings
        self.stop = stop

    def _command(self, *arguments: str) -> list:
        return [str(self.settings.impact_python), "-u", "-m", _IMPACT_DRIVER, *arguments]

    def _environment(self) -> dict:
        # Public data in, a file out, and the Atlas's key, which nothing else
        # the worker starts is given. Not the database or the storage key.
        settings = self.settings
        environment = {key: value for key, value in os.environ.items()
                       if key not in _CREDENTIALS}
        if settings.uv.is_absolute():
            environment["PATH"] = os.pathsep.join(
                [str(settings.uv.parent), environment.get("PATH", os.defpath)])
        if settings.alphagenome_key:
            environment["ALPHAGENOME_API_KEY"] = settings.alphagenome_key
        environment.update(PYTHONUNBUFFERED="1")
        return environment

    def unready(self) -> Optional[str]:
        """Why no AVI track could be made here now. None when one could.

        Asked before a bake is claimed, as the scorer is, and at the cost of
        one request to the Atlas: a key revoked, or a skill gone, would break
        every protein alike, three times each, and leave each one's AVI and
        ClinVar refused.
        """
        try:
            done = subprocess.run(
                self._command("--check"), cwd=str(BACKEND), env=self._environment(),
                stdin=subprocess.DEVNULL, capture_output=True, text=True,
                timeout=_CHECK_TIMEOUT)
        except subprocess.TimeoutExpired:
            return f"The AVI bake's check had not finished after {_CHECK_TIMEOUT} s."
        except OSError as exc:
            return f"The AVI bake could not be started: {exc}"
        return None if done.returncode == 0 else _last_words(
            done.returncode, done.stderr, "The AVI bake")

    def _wait(self, child: subprocess.Popen, started: float) -> Optional[str]:
        """Why the bake was cut short, or None once it has ended by itself."""
        while True:
            try:
                child.wait(timeout=_CHILD_POLL)
                return None
            except subprocess.TimeoutExpired:
                pass
            if self.stop.is_set():
                return STOPPED_PLACING
            if time.monotonic() - started > EVIDENCE_TIMEOUT:
                return (f"AlphaGenome's scores were still being fetched after "
                        f"{EVIDENCE_TIMEOUT / 60:.0f} minutes, and the bake was stopped.")

    def __call__(self, target, record: bytes) -> bytes:
        from pipeline.resolver import impact_local, worker

        settings = self.settings
        if self.stop.is_set():
            raise RuntimeError(STOPPED_PLACING)
        settings.fetching.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=settings.fetching,
                                         prefix=f"{target.slug}-") as folder:
            work = Path(folder)
            target_file, record_file, out = work / "target.pkl", work / "record.json", work / "out"
            printed_file, errors_file = work / "stdout.txt", work / "stderr.txt"
            target_file.write_bytes(pickle.dumps(target))
            record_file.write_bytes(record)

            started = time.monotonic()
            with open(printed_file, "wb") as printed, open(errors_file, "wb") as errors:
                # Its own session: the GENCODE lookup is its child, under `uv`.
                child = subprocess.Popen(
                    self._command(str(target_file), str(record_file), str(out),
                                  "--state", str(settings.impact_state)),
                    cwd=str(BACKEND), env=self._environment(), stdin=subprocess.DEVNULL,
                    stdout=printed, stderr=errors, start_new_session=True)
                try:
                    cut = self._wait(child, started)
                finally:
                    _end_all(child)
            seconds = time.monotonic() - started
            if cut:
                raise RuntimeError(cut)

            code = child.returncode
            refusal, track = out / "refusal.txt", out / "impact.json"
            if code == impact_local.REFUSED and refusal.exists():
                raise worker.Unplaced(refusal.read_text(encoding="utf-8"))
            if code == 0 and track.exists():
                log.info("%s: %s", target.slug, _avi_facts(_read(printed_file), track, seconds))
                return track.read_bytes()

            said = _read(errors_file)
            if said.strip():
                log.error("%s: the AVI bake's last words:\n%s", target.slug,
                          said[-2000:].rstrip())
            if code == 0:
                raise RuntimeError("The AVI bake exited with status 0 and wrote no track.")
            if code == impact_local.REFUSED:
                raise RuntimeError("The AVI bake exited with status 3 and gave no reason.")
            raise RuntimeError(_last_words(code, said, "The AVI bake"))


def _avi_facts(printed: str, track: Path, seconds: float) -> str:
    """What the AVI bake did, on one line: how much it asked the Atlas for,
    how the record's bases agreed with GRCh38's, how large, how long."""
    windows, gate = _WINDOWS.findall(printed), _GATE.search(printed)
    asked = sum(int(count) for count, _ in windows)
    line = f"AVI, {asked} Atlas window(s)" if windows else "AVI"
    if gate:
        line += f", sequence gate {gate.group(1)}"
    return line + f", {track.stat().st_size:,} B, {seconds:.1f} s"


class ClinVarBaker:
    """`worker.ClinVar` on the Mac: the twenty's baker in this process
    (`worker.clinvar_in_process`), its raw responses in the worker's own
    folder only while it runs, a stop heard between two of its batches, and
    what it found said in the log."""

    def __init__(self, settings: Settings, stop: threading.Event):
        self.settings = settings
        self.stop = stop

    def __call__(self, target) -> bytes:
        from pipeline.resolver import worker

        if self.stop.is_set():
            raise worker.Unfetched(worker.STOPPED_FETCHING)
        started = time.monotonic()
        payload = worker.clinvar_in_process(target, self.settings.fetching / "clinvar",
                                            stopped=self.stop.is_set)
        said = json.loads(payload)
        log.info("%s: ClinVar, %s of %s record(s) mapped, %d batch(es), %s B, %.1f s",
                 target.slug, f"{len(said['variants']):,}", f"{said['searched_records']:,}",
                 len(said["source_batches"]), f"{len(payload):,}", time.monotonic() - started)
        return payload


# ------------------------------------------------------------ one cycle


class Awake:
    """The Mac kept from idle sleep, from the first `hold()` until `release()`."""

    def __init__(self):
        self._child: Optional[subprocess.Popen] = None

    def hold(self) -> None:
        if self._child is not None:
            return
        try:
            self._child = subprocess.Popen(
                [*CAFFEINATE, "-w", str(os.getpid())], stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except OSError as exc:
            log.warning("caffeinate could not be started (%s); the Mac may sleep mid-build", exc)

    def release(self) -> None:
        if self._child is not None:
            _end(self._child)
            self._child = None


def clear_workspace(settings: Settings) -> None:
    """Empty the bakers' two directories and the hand-overs: the scorer's, the
    modeller's, and the AVI and ClinVar bakes'.

    On Modal they are `/tmp` and go with the container, so no protein meets
    what another left. Here they would stay, and the record builder reads its
    flat-file cache before it asks NCBI: a reply garbled once, which `resolve`
    takes for a refusal, would be read again when that gene is next asked for.
    What an AVI bake keeps between tries (`impact_state`) is not emptied.
    """
    for folder in (settings.data, settings.genbank, settings.scoring, settings.modelling,
                   settings.fetching):
        if BACKEND == folder or BACKEND in folder.parents:
            raise RuntimeError(f"{folder} is inside the repository; the worker empties only "
                               f"directories of its own.")
        shutil.rmtree(folder, ignore_errors=True)
        folder.mkdir(parents=True, exist_ok=True)


def _conninfo(url: str) -> str:
    """`DATABASE_URL` with the worker's connection settings under its own."""
    from psycopg.conninfo import conninfo_to_dict, make_conninfo

    return make_conninfo(**{**_CONNECTION, **conninfo_to_dict(url)})


def _say(name: str, what: str, outcome: dict, seconds: float) -> None:
    """One line for one protein: who, what became of it, how long, and why."""
    state = outcome["state"]
    told = {"done": f"done as {outcome.get('slug')}",
            "queued": "back on the queue"}.get(state, state)
    line = f"{name}: {what} {told} in {seconds:.1f} s"
    if outcome.get("reason"):
        line += f": {outcome['reason']}"
    log.log(logging.WARNING if state in ("queued", "failed") else logging.INFO, "%s", line)


def _resolve(conn, storage, stop: threading.Event, resolved: list,
             structures: bool = False, evidence: bool = False) -> int:
    """Reap, then resolve what is queued. Returns how many bakes wait.

    `structures` has each protein resolved queue its model too, and
    `evidence` its AVI and ClinVar bakes (`worker.sweep`): each said only by a
    cycle that can make them."""
    from pipeline.resolver import worker

    asked = {}
    if structures:
        asked["structures"] = True
    if evidence:
        asked["evidence"] = True
    waiting = 0
    for _ in range(SWEEP_LIMIT):
        started = time.monotonic()
        summary = worker.sweep(conn, storage, limit=1, **asked)
        waiting = summary["constraint_queued"]
        for outcome in summary["resolved"]:
            resolved.append(outcome)
            _say(outcome["gene"], "request", outcome, time.monotonic() - started)
        if (not summary["resolved"] or stop.is_set()
                or summary["resolved"][-1]["state"] == "queued"):
            break
    return waiting


def _score(conn, storage, scorer: Scorer, stop: threading.Event, scored: list) -> None:
    from pipeline.resolver import worker

    for _ in range(SCORE_LIMIT):
        if stop.is_set():
            break
        started = time.monotonic()
        done = worker.score_all(conn, storage, limit=1, score=scorer)
        for outcome in done:
            scored.append(outcome)
            _say(outcome["slug"], "constraint", outcome, time.monotonic() - started)
        if not done or done[-1]["state"] == "queued":
            break


def _model(conn, storage, modeller: Modeller, stop: threading.Event, modelled: list) -> None:
    from pipeline.resolver import worker

    for _ in range(SCORE_LIMIT):
        if stop.is_set():
            break
        started = time.monotonic()
        done = worker.structure_all(conn, storage, limit=1, model=modeller)
        for outcome in done:
            modelled.append(outcome)
            _say(outcome["slug"], "structure", outcome, time.monotonic() - started)
        if not done or done[-1]["state"] == "queued":
            break


def _place(conn, storage, evidencer: Evidencer, stop: threading.Event, placed: list) -> None:
    from pipeline.resolver import worker

    for _ in range(SCORE_LIMIT):
        if stop.is_set():
            break
        started = time.monotonic()
        done = worker.impact_all(conn, storage, limit=1, avi=evidencer)
        for outcome in done:
            placed.append(outcome)
            _say(outcome["slug"], "impact", outcome, time.monotonic() - started)
        if not done or done[-1]["state"] == "queued":
            break


def _fetch(conn, storage, clinvar: ClinVarBaker, stop: threading.Event, fetched: list) -> None:
    from pipeline.resolver import worker

    for _ in range(SCORE_LIMIT):
        if stop.is_set():
            break
        started = time.monotonic()
        done = worker.clinvar_all(conn, storage, limit=1, clinvar=clinvar)
        for outcome in done:
            fetched.append(outcome)
            _say(outcome["slug"], "clinvar", outcome, time.monotonic() - started)
        if not done or done[-1]["state"] == "queued":
            break


def cycle(settings: Settings, scorer: Scorer, stop: threading.Event,
          modeller: Optional[Modeller] = None, evidencer: Optional[Evidencer] = None,
          clinvar: Optional[ClinVarBaker] = None) -> dict:
    """One pass over the queue: `modal_app`'s `sweep` and then its `score`.

    Where the settings say this worker makes models (`structures`, which the
    Mac's does unless `.env` says `RESOLVER_STRUCTURES=0`), each protein
    resolved also has its model queued, and the models are made between the
    two: seconds each, so no fold page waits on an hour of ESM-2. `modeller`
    is what makes them, a `Modeller` unless one is handed in. Otherwise no
    model is queued and none is made.

    Where they say it makes the variant evidence (`evidence`, unless `.env`
    says `RESOLVER_EVIDENCE=0`), each protein resolved also has its AVI and
    ClinVar bakes queued, and they are made after the models and before the
    scoring, which waits for them (`store.WAITS_FOR`): AVI by `evidencer`,
    ClinVar by `clinvar`, an `Evidencer` and a `ClinVarBaker` unless handed in.

    Returns what was resolved and scored (and `modelled`, where models are
    made, and `impact` and `clinvar`, where the evidence is), and whether there
    was nothing to do. Raises what kept it from the queue, or took the queue
    away part-way: that is no protein's error, and the loop waits it out.
    """
    from pipeline.resolver import store

    if modeller is None and settings.structures:
        modeller = Modeller(settings, stop)
    if evidencer is None and settings.evidence:
        evidencer = Evidencer(settings, stop)
    if clinvar is None and settings.evidence:
        clinvar = ClinVarBaker(settings, stop)
    making_evidence = evidencer is not None and clinvar is not None
    resolved: list = []
    scored: list = []
    modelled: list = []
    placed: list = []
    fetched: list = []
    queued = waiting = models = evidence = 0
    awake = Awake()
    conn = store.connect(_conninfo(settings.database_url))
    try:
        with conn:
            storage = store.TrackStorage(settings.supabase_url, settings.service_key)
            queued = store.queued_requests(conn) + store.queued_bakes(conn, "constraint")
            if modeller is not None:
                queued += store.queued_bakes(conn, "structure")
            if making_evidence:
                queued += sum(store.queued_bakes(conn, kind) for kind in store.EVIDENCE)
            if queued:
                awake.hold()
            # Every cycle, idle or not: the reaper is in the sweep.
            waiting = _resolve(conn, storage, stop, resolved, structures=modeller is not None,
                               evidence=making_evidence)
            if modeller is not None and not stop.is_set():
                models = store.queued_bakes(conn, "structure")
                if models:
                    awake.hold()
                    unready = modeller.unready()
                    if unready:
                        log.error("%d model(s) left on the queue: no model can be made here "
                                  "now. %s", models, unready)
                    else:
                        _model(conn, storage, modeller, stop, modelled)
            if making_evidence and not stop.is_set():
                placing = store.queued_bakes(conn, "impact")
                evidence = placing + store.queued_bakes(conn, "clinvar")
                if evidence:
                    awake.hold()
                if placing:
                    unready = evidencer.unready()
                    if unready:
                        log.error("%d AVI bake(s) left on the queue, and their ClinVar bakes "
                                  "behind them: no AVI track can be made here now. %s",
                                  placing, unready)
                    else:
                        _place(conn, storage, evidencer, stop, placed)
                # ClinVar needs NCBI alone; it is claimed only once its AVI
                # bake is over (`store.WAITS_FOR`).
                if evidence and not stop.is_set():
                    _fetch(conn, storage, clinvar, stop, fetched)
            if waiting and not stop.is_set():
                awake.hold()
                unready = scorer.unready()
                if unready:
                    log.error("%d bake(s) left on the queue: no protein can be scored here "
                              "now. %s", waiting, unready)
                else:
                    _score(conn, storage, scorer, stop, scored)
    finally:
        awake.release()
        worked = (queued or waiting or models or evidence or resolved or scored or modelled
                  or placed or fetched)
        if worked:
            clear_workspace(settings)
    done = {"resolved": resolved, "scored": scored, "idle": not worked}
    if modeller is not None:
        done["modelled"] = modelled
    if making_evidence:
        done["impact"] = placed
        done["clinvar"] = fetched
    return done


# ------------------------------------------------------------ the loop


class Stop(threading.Event):
    """Set once the worker has been asked to leave. `by` is the signal that asked."""

    def __init__(self):
        super().__init__()
        self.by: Optional[str] = None


def _listen(stop: Stop) -> None:
    def asked(number, frame):
        # Only the flag: a handler that logged could find the log mid-line.
        stop.by = signal.Signals(number).name
        stop.set()

    for number in (signal.SIGTERM, signal.SIGINT):
        signal.signal(number, asked)


def _backoff(poll: float, failures: int) -> float:
    """The poll interval, doubled for each failure after the first, to ten minutes."""
    return min(poll * 2 ** min(failures - 1, 16), max(MAX_BACKOFF, poll))


def _transient(error: BaseException) -> bool:
    """Whether an error is the network's or the database's: one line's worth."""
    import psycopg

    return isinstance(error, (OSError, psycopg.OperationalError))


def serve(settings: Settings, stop: threading.Event, once: bool = False) -> int:
    """Cycle until stopped, or once. Returns the exit status.

    Nothing a cycle raises ends the loop. The database out of reach, DNS not
    back after a wake: each is logged and waited out, longer each time. A
    crash is something else, and launchd restarts it.
    """
    scorer = Scorer(settings, stop)
    failures = 0
    while not stop.is_set():
        try:
            done = cycle(settings, scorer, stop)
        except Exception as exc:  # noqa: BLE001 -- see the docstring
            failures += 1
            wait = _backoff(settings.poll, failures)
            log.error("the cycle stopped short: %s: %s%s", type(exc).__name__,
                      " ".join(str(exc).split()),
                      "" if once else f". Next try in {wait:.0f} s",
                      exc_info=not _transient(exc))
            if once:
                return 1
        else:
            failures = 0
            wait = settings.poll
            if done["idle"]:
                log.info("nothing queued")
            if once:
                return 0
        stop.wait(wait)
    return 0


# ------------------------------------------------------------ the process


def acquire(path: Path):
    """The lock that makes this the machine's one worker, or None where another
    holds it. `flock`, so it goes when the process does, however that ends."""
    import fcntl

    path.parent.mkdir(parents=True, exist_ok=True)
    handle = open(path, "a+", encoding="utf-8")
    try:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        handle.close()
        return None
    handle.seek(0)
    handle.truncate()
    handle.write(f"{os.getpid()}\n")
    handle.flush()
    return handle


def holder(path: Path) -> Optional[int]:
    """The pid the lock's holder wrote into it, where it can be read."""
    try:
        return int(path.read_text(encoding="utf-8").split()[0])
    except (OSError, ValueError, IndexError):
        return None


class _Formatter(logging.Formatter):
    """One timestamped line an event, with the credentials kept out of it.

    Whoever wrote the line: libpq quotes a connection string it cannot parse,
    password and all, in the error it raises, and that error is what a failed
    cycle logs. So this is done to the finished line, traceback included.
    """

    def __init__(self, secrets: tuple):
        super().__init__("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
        self._secrets = secrets

    def format(self, record: logging.LogRecord) -> str:
        return _scrub(super().format(record), self._secrets)


def _scrub(text: str, secrets: tuple) -> str:
    for secret in secrets:
        text = text.replace(secret, "[redacted]")
    return text


def start_logging(settings: Settings) -> list:
    """Open the log: a file the user reads, 1 MB and five older ones.

    A terminal gets the lines as well. launchd's copy of stderr does not: it is
    a file nothing rotates, kept for what Python says before this has run.
    Returns the handlers added.
    """
    settings.log_file.parent.mkdir(parents=True, exist_ok=True)
    handlers: list = [logging.handlers.RotatingFileHandler(
        settings.log_file, maxBytes=1_000_000, backupCount=5, encoding="utf-8")]
    if sys.stderr.isatty():
        handlers.append(logging.StreamHandler())
    formatter = _Formatter(settings.secrets)
    for handler in handlers:
        handler.setFormatter(formatter)
        log.addHandler(handler)
    log.setLevel(logging.INFO)
    return handlers


def take_directories(settings: Settings) -> None:
    """Name the bakers' directories before any baker is imported, and hold them to it.

    Named any later, a record would be scored into `pipeline/data/`, which
    holds the twenty's stored tracks: the bytes every sha256 proof is read
    against.
    """
    os.environ["NCBI_EMAIL"] = settings.ncbi_email
    os.environ["HELIXPEEK_DATA"] = str(settings.data)
    os.environ["HELIXPEEK_GB_CACHE"] = str(settings.genbank)

    from pipeline import paths
    from pipeline.mock import build_gene_record as builder

    if paths.DATA != settings.data or builder.CACHE != settings.genbank:
        raise Misconfigured(
            f"The bakers were imported before the worker named its directories: they write "
            f"to {paths.DATA} and {builder.CACHE}, not {settings.data} and {settings.genbank}.")


def queue(conn) -> str:
    """What waits, and what was last asked for, as lines to print. Reads only."""
    with conn.transaction():
        conn.execute("set transaction read only")
        requests = dict(conn.execute(
            "select state, count(*) from resolve_request group by state").fetchall())
        bakes = dict(conn.execute(
            "select state, count(*) from bake_job group by state").fetchall())
        held = conn.execute(
            "select slug, kind, state, attempts, error from bake_job "
            "where state in ('queued', 'running') order by requested_at, id").fetchall()
        last = conn.execute(
            "select requested_at, gene, state, reason from resolve_request "
            "order by requested_at desc, id desc limit 5").fetchall()

    def counted(found: dict, finished: tuple) -> str:
        return (f"{found.get('queued', 0)} queued, {found.get('running', 0)} running ("
                + ", ".join(f"{found.get(state, 0)} {state}" for state in finished) + ")")

    lines = ["requests  " + counted(requests, ("done", "refused", "failed")),
             "bakes     " + counted(bakes, ("done", "failed"))]
    for slug, kind, state, attempts, error in held:
        lines.append(f"  {slug} {kind}: {state}, claimed {attempts} time(s)"
                     + (f"; last: {error}" if error else ""))
    lines.append("last requests" if last else "no request yet")
    for requested, gene, state, reason in last:
        lines.append(f"  {requested.astimezone():%Y-%m-%d %H:%M}  {gene:<10} {state}"
                     + (f": {reason}" if reason else ""))
    return "\n".join(lines)


def _print_queue(settings: Settings) -> int:
    try:
        take_directories(settings)
        from pipeline.resolver import store

        with store.connect(_conninfo(settings.database_url)) as conn:
            lines = queue(conn)
    except Exception as exc:  # noqa: BLE001 -- said in a line, for `status` to show
        print("resolver worker: the queue could not be read: "
              + _scrub(" ".join(f"{type(exc).__name__}: {exc}".split()), settings.secrets),
              file=sys.stderr)
        return 1
    try:
        print(lines, flush=True)
    except BrokenPipeError:
        # `status | head`: the reader has what it wanted and has gone. Not a
        # queue that could not be read, and not for Python to report at exit.
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
    return 0


def main(arguments: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m pipeline.resolver.local_worker",
        description="The resolver's worker on a Mac: resolves queued requests, makes "
                    "their models, fetches their AVI and ClinVar tracks and scores their "
                    "ESM-2 tracks, until stopped.")
    how = parser.add_mutually_exclusive_group()
    how.add_argument("--once", action="store_true",
                     help="one cycle, then exit: 0 if it reached the queue, 1 if not")
    how.add_argument("--queue", action="store_true",
                     help="print what waits and what was last asked for; reads only")
    asked = parser.parse_args(arguments)

    try:
        settings = load_settings()
    except Misconfigured as exc:
        print(f"resolver worker: {exc}", file=sys.stderr)
        return _EX_CONFIG
    if asked.queue:
        return _print_queue(settings)

    held = acquire(settings.lock)
    if held is None:
        other = holder(settings.lock)
        print("resolver worker: another one is running" + (f" (pid {other})" if other else "")
              + "; this one is not needed.", file=sys.stderr)
        return 0

    start_logging(settings)
    try:
        take_directories(settings)
    except Misconfigured as exc:
        log.critical("%s", exc)
        return _EX_CONFIG
    socket.setdefaulttimeout(NETWORK_TIMEOUT)
    stop = Stop()
    _listen(stop)

    log.info("started (pid %d): asking the queue every %g s, scoring with %s, %s, %s, "
             "state in %s",
             os.getpid(), settings.poll, settings.esm_python,
             f"making models with {settings.structure_python}" if settings.structures
             else "making no models",
             f"fetching AVI with {settings.impact_python} and ClinVar here" if settings.evidence
             else "fetching no AVI or ClinVar", settings.state)
    for note in settings.notes:
        log.warning("%s", note)
    try:
        clear_workspace(settings)
        status = serve(settings, stop, once=asked.once)
    except BaseException:
        log.critical("crashed", exc_info=True)
        raise
    log.info("stopped%s", f" ({stop.by})" if stop.by else "")
    return status


if __name__ == "__main__":
    sys.exit(main())
