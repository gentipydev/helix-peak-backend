"""The resolver's worker on a Mac: what it is set up with, and where it runs ESM-2.

`worker.score_all` scores with whatever it is handed (`worker.Score`). On Modal
that is `scoring.score_with_esm` in the worker's own process, because one image
holds torch and psycopg both. On the dev Mac no environment does: the scorer's,
`pipeline/.esm-venv`, baked the twenty and is given nothing new. So a `Scorer`
is handed over instead, which runs `score_local.py` there as a subprocess, one
protein at a time, and comes back with what `score_with_esm` would have: the
track's bytes, or the scorer's own `ValueError`.

Nothing under `pipeline` is imported when this module is. `paths.DATA` and the
record builder's `CACHE` are fixed the first time they are imported, from the
environment as it is then, and the worker has to name its own directories
before that: never `pipeline/data/`, which holds the twenty's stored tracks.
"""

from __future__ import annotations

import logging
import math
import os
import pickle
import re
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Optional

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

log = logging.getLogger("resolver-worker")

# The uploader's three credentials and NCBI's contact address: the four keys of
# Modal's `helix-peak-resolver` secret (`modal_app.py`).
KEYS = ("DATABASE_URL", "SUPABASE_URL", "SUPABASE_SERVICE_KEY", "NCBI_EMAIL")

# What the scorer is started without. It reads a record and writes a track,
# both as files; it has no use for the database or the storage key.
_CREDENTIALS = ("DATABASE_URL", "SUPABASE_SERVICE_KEY")

# `store.STALE_BAKE` is 90 minutes: a sweep takes a bake that has run that long
# for a dead worker's and queues it again. A scoring let run past it could be
# claimed twice, so the timeout stays under it whatever is asked for.
SCORE_TIMEOUT_CEILING = 85 * 60

_DRIVER = "pipeline.resolver.score_local"
_CHECK_TIMEOUT = 120
_CHILD_POLL = 0.5
_PROGRESS_EVERY = 60
_LAST_WORDS = 220

STOPPED = "The worker was stopped while ESM-2 was scoring."

_DEVICE = re.compile(r"^Loading \S+ on (\w+) ", re.MULTILINE)
_BEFORE = re.compile(r"over the residue before:\s+([\d.]+%)")
_AFTER = re.compile(r"over the residue after:\s+([\d.]+%)")
_PROGRESS = re.compile(r"^Scored \d+/\d+ .*$", re.MULTILINE)


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
    def lock(self) -> Path:
        return self.state / "worker.lock"


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
        notes=tuple(notes),
    )


# ------------------------------------------------------------ the scorer


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""


def _last_words(code: int, errors: str) -> str:
    """Why the scorer broke, inside `worker._said`'s 300 characters, keeping the
    end of what it said: a traceback names the fault on its last line."""
    how = f"was ended by signal {-code}" if code < 0 else f"exited with status {code}"
    said = " ".join(errors.split())
    if not said:
        return f"The scorer {how} and said nothing."
    if len(said) > _LAST_WORDS:
        said = "..." + said[-_LAST_WORDS:]
    return f"The scorer {how}: {said}"


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
    """Make sure the scorer is gone: asked first, then not asked."""
    if child.poll() is not None:
        return
    child.terminate()
    try:
        child.wait(timeout=5)
    except subprocess.TimeoutExpired:
        child.kill()
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
              started: float) -> Optional[str]:
        """Why the scorer was cut short, or None once it has ended by itself."""
        report_at = started + _PROGRESS_EVERY
        reported = None
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
            if now >= report_at:
                report_at = now + _PROGRESS_EVERY
                progress = _PROGRESS.findall(_read(printed))
                if progress and progress[-1] != reported:
                    reported = progress[-1]
                    log.info("%s: %s", slug, reported)

    def __call__(self, target, record: bytes) -> bytes:
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
                    cut = self._wait(child, target.slug, printed_file, started)
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
