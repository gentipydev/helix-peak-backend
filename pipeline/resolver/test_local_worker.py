"""The Mac's worker, offline: its settings, its bridge to the scorer, the driver,
one cycle, the loop, and the process around them.

Nothing here has torch, a database or the network. The scorer's environment is
stood in for by `fake_python`: a script this interpreter runs in the place of
`pipeline/.esm-venv/bin/python`, started with the command line the bridge gives
the real one. The queue is stood in for by `Queue`, and `caffeinate` by a
process that only waits. `test_worker_pg.py` runs a cycle on a real Postgres.
"""

from __future__ import annotations

import dataclasses
import json
import logging
import os
import re
import subprocess
import sys
import textwrap
import threading
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.resolver import local_worker, store, worker  # noqa: E402
from pipeline.targets import Region, Source, Target  # noqa: E402

pytestmark = pytest.mark.skipif(
    sys.platform == "win32", reason="the worker is a Mac's: launchd, flock, POSIX signals")

INS = Target(
    slug="ins", gene="INS", uniprot="P01308", display="Insulin",
    source=Source("NG_007114.1", protein_id="NP_000198.1", transcript_id="NM_000207.3"),
    structure=None, cleaved=True, chain_label="Insulin",
    regions=(Region("Signal peptide", "S", 1, 24, kept=False),
             Region("Insulin B chain", "B", 25, 54)),
    disulfides=((31, 96),), aa=110, mature_peptides=False,
)
RECORD = b'{"gene": "INS"}\n'

GATE_FAILED = "Alignment gate FAILED: 40.0% before, 45.0% after. No JSON written."


def fake_python(folder: Path, body: str, check: int = 0, check_says: str = "") -> Path:
    """A stand-in for the scorer's interpreter: this one, running `body`.

    The bridge starts it as it starts the real one, `-u -m <driver>` and then
    the four files or `--check`. `body` runs where the driver would, with the
    four paths under the driver's own names and the repository importable.
    `--check` exits `check`, having said `check_says` on stderr.
    """
    script = folder / "fake-python"
    script.write_text(
        f"#!{sys.executable}\n"
        "import os, sys\n"
        "from pathlib import Path\n"
        "sys.path.insert(0, os.getcwd())\n"
        "if sys.argv[-1] == '--check':\n"
        f"    sys.stderr.write({check_says!r})\n"
        f"    sys.exit({check})\n"
        "target_file, record_file, track_file, refusal_file = (Path(p) for p in sys.argv[-4:])\n"
        + textwrap.dedent(body), encoding="utf-8")
    script.chmod(0o755)
    return script


def real_driver(folder: Path, esm: str) -> Path:
    """`score_local.py` itself, in this interpreter, with ESM-2 replaced by `esm`:
    statements that leave `scoring.score_with_esm` bound to a stand-in."""
    return fake_python(folder, "\n".join([
        "import runpy",
        "from pipeline.resolver import scoring",
        textwrap.dedent(esm),
        "sys.argv = [sys.argv[0]] + sys.argv[-4:]",
        "runpy.run_module('pipeline.resolver.score_local', run_name='__main__')",
    ]))


SCORES = "scoring.score_with_esm = lambda target, record: f'{target.slug} {len(record)}'.encode()"
REFUSES = f"""
    def refuse(target, record):
        raise ValueError({GATE_FAILED!r})
    scoring.score_with_esm = refuse
"""
BREAKS = """
    def broken(target, record):
        raise RuntimeError("MPS backend out of memory")
    scoring.score_with_esm = broken
"""

# A `fake_python` body for the tests that go as far as the upload gate: what
# the scorer writes, as far as the gate reads it.
SCORES_A_TRACK = """
    import json, pickle
    target = pickle.loads(target_file.read_bytes())
    track_file.write_text(json.dumps({
        "gene": target.gene, "uniprot": target.uniprot,
        "model": "facebook/esm2_t33_650M_UR50D", "revision": "test",
        "method": "masked_marginals", "normalization": "minmax",
        "score_units": "natural_log_ratio_to_wildtype", "entropy_vocabulary": "full",
        "vocabulary_size": 33, "context": {"mode": "full", "residues": target.aa},
    }) + "\\n")
"""


def make_settings(tmp_path: Path, python: Path, **said) -> local_worker.Settings:
    return local_worker.Settings(
        database_url="postgresql://worker:hunter2-secret@db.invalid:5432/postgres",
        supabase_url="https://project.invalid", service_key="sb_secret_not_a_real_key",
        ncbi_email="tests@example.com", state=tmp_path / "state",
        log_file=tmp_path / "logs" / "resolver-worker.log", esm_python=python, **said)


@pytest.fixture
def spawned(monkeypatch):
    """Every process the worker starts, and a child polled often enough to test."""
    started = []
    popen = subprocess.Popen

    def recording(*args, **kwargs):
        started.append(popen(*args, **kwargs))
        return started[-1]

    monkeypatch.setattr(subprocess, "Popen", recording)
    monkeypatch.setattr(local_worker, "_CHILD_POLL", 0.02)
    return started


# ------------------------------------------------------------ settings

_FOUR = {"DATABASE_URL": "postgresql://worker:hunter2-secret@db.invalid/postgres",
         "SUPABASE_URL": "https://project.invalid",
         "SUPABASE_SERVICE_KEY": "sb_secret_not_a_real_key", "NCBI_EMAIL": "tests@example.com"}


def test_settings_are_the_environments_word_first_and_then_the_env_files(tmp_path):
    env_file = tmp_path / ".env"
    env_file.write_text(
        "NCBI_EMAIL=file@example.com\nDATABASE_URL=postgresql://from-the-file\n"
        "SUPABASE_URL=https://file.invalid\nSUPABASE_SERVICE_KEY=file-key\n"
        "RESOLVER_POLL_SECONDS=30\n")
    settings = local_worker.load_settings(
        {"DATABASE_URL": "postgresql://from-the-environment", "PATH": "/usr/bin"},
        env_file, tmp_path)

    assert settings.database_url == "postgresql://from-the-environment"
    assert (settings.supabase_url, settings.service_key, settings.ncbi_email) == \
        ("https://file.invalid", "file-key", "file@example.com")
    assert (settings.poll, settings.score_timeout, settings.notes) == (30, 3600, ())
    assert settings.esm_python == BACKEND / "pipeline" / ".esm-venv" / "bin" / "python"

    state = tmp_path / "Library" / "Application Support" / "HelixPeek" / "resolver-worker"
    assert (settings.state, settings.data, settings.genbank, settings.scoring, settings.lock) == \
        (state, state / "data", state / "genbank", state / "scoring", state / "worker.lock")
    assert settings.log_file == tmp_path / "Library" / "Logs" / "HelixPeek" / "resolver-worker.log"
    # Never where the twenty's stored tracks are.
    assert BACKEND not in settings.data.parents and BACKEND not in settings.genbank.parents


def test_a_missing_key_is_named_and_no_value_is_said(tmp_path):
    with pytest.raises(local_worker.Misconfigured) as raised:
        local_worker.load_settings(
            {"DATABASE_URL": _FOUR["DATABASE_URL"], "NCBI_EMAIL": "tests@example.com"},
            tmp_path / "no.env", tmp_path)
    said = str(raised.value)
    assert said.startswith("SUPABASE_URL, SUPABASE_SERVICE_KEY are not set in ")
    assert "hunter2" not in said and "DATABASE_URL" not in said


def test_a_timeout_past_the_reapers_patience_is_cut_short_and_said(tmp_path):
    settings = local_worker.load_settings(
        {**_FOUR, "RESOLVER_SCORE_TIMEOUT": "7200", "RESOLVER_ESM_PYTHON": "/opt/esm/python"},
        tmp_path / "no.env", tmp_path)
    assert settings.score_timeout == local_worker.SCORE_TIMEOUT_CEILING
    assert len(settings.notes) == 1 and "7200 s" in settings.notes[0]
    assert settings.esm_python == Path("/opt/esm/python")

    # The ceiling has to stay under what `reap` waits before it takes a
    # running bake for a dead worker's.
    minutes, unit = store.STALE_BAKE.split()
    assert unit == "minutes" and local_worker.SCORE_TIMEOUT_CEILING < int(minutes) * 60


@pytest.mark.parametrize("name", ["RESOLVER_POLL_SECONDS", "RESOLVER_SCORE_TIMEOUT"])
@pytest.mark.parametrize("value", ["soon", "0", "-5", "nan", "inf"])
def test_an_override_that_is_not_a_length_of_time_is_refused(tmp_path, name, value):
    with pytest.raises(local_worker.Misconfigured, match=name):
        local_worker.load_settings({**_FOUR, name: value}, tmp_path / "no.env", tmp_path)


def test_importing_the_worker_imports_no_baker_and_sets_nothing():
    # `paths.DATA` is fixed when `pipeline.paths` is first imported, so the
    # module must be importable without that having happened.
    code = ("import os, sys; before = dict(os.environ); import pipeline.resolver.local_worker; "
            "print(sorted(n for n in sys.modules if n.startswith('pipeline')), "
            "dict(os.environ) == before)")
    done = subprocess.run([sys.executable, "-c", code], cwd=BACKEND, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == \
        "['pipeline', 'pipeline.resolver', 'pipeline.resolver.local_worker'] True"


# ------------------------------------------------------------ the bridge


def test_the_bridge_hands_over_the_protein_and_returns_the_track(tmp_path, monkeypatch, caplog):
    monkeypatch.setenv("DATABASE_URL", _FOUR["DATABASE_URL"])
    monkeypatch.setenv("SUPABASE_SERVICE_KEY", _FOUR["SUPABASE_SERVICE_KEY"])
    python = fake_python(tmp_path, """
        import hashlib, json, pickle
        print("Loading facebook/esm2_t33_650M_UR50D@08e4846e on mps (float32).")
        print("  over the residue before: 95.6%")
        print("  over the residue after:  94.7%")
        track_file.write_text(json.dumps({
            "target": repr(pickle.loads(target_file.read_bytes())),
            "record": hashlib.sha256(record_file.read_bytes()).hexdigest(),
            "started_with": sys.argv[1:4], "cwd": os.getcwd(),
            "environment": {key: os.environ.get(key) for key in (
                "HELIXPEEK_DATA", "HELIXPEEK_GB_CACHE", "HF_HUB_OFFLINE",
                "DATABASE_URL", "SUPABASE_SERVICE_KEY")},
        }))
    """)
    settings = make_settings(tmp_path, python)
    with caplog.at_level(logging.INFO, logger="resolver-worker"):
        track = json.loads(local_worker.Scorer(settings, threading.Event())(INS, RECORD))

    import hashlib
    assert track["target"] == repr(INS)
    assert track["record"] == hashlib.sha256(RECORD).hexdigest()
    assert track["started_with"] == ["-u", "-m", "pipeline.resolver.score_local"]
    assert Path(track["cwd"]) == BACKEND
    # Where the bakers write, the weights offline, and no credential.
    assert track["environment"] == {
        "HELIXPEEK_DATA": str(settings.data), "HELIXPEEK_GB_CACHE": str(settings.genbank),
        "HF_HUB_OFFLINE": "1", "DATABASE_URL": None, "SUPABASE_SERVICE_KEY": None}
    assert list(settings.scoring.iterdir()) == []
    assert [r.getMessage() for r in caplog.records][0].startswith(
        "ins: ESM-2 on mps, 110 residues, ")
    assert caplog.records[0].getMessage().endswith("; alignment 95.6% before, 94.7% after")


def test_the_scorers_refusal_comes_back_as_the_refusal_it_was(tmp_path):
    python = fake_python(tmp_path, f"""
        refusal_file.write_text({GATE_FAILED!r})
        sys.exit(3)
    """)
    settings = make_settings(tmp_path, python)
    with pytest.raises(ValueError) as raised:
        local_worker.Scorer(settings, threading.Event())(INS, RECORD)

    assert type(raised.value) is ValueError and str(raised.value) == GATE_FAILED
    # So `score_next` tells it as it tells the scorer's own, in-process.
    assert worker._refusal(raised.value).startswith(
        "ESM-2 prefers the residue that is there to the one before it 40.0% of the time, "
        "and to the one after it 45.0%.")
    assert list(settings.scoring.iterdir()) == []


def test_a_scorer_that_breaks_is_an_error_carrying_its_last_words(tmp_path, caplog):
    python = fake_python(tmp_path, r"""
        sys.stderr.write("Traceback (most recent call last):\n  File \"score_protein.py\", "
                         "line 291, in score_protein\n    " + "logits " * 80 + "\n"
                         "RuntimeError: MPS backend out of memory\n")
        sys.exit(1)
    """)
    with caplog.at_level(logging.INFO, logger="resolver-worker"):
        with pytest.raises(RuntimeError) as raised:
            local_worker.Scorer(make_settings(tmp_path, python), threading.Event())(INS, RECORD)

    said = str(raised.value)
    assert not isinstance(raised.value, ValueError)  # tried again, not refused
    assert said.startswith("The scorer exited with status 1: ...")
    assert said.endswith("RuntimeError: MPS backend out of memory")
    # The line that names the fault survives what a track's reason may hold.
    assert worker._said(raised.value) == said
    assert "Traceback (most recent call last):" in caplog.text


def test_a_scorer_that_runs_past_its_time_is_killed(tmp_path, spawned):
    python = fake_python(tmp_path, "import time; time.sleep(120)")
    settings = make_settings(tmp_path, python, score_timeout=0.2)
    with pytest.raises(RuntimeError, match="ESM-2 was still scoring after .* and was stopped"):
        local_worker.Scorer(settings, threading.Event())(INS, RECORD)
    assert len(spawned) == 1 and spawned[0].poll() is not None
    assert list(settings.scoring.iterdir()) == []


def test_a_stop_ends_the_scorer_and_starts_no_other(tmp_path, spawned):
    python = fake_python(tmp_path, "import time; time.sleep(120)")
    stop = threading.Event()
    scorer = local_worker.Scorer(make_settings(tmp_path, python), stop)
    threading.Timer(0.2, stop.set).start()
    with pytest.raises(RuntimeError) as raised:
        scorer(INS, RECORD)
    assert str(raised.value) == local_worker.STOPPED
    assert not isinstance(raised.value, ValueError)  # back on the queue, not refused
    assert len(spawned) == 1 and spawned[0].poll() is not None

    with pytest.raises(RuntimeError, match="stopped"):
        scorer(INS, RECORD)
    assert len(spawned) == 1


def test_a_long_scoring_says_how_far_it_is(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(local_worker, "_CHILD_POLL", 0.02)
    monkeypatch.setattr(local_worker, "_PROGRESS_EVERY", 0.05)
    python = fake_python(tmp_path, """
        import time
        print("Scored 25/110 (1.0s, 0.1 min left)", flush=True)
        time.sleep(0.6)
        track_file.write_bytes(b"{}")
    """)
    with caplog.at_level(logging.INFO, logger="resolver-worker"):
        local_worker.Scorer(make_settings(tmp_path, python), threading.Event())(INS, RECORD)
    said = [r.getMessage() for r in caplog.records]
    assert said.count("ins: Scored 25/110 (1.0s, 0.1 min left)") == 1


@pytest.mark.parametrize("body, said", [
    ("sys.exit(0)", "The scorer exited with status 0 and wrote no track."),
    ("sys.exit(3)", "The scorer exited with status 3 and gave no reason."),
    ("os.kill(os.getpid(), 9)", "The scorer was ended by signal 9 and said nothing."),
])
def test_an_ending_outside_the_protocol_is_an_error(tmp_path, body, said):
    with pytest.raises(RuntimeError) as raised:
        local_worker.Scorer(make_settings(tmp_path, fake_python(tmp_path, body)),
                            threading.Event())(INS, RECORD)
    assert str(raised.value) == said


def test_a_scorer_that_can_start_is_ready_and_one_that_cannot_says_why(tmp_path):
    stop = threading.Event()
    ready = fake_python(tmp_path, "sys.exit(1)")
    assert local_worker.Scorer(make_settings(tmp_path, ready), stop).unready() is None

    no_weights = fake_python(tmp_path, "sys.exit(1)", check=1,
                             check_says="esm2_t33_650M_UR50D is not whole in the cache.\n")
    assert local_worker.Scorer(make_settings(tmp_path, no_weights), stop).unready() == \
        "The scorer exited with status 1: esm2_t33_650M_UR50D is not whole in the cache."

    absent = local_worker.Scorer(make_settings(tmp_path, tmp_path / "no-venv" / "python"), stop)
    assert absent.unready().startswith("The scorer could not be started: ")


# ------------------------------------------------------------ the driver


def _hand_over(tmp_path: Path) -> list:
    import pickle

    work = tmp_path / "hand-over"
    work.mkdir()
    (work / "target.pkl").write_bytes(pickle.dumps(INS))
    (work / "record.json").write_bytes(RECORD)
    return [str(work / name) for name in ("target.pkl", "record.json", "track.json", "refusal.txt")]


def _drive(python: Path, files: list, tmp_path: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [str(python), "-u", "-m", "pipeline.resolver.score_local", *files], cwd=BACKEND,
        env={**os.environ, "HELIXPEEK_DATA": str(tmp_path / "data")},
        capture_output=True, text=True)


def test_the_driver_writes_the_track_and_exits_0(tmp_path):
    files = _hand_over(tmp_path)
    done = _drive(real_driver(tmp_path, SCORES), files, tmp_path)
    assert done.returncode == 0, done.stderr
    assert Path(files[2]).read_bytes() == b"ins 16" and not Path(files[3]).exists()


def test_the_driver_writes_the_scorers_refusal_and_exits_3(tmp_path):
    files = _hand_over(tmp_path)
    done = _drive(real_driver(tmp_path, REFUSES), files, tmp_path)
    assert done.returncode == 3, done.stderr
    assert Path(files[3]).read_text() == GATE_FAILED and not Path(files[2]).exists()


def test_the_driver_breaking_is_neither(tmp_path):
    files = _hand_over(tmp_path)
    done = _drive(real_driver(tmp_path, BREAKS), files, tmp_path)
    assert done.returncode == 1
    assert done.stderr.rstrip().endswith("RuntimeError: MPS backend out of memory")
    assert not Path(files[2]).exists() and not Path(files[3]).exists()


def test_a_hand_over_it_cannot_read_is_the_driver_breaking_not_a_refusal(tmp_path):
    files = _hand_over(tmp_path)
    Path(files[0]).write_bytes(b"not a pickle")
    done = _drive(real_driver(tmp_path, SCORES), files, tmp_path)
    assert done.returncode not in (0, 3)
    assert not Path(files[2]).exists() and not Path(files[3]).exists()


def test_the_driver_scores_nothing_without_a_data_directory_of_its_own(tmp_path):
    # Unset, `paths.DATA` is `pipeline/data/`: the twenty's stored tracks.
    files = _hand_over(tmp_path)
    environment = {key: value for key, value in os.environ.items() if key != "HELIXPEEK_DATA"}
    for arguments in (files, files[:1]):
        done = subprocess.run(
            [sys.executable, "-m", "pipeline.resolver.score_local", *arguments], cwd=BACKEND,
            env=environment, capture_output=True, text=True)
        assert done.returncode == 2
    assert not Path(files[2]).exists() and not Path(files[3]).exists()


def test_the_bridge_and_the_driver_agree(tmp_path):
    stop = threading.Event()
    scores = local_worker.Scorer(make_settings(tmp_path, real_driver(tmp_path, SCORES)), stop)
    assert scores(INS, RECORD) == b"ins 16"

    refuses = local_worker.Scorer(make_settings(tmp_path, real_driver(tmp_path, REFUSES)), stop)
    with pytest.raises(ValueError) as raised:
        refuses(INS, RECORD)
    assert str(raised.value) == GATE_FAILED

    breaks = local_worker.Scorer(make_settings(tmp_path, real_driver(tmp_path, BREAKS)), stop)
    with pytest.raises(RuntimeError, match="MPS backend out of memory$") as raised:
        breaks(INS, RECORD)
    assert not isinstance(raised.value, ValueError)


def test_the_driver_imports_where_only_the_scorer_is():
    # `pipeline/.esm-venv` has no psycopg and no Biopython, and gets neither.
    code = textwrap.dedent('''
        import sys

        class Absent:
            def find_spec(self, name, path=None, target=None):
                if name.split(".")[0] in {"psycopg", "psycopg_pool", "Bio", "dotenv", "fastapi",
                                          "pydantic", "pydantic_settings", "modal", "app"}:
                    raise ImportError(name + " is not in the scorer's environment")

        sys.meta_path.insert(0, Absent())
        import pipeline.resolver.score_local
        import pipeline.constraint.score_protein
        print(sorted(name for name in sys.modules if name.startswith("pipeline")))
    ''')
    done = subprocess.run([sys.executable, "-c", code], cwd=BACKEND, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == str([
        "pipeline", "pipeline.constraint", "pipeline.constraint.score_protein", "pipeline.paths",
        "pipeline.resolver", "pipeline.resolver.score_local", "pipeline.resolver.scoring",
        "pipeline.targets"])


def test_score_with_esm_puts_the_record_where_the_scorer_reads_it(tmp_path, monkeypatch):
    from pipeline.constraint import score_protein
    from pipeline.resolver import scoring

    calls = []

    def scored(target, destination, context):
        calls.append((target, destination, context))
        assert (tmp_path / "assets" / "mock" / "gene_ins.json").read_bytes() == RECORD
        destination.parent.mkdir(parents=True)
        destination.write_bytes(b"track\n")

    monkeypatch.setattr(score_protein, "DATA", tmp_path)
    monkeypatch.setattr(score_protein, "score_protein", scored)
    assert scoring.score_with_esm(INS, RECORD) == b"track\n"
    assert calls == [(INS, tmp_path / "assets" / "constraint" / "ins_esm_constraint.json", 1022)]


# ------------------------------------------------------------ one cycle


class Connection:
    """All the cycle asks of a connection here: to be closed when it is done."""

    def __init__(self):
        self.closed = False

    def __enter__(self):
        return self

    def __exit__(self, *raised):
        self.closed = True


class Queue:
    """`store` and `worker` as the cycle calls them, over two lists.

    `requests` and `bakes` are the outcomes still to hand out, oldest first.
    One that goes back on the queue ("queued") stays at the front, as a
    requeued row does, so a cycle that asked again at once would get it again.
    """

    def __init__(self, monkeypatch, requests=(), bakes=()):
        self.requests, self.bakes = list(requests), list(bakes)
        self.connections, self.sweeps, self.scorings = [], [], []
        self.conninfo = None
        self.during = lambda: None  # run while a protein is being worked
        self.seen = True            # whether the first look sees what is queued
        monkeypatch.setattr(store, "connect", self._connect)
        monkeypatch.setattr(store, "TrackStorage", lambda url, key: ("storage", url, key))
        monkeypatch.setattr(store, "queued_requests",
                            lambda conn: len(self.requests) if self.seen else 0)
        monkeypatch.setattr(store, "queued_bakes",
                            lambda conn, kind: len(self.bakes) if self.seen else 0)
        monkeypatch.setattr(worker, "sweep", self._sweep)
        monkeypatch.setattr(worker, "score_all", self._score_all)

    def _connect(self, conninfo):
        self.conninfo = conninfo
        self.connections.append(Connection())
        return self.connections[-1]

    def _take(self, waiting):
        if not waiting:
            return []
        self.during()
        if waiting[0]["state"] == "queued":
            return [waiting[0]]
        return [waiting.pop(0)]

    def _sweep(self, conn, storage, *, limit):
        self.sweeps.append(limit)
        return {"resolved": self._take(self.requests), "constraint_queued": len(self.bakes)}

    def _score_all(self, conn, storage, *, limit, score):
        self.scorings.append((limit, score))
        return self._take(self.bakes)


class Ready:
    """A scorer that could start, or that says why it could not."""

    def __init__(self, unready=None):
        self.why, self.asked = unready, 0

    def unready(self):
        self.asked += 1
        return self.why


def _request(gene="INS", state="done", **said):
    return {"id": 1, "gene": gene, "state": state, "slug": gene.lower(), **said}


def _bake(slug="ins", state="ready", **said):
    return {"id": 1, "slug": slug, "state": state, **said}


@pytest.fixture
def caffeinate(monkeypatch, spawned):
    """`caffeinate` as a process that only waits, and every one that was started."""
    monkeypatch.setattr(local_worker, "CAFFEINATE",
                        (sys.executable, "-c", "import time; time.sleep(120)"))
    return spawned


@pytest.fixture
def settings(tmp_path):
    return make_settings(tmp_path, tmp_path / "no-scorer")


def test_an_idle_cycle_reaps_and_keeps_nothing_awake(settings, monkeypatch, caffeinate):
    queue, scorer = Queue(monkeypatch), Ready()
    assert local_worker.cycle(settings, scorer, threading.Event()) == \
        {"resolved": [], "scored": [], "idle": True}
    # The sweep runs all the same: the reaper is in it.
    assert queue.sweeps == [1]
    assert queue.scorings == [] and scorer.asked == 0
    assert caffeinate == []
    assert [connection.closed for connection in queue.connections] == [True]


def test_a_cycle_with_work_keeps_the_mac_awake_until_it_is_done(
        settings, monkeypatch, caffeinate, caplog):
    from psycopg.conninfo import conninfo_to_dict

    queue, scorer = Queue(monkeypatch, requests=[_request()], bakes=[_bake()]), Ready()
    left = settings.genbank / "NG_007114.1.gb"
    left.parent.mkdir(parents=True)
    left.write_text("LOCUS")
    awake = []
    queue.during = lambda: awake.append([child.poll() for child in caffeinate])
    with caplog.at_level(logging.INFO, logger="resolver-worker"):
        done = local_worker.cycle(settings, scorer, threading.Event())

    assert [o["state"] for o in done["resolved"]] == ["done"]
    assert [o["state"] for o in done["scored"]] == ["ready"]
    assert done["idle"] is False
    # One protein a call, so a stop can land between two; then none is left.
    assert queue.sweeps == [1, 1]
    assert queue.scorings == [(1, scorer), (1, scorer)]
    # One caffeinate, tied to this process, there through both and gone after.
    assert len(caffeinate) == 1
    assert caffeinate[0].args[-2:] == ["-w", str(os.getpid())]
    assert awake == [[None], [None]]
    assert caffeinate[0].poll() is not None
    # What the bakers left is not there for the next protein to meet.
    assert not left.exists() and settings.genbank.is_dir()

    said = [record.getMessage() for record in caplog.records]
    assert re.fullmatch(r"INS: request done as ins in \d+\.\d s", said[0])
    assert re.fullmatch(r"ins: constraint ready in \d+\.\d s", said[1])

    # The connection is the URL's, probed while a protein is scored.
    opened = conninfo_to_dict(queue.conninfo)
    assert (opened["host"], opened["user"], opened["password"]) == \
        ("db.invalid", "worker", "hunter2-secret")
    assert (opened["keepalives_idle"], opened["connect_timeout"]) == ("60", "20")


def test_a_url_that_sets_its_own_connection_settings_keeps_them():
    from psycopg.conninfo import conninfo_to_dict

    opened = conninfo_to_dict(local_worker._conninfo(
        "postgresql://worker:pw@db.invalid/postgres?sslmode=require&keepalives_idle=5"))
    assert (opened["sslmode"], opened["keepalives_idle"], opened["keepalives_count"]) == \
        ("require", "5", "10")


def test_a_bake_with_no_request_is_work_and_so_is_one_the_reaper_puts_back(
        settings, monkeypatch, caffeinate):
    queue = Queue(monkeypatch, bakes=[_bake()])
    done = local_worker.cycle(settings, Ready(), threading.Event())
    assert [o["state"] for o in done["scored"]] == ["ready"] and queue.sweeps == [1]
    assert len(caffeinate) == 1 and caffeinate[0].poll() is not None

    # Nothing was queued at the first look; the sweep's reaper then put a dead
    # worker's bake back. It is scored awake like any other.
    queue = Queue(monkeypatch, bakes=[_bake()])
    queue.seen = False
    done = local_worker.cycle(settings, Ready(), threading.Event())
    assert [o["state"] for o in done["scored"]] == ["ready"] and done["idle"] is False
    assert len(caffeinate) == 2 and caffeinate[1].poll() is not None


def test_no_bake_is_claimed_while_the_scorer_cannot_start(
        settings, monkeypatch, caffeinate, caplog):
    queue = Queue(monkeypatch, bakes=[_bake(), _bake("b2m")])
    scorer = Ready("The scorer exited with status 1: No module named 'torch'")
    done = local_worker.cycle(settings, scorer, threading.Event())

    assert done == {"resolved": [], "scored": [], "idle": False}
    assert queue.scorings == [] and len(queue.bakes) == 2
    assert ("2 bake(s) left on the queue: no protein can be scored here now. "
            "The scorer exited with status 1: No module named 'torch'") in caplog.text
    assert caffeinate[0].poll() is not None


def test_a_stop_lands_between_two_proteins(settings, monkeypatch, caffeinate):
    stop = threading.Event()
    queue = Queue(monkeypatch, bakes=[_bake(), _bake("b2m")])
    queue.during = stop.set
    done = local_worker.cycle(settings, Ready(), stop)
    assert [o["slug"] for o in done["scored"]] == ["ins"]
    assert len(queue.scorings) == 1 and [b["slug"] for b in queue.bakes] == ["b2m"]

    # While a request is resolved: the next is left, and no scoring is started.
    stop = threading.Event()
    queue = Queue(monkeypatch, requests=[_request(), _request("B2M")], bakes=[_bake()])
    queue.during = stop.set
    scorer = Ready()
    done = local_worker.cycle(settings, scorer, stop)
    assert [o["gene"] for o in done["resolved"]] == ["INS"]
    assert queue.sweeps == [1] and queue.scorings == [] and scorer.asked == 0


def test_a_protein_put_back_on_the_queue_waits_for_the_next_cycle(
        settings, monkeypatch, caffeinate, caplog):
    # `sweep` and `score_all` with their own limits would claim it again at
    # once, and a third time, and that is all the tries a protein gets.
    queue = Queue(monkeypatch, requests=[
        _request(state="queued", reason="rest.uniprot.org did not answer"), _request("B2M")])
    local_worker.cycle(settings, Ready(), threading.Event())
    assert queue.sweeps == [1] and len(queue.requests) == 2
    assert re.search(r"INS: request back on the queue in \d+\.\d s: rest.uniprot.org did not "
                     r"answer", caplog.text)

    queue = Queue(monkeypatch, bakes=[
        _bake(state="queued", reason=local_worker.STOPPED), _bake("b2m")])
    local_worker.cycle(settings, Ready(), threading.Event())
    assert len(queue.scorings) == 1 and len(queue.bakes) == 2


def test_a_cycle_that_breaks_still_lets_the_mac_sleep(settings, monkeypatch, caffeinate):
    queue = Queue(monkeypatch, bakes=[_bake()])

    def gone():
        raise OSError("the network went with the lid")

    queue.during = gone
    with pytest.raises(OSError, match="went with the lid"):
        local_worker.cycle(settings, Ready(), threading.Event())
    assert len(caffeinate) == 1 and caffeinate[0].poll() is not None
    assert queue.connections[0].closed


# ------------------------------------------------------------ the loop


class Waits(local_worker.Stop):
    """A stop that comes after so many waits, remembering how long each was."""

    def __init__(self, limit):
        super().__init__()
        self.limit, self.asked = limit, []

    def wait(self, timeout=None):
        self.asked.append(timeout)
        if len(self.asked) >= self.limit:
            self.set()
        return self.is_set()


IDLE = {"resolved": [], "scored": [], "idle": True}


def _cycles(monkeypatch, outcomes):
    """`cycle` as a list of what each one does: returns it, or raises it."""
    outcomes = list(outcomes)

    def cycle(settings, scorer, stop):
        outcome = outcomes.pop(0)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    monkeypatch.setattr(local_worker, "cycle", cycle)


def test_a_database_out_of_reach_is_waited_out_longer_each_time(settings, monkeypatch, caplog):
    import psycopg

    unreachable = psycopg.OperationalError("connection to server at db.invalid failed")
    _cycles(monkeypatch, [unreachable] * 6 + [IDLE, IDLE])
    stop = Waits(8)
    with caplog.at_level(logging.INFO, logger="resolver-worker"):
        assert local_worker.serve(settings, stop) == 0

    assert stop.asked == [60, 120, 240, 480, 600, 600, 60, 60]
    said = [record.getMessage() for record in caplog.records]
    assert said[0] == ("the cycle stopped short: OperationalError: connection to server at "
                       "db.invalid failed. Next try in 60 s")
    assert said[5].endswith("Next try in 600 s") and said[6:] == ["nothing queued"] * 2
    # The network's fault is one line. Nobody needs its traceback.
    assert not any(record.exc_info for record in caplog.records)


def test_an_error_nobody_expected_is_logged_whole_and_the_loop_goes_on(
        settings, monkeypatch, caplog):
    _cycles(monkeypatch, [KeyError("constraint_queued"), IDLE])
    stop = Waits(2)
    assert local_worker.serve(settings, stop) == 0
    assert stop.asked == [60, 60]
    assert "KeyError: 'constraint_queued'" in caplog.text and "Traceback" in caplog.text


def test_once_is_one_cycle_and_its_status_says_whether_the_queue_was_reached(
        settings, monkeypatch):
    stop = Waits(1)
    _cycles(monkeypatch, [IDLE])
    assert local_worker.serve(settings, stop, once=True) == 0
    _cycles(monkeypatch, [OSError("no route to host")])
    assert local_worker.serve(settings, stop, once=True) == 1
    assert stop.asked == []


def test_a_stop_ends_the_loop_without_another_cycle(settings, monkeypatch):
    stop = local_worker.Stop()
    ran = []

    def cycle(settings, scorer, stopping):
        ran.append(True)
        stop.set()
        return IDLE

    monkeypatch.setattr(local_worker, "cycle", cycle)
    assert local_worker.serve(settings, stop) == 0
    assert ran == [True]


def test_the_wait_never_grows_past_ten_minutes_or_shrinks_below_the_poll():
    assert [local_worker._backoff(60, n) for n in (1, 2, 3, 4, 5, 6, 5000)] == \
        [60, 120, 240, 480, 600, 600, 600]
    assert [local_worker._backoff(900, n) for n in (1, 2, 3)] == [900, 900, 900]


# ------------------------------------------------------------ the process


def test_one_worker_holds_the_lock(tmp_path):
    path = tmp_path / "state" / "worker.lock"
    first = local_worker.acquire(path)
    assert first is not None and local_worker.holder(path) == os.getpid()
    assert local_worker.acquire(path) is None
    first.close()
    second = local_worker.acquire(path)
    assert second is not None
    second.close()


def test_the_workspace_is_emptied_and_is_only_ever_the_workers_own(
        settings, tmp_path, monkeypatch):
    for folder in (settings.data, settings.genbank, settings.scoring):
        (folder / "assets").mkdir(parents=True)
        (folder / "assets" / "left.json").write_text("{}")
    settings.lock.write_text("1\n")
    local_worker.clear_workspace(settings)
    assert [list(folder.iterdir()) for folder in
            (settings.data, settings.genbank, settings.scoring)] == [[], [], []]
    assert settings.lock.exists()

    # Never a directory of the repository's: `pipeline/data/` holds the
    # twenty's stored tracks. A stand-in repository, so nothing real is at risk.
    repository = tmp_path / "repository"
    stored = repository / "pipeline" / "data" / "assets" / "insulin.json"
    stored.parent.mkdir(parents=True)
    stored.write_text("{}")
    monkeypatch.setattr(local_worker, "BACKEND", repository)
    with pytest.raises(RuntimeError, match="inside the repository"):
        local_worker.clear_workspace(
            dataclasses.replace(settings, state=repository / "pipeline"))
    assert stored.exists()


def test_the_bakers_are_held_to_the_workers_directories(settings, monkeypatch):
    from pipeline import paths
    from pipeline.mock import build_gene_record as builder

    for key in ("NCBI_EMAIL", "HELIXPEEK_DATA", "HELIXPEEK_GB_CACHE"):
        monkeypatch.delenv(key, raising=False)  # so each is put back afterwards
    # Both were imported long ago, by this test run: what the check is for.
    with pytest.raises(local_worker.Misconfigured, match="before the worker named"):
        local_worker.take_directories(settings)
    assert os.environ["HELIXPEEK_DATA"] == str(settings.data)
    assert os.environ["HELIXPEEK_GB_CACHE"] == str(settings.genbank)

    monkeypatch.setattr(paths, "DATA", settings.data)
    monkeypatch.setattr(builder, "CACHE", settings.genbank)
    local_worker.take_directories(settings)


def test_the_secrets_are_the_key_and_the_url_and_its_password_however_written(settings):
    assert settings.secrets == (settings.database_url, "sb_secret_not_a_real_key",
                                "hunter2-secret")
    for url, password in [
        ("postgresql://postgres.ref:p%40ss%2Fword@pooler.invalid:5432/postgres", "p@ss/word"),
        ("host=db.invalid user=worker password=hunter2-secret dbname=postgres", "hunter2-secret"),
        ("postgres//worker:hunter2-secret@db.invalid/postgres", "hunter2-secret"),
    ]:
        assert password in dataclasses.replace(settings, database_url=url).secrets


def test_no_line_of_the_log_holds_a_credential(settings):
    handlers = local_worker.start_logging(settings)
    try:
        local_worker.log.info("connecting to %s", settings.database_url)
        local_worker.log.warning("the key is %s", settings.service_key)
        try:
            # What libpq says of a connection string it cannot parse.
            raise ValueError(f'missing "=" after "{settings.database_url}" in connection info')
        except ValueError:
            local_worker.log.error("the cycle stopped short", exc_info=True)
    finally:
        for handler in handlers:
            local_worker.log.removeHandler(handler)
            handler.close()

    text = settings.log_file.read_text()
    assert "hunter2-secret" not in text and settings.service_key not in text
    lines = text.splitlines()
    assert re.fullmatch(r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\d INFO    connecting to \[redacted\]",
                        lines[0])
    assert lines[1].endswith("WARNING the key is [redacted]")
    assert lines[-1] == 'ValueError: missing "=" after "[redacted]" in connection info'


_NOWHERE = "postgresql://worker:hunter2-secret@/postgres?host=/nonexistent-helixpeek"
_MISTYPED = "postgres//worker:hunter2-secret@db.invalid/postgres"


def _run(tmp_path: Path, *arguments: str, database_url: str):
    """The worker as launchd starts it, with a home of its own and no database:
    every key is named here, so the repository's `.env` is never what is read."""
    home = tmp_path / "home"
    home.mkdir(parents=True, exist_ok=True)
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith(("RESOLVER_", "HELIXPEEK_"))}
    environment.update(_FOUR, DATABASE_URL=database_url, HOME=str(home),
                       RESOLVER_ESM_PYTHON=str(tmp_path / "no-scorer"))
    done = subprocess.run(
        [sys.executable, "-m", "pipeline.resolver.local_worker", *arguments], cwd=BACKEND,
        env=environment, capture_output=True, text=True, timeout=120)
    return home, done


def test_run_once_with_no_database_it_says_so_and_writes_only_under_its_own_home(tmp_path):
    home, done = _run(tmp_path, "--once", database_url=_NOWHERE)
    assert done.returncode == 1, done.stderr

    state = home / "Library" / "Application Support" / "HelixPeek" / "resolver-worker"
    lines = (home / "Library" / "Logs" / "HelixPeek" / "resolver-worker.log").read_text() \
        .splitlines()
    # It started, which it does only once the bakers write where it says.
    assert "started (pid " in lines[0] and lines[0].endswith(f"state in {state}")
    assert "the cycle stopped short: OperationalError: " in lines[1]
    assert lines[-1].endswith("stopped") and len(lines) == 3
    assert sorted(path.name for path in state.iterdir()) == \
        ["data", "genbank", "scoring", "worker.lock"]
    # Not a terminal, as under launchd: the lines are the file's alone.
    assert "started (pid " not in done.stderr
    assert "hunter2-secret" not in "\n".join(lines) + done.stderr + done.stdout


def test_a_mistyped_database_url_is_not_copied_into_the_log(tmp_path):
    home, done = _run(tmp_path, "--once", database_url=_MISTYPED)
    assert done.returncode == 1, done.stderr
    text = (home / "Library" / "Logs" / "HelixPeek" / "resolver-worker.log").read_text()
    # libpq quotes the string back, and the traceback carries it.
    assert "ProgrammingError" in text and '"[redacted]"' in text
    assert "hunter2-secret" not in text + done.stderr + done.stdout


def test_a_second_worker_says_it_is_not_needed_and_leaves(tmp_path):
    home = tmp_path / "home"
    held = local_worker.acquire(
        home / "Library" / "Application Support" / "HelixPeek" / "resolver-worker" / "worker.lock")
    try:
        _, done = _run(tmp_path, "--once", database_url=_NOWHERE)
    finally:
        held.close()
    assert done.returncode == 0
    assert f"another one is running (pid {os.getpid()}); this one is not needed." in done.stderr
    assert not (home / "Library" / "Logs").exists()


def test_a_queue_that_cannot_be_read_is_one_line_with_no_credential(tmp_path):
    home, done = _run(tmp_path, "--queue", database_url=_MISTYPED)
    assert done.returncode == 1
    assert done.stderr.startswith("resolver worker: the queue could not be read: ")
    assert "hunter2-secret" not in done.stderr + done.stdout
    # Reading the queue opens no log and makes no directory of the worker's.
    assert not (home / "Library" / "Logs").exists()
    assert not (home / "Library" / "Application Support").exists()
