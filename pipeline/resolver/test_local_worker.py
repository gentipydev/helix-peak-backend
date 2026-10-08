"""The Mac's worker, offline: its settings, its bridge to the scorer, the driver.

Nothing here has torch, a database or the network. The scorer's environment is
stood in for by `fake_python`: a script this interpreter runs in the place of
`pipeline/.esm-venv/bin/python`, started with the command line the bridge gives
the real one.
"""

from __future__ import annotations

import json
import logging
import os
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
