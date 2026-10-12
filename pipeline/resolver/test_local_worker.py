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


def fake_model_python(folder: Path, body: str, check: int = 0, check_says: str = "") -> Path:
    """A stand-in for the structure bake's interpreter: this one, running `body`.

    The modeller starts it as it starts the real one: `-u -m <driver> --app
    <checkout> --dart <dart>`, and then `--check`, or the job's file and the
    folder the model goes into. `body` runs where the driver would, with
    those two as `job_file` and `out` and the repository importable.
    `--check` exits `check`, having said `check_says` on stderr.
    """
    script = folder / "fake-structure-python"
    script.write_text(
        f"#!{sys.executable}\n"
        "import json, os, sys\n"
        "from pathlib import Path\n"
        "sys.path.insert(0, os.getcwd())\n"
        "if sys.argv[-1] == '--check':\n"
        f"    sys.stderr.write({check_says!r})\n"
        f"    sys.exit({check})\n"
        "job_file, out = (Path(p) for p in sys.argv[-2:])\n"
        + textwrap.dedent(body), encoding="utf-8")
    script.chmod(0o755)
    return script


# A `fake_model_python` body for the tests that go as far as the worker's gate:
# what `model_local.py` leaves, as far as the gate reads it.
MAKES_A_MODEL = """
    import struct
    job = json.loads(job_file.read_text())
    record = (job_file.parent / job["record"]).read_bytes()
    nodes = ["plddtVeryHigh", "bonds"]
    document = json.dumps({
        "asset": {"version": "2.0"}, "scenes": [{"nodes": [0, 1]}],
        "nodes": [{"name": name, "mesh": index} for index, name in enumerate(nodes)],
        "meshes": [{"primitives": []}, {"primitives": []}]}).encode()
    document += b" " * (-len(document) % 4)
    entry = "AF-" + job["accession"] + "-F1"
    out.mkdir(parents=True, exist_ok=True)
    (out / "model.glb").write_bytes(
        struct.pack("<4sII", b"glTF", 2, 20 + len(document))
        + struct.pack("<I4s", len(document), b"JSON") + document)
    (out / "model.fsceneb").write_bytes(b"a compiled scene")
    (out / "described.json").write_text(json.dumps({
        "chrome": {"pdb": entry, "modelled": None, "label": "the predicted fold",
                   "count": len(job["protein"]), "unit": "residues",
                   "sentence": "AlphaFold prediction, mean pLDDT 90.0.",
                   "semantics": "AlphaFold's predicted fold. Drag to turn it."},
        "chains": [{"node": "plddtVeryHigh", "tint": "plddtVeryHigh"},
                   {"node": "bonds", "tint": "cysteine"}],
        "provenance": {"pdb": entry, "nodes": nodes, "entry": entry, "model_version": 6,
                       "span": [1, len(job["protein"])], "mean_plddt": 90.0, "sampling": 8},
    }))
    # The whole search found no experimental entry: AlphaFold's, said so.
    (out / "outcome.json").write_text(json.dumps({"made": "alphafold", "selection": {
        "policy": 1, "inputs": "the-inputs-digest", "outcome": "none", "refused": [],
        "search": {"query": "UniProt " + job["accession"], "release": "2026-10-07",
                   "entities": 0, "sha256": "e3b0c442"},
        "record_bytes": len(record)}}))
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


def test_a_reader_is_told_how_far_esm2_has_got_sooner_than_the_log(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(local_worker, "_CHILD_POLL", 0.02)
    monkeypatch.setattr(local_worker, "_REPORT_EVERY", 0.05)
    python = fake_python(tmp_path, """
        import time
        print("Scored 25/110 (1.0s, 0.1 min left)", flush=True)
        time.sleep(0.4)
        print("Scored 50/110 (2.0s, 0.0 min left)", flush=True)
        time.sleep(0.4)
        track_file.write_bytes(b"{}")
    """)
    told = []
    with caplog.at_level(logging.INFO, logger="resolver-worker"):
        local_worker.Scorer(make_settings(tmp_path, python), threading.Event())(
            INS, RECORD, lambda *said: told.append(said))
    # Each line once, however often the file was read.
    assert [(done, total) for done, total, _ in told] == [(25, 110), (50, 110)]
    assert told[0][2] == pytest.approx(1.0 / 25 * 85)
    # The log keeps its own minute.
    assert not [r for r in caplog.records if "Scored" in r.getMessage()]


def test_a_stop_heard_between_lines_ends_the_scorer(tmp_path, monkeypatch, spawned):
    monkeypatch.setattr(local_worker, "_REPORT_EVERY", 0.05)
    # The model loading: nothing printed for a long while.
    python = fake_python(tmp_path, "import time; time.sleep(120)")

    class Report:
        checks = 0

        def __call__(self, done, total, left):
            raise AssertionError("no line was printed")

        def check(self):
            Report.checks += 1
            if Report.checks >= 3:
                raise worker.Stopped()

    with pytest.raises(worker.Stopped):
        local_worker.Scorer(make_settings(tmp_path, python), threading.Event())(
            INS, RECORD, Report())
    assert Report.checks == 3
    assert len(spawned) == 1 and spawned[0].poll() is not None


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

    def __init__(self, monkeypatch, requests=(), bakes=(), impacts=(), clinvars=()):
        self.requests, self.bakes = list(requests), list(bakes)
        self.impacts, self.clinvars = list(impacts), list(clinvars)
        self.connections, self.sweeps, self.scorings = [], [], []
        self.order = []             # what was worked, in the order it was
        self.conninfo = None
        self.during = lambda: None  # run while a protein is being worked
        self.seen = True            # whether the first look sees what is queued
        monkeypatch.setattr(store, "connect", self._connect)
        monkeypatch.setattr(store, "TrackStorage", lambda url, key: ("storage", url, key))
        monkeypatch.setattr(store, "queued_requests",
                            lambda conn: len(self.requests) if self.seen else 0)
        monkeypatch.setattr(store, "queued_bakes",
                            lambda conn, kind: len(self._waiting(kind)) if self.seen else 0)
        monkeypatch.setattr(worker, "sweep", self._sweep)
        monkeypatch.setattr(worker, "score_all", self._score_all)
        monkeypatch.setattr(worker, "impact_all", self._impact_all)
        monkeypatch.setattr(worker, "clinvar_all", self._clinvar_all)

    def _waiting(self, kind):
        return {"impact": self.impacts, "clinvar": self.clinvars}.get(kind, self.bakes)

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

    def _sweep(self, conn, storage, *, limit, **asked):
        self.sweeps.append(limit)
        self.order.append(("sweep", asked))
        return {"resolved": self._take(self.requests), "constraint_queued": len(self.bakes)}

    def _score_all(self, conn, storage, *, limit, score):
        self.scorings.append((limit, score))
        self.order.append(("constraint", limit))
        return self._take(self.bakes)

    def _impact_all(self, conn, storage, *, limit, avi):
        self.order.append(("impact", avi))
        return self._take(self.impacts)

    def _clinvar_all(self, conn, storage, *, limit, clinvar):
        self.order.append(("clinvar", clinvar))
        return self._take(self.clinvars)


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
        ["data", "fetching", "genbank", "modelling", "scoring", "worker.lock"]
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


# ------------------------------------------------------------ the modeller

A_RECORD = json.dumps({"gene": "INS", "protein": {"translation": "M" * 110}}).encode()
# INS's row as the resolver writes it, as far as the app's check reads it.
_A_ROW = {"slug": "ins", "display": "Insulin", "gene": "INS", "uniprot": "P01308",
          "accession": "NG_007114.1", "summary": "A hormone.", "residues": 110, "exons": 3,
          "chains": 2, "bridges": 1, "catalog_order": None, "chain_name": None,
          "structure": None}


def _modeller(tmp_path, body, stop=None, **check):
    folder = tmp_path / "structure"
    folder.mkdir(exist_ok=True)
    settings = make_settings(
        tmp_path, tmp_path / "no-esm-python",
        structure_python=fake_model_python(folder, body, **check),
        dart=Path("/opt/dart-sdk/bin/dart"), flutter=Path("/opt/flutter/bin/flutter"),
        app=tmp_path / "helix-peek")
    return local_worker.Modeller(settings, stop or threading.Event())


def test_this_worker_makes_models_unless_its_env_says_not(tmp_path):
    settings = local_worker.load_settings(_FOUR, tmp_path / "no.env", tmp_path)
    assert settings.structures is True
    assert settings.structure_python == \
        BACKEND / "pipeline" / "structure" / "venv" / "bin" / "python"
    assert settings.app == BACKEND.parent / "helix-peek"
    assert settings.modelling == settings.state / "modelling"
    assert BACKEND not in settings.modelling.parents

    for word in ("0", "false", "No", "off"):
        said = local_worker.load_settings(
            {**_FOUR, "RESOLVER_STRUCTURES": word}, tmp_path / "no.env", tmp_path)
        assert said.structures is False
    told = local_worker.load_settings(
        {**_FOUR, "RESOLVER_STRUCTURE_PYTHON": "/opt/bake/python", "RESOLVER_DART": "/opt/dart",
         "RESOLVER_FLUTTER": "/opt/flutter/bin/flutter", "RESOLVER_APP_DIR": "/src/app"},
        tmp_path / "no.env", tmp_path)
    assert (told.structure_python, told.dart, told.flutter, told.app) == \
        (Path("/opt/bake/python"), Path("/opt/dart"), Path("/opt/flutter/bin/flutter"),
         Path("/src/app"))
    # A worker built by hand, as every test before this one builds it, makes none.
    assert make_settings(tmp_path, tmp_path / "python").structures is False


def test_dart_is_the_sdks_own_binary_where_a_flutter_checkout_has_one(tmp_path, monkeypatch):
    # `flutter/bin/dart` is a script that takes Flutter's start-up lock first.
    sdk = tmp_path / "flutter" / "bin" / "cache" / "dart-sdk" / "bin" / "dart"
    sdk.parent.mkdir(parents=True)
    sdk.write_text("")
    assert local_worker.default_dart(tmp_path) == sdk
    monkeypatch.setattr(local_worker.shutil, "which", lambda name: "/usr/local/bin/dart")
    assert local_worker.default_dart(tmp_path / "elsewhere") == Path("/usr/local/bin/dart")
    monkeypatch.setattr(local_worker.shutil, "which", lambda name: None)
    assert local_worker.default_dart(tmp_path / "elsewhere") == Path("dart")


def test_flutter_is_the_home_checkouts_where_there_is_one(tmp_path, monkeypatch):
    flutter = tmp_path / "flutter" / "bin" / "flutter"
    flutter.parent.mkdir(parents=True)
    flutter.write_text("")
    assert local_worker.default_flutter(tmp_path) == flutter
    monkeypatch.setattr(local_worker.shutil, "which", lambda name: "/usr/local/bin/flutter")
    assert local_worker.default_flutter(tmp_path / "elsewhere") == Path("/usr/local/bin/flutter")
    monkeypatch.setattr(local_worker.shutil, "which", lambda name: None)
    assert local_worker.default_flutter(tmp_path / "elsewhere") == Path("flutter")


def test_a_models_time_is_inside_the_reapers_patience_for_one():
    minutes, unit = store.STALE_STRUCTURE.split()
    assert unit == "minutes" and local_worker.MODEL_TIMEOUT < int(minutes) * 60


def test_the_modeller_hands_over_the_job_and_returns_the_model(tmp_path, monkeypatch, caplog):
    monkeypatch.setenv("DATABASE_URL", _FOUR["DATABASE_URL"])
    monkeypatch.setenv("SUPABASE_SERVICE_KEY", _FOUR["SUPABASE_SERVICE_KEY"])
    modeller = _modeller(tmp_path, MAKES_A_MODEL + """
    (out / "seen.json").write_text(json.dumps({
        "started_with": sys.argv[1:10], "cwd": os.getcwd(), "job": job,
        "record": record.decode(),
        "environment": {key: os.environ.get(key) for key in (
            "DATABASE_URL", "SUPABASE_SERVICE_KEY", "PYTHONUNBUFFERED")}}))
    (Path(os.environ["SEEN"])).write_text((out / "seen.json").read_text())
""")
    monkeypatch.setenv("SEEN", str(tmp_path / "seen.json"))
    ask = worker.FoldAsk(_A_ROW, {"none": None, "progress": None, "active": None,
                                  "refresh": False})
    with caplog.at_level(logging.INFO, logger="resolver-worker"):
        folded = modeller(INS, A_RECORD, ask)

    assert isinstance(folded, worker.Folded) and folded.source == "alphafold"
    assert (folded.folding, folded.audit) == (None, None)
    assert folded.selection["outcome"] == "none" and folded.selection["record_bytes"] == \
        len(A_RECORD)
    built = folded.model
    assert worker.glb_nodes(built.glb) == ["plddtVeryHigh", "bonds"]
    assert built.scene == b"a compiled scene"
    assert built.chrome["pdb"] == "AF-P01308-F1" and built.chrome["count"] == 110
    assert [chain["tint"] for chain in built.chains] == ["plddtVeryHigh", "cysteine"]

    seen = json.loads((tmp_path / "seen.json").read_text())
    assert seen["started_with"] == [
        "-u", "-m", "pipeline.resolver.model_local",
        "--app", str(tmp_path / "helix-peek"), "--dart", "/opt/dart-sdk/bin/dart",
        "--flutter", "/opt/flutter/bin/flutter"]
    assert Path(seen["cwd"]) == BACKEND
    # The job is the protein as its record and row have it (`worker.model_job`),
    # with what earlier searches decided, and the record beside it.
    assert seen["job"] == worker.model_job(INS, A_RECORD, ask)
    assert seen["job"]["prior"] == ask.prior and seen["record"] == A_RECORD.decode()
    assert seen["job"]["row"]["slug"] == "ins" and seen["job"]["target"]["slug"] == "ins"
    assert seen["job"]["kept"] == [[25, 54]] and seen["job"]["disulfides"] == [[31, 96]]
    # Public data in, files out: the baker is given neither credential.
    assert seen["environment"] == {
        "DATABASE_URL": None, "SUPABASE_SERVICE_KEY": None, "PYTHONUNBUFFERED": "1"}
    said = [record.getMessage() for record in caplog.records]
    assert len(said) == 1 and said[0].startswith(
        "ins: AF-P01308-F1 v6, residues 1-110, mean pLDDT 90.0, sampling 8, scene 16 B, ")
    # Nothing of the hand-over is left behind.
    assert list(modeller.settings.modelling.iterdir()) == []


def test_no_model_for_a_protein_is_told_as_that_and_says_why(tmp_path):
    modeller = _modeller(tmp_path, """
        out.mkdir(parents=True, exist_ok=True)
        (out / "refusal.txt").write_text(
            "AlphaFold DB has no model of proteins over 2,700 residues; this one has 4,834.")
        sys.exit(3)
    """)
    with pytest.raises(worker.Unmodelled) as none:
        modeller(INS, A_RECORD)
    assert str(none.value) == ("AlphaFold DB has no model of proteins over 2,700 residues; "
                               "this one has 4,834.")
    assert none.value.selection is None
    # And with what the search decided, which is kept beside the refusal.
    said = {"policy": 1, "inputs": "d", "outcome": "none", "refused": ["1ABC:A@1"]}
    modeller = _modeller(tmp_path, f"""
        out.mkdir(parents=True, exist_ok=True)
        (out / "refusal.txt").write_text("No model.")
        (out / "outcome.json").write_text(json.dumps({{"made": None, "selection": {said!r}}}))
        sys.exit(3)
    """)
    with pytest.raises(worker.Unmodelled) as none:
        modeller(INS, A_RECORD)
    assert (str(none.value), none.value.selection) == ("No model.", said)


def test_a_search_cut_short_is_told_as_unfinished_with_how_far_it_got(tmp_path, caplog):
    said = {"policy": 1, "inputs": "d", "outcome": "incomplete",
            "search": {"sha256": "abc"}, "refused": ["1ABC:A@1", "2DEF:B@3"]}
    modeller = _modeller(tmp_path, f"""
        out.mkdir(parents=True, exist_ok=True)
        (out / "outcome.json").write_text(json.dumps(
            {{"made": None, "unfinished": True, "selection": {said!r}}}))
        sys.exit(4)
    """)
    with caplog.at_level(logging.INFO, logger="resolver-worker"):
        with pytest.raises(worker.Unfinished) as unfinished:
            modeller(INS, A_RECORD)
    assert unfinished.value.selection == said
    assert "goes on next cycle" in caplog.records[-1].getMessage()


A_FOLDING = {"slug": "ins", "gene": "INS", "uniprot": "P01308", "pdb": "1ABC",
             "schema_version": 2, "chains": [{"node": "chainA"}], "bridges": []}


def test_an_experimental_fold_is_handed_over_with_its_track_and_audit(tmp_path, caplog):
    modeller = _modeller(tmp_path, f"""
        import struct
        document = json.dumps({{"asset": {{"version": "2.0"}}, "scenes": [{{"nodes": [0]}}],
                               "nodes": [{{"name": "chainA", "mesh": 0}}],
                               "meshes": [{{"primitives": []}}]}}).encode()
        document += b" " * (-len(document) % 4)
        out.mkdir(parents=True, exist_ok=True)
        (out / "model.glb").write_bytes(struct.pack("<4sII", b"glTF", 2, 20 + len(document))
                                        + struct.pack("<I4s", len(document), b"JSON") + document)
        (out / "model.fsceneb").write_bytes(b"a compiled scene")
        (out / "folding.json").write_text(json.dumps({A_FOLDING!r}))
        (out / "audit.json").write_text(json.dumps({{"entry": "1ABC"}}))
        (out / "described.json").write_text(json.dumps({{
            "chrome": {{"pdb": "1ABC"}}, "chains": [{{"node": "chainA", "tint": "mature1"}}],
            "provenance": {{"pdb": "1ABC", "entry": "1ABC", "source": "PDB",
                            "method": ["X-RAY DIFFRACTION"], "resolution": 1.5,
                            "chains": [{{"placed": 51}}], "sampling": 8}}}}))
        (out / "outcome.json").write_text(json.dumps({{"made": "pdb", "selection": {{
            "outcome": "chosen", "candidate": "1ABC:A@1"}}}}))
    """)
    with caplog.at_level(logging.INFO, logger="resolver-worker"):
        folded = modeller(INS, A_RECORD)
    assert folded.source == "pdb" and json.loads(folded.folding) == A_FOLDING
    assert json.loads(folded.audit) == {"entry": "1ABC"}
    assert folded.selection == {"outcome": "chosen", "candidate": "1ABC:A@1"}
    assert caplog.records[-1].getMessage().startswith(
        "ins: PDB 1ABC (X-RAY DIFFRACTION, 1.5 A), 1 chain(s), 51 residues placed, sampling 8, "
        "scene 16 B, with its fold, ")
    # An experimental fold with no track to animate it is no fold.
    modeller = _modeller(tmp_path, """
        out.mkdir(parents=True, exist_ok=True)
        for name in ("model.glb", "model.fsceneb"):
            (out / name).write_bytes(b"glTF")
        (out / "described.json").write_text(json.dumps({"chrome": {}, "chains": [],
                                                        "provenance": {}}))
        (out / "outcome.json").write_text(json.dumps({"made": "pdb", "selection": {}}))
    """)
    with pytest.raises(RuntimeError, match="left no folding track or audit"):
        modeller(INS, A_RECORD)


@pytest.mark.parametrize("body, said", [
    ("sys.exit(0)", "The modeller exited with status 0 and wrote no model."),
    ("out.mkdir(parents=True); (out / 'outcome.json').write_text('{\"made\": \"alphafold\"}')",
     "The modeller exited with status 0 and wrote no model."),
    ("sys.exit(3)", "The modeller exited with status 3 and gave no reason."),
    ("sys.exit(4)", "The modeller exited with status 4 and said nothing of its search."),
    ("sys.stderr.write('pymol failed:\\nSegmentation fault'); sys.exit(1)",
     "The modeller exited with status 1: pymol failed: Segmentation fault"),
    ("os.kill(os.getpid(), 9)", "The modeller was ended by signal 9 and said nothing."),
])
def test_a_modeller_that_ends_any_other_way_broke(tmp_path, body, said, spawned):
    modeller = _modeller(tmp_path, body)
    with pytest.raises(RuntimeError) as broke:
        modeller(INS, A_RECORD)
    assert not isinstance(broke.value, (worker.Unmodelled, worker.Unfinished))
    assert str(broke.value) == said


def test_a_stop_ends_the_modeller_and_its_bake_can_go_back(tmp_path, spawned):
    stop = threading.Event()
    modeller = _modeller(tmp_path, """
        import time
        Path(os.environ["MODELLING_STARTED"]).write_text("")
        time.sleep(120)
    """, stop=stop)
    started = tmp_path / "modelling-has-started"
    os.environ["MODELLING_STARTED"] = str(started)

    def stop_once_started():
        import time
        while not started.exists():
            time.sleep(0.01)
        stop.set()

    try:
        threading.Thread(target=stop_once_started, daemon=True).start()
        with pytest.raises(RuntimeError) as cut:
            modeller(INS, A_RECORD)
    finally:
        del os.environ["MODELLING_STARTED"]
    assert str(cut.value) == local_worker.STOPPED_MODELLING
    assert spawned[-1].poll() is not None            # ended, not left running
    with pytest.raises(RuntimeError, match="stopped while a model was being made"):
        modeller(INS, A_RECORD)                      # and none is started once stopping


def test_a_stop_ends_what_the_modeller_started_too(tmp_path, spawned):
    # PyMOL and the scene importer are the modeller's own children. Ended alone,
    # the modeller would leave either running, with nobody to read it.
    stop = threading.Event()
    grandchild = tmp_path / "grandchild.pid"
    modeller = _modeller(tmp_path, """
        import subprocess, time
        child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
        Path(os.environ["GRANDCHILD"]).write_text(str(child.pid))
        time.sleep(120)
    """, stop=stop)
    os.environ["GRANDCHILD"] = str(grandchild)

    def stop_once_started():
        import time
        while not grandchild.exists() or not grandchild.read_text():
            time.sleep(0.01)
        stop.set()

    try:
        threading.Thread(target=stop_once_started, daemon=True).start()
        with pytest.raises(RuntimeError, match="stopped while a model was being made"):
            modeller(INS, A_RECORD)
    finally:
        del os.environ["GRANDCHILD"]

    import time
    pid = int(grandchild.read_text())
    for _ in range(100):
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            break
        time.sleep(0.05)
    else:
        os.kill(pid, 9)
        raise AssertionError("what the modeller started outlived it")


def test_a_modeller_that_can_start_is_ready_and_one_that_cannot_says_why(tmp_path):
    assert _modeller(tmp_path, "pass").unready() is None
    pinned = ("The app's checkout pins flutter_scene 0.24.0; every stored scene is "
              "compiled by 0.23.0.")
    assert _modeller(tmp_path, "pass", check=1, check_says=pinned + "\n").unready() == \
        "The modeller exited with status 1: " + pinned
    settings = make_settings(tmp_path, tmp_path / "python",
                             structure_python=tmp_path / "no-such-venv" / "python")
    assert local_worker.Modeller(settings, threading.Event()).unready().startswith(
        "The modeller could not be started: ")


def _no_entries():
    """RCSB and UniProt as they answer for a protein the PDB holds nothing of."""
    from types import SimpleNamespace

    from pipeline.structure.test_fold_job import NOW, services

    post, get, asked = services(ids=[])
    entry = SimpleNamespace(body={"sequence": {"value": "M" * 110}, "features": []},
                            release="2026_03")
    return {"post": post, "get": get, "now": lambda: NOW,
            "fetch_entry": lambda accession: entry, "pfam": lambda accession: ()}, asked


def _driver_job(prior=None) -> dict:
    """A job as the Modeller writes one (`worker.model_job`)."""
    return worker.model_job(INS, A_RECORD, worker.FoldAsk(_A_ROW, prior or {"refresh": False}))


def _alphafold(monkeypatch, *, refuse=None, broken=None):
    alphafold = pytest.importorskip("pipeline.structure.alphafold")   # needs numpy
    asked = []

    def build(job, workspace, app, dart, **said):
        asked.append((job, workspace, app, dart))
        if refuse:
            raise alphafold.Refused(refuse)
        if broken:
            raise RuntimeError(broken)
        return alphafold.Built(b"glTF", b"scene", {"pdb": "AF-P01308-F1"},
                               [{"node": "plddtLow", "tint": "plddtLow"}], {"sampling": 8})

    monkeypatch.setattr(alphafold, "build", build)
    return alphafold, asked


def test_the_model_driver_draws_alphafolds_model_where_the_pdb_has_none(tmp_path, monkeypatch):
    from pipeline.resolver import model_local

    alphafold, asked = _alphafold(monkeypatch)
    services, searched = _no_entries()
    out = tmp_path / "out"
    assert model_local.make(_driver_job(), A_RECORD, out, Path("/src/app"), "/opt/dart",
                            "/opt/flutter", **services) == 0
    job, workspace, app, dart = asked[0]
    assert job == alphafold.Job(slug="ins", accession="P01308", display="Insulin",
                                protein="M" * 110, kept=((25, 54),), disulfides=((31, 96),),
                                allowed=3)
    assert (workspace, app, dart) == (out / "alphafold", Path("/src/app"), "/opt/dart")
    assert (out / "model.glb").read_bytes() == b"glTF"
    assert (out / "model.fsceneb").read_bytes() == b"scene"
    assert json.loads((out / "described.json").read_text()) == {
        "chrome": {"pdb": "AF-P01308-F1"},
        "chains": [{"node": "plddtLow", "tint": "plddtLow"}], "provenance": {"sampling": 8}}
    said = json.loads((out / "outcome.json").read_text())
    # The PDB was searched, whole, first.
    assert said["made"] == "alphafold" and searched.search == 1
    assert (said["selection"]["outcome"], said["selection"]["search"]["entities"],
            said["selection"]["policy"]) == ("none", 0, 1)
    assert not (out / "folding.json").exists()


def test_a_search_that_found_nothing_lately_is_not_made_again(tmp_path, monkeypatch):
    from pipeline.resolver import model_local
    from pipeline.structure import experimental, fold_report

    _alphafold(monkeypatch)
    services, searched = _no_entries()
    inputs, _ = experimental.inputs_for(INS, A_RECORD, fetch_entry=services["fetch_entry"],
                                        pfam=services["pfam"])
    none = {"id": 7, "policy": 1, "inputs": fold_report.inputs_digest(inputs), "entities": 0}
    out = tmp_path / "cached"
    assert model_local.make(_driver_job({"none": none, "refresh": False}), A_RECORD, out,
                            Path("/src/app"), "/opt/dart", "/opt/flutter", **services) == 0
    assert searched.search == 0
    assert json.loads((out / "outcome.json").read_text()) == {
        "made": "alphafold", "selection": {"cached": 7}}
    # Asked to start afresh, or for other inputs, or under another rule, it searches.
    for prior in ({"none": none, "refresh": True},
                  {"none": {**none, "inputs": "other"}, "refresh": False},
                  {"none": {**none, "policy": 2}, "refresh": False}):
        services, searched = _no_entries()
        assert model_local.make(_driver_job(prior), A_RECORD, tmp_path / "again",
                                Path("/src/app"), "/opt/dart", "/opt/flutter", **services) == 0
        assert searched.search == 1


def test_the_model_driver_says_why_neither_source_has_a_model_and_exits_3(tmp_path, monkeypatch):
    from pipeline.resolver import model_local

    _alphafold(monkeypatch, refuse="AlphaFold DB holds no model of UniProt P01308.")
    services, _ = _no_entries()
    out = tmp_path / "out"
    assert model_local.make(_driver_job(), A_RECORD, out, Path("/src/app"), "/opt/dart",
                            "/opt/flutter", **services) == model_local.REFUSED == 3
    assert (out / "refusal.txt").read_text() == (
        "The PDB holds no experimental structure of it. AlphaFold DB holds no model of "
        "UniProt P01308.")
    assert json.loads((out / "outcome.json").read_text())["selection"]["outcome"] == "none"
    assert not (out / "model.glb").exists()
    assert model_local.no_model(4, "AlphaFold's model is not confident here.") == (
        "None of its experimental structures in the PDB passes the checks for drawing. "
        "AlphaFold's model is not confident here.")
    # Anything else is the bake breaking, and is not a verdict on the protein.
    _alphafold(monkeypatch, broken="the scene importer failed (255)")
    services, _ = _no_entries()
    with pytest.raises(RuntimeError, match="scene importer"):
        model_local.make(_driver_job(), A_RECORD, tmp_path / "broke", Path("/src/app"),
                         "/opt/dart", "/opt/flutter", **services)


def _relaxin_job(prior=None) -> tuple:
    """Relaxin, whose entries the fixtures hold: its job, record and services."""
    from types import SimpleNamespace

    from pipeline.structure import test_experimental as fixtures
    from pipeline.structure.test_fold_job import NOW

    job = {**worker.model_job(fixtures.RELAXIN, fixtures.RECORD), "row": fixtures.ROW,
           "prior": prior or {"refresh": False}}
    entry = SimpleNamespace(body=fixtures.UNIPROT, release="2026_03")
    return job, fixtures.RECORD, {"now": lambda: NOW, "fetch_entry": lambda accession: entry,
                                  "pfam": lambda accession: ()}


def _baked(monkeypatch):
    """The bake and the app's check, stood in for: each candidate baked is
    said, and comes back as a fold."""
    from pipeline.structure import experimental

    drawn = []

    def draw(self, candidate, examined, entry, entity):
        drawn.append(candidate.entry)
        return experimental.Built(
            glb=b"glTF", scene=b"scene", folding=b'{"pdb": "' + examined.entry.encode() + b'"}',
            chrome={"pdb": examined.entry}, chains=[{"node": "chainA", "tint": "mature1"}],
            provenance={"pdb": examined.entry, "source": "PDB"}, folding_provenance={"pdb": 1},
            working_copy="", audit={"entry": examined.entry}, measured={"scene_bytes": 5})

    monkeypatch.setattr(experimental.Drawer, "__call__", draw)
    return drawn


def test_the_model_driver_writes_an_experimental_fold_and_exits_0(tmp_path, monkeypatch):
    from pipeline.resolver import model_local
    from pipeline.structure.test_fold_job import MV1, RLX, services

    drawn = _baked(monkeypatch)
    job, record, said = _relaxin_job()
    post, get, _ = services(files={"6RLX": RLX, "2MV1": MV1})
    out = tmp_path / "out"
    assert model_local.make(job, record, out, Path("/src/app"), "/opt/dart", "/opt/flutter",
                            post=post, get=get, **said) == 0
    assert drawn == ["6RLX"]
    assert (out / "model.glb").read_bytes() == b"glTF"
    assert (out / "model.fsceneb").read_bytes() == b"scene"
    assert json.loads((out / "folding.json").read_text()) == {"pdb": "6RLX"}
    described = json.loads((out / "described.json").read_text())
    assert described == {"chrome": {"pdb": "6RLX"},
                         "chains": [{"node": "chainA", "tint": "mature1"}],
                         "provenance": {"pdb": "6RLX", "source": "PDB"},
                         "folding_provenance": {"pdb": 1}}
    outcome = json.loads((out / "outcome.json").read_text())
    assert outcome["made"] == "pdb" and outcome["measured"] == {"scene_bytes": 5}
    assert (outcome["selection"]["outcome"], outcome["selection"]["candidate"]) == (
        "chosen", "6RLX:B@26,A@162")
    audit = json.loads((out / "audit.json").read_text())
    assert audit["entry"] == "6RLX" and audit["selection"] == outcome["selection"]


def test_a_search_cut_short_exits_4_and_the_next_goes_on_from_there(tmp_path, monkeypatch):
    from pipeline.resolver import model_local
    from pipeline.structure.test_fold_job import MV1, RLX_ELSEWHERE, services

    drawn = _baked(monkeypatch)
    job, record, said = _relaxin_job()
    post, get, asked = services(files={"6RLX": RLX_ELSEWHERE, "2MV1": MV1})
    out = tmp_path / "first"
    assert model_local.make(job, record, out, Path("/src/app"), "/opt/dart", "/opt/flutter",
                            post=post, get=get, limit=1, **said) == model_local.UNFINISHED == 4
    selection = json.loads((out / "outcome.json").read_text())["selection"]
    assert selection["outcome"] == "incomplete" and asked.files == ["6RLX"]
    assert selection["refused"] and all(c.startswith("6RLX:") for c in selection["refused"])
    assert not (out / "model.glb").exists() and drawn == []

    # The next attempt skips what was refused, where the search finds the same entries.
    progress = {"policy": 1, "inputs": selection["inputs"],
                "search": selection["search"]["sha256"], "refused": selection["refused"]}
    job, record, said = _relaxin_job({"progress": progress, "refresh": False})
    post, get, asked = services(files={"6RLX": RLX_ELSEWHERE, "2MV1": MV1})
    out = tmp_path / "second"
    assert model_local.make(job, record, out, Path("/src/app"), "/opt/dart", "/opt/flutter",
                            post=post, get=get, limit=1, **said) == 0
    assert asked.files == ["2MV1"] and drawn == ["2MV1"]
    outcome = json.loads((out / "outcome.json").read_text())
    assert outcome["selection"]["resumed"] == len(progress["refused"])
    # Where the PDB has moved on, it starts over.
    job, record, said = _relaxin_job({"progress": {**progress, "search": "another"},
                                      "refresh": False})
    post, get, asked = services(files={"6RLX": RLX_ELSEWHERE, "2MV1": MV1})
    assert model_local.make(job, record, tmp_path / "third", Path("/src/app"), "/opt/dart",
                            "/opt/flutter", post=post, get=get, limit=1, **said) == 4
    assert asked.files == ["6RLX"]


def test_the_model_driver_reads_its_job_and_the_record_beside_it(tmp_path, monkeypatch):
    from pipeline.resolver import model_local

    seen = []
    monkeypatch.setattr(model_local, "make", lambda job, record, out, app, dart, flutter:
                        seen.append((job, record, out, app, dart, flutter)) or 0)
    job = _driver_job()
    (tmp_path / job["record"]).write_bytes(A_RECORD)
    (tmp_path / "job.json").write_text(json.dumps(job))
    assert model_local.main([str(tmp_path / "job.json"), str(tmp_path / "out"), "--app",
                             "/src/app", "--dart", "/opt/dart", "--flutter", "/opt/flutter"]) == 0
    assert seen == [(job, A_RECORD, tmp_path / "out", Path("/src/app"), "/opt/dart",
                     "/opt/flutter")]
    assert model_local.main(["--app", "/src/app", "--dart", "/opt/dart",
                             "--flutter", "/opt/flutter"]) == 2


def test_the_model_driver_imports_nothing_heavy_until_it_is_asked():
    # The worker imports it for its exit status, in an environment with no numpy.
    code = ("import sys; import pipeline.resolver.model_local as m; "
            "print(m.REFUSED, sorted(n for n in sys.modules if n.split('.')[0] in "
            "('numpy', 'scipy', 'trimesh')))")
    done = subprocess.run([sys.executable, "-c", code], cwd=BACKEND, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == "3 []"
    assert subprocess.run([sys.executable, "-c", code.replace("m.REFUSED", "m.UNFINISHED")],
                          cwd=BACKEND, capture_output=True, text=True).stdout.strip() == "4 []"


def test_a_model_can_be_made_here_or_the_check_says_why_not(tmp_path):
    pytest.importorskip("trimesh")
    pytest.importorskip("scipy")
    from pipeline.paths import CLIENT
    from pipeline.resolver import model_local
    from pipeline.structure import alphafold, bake

    dart = local_worker.default_dart(Path.home())
    flutter = local_worker.default_flutter(Path.home())
    if (not Path(bake.PYMOL).exists() or not dart.exists() or not flutter.exists()
            or alphafold.importer_unready(CLIENT)):
        pytest.skip("needs PyMOL, a Dart SDK, Flutter and the app's checkout")
    assert model_local.check(CLIENT, str(dart), str(flutter)) == 0
    # A checkout that is not the app's is said, before anything is run.
    assert model_local.check(tmp_path, str(dart), str(flutter)) == 1
    # And a Flutter that is not there to hold a fold to the app's rules.
    assert model_local.check(CLIENT, str(dart), str(tmp_path / "no-flutter")) == 1


# ------------------------------------------------------------ the evidence

ATLAS_KEY = "AIzaSy-not-a-real-atlas-key"


def fake_impact_python(folder: Path, body: str, check: int = 0, check_says: str = "") -> Path:
    """A stand-in for the AVI bake's interpreter: this one, running `body`.

    The evidencer starts it as it starts the real one: `-u -m <driver>` and
    then `--check`, or the target's file, the record's, the folder the track
    goes into and `--state <folder>`. `body` runs where the driver would, with
    those as `target_file`, `record_file`, `out` and `state`, and the
    repository importable. `--check` exits `check`, having said `check_says`.
    """
    script = folder / "fake-impact-python"
    script.write_text(
        f"#!{sys.executable}\n"
        "import json, os, sys\n"
        "from pathlib import Path\n"
        "sys.path.insert(0, os.getcwd())\n"
        "if sys.argv[-1] == '--check':\n"
        f"    sys.stderr.write({check_says!r})\n"
        f"    sys.exit({check})\n"
        "target_file, record_file, out = (Path(p) for p in sys.argv[-5:-2])\n"
        "state = Path(sys.argv[-1])\n"
        + textwrap.dedent(body), encoding="utf-8")
    script.chmod(0o755)
    return script


# What `impact_local.py` prints and leaves, as far as the worker reads them.
PLACES_A_TRACK = """
    print("    chr11 - ENST00000381330.5 (3 exons)")
    print("    1 window(s), 1,559 bp")
    print("    sequence gate: 1,431 of 1,431 bases agree")
    out.mkdir(parents=True, exist_ok=True)
    (out / "impact.json").write_text('{"gene": "INS"}')
"""


def _evidencer(tmp_path, body, stop=None, **check):
    folder = tmp_path / "impact"
    folder.mkdir(exist_ok=True)
    settings = make_settings(
        tmp_path, tmp_path / "no-esm-python", evidence=True,
        impact_python=fake_impact_python(folder, body, **check),
        uv=Path("/opt/uv/bin/uv"), alphagenome_key=ATLAS_KEY)
    return local_worker.Evidencer(settings, stop or threading.Event())


def test_this_worker_fetches_evidence_unless_its_env_says_not(tmp_path):
    (tmp_path / ".env").write_text(f"ALPHAGENOME_API_KEY={ATLAS_KEY}\n")
    settings = local_worker.load_settings(_FOUR, tmp_path / ".env", tmp_path)
    assert settings.evidence is True and settings.alphagenome_key == ATLAS_KEY
    assert settings.impact_python == BACKEND / "pipeline" / "impact" / "venv" / "bin" / "python"
    assert (settings.fetching, settings.impact_state) == \
        (settings.state / "fetching", settings.state / "impact")
    assert BACKEND not in settings.fetching.parents and BACKEND not in settings.impact_state.parents

    for word in ("0", "false", "No", "off"):
        said = local_worker.load_settings(
            {**_FOUR, "RESOLVER_EVIDENCE": word}, tmp_path / "no.env", tmp_path)
        assert said.evidence is False and said.alphagenome_key is None
    told = local_worker.load_settings(
        {**_FOUR, "RESOLVER_IMPACT_PYTHON": "/opt/avi/python", "RESOLVER_UV": "/opt/uv"},
        tmp_path / "no.env", tmp_path)
    assert (told.impact_python, told.uv) == (Path("/opt/avi/python"), Path("/opt/uv"))
    # A worker built by hand, as every test before these builds it, fetches none.
    assert make_settings(tmp_path, tmp_path / "python").evidence is False


def test_uv_is_where_its_installer_puts_it_or_on_the_path(tmp_path, monkeypatch):
    installed = tmp_path / ".local" / "bin" / "uv"
    installed.parent.mkdir(parents=True)
    installed.write_text("")
    assert local_worker.default_uv(tmp_path) == installed
    monkeypatch.setattr(local_worker.shutil, "which", lambda name: "/opt/homebrew/bin/uv")
    assert local_worker.default_uv(tmp_path / "elsewhere") == Path("/opt/homebrew/bin/uv")
    monkeypatch.setattr(local_worker.shutil, "which", lambda name: None)
    assert local_worker.default_uv(tmp_path / "elsewhere") == Path("uv")


def test_an_avi_bakes_time_is_inside_the_reapers_patience_for_one():
    minutes, unit = store.STALE_EVIDENCE.split()
    assert unit == "minutes" and local_worker.EVIDENCE_TIMEOUT < int(minutes) * 60


def test_the_atlas_key_is_a_secret_the_log_never_holds(settings):
    told = dataclasses.replace(settings, alphagenome_key=ATLAS_KEY)
    assert ATLAS_KEY in told.secrets
    handlers = local_worker.start_logging(told)
    try:
        local_worker.log.error("The AVI bake exited with status 1: key %s rejected", ATLAS_KEY)
    finally:
        for handler in handlers:
            local_worker.log.removeHandler(handler)
            handler.close()
    text = told.log_file.read_text()
    assert ATLAS_KEY not in text and "key [redacted] rejected" in text


def test_the_evidencer_hands_over_the_protein_and_returns_the_track(
        tmp_path, monkeypatch, caplog):
    monkeypatch.setenv("DATABASE_URL", _FOUR["DATABASE_URL"])
    monkeypatch.setenv("SUPABASE_SERVICE_KEY", _FOUR["SUPABASE_SERVICE_KEY"])
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "a-key-from-elsewhere")
    monkeypatch.setenv("SEEN", str(tmp_path / "seen.json"))
    evidencer = _evidencer(tmp_path, PLACES_A_TRACK + """
    import pickle
    Path(os.environ["SEEN"]).write_text(json.dumps({
        "started_with": sys.argv[1:4], "cwd": os.getcwd(), "state": str(state),
        "target": pickle.loads(target_file.read_bytes()).slug,
        "record": record_file.read_text(),
        "path": os.environ["PATH"].split(os.pathsep)[0],
        "environment": {key: os.environ.get(key) for key in (
            "DATABASE_URL", "SUPABASE_SERVICE_KEY", "ALPHAGENOME_API_KEY", "PYTHONUNBUFFERED")}}))
""")
    with caplog.at_level(logging.INFO, logger="resolver-worker"):
        assert evidencer(INS, RECORD) == b'{"gene": "INS"}'

    seen = json.loads((tmp_path / "seen.json").read_text())
    assert seen["started_with"] == ["-u", "-m", "pipeline.resolver.impact_local"]
    assert Path(seen["cwd"]) == BACKEND
    assert (seen["target"], seen["record"]) == ("ins", RECORD.decode())
    # What a bake keeps between tries is the worker's, and is kept.
    assert seen["state"] == str(evidencer.settings.impact_state)
    # uv's folder first, where launchd's PATH names none of it.
    assert seen["path"] == "/opt/uv/bin"
    # The Atlas's key from the worker's settings, and neither database credential.
    assert seen["environment"] == {"DATABASE_URL": None, "SUPABASE_SERVICE_KEY": None,
                                   "ALPHAGENOME_API_KEY": ATLAS_KEY, "PYTHONUNBUFFERED": "1"}
    said = [record.getMessage() for record in caplog.records]
    assert len(said) == 1 and re.fullmatch(
        r"ins: AVI, 1 Atlas window\(s\), sequence gate 1,431 of 1,431, 15 B, \d+\.\d s", said[0])
    assert list(evidencer.settings.fetching.iterdir()) == []


def test_a_gate_is_told_as_a_verdict_in_its_own_words(tmp_path):
    evidencer = _evidencer(tmp_path, """
        out.mkdir(parents=True, exist_ok=True)
        (out / "refusal.txt").write_text("INS: record has 3 exons, MANE Select has 4")
        sys.exit(3)
    """)
    with pytest.raises(worker.Unplaced) as unplaced:
        evidencer(INS, RECORD)
    assert str(unplaced.value) == "INS: record has 3 exons, MANE Select has 4"


@pytest.mark.parametrize("body, said", [
    ("sys.exit(0)", "The AVI bake exited with status 0 and wrote no track."),
    ("sys.exit(3)", "The AVI bake exited with status 3 and gave no reason."),
    ("sys.stderr.write('Unreachable: chr11:1-2 failed after 7 attempts:\\nUNAVAILABLE'); "
     "sys.exit(1)",
     "The AVI bake exited with status 1: Unreachable: chr11:1-2 failed after 7 attempts: "
     "UNAVAILABLE"),
    ("os.kill(os.getpid(), 9)", "The AVI bake was ended by signal 9 and said nothing."),
])
def test_an_avi_bake_that_ends_any_other_way_broke(tmp_path, body, said, spawned):
    evidencer = _evidencer(tmp_path, body)
    with pytest.raises(RuntimeError) as broke:
        evidencer(INS, RECORD)
    assert not isinstance(broke.value, worker.Unplaced) and str(broke.value) == said


def test_a_stop_ends_the_avi_bake_and_what_it_started(tmp_path, spawned):
    # GENCODE's lookup runs under uv, the bake's own child.
    stop = threading.Event()
    grandchild = tmp_path / "grandchild.pid"
    evidencer = _evidencer(tmp_path, """
        import subprocess, time
        child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
        Path(os.environ["GRANDCHILD"]).write_text(str(child.pid))
        time.sleep(120)
    """, stop=stop)
    os.environ["GRANDCHILD"] = str(grandchild)

    def stop_once_started():
        import time
        while not grandchild.exists() or not grandchild.read_text():
            time.sleep(0.01)
        stop.set()

    try:
        threading.Thread(target=stop_once_started, daemon=True).start()
        with pytest.raises(RuntimeError) as cut:
            evidencer(INS, RECORD)
    finally:
        del os.environ["GRANDCHILD"]
    assert str(cut.value) == local_worker.STOPPED_PLACING

    import time
    pid = int(grandchild.read_text())
    for _ in range(100):
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            break
        time.sleep(0.05)
    else:
        os.kill(pid, 9)
        raise AssertionError("what the AVI bake started outlived it")
    with pytest.raises(RuntimeError, match="AlphaGenome's scores were being fetched"):
        evidencer(INS, RECORD)                       # and none is started once stopping


def test_an_avi_bake_running_past_its_time_is_ended(tmp_path, monkeypatch, spawned):
    monkeypatch.setattr(local_worker, "EVIDENCE_TIMEOUT", 0.2)
    evidencer = _evidencer(tmp_path, "import time; time.sleep(120)")
    with pytest.raises(RuntimeError, match="still being fetched after 0 minutes"):
        evidencer(INS, RECORD)
    assert spawned[-1].poll() is not None


def test_an_evidencer_that_can_ask_the_atlas_is_ready_and_one_that_cannot_says_why(tmp_path):
    assert _evidencer(tmp_path, "pass").unready() is None
    said = "The Atlas did not answer: API key not valid. Please pass a valid API key."
    assert _evidencer(tmp_path, "pass", check=1, check_says=said + "\n").unready() == \
        "The AVI bake exited with status 1: " + said
    settings = make_settings(tmp_path, tmp_path / "python", evidence=True,
                             impact_python=tmp_path / "no-such-venv" / "python")
    assert local_worker.Evidencer(settings, threading.Event()).unready().startswith(
        "The AVI bake could not be started: ")


def test_clinvar_is_fetched_here_into_the_workers_own_folder(settings, monkeypatch, caplog):
    asked = []
    track = json.dumps({"variants": [{}, {}], "searched_records": 503,
                        "source_batches": [{}] * 6}).encode()

    def fetched(target, cache, stopped):
        asked.append((target, cache, stopped()))
        return track

    monkeypatch.setattr(worker, "clinvar_in_process", fetched)
    stop = threading.Event()
    baker = local_worker.ClinVarBaker(settings, stop)
    with caplog.at_level(logging.INFO, logger="resolver-worker"):
        assert baker(INS) == track
    assert asked == [(INS, settings.fetching / "clinvar", False)]
    assert re.fullmatch(r"ins: ClinVar, 2 of 503 record\(s\) mapped, 6 batch\(es\), \d+ B, "
                        r"\d+\.\d s", caplog.records[-1].getMessage())
    stop.set()
    with pytest.raises(worker.Unfetched, match="stopped while ClinVar"):
        baker(INS)


class Fake:
    """An evidencer or a ClinVar baker the cycle only hands on, and may ask
    whether it is ready."""

    def __init__(self, unready=None):
        self.why, self.asked = unready, 0

    def unready(self):
        self.asked += 1
        return self.why


def test_the_evidence_is_fetched_after_the_resolving_and_before_the_scoring(
        settings, monkeypatch, caffeinate, caplog):
    queue = Queue(monkeypatch, requests=[_request()], bakes=[_bake()],
                  impacts=[_bake()], clinvars=[_bake()])
    evidencer, clinvar = Fake(), Fake()
    with caplog.at_level(logging.INFO, logger="resolver-worker"):
        done = local_worker.cycle(settings, Ready(), threading.Event(), evidencer=evidencer,
                                  clinvar=clinvar)

    assert [kind for kind, _ in queue.order] == \
        ["sweep", "sweep", "impact", "impact", "clinvar", "clinvar", "constraint", "constraint"]
    # The protein resolved has its evidence queued; the bakers are the ones handed in.
    assert queue.order[0] == ("sweep", {"evidence": True})
    assert (queue.order[2][1], queue.order[4][1]) == (evidencer, clinvar)
    assert [o["state"] for o in done["impact"]] == ["ready"]
    assert [o["state"] for o in done["clinvar"]] == ["ready"]
    assert done["idle"] is False and evidencer.asked == 1
    said = [record.getMessage() for record in caplog.records]
    assert re.fullmatch(r"ins: impact ready in \d+\.\d s", said[1])
    assert re.fullmatch(r"ins: clinvar ready in \d+\.\d s", said[2])
    assert len(caffeinate) == 1 and caffeinate[0].poll() is not None

    # A worker that fetches none says nothing of evidence to the sweep.
    queue = Queue(monkeypatch, requests=[_request()])
    assert "impact" not in local_worker.cycle(settings, Ready(), threading.Event())
    assert queue.order[0] == ("sweep", {})


def test_no_avi_bake_is_claimed_while_the_atlas_cannot_be_asked(
        settings, monkeypatch, caffeinate, caplog):
    queue = Queue(monkeypatch, bakes=[_bake()], impacts=[_bake(), _bake("b2m")])
    evidencer = Fake("The AVI bake exited with status 1: ALPHAGENOME_API_KEY is not set in "
                     "the worker's .env.")
    done = local_worker.cycle(settings, Ready(), threading.Event(), evidencer=evidencer,
                              clinvar=Fake())
    assert done["impact"] == [] and len(queue.impacts) == 2
    assert ("2 AVI bake(s) left on the queue, and their ClinVar bakes behind them: no AVI "
            "track can be made here now. The AVI bake exited with status 1: "
            "ALPHAGENOME_API_KEY is not set in the worker's .env.") in caplog.text
    # ClinVar needs only NCBI, and the scoring goes on for whatever the claim
    # lets through: here all of it, where the real claim holds back a protein
    # whose evidence waits (`store.WAITS_FOR`, held on Postgres).
    assert [kind for kind, _ in queue.order] == ["sweep", "clinvar", "constraint", "constraint"]


def test_a_protein_whose_avi_goes_back_on_the_queue_waits_for_the_next_cycle(
        settings, monkeypatch, caffeinate):
    queue = Queue(monkeypatch, impacts=[
        _bake(state="queued", reason="RESOURCE_EXHAUSTED"), _bake("b2m")])
    local_worker.cycle(settings, Ready(), threading.Event(), evidencer=Fake(), clinvar=Fake())
    assert [kind for kind, _ in queue.order] == ["sweep", "impact", "clinvar"]
    assert len(queue.impacts) == 2


def test_the_avi_hand_over_is_emptied_and_what_a_bake_keeps_is_not(settings):
    for folder in (settings.fetching, settings.impact_state / "checkpoints"):
        folder.mkdir(parents=True)
        (folder / "left.json").write_text("{}")
    local_worker.clear_workspace(settings)
    assert list(settings.fetching.iterdir()) == []
    assert (settings.impact_state / "checkpoints" / "left.json").exists()


def test_no_child_but_the_avi_bake_is_given_the_atlas_key(settings, monkeypatch):
    # Started by hand from a shell that read ~/.env, the worker's own
    # environment may hold it: the scorer and the modeller are never handed it.
    monkeypatch.setenv("ALPHAGENOME_API_KEY", ATLAS_KEY)
    told = dataclasses.replace(settings, alphagenome_key=ATLAS_KEY)
    stop = threading.Event()
    for child in (local_worker.Scorer(told, stop), local_worker.Modeller(told, stop)):
        assert "ALPHAGENOME_API_KEY" not in child._environment()
    assert local_worker.Evidencer(told, stop)._environment()["ALPHAGENOME_API_KEY"] == ATLAS_KEY
