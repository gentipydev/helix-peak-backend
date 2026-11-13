"""The Modal app as deployed: its functions, and what it uploads.

Skipped where the Modal client is not installed (it is not a service
dependency; `pipeline/resolver/requirements.txt` installs it to deploy).
Importing the app defines its images and functions without contacting Modal.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

pytest.importorskip("modal")
os.environ.setdefault("NCBI_EMAIL", "tests@example.com")

from pipeline.resolver import modal_app  # noqa: E402

PIPELINE = BACKEND / "pipeline"


def _uploaded(package: Path) -> set[str]:
    """Every file `add_local_python_source` would send, as Modal walks it."""
    found = set()
    for root, _, files in os.walk(package):
        for name in files:
            relative = (Path(root) / name).relative_to(package)
            if any(part.startswith(".") for part in relative.parts):
                continue  # Modal leaves dot-prefixed paths out itself
            if not modal_app.not_code(relative):
                found.add(relative.as_posix())
    return found


def test_the_app_defines_its_three_functions():
    import modal

    assert modal_app.app.name == "helix-peak-resolver"
    for function in (modal_app.sweep, modal_app.score, modal_app.wake):
        assert isinstance(function, modal.Function)


def test_the_scorer_image_is_built_from_the_scorers_own_lock():
    assert modal_app.LOCK == PIPELINE / "constraint" / "requirements-lock.txt"
    assert modal_app.LOCK.exists()
    assert modal_app.MODEL == "facebook/esm2_t33_650M_UR50D"


def test_only_code_goes_up():
    sent = _uploaded(PIPELINE)
    assert {"resolver/resolve.py", "resolver/worker.py", "resolver/store.py",
            "resolver/scoring.py",
            "mock/build_gene_record.py", "constraint/score_protein.py",
            "upload_tracks.py", "seed_catalog.py", "targets.py"} <= sent
    assert all(path.endswith(".py") for path in sent)
    assert not [path for path in sent if path.startswith(("data/", "resolver/fixtures/"))]


def test_a_venv_or_baked_data_never_goes_up():
    for path in ("structure/venv/lib/python3.12/site-packages/torch/__init__.py",
                 "impact/venv/bin/activate_this.py",
                 "data/assets/mock/gene_ins.json", "data/notes.py",
                 "locus/fixtures/chromalias_hg38_five.json", "resolver/__pycache__/resolve.cpython-312.pyc"):
        assert modal_app.not_code(Path(path)), path
    assert not modal_app.not_code(Path("resolver/worker.py"))


# ------------------------------------------------- what starts the GPU


def _minutes(interval: str) -> int:
    number, unit = interval.split()
    assert unit == "minutes"
    return int(number)


def test_no_call_outlives_the_reapers_patience_with_its_claim():
    from pipeline.resolver import store

    assert modal_app.SCORE_TIMEOUT < _minutes(store.STALE_BAKE) * 60
    assert modal_app.SWEEP_TIMEOUT < _minutes(store.STALE_RESOLVE) * 60


class _Sweeps:
    """`worker.sweep` as a sequence of what each call resolved."""

    def __init__(self, *rounds):
        self.rounds = list(rounds)
        self.limits = []

    def __call__(self, conn, storage, *, limit):
        self.limits.append(limit)
        resolved = self.rounds.pop(0) if self.rounds else []
        return {"resolved": resolved, "constraint_queued": 0, "structure_queued": 0}


def test_requests_are_resolved_one_at_a_time_until_none_is_left():
    sweeps = _Sweeps([{"gene": "B2M", "state": "done"}], [{"gene": "SLN", "state": "refused"}])
    summary = modal_app.resolve_queued(None, None, sweep=sweeps)
    assert [r["gene"] for r in summary["resolved"]] == ["B2M", "SLN"]
    assert sweeps.limits == [1, 1, 1]


def test_a_request_put_back_on_the_queue_ends_the_sweep():
    sweeps = _Sweeps([{"gene": "B2M", "state": "queued"}], [{"gene": "B2M", "state": "done"}])
    summary = modal_app.resolve_queued(None, None, sweep=sweeps)
    assert summary["resolved"] == [{"gene": "B2M", "state": "queued"}]
    assert sweeps.limits == [1]


def test_a_sweep_takes_at_most_its_limit():
    sweeps = _Sweeps(*[[{"gene": f"G{i}", "state": "done"}] for i in range(20)])
    summary = modal_app.resolve_queued(None, None, sweep=sweeps)
    assert len(summary["resolved"]) == modal_app.SWEEP_LIMIT


def test_score_goes_on_only_after_a_bake_it_did_not_put_back():
    assert not modal_app.goes_on([])
    assert not modal_app.goes_on([{"slug": "b2m", "state": "queued"}])
    for state in ("ready", "refused", "stopped"):
        assert modal_app.goes_on([{"slug": "b2m", "state": state}])


# The functions as defined, kept before a test stands a fake in for `score`.
SWEEP, SCORE = modal_app.sweep, modal_app.score


class _Conn:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class _Spawns:
    def __init__(self):
        self.count = 0

    def spawn(self):
        self.count += 1


@pytest.fixture
def body(monkeypatch):
    """The functions' bodies run here, on fakes: no Modal, no database."""
    from pipeline.resolver import store, worker

    spawned = _Spawns()
    monkeypatch.setattr(modal_app, "_connect", lambda: (_Conn(), object()))
    monkeypatch.setattr(modal_app, "score", spawned)
    monkeypatch.setattr(modal_app.socket, "setdefaulttimeout", lambda seconds: None)
    monkeypatch.setattr(worker, "sweep", _Sweeps())
    state = {"claimable": 0, "scored": []}
    monkeypatch.setattr(store, "claimable_bakes", lambda conn, kind: state["claimable"])
    monkeypatch.setattr(worker, "score_all", lambda conn, storage, limit: state["scored"])
    state["spawned"] = spawned
    return state


@pytest.mark.parametrize("claimable, spawns", [(0, 0), (2, 1)])
def test_a_sweep_starts_the_gpu_only_for_a_bake_it_can_take(body, claimable, spawns):
    body["claimable"] = claimable
    summary = SWEEP.local()
    assert summary["constraint_claimable"] == claimable
    assert body["spawned"].count == spawns


@pytest.mark.parametrize("scored, claimable, spawns", [
    ([], 3, 0),                                        # nothing it could take
    ([{"slug": "b2m", "state": "queued"}], 3, 0),       # put back: the next sweep tries again
    ([{"slug": "b2m", "state": "ready"}], 0, 0),        # nothing more waits
    ([{"slug": "b2m", "state": "ready"}], 1, 1),        # the next, in a call of its own
])
def test_score_starts_itself_again_only_while_a_bake_waits(body, scored, claimable, spawns):
    body["scored"], body["claimable"] = scored, claimable
    assert SCORE.local() == scored
    assert body["spawned"].count == spawns
