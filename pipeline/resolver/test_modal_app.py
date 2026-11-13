"""The Modal app as deployed: its functions, what starts each, and what it uploads.

Skipped where the Modal client is not installed (it is not a service
dependency; `pipeline/resolver/requirements.txt` installs it to deploy).
Importing the app defines its images and functions without contacting Modal.
"""

from __future__ import annotations

import os
import re
import subprocess
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


def test_the_app_defines_its_five_functions():
    import modal

    assert modal_app.app.name == "helix-peak-resolver"
    for function in (modal_app.sweep, modal_app.score, modal_app.wake, modal_app.fold,
                     modal_app.evidence):
        assert isinstance(function, modal.Function)


def test_the_secret_holds_the_macs_four_keys_and_the_atlases():
    from pipeline.resolver import local_worker

    assert modal_app.KEYS == (*local_worker.KEYS, "ALPHAGENOME_API_KEY")


def test_the_scorer_image_is_built_from_the_scorers_own_lock():
    assert modal_app.LOCK == PIPELINE / "constraint" / "requirements-lock.txt"
    assert modal_app.LOCK.exists()
    assert modal_app.MODEL == "facebook/esm2_t33_650M_UR50D"


def _pins(path: Path) -> dict:
    found = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        name, _, version = line.partition("==")
        assert version and re.fullmatch(r"[A-Za-z0-9_.-]+", name), f"not a pin: {line}"
        found[name.lower().replace("_", "-")] = version
    return found


def test_the_bakes_environments_are_pinned_as_the_dev_macs():
    impact = _pins(modal_app.IMPACT_LOCK)
    assert impact["alphagenome"] == "0.9.0"
    assert _pins(modal_app.STRUCTURE_REQUIREMENTS) == {
        "numpy": "2.5.3", "scipy": "1.18.1", "trimesh": "5.1.0"}
    # The worker's own numpy is the one both bakes run.
    assert modal_app._WORKER_NUMPY == f"numpy=={impact['numpy']}"


def test_the_images_tools_are_the_ones_the_bakes_hold_a_fold_and_a_lookup_to():
    pytest.importorskip("numpy")
    from pipeline.impact import bake_impact
    from pipeline.structure import experimental

    assert modal_app.PYMOL_VERSION == experimental.PYMOL
    assert modal_app.SKILL.name == bake_impact.DEFAULT_SKILL.name
    assert modal_app.DART.endswith("/bin/cache/dart-sdk/bin/dart")


def test_only_code_goes_up():
    sent = _uploaded(PIPELINE)
    assert {"resolver/resolve.py", "resolver/worker.py", "resolver/store.py",
            "resolver/scoring.py", "resolver/local_worker.py", "resolver/model_local.py",
            "resolver/impact_local.py",
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


_NOT_A_FOLDS = ("build", ".dart_tool", "android", "ios", "macos", "linux", "web", "windows")


def test_only_what_a_fold_reads_of_the_app_goes_up_and_never_its_env():
    parts = modal_app.APP_FILES + modal_app.APP_DIRS
    assert ".env" not in parts
    assert not [part for part in parts if part.split("/")[0] in _NOT_A_FOLDS]
    assert {"pubspec.lock", "test/fold_bundle_check.dart", "hook"} <= set(parts)
    if not modal_app.CLIENT.is_dir():
        pytest.skip("no app checkout beside this repository")
    for part in parts:
        assert (modal_app.CLIENT / part).exists(), part
    # What the check imports from beside it goes up with it.
    check = (modal_app.CLIENT / "test/fold_bundle_check.dart").read_text(encoding="utf-8")
    for imported in re.findall(r"^import '([^':]+)';", check, re.M):
        assert any(f"test/{imported}".startswith(f"{part}/") for part in modal_app.APP_DIRS), \
            imported


def test_the_app_imports_with_neither_the_app_nor_the_skill_beside_it(tmp_path):
    environment = {**os.environ, "HELIXPEEK_CLIENT": str(tmp_path / "no-app"),
                   "ALPHAGENOME_SKILLS": str(tmp_path / "no-skills")}
    done = subprocess.run(
        [sys.executable, "-c", "from pipeline.resolver import modal_app; "
                               "print(modal_app.CLIENT); print(modal_app.SKILL)"],
        cwd=str(BACKEND), env=environment, capture_output=True, text=True, timeout=300)
    assert done.returncode == 0, done.stderr
    assert done.stdout.split() == [str(tmp_path / "no-app"),
                                   str(tmp_path / "no-skills" / modal_app.SKILL.name)]


def test_a_containers_settings_are_the_macs_with_the_images_tools(monkeypatch):
    for key in modal_app.KEYS:
        monkeypatch.setenv(key, f"{key.lower()}-said")
    settings = modal_app._settings()
    assert (settings.database_url, settings.supabase_url, settings.service_key,
            settings.ncbi_email, settings.alphagenome_key) == (
        "database_url-said", "supabase_url-said", "supabase_service_key-said",
        "ncbi_email-said", "alphagenome_api_key-said")
    assert settings.structures and settings.evidence
    # The worker's directories are the ones the bakers are pointed at.
    assert settings.data == Path(modal_app._ENV["HELIXPEEK_DATA"])
    assert settings.genbank == Path(modal_app._ENV["HELIXPEEK_GB_CACHE"])
    assert (settings.structure_python, settings.dart, settings.flutter, settings.app,
            settings.impact_python, settings.uv) == tuple(Path(path) for path in (
                modal_app.STRUCTURE_PYTHON, modal_app.DART, modal_app.FLUTTER,
                modal_app.APP_DIR, modal_app.IMPACT_PYTHON, modal_app.UV))


@pytest.mark.parametrize("said, on", [
    (None, True), ("", True), ("1", True), ("yes", True),
    ("0", False), ("false", False), (" No ", False), ("OFF", False),
])
def test_a_switch_reads_as_the_macs_env_says_it(said, on):
    environ = {} if said is None else {"RESOLVER_STRUCTURES": said}
    assert modal_app.switched_on("RESOLVER_STRUCTURES", environ) is on


def test_the_sweeps_containers_carry_the_switches_as_deployed():
    assert modal_app._SWITCHES == {
        "RESOLVER_STRUCTURES": "1" if modal_app.STRUCTURES else "0",
        "RESOLVER_EVIDENCE": "1" if modal_app.EVIDENCE else "0"}


# ------------------------------------------------- what starts what


def _minutes(interval: str) -> int:
    number, unit = interval.split()
    assert unit == "minutes"
    return int(number)


def test_no_call_outlives_the_reapers_patience_with_its_claim():
    from pipeline.resolver import store

    assert modal_app.SCORE_TIMEOUT < _minutes(store.STALE_BAKE) * 60
    assert modal_app.SWEEP_TIMEOUT < _minutes(store.STALE_RESOLVE) * 60
    assert modal_app.FOLD_TIMEOUT < _minutes(store.STALE_STRUCTURE) * 60
    assert modal_app.EVIDENCE_TIMEOUT < _minutes(store.STALE_EVIDENCE) * 60


def test_a_drivers_own_limits_end_its_bake_before_modal_ends_the_call():
    from pipeline.resolver import local_worker

    assert local_worker._CHECK_TIMEOUT + local_worker.MODEL_TIMEOUT < modal_app.FOLD_TIMEOUT
    assert (local_worker._CHECK_TIMEOUT + local_worker.EVIDENCE_TIMEOUT
            < modal_app.EVIDENCE_TIMEOUT)


class _Sweeps:
    """`worker.sweep` as a sequence of what each call resolved."""

    def __init__(self, *rounds):
        self.rounds = list(rounds)
        self.limits = []
        self.asked = []

    def __call__(self, conn, storage, *, limit, **asked):
        self.limits.append(limit)
        self.asked.append(asked)
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


@pytest.mark.parametrize("structures, evidence, asked", [
    (False, False, {}),
    (True, False, {"structures": True}),
    (False, True, {"evidence": True}),
    (True, True, {"structures": True, "evidence": True}),
])
def test_a_protein_resolved_has_its_fold_and_evidence_queued_only_where_they_are_made(
        structures, evidence, asked):
    sweeps = _Sweeps([{"gene": "B2M", "state": "done"}])
    modal_app.resolve_queued(None, None, sweep=sweeps, structures=structures, evidence=evidence)
    assert sweeps.asked == [asked, asked]


def test_a_baker_goes_on_only_after_work_it_did_not_put_back():
    assert not modal_app.goes_on([])
    assert not modal_app.goes_on([{"slug": "b2m", "state": "queued"}])
    for state in ("ready", "refused", "stopped", "kept", "superseded"):
        assert modal_app.goes_on([{"slug": "b2m", "state": state}])
    # An AVI bake done and its ClinVar bake put back: the next sweep tries again.
    assert not modal_app.goes_on([{"slug": "ttr", "state": "ready"},
                                  {"slug": "ttr", "state": "queued"}])
    assert modal_app.goes_on([{"slug": "ttr", "state": "ready"},
                              {"slug": "ttr", "state": "refused"}])


@pytest.mark.parametrize("claimable, structures, evidence, started", [
    ({}, True, True, []),
    ({"constraint": 1}, False, False, ["score"]),
    ({"structure": 1, "impact": 1}, False, False, []),
    ({"structure": 1, "impact": 1}, True, True, ["fold", "evidence"]),
    ({"structure": 1, "impact": 1}, True, False, ["fold"]),
    ({"clinvar": 2}, True, True, ["evidence"]),
    ({"constraint": 1, "structure": 2, "impact": 3, "clinvar": 4}, True, True,
     ["score", "fold", "evidence"]),
])
def test_a_sweep_starts_each_baker_only_for_a_bake_it_can_take(
        claimable, structures, evidence, started):
    assert modal_app.starts(claimable, structures, evidence) == started


# The functions as defined, kept before a test stands fakes in for them.
SWEEP, SCORE, FOLD, EVIDENCE = modal_app.sweep, modal_app.score, modal_app.fold, modal_app.evidence


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


class _Driver:
    """A driver of `local_worker` whose check says what the test says."""

    def __init__(self, said):
        self.said = said

    def unready(self):
        return self.said


@pytest.fixture
def body(monkeypatch):
    """The functions' bodies run here, on fakes: no Modal, no database, no bake."""
    from pipeline.resolver import local_worker, store, worker

    state = {"claimable": {}, "folds_ready": True, "unready": {"model": None, "avi": None},
             "scored": [], "modelled": [], "placed": [], "fetched": [], "calls": []}
    spawned = {name: _Spawns() for name in ("score", "fold", "evidence")}
    for name, fake in spawned.items():
        monkeypatch.setattr(modal_app, name, fake)
    monkeypatch.setattr(modal_app, "_connect", lambda: (_Conn(), object()))
    monkeypatch.setattr(modal_app, "_settings", lambda: "settings")
    monkeypatch.setattr(modal_app, "_work_here", lambda settings: None)
    monkeypatch.setattr(modal_app.socket, "setdefaulttimeout", lambda seconds: None)
    monkeypatch.setattr(modal_app, "STRUCTURES", True)
    monkeypatch.setattr(modal_app, "EVIDENCE", True)
    monkeypatch.setattr(worker, "sweep", _Sweeps())
    monkeypatch.setattr(store, "claimable_bakes",
                        lambda conn, kind: state["claimable"].get(kind, 0))
    monkeypatch.setattr(store, "folds_ready", lambda conn: state["folds_ready"])
    monkeypatch.setattr(local_worker, "Modeller",
                        lambda settings, stop: _Driver(state["unready"]["model"]))
    monkeypatch.setattr(local_worker, "Evidencer",
                        lambda settings, stop: _Driver(state["unready"]["avi"]))
    monkeypatch.setattr(local_worker, "ClinVarBaker", lambda settings, stop: "ClinVar's baker")
    monkeypatch.setattr(local_worker, "clear_workspace",
                        lambda settings: state["calls"].append("cleared"))

    def called(what, result):
        def call(conn, storage, *, limit, **given):
            state["calls"].append((what, limit, sorted(given)))
            return state[result]
        return call

    monkeypatch.setattr(worker, "score_all", called("score_all", "scored"))
    monkeypatch.setattr(worker, "structure_all", called("structure_all", "modelled"))
    monkeypatch.setattr(worker, "impact_all", called("impact_all", "placed"))
    monkeypatch.setattr(worker, "clinvar_all", called("clinvar_all", "fetched"))
    state["spawned"] = spawned
    return state


def _spawns(body) -> dict:
    return {name: fake.count for name, fake in body["spawned"].items()}


def _work(body) -> list:
    return [call for call in body["calls"] if call != "cleared"]


@pytest.mark.parametrize("switches, claimable, spawns", [
    ((True, True), {}, {"score": 0, "fold": 0, "evidence": 0}),
    ((True, True), {"constraint": 2}, {"score": 1, "fold": 0, "evidence": 0}),
    ((True, True), {"structure": 1, "impact": 1}, {"score": 0, "fold": 1, "evidence": 1}),
    ((False, False), {"constraint": 1, "structure": 1, "impact": 1},
     {"score": 1, "fold": 0, "evidence": 0}),
])
def test_a_sweep_starts_the_bakers_it_counts_bakes_for(body, monkeypatch, switches, claimable,
                                                       spawns):
    monkeypatch.setattr(modal_app, "STRUCTURES", switches[0])
    monkeypatch.setattr(modal_app, "EVIDENCE", switches[1])
    body["claimable"] = claimable
    summary = SWEEP.local()
    assert _spawns(body) == spawns
    assert summary["constraint_claimable"] == claimable.get("constraint", 0)
    assert (summary["structures"], summary["evidence"]) == switches
    assert summary["started"] == [name for name, count in spawns.items() if count]


@pytest.mark.parametrize("switches, asked", [
    ((True, True), {"structures": True, "evidence": True}),
    ((False, False), {}),
])
def test_a_sweep_asks_for_folds_and_evidence_only_when_switched_on(body, monkeypatch, switches,
                                                                   asked):
    from pipeline.resolver import worker

    monkeypatch.setattr(modal_app, "STRUCTURES", switches[0])
    monkeypatch.setattr(modal_app, "EVIDENCE", switches[1])
    SWEEP.local()
    assert worker.sweep.asked == [asked]


@pytest.mark.parametrize("scored, claimable, spawns", [
    ([], 3, 0),                                        # nothing it could take
    ([{"slug": "b2m", "state": "queued"}], 3, 0),       # put back: the next sweep tries again
    ([{"slug": "b2m", "state": "ready"}], 0, 0),        # nothing more waits
    ([{"slug": "b2m", "state": "ready"}], 1, 1),        # the next, in a call of its own
])
def test_score_starts_itself_again_only_while_a_bake_waits(body, scored, claimable, spawns):
    body["scored"], body["claimable"] = scored, {"constraint": claimable}
    assert SCORE.local() == scored
    assert body["spawned"]["score"].count == spawns


@pytest.mark.parametrize("folds_ready, unready, said", [
    (True, "PyMOL (pymol) is 2.5.0; every model is drawn with 3.1.0.",
     "PyMOL (pymol) is 2.5.0; every model is drawn with 3.1.0."),
    (False, None, modal_app.NO_FOLDS),
])
def test_a_fold_claims_nothing_while_none_could_be_made_here(body, folds_ready, unready, said):
    body["folds_ready"], body["unready"]["model"] = folds_ready, unready
    body["claimable"] = {"structure": 3}
    assert FOLD.local() == {"unready": said, "modelled": []}
    assert _work(body) == []
    assert _spawns(body) == {"score": 0, "fold": 0, "evidence": 0}
    assert body["calls"] == ["cleared"]


@pytest.mark.parametrize("modelled, claimable, spawns", [
    ([], 3, 0),                                                # another worker took it
    ([{"slug": "b2m", "state": "queued"}], 3, 0),               # put back: the next sweep
    ([{"slug": "b2m", "state": "ready", "source": "pdb"}], 0, 0),
    ([{"slug": "b2m", "state": "refused"}], 1, 1),
    ([{"slug": "b2m", "state": "kept"}], 1, 1),
])
def test_a_fold_is_made_one_protein_a_call_and_starts_itself_again_while_one_waits(
        body, modelled, claimable, spawns):
    body["modelled"], body["claimable"] = modelled, {"structure": claimable}
    assert FOLD.local() == {"unready": None, "modelled": modelled}
    assert _work(body) == [("structure_all", 1, ["model"])]
    assert body["calls"][-1] == "cleared"
    assert _spawns(body) == {"score": 0, "fold": spawns, "evidence": 0}


def test_no_avi_bake_is_claimed_while_the_atlas_cannot_be_asked_and_clinvar_goes_on(body):
    body["unready"]["avi"] = "The Atlas did not answer: StatusCode.UNAVAILABLE"
    body["fetched"] = [{"slug": "ttr", "state": "ready"}]
    body["claimable"] = {"impact": 2, "constraint": 1}
    outcome = EVIDENCE.local()
    assert outcome == {"unready": "The Atlas did not answer: StatusCode.UNAVAILABLE",
                       "impact": [], "clinvar": body["fetched"], "scoring": True}
    assert _work(body) == [("clinvar_all", 1, ["clinvar"])]
    # Nothing it could take waits (the AVI bakes, while the Atlas is away); TTR's
    # scoring waits no longer.
    assert _spawns(body) == {"score": 1, "fold": 0, "evidence": 0}
    assert body["calls"][-1] == "cleared"


@pytest.mark.parametrize("placed, fetched, claimable, spawns", [
    ([], [], {"impact": 1}, {"evidence": 0, "score": 0}),        # another worker took it
    ([{"slug": "ttr", "state": "ready"}], [{"slug": "ttr", "state": "ready"}],
     {"impact": 1, "constraint": 1}, {"evidence": 1, "score": 1}),
    ([{"slug": "ttr", "state": "queued"}], [], {"impact": 1}, {"evidence": 0, "score": 0}),
    ([{"slug": "ttr", "state": "ready"}], [{"slug": "prl", "state": "queued"}],
     {"clinvar": 2}, {"evidence": 0, "score": 0}),
    ([{"slug": "ttr", "state": "refused"}], [{"slug": "ttr", "state": "refused"}],
     {"constraint": 1}, {"evidence": 0, "score": 1}),
    ([{"slug": "ttr", "state": "ready"}], [], {"clinvar": 1, "constraint": 1},
     {"evidence": 1, "score": 0}),
])
def test_evidence_starts_itself_again_while_some_waits_and_the_gpu_once_a_scoring_does_not(
        body, placed, fetched, claimable, spawns):
    body["placed"], body["fetched"], body["claimable"] = placed, fetched, claimable
    outcome = EVIDENCE.local()
    assert (outcome["impact"], outcome["clinvar"]) == (placed, fetched)
    assert _work(body) == [("impact_all", 1, ["avi"]), ("clinvar_all", 1, ["clinvar"])]
    assert _spawns(body) == {"score": spawns["score"], "fold": 0, "evidence": spawns["evidence"]}
    assert body["calls"][-1] == "cleared"
