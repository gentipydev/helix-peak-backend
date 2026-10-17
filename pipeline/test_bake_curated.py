"""`bake_curated.py`'s plan, offline: which steps a row gets, in which
environment, and what each is handed. Nothing is baked or fetched.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import bake_curated  # noqa: E402
from pipeline.structure import pdb  # noqa: E402
from pipeline.targets import BY_SLUG  # noqa: E402

ENVIRON = {"PATH": "/usr/bin", "NCBI_EMAIL": "someone@example.com",
           bake_curated.SECRET: "the key", "HOME": "/home/x"}


def _plan(target, working_copy=False, environ=ENVIRON):
    return bake_curated.plan(target, dict(environ), Path("/home/x"), working_copy)


def test_every_step_is_handed_this_row_and_no_other():
    steps = _plan(BY_SLUG["insulin"])
    assert [s["step"] for s in steps] == ["record", "structure", "folding", "constraint",
                                          "impact", "clinvar", "locus", "check", "check",
                                          "check"]
    for s in steps:
        assert s["command"][-2:] == ["--target", "insulin"], s["step"]
        assert s["command"][1] == "-u"


def test_each_step_runs_in_its_own_environment():
    python = {s["step"]: Path(s["command"][0]) for s in _plan(BY_SLUG["insulin"], True)}
    assert python["record"] == python["clinvar"] == python["locus"] == python["check"] \
        == bake_curated.VENV
    assert python["structure"] == python["folding"] == python["bundle"] \
        == bake_curated.STRUCTURE
    assert python["constraint"] == bake_curated.ESM
    assert python["impact"] == bake_curated.IMPACT


def test_the_key_goes_to_the_impact_step_alone():
    from pipeline.resolver.local_worker import default_uv

    for s in _plan(BY_SLUG["insulin"]):
        if s["step"] == "impact":
            assert s["env"][bake_curated.SECRET] == "the key"
            # `uv` where the worker finds it, for the bake's GENCODE lookup.
            assert s["env"]["PATH"].startswith(str(default_uv(Path("/home/x")).parent))
        else:
            assert bake_curated.SECRET not in s["env"], s["step"]
    ncbi = {s["step"]: s["env"].get("NCBI_EMAIL") for s in _plan(BY_SLUG["insulin"])}
    assert ncbi["record"] == ncbi["clinvar"] == "someone@example.com"


def test_a_row_cut_from_a_working_copy_is_held_to_the_apps_check():
    names = [s["step"] for s in _plan(BY_SLUG["insulin"], working_copy=True)]
    assert names.index("folding") + 1 == names.index("bundle")
    assert "bundle" not in [s["step"] for s in _plan(BY_SLUG["insulin"])]


def test_a_row_without_a_track_skips_its_step():
    bare = dataclasses.replace(BY_SLUG["insulin"], scored=False, impact_scored=False,
                               clinvar_available=False)
    assert [s["step"] for s in _plan(bare)] == ["record", "structure", "folding", "locus",
                                                "check", "check", "check"]


def test_a_run_resumes_from_a_step_or_takes_only_some():
    steps = _plan(BY_SLUG["insulin"])
    assert [s["step"] for s in bake_curated.chosen(steps, [], "impact")] == [
        "impact", "clinvar", "locus", "check", "check", "check"]
    assert [s["step"] for s in bake_curated.chosen(steps, ["structure", "folding"], None)] == [
        "structure", "folding"]
    with pytest.raises(SystemExit, match="no 'bundle' step"):
        bake_curated.chosen(steps, [], "bundle")


def test_the_key_is_found_where_the_bakes_keep_it(tmp_path, monkeypatch):
    monkeypatch.setattr(bake_curated, "BACKEND", tmp_path / "backend")
    (tmp_path / ".env").write_text(f"# keys\nexport {bake_curated.SECRET}='from home'\n")
    assert bake_curated.api_key({}, tmp_path) == "from home"
    assert bake_curated.api_key({bake_curated.SECRET: "here"}, tmp_path) == "here"
    assert bake_curated.api_key({}, tmp_path / "nowhere") is None


def test_a_dry_run_names_every_command_and_no_key(monkeypatch, capsys):
    monkeypatch.setattr(bake_curated.os, "environ", dict(ENVIRON))
    assert bake_curated.main(["--target", "insulin", "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "[impact] (with the AlphaGenome key)" in out and "the key" not in out
    assert "pipeline/mock/build_gene_record.py --target insulin" in out


def test_a_record_is_held_to_the_protein_its_fold_was_chosen_for(tmp_path, monkeypatch):
    target = dataclasses.replace(BY_SLUG["insulin"], slug="insulin-copy")
    monkeypatch.setattr(pdb, "WORKING", tmp_path / "working")
    (tmp_path / "working").mkdir()
    (tmp_path / "working" / "insulin-copy.pdb").write_text(
        "REMARK   1 3I40 MODEL 1, AS DRAWN: THE RECORD'S NUMBERING\n")
    data = tmp_path / "data"
    record = data / target.mock_asset
    record.parent.mkdir(parents=True)
    record.write_text(json.dumps({"protein": {"translation": "MALWMR"}}))
    monkeypatch.setattr("pipeline.paths.DATA", data)

    def sidecar(translation: str) -> None:
        (tmp_path / "working" / "insulin-copy.json").write_text(json.dumps(
            {"sha256": {"translation": hashlib.sha256(translation.encode()).hexdigest()}}))

    sidecar("MALWMR")
    assert bake_curated.translation_held(target) is None
    sidecar("MALWMRL")
    assert "another protein" in bake_curated.translation_held(target)
    assert bake_curated.translation_held(BY_SLUG["insulin"]) is None
