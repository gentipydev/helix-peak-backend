"""The uploader's own decisions, offline: nothing is written and no network is
reached. A structure goes up with the scene the app's importer compiles from
it; a target that cannot be uploaded is refused out loud.
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import shutil
import sys
import urllib.error
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import upload_tracks  # noqa: E402
from pipeline.paths import CLIENT, DATA  # noqa: E402
from pipeline.targets import BY_SLUG  # noqa: E402

INSULIN_GLB = DATA / BY_SLUG["insulin"].structure_asset
needs_model = pytest.mark.skipif(not INSULIN_GLB.exists(),
                                 reason="run fetch_tracks.py --kind structure first")

# Insulin's stored scene, as its row names it: the `.fsceneb` the build hook
# compiled from its `.glb` with flutter_scene 0.23.0.
INSULIN_SCENE_SHA256 = "2011fd6cb8ac1bc11212e350e12e2b6ca77b9a8f55c46c0de11c08e0c0e46700"


def _run(monkeypatch, *arguments: str) -> int:
    monkeypatch.setattr(sys, "argv", ["upload_tracks.py", *arguments])
    return upload_tracks.main()


@needs_model
def test_a_curated_structure_goes_up_with_the_scene_compiled_from_it(monkeypatch, capsys):
    compiled = []

    def compile_scene(glb: Path, app: Path, dart: Path) -> bytes:
        compiled.append((glb, app, dart))
        return b"FSCB" + glb.read_bytes()[:64]

    monkeypatch.setattr(upload_tracks, "compiled_scene", compile_scene)
    assert _run(monkeypatch, "--kind", "structure", "--target", "insulin", "--dry-run",
                "--dart", "/nowhere/dart") == 0
    glb = INSULIN_GLB.read_bytes()
    scene = b"FSCB" + glb[:64]
    out = capsys.readouterr().out
    assert f"structure/insulin.{hashlib.sha256(glb).hexdigest()[:12]}.glb" in out
    assert f"structure/insulin.{hashlib.sha256(scene).hexdigest()[:12]}.fsceneb" in out
    assert compiled == [(INSULIN_GLB, CLIENT, Path("/nowhere/dart"))]


@needs_model
def test_a_structure_with_no_scene_is_refused_and_the_run_says_so(monkeypatch, capsys):
    def no_scene(*_):
        raise RuntimeError("the importer failed")

    monkeypatch.setattr(upload_tracks, "compiled_scene", no_scene)
    assert _run(monkeypatch, "--kind", "structure", "--target", "insulin", "--dry-run") == 1
    said = capsys.readouterr().err
    assert "insulin: REFUSED -- no scene could be compiled: the importer failed" in said
    assert "1 refused: insulin" in said


def test_an_unknown_target_is_no_target(monkeypatch):
    with pytest.raises(SystemExit, match="no such target"):
        _run(monkeypatch, "--kind", "record", "--target", "no-such-protein", "--dry-run")


def test_against_is_for_a_dry_run(monkeypatch):
    with pytest.raises(SystemExit):
        _run(monkeypatch, "--kind", "record", "--against", "https://example.invalid")


def test_a_checkout_pinning_another_importer_compiles_nothing(tmp_path):
    (tmp_path / "pubspec.lock").write_text(
        'packages:\n  flutter_scene:\n    dependency: "direct main"\n    version: "0.9.0"\n')
    with pytest.raises(RuntimeError, match="0.9.0"):
        upload_tracks.compiled_scene(tmp_path / "model.glb", tmp_path, Path("dart"))


def test_the_dart_is_the_workers_unless_one_is_named(monkeypatch):
    monkeypatch.setenv("RESOLVER_DART", "/opt/dart/bin/dart")
    assert upload_tracks.default_dart() == Path("/opt/dart/bin/dart")
    monkeypatch.delenv("RESOLVER_DART")
    from pipeline.resolver.local_worker import default_dart
    assert upload_tracks.default_dart() == default_dart(Path.home())


class _Answer(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False


def test_a_dry_run_says_whether_the_service_serves_those_bytes(monkeypatch):
    row = {"slug": "insulin", "kind": "structure", "sha256": "scene",
           "provenance": {"pdb": "3I40", "glb": {"sha256": "model"}}}
    served = {"structure": {"state": "ready", "sha256": "scene",
                            "provenance": {"glb": {"sha256": "model"}}}}
    monkeypatch.setattr(upload_tracks.urllib.request, "urlopen",
                        lambda *_, **__: _Answer(json.dumps(served).encode()))
    assert upload_tracks.served_as("https://example.invalid", row) == "same"
    served["structure"]["provenance"]["glb"]["sha256"] = "another model"
    assert upload_tracks.served_as("https://example.invalid", row) == "differs"
    served["structure"]["state"] = "absent"
    assert upload_tracks.served_as("https://example.invalid", row) == "not ready (absent)"

    def missing(*_, **__):
        raise urllib.error.HTTPError("https://example.invalid", 404, "Not Found", None, None)

    monkeypatch.setattr(upload_tracks.urllib.request, "urlopen", missing)
    assert upload_tracks.served_as("https://example.invalid", row) == "not served (404)"


@needs_model
@pytest.mark.skipif(os.environ.get("HELIXPEEK_COMPILE") != "1"
                    or shutil.which(str(upload_tracks.default_dart())) is None,
                    reason="set HELIXPEEK_COMPILE=1 to run the app's importer (needs Dart and the app)")
def test_the_importer_compiles_insulins_stored_scene_byte_for_byte():
    scene = upload_tracks.compiled_scene(INSULIN_GLB, CLIENT, upload_tracks.default_dart())
    assert hashlib.sha256(scene).hexdigest() == INSULIN_SCENE_SHA256
