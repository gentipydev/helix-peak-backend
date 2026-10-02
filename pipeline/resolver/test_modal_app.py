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
            "mock/build_gene_record.py", "constraint/score_protein.py",
            "upload_tracks.py", "seed_catalog.py", "targets.py"} <= sent
    assert all(path.endswith(".py") for path in sent)
    assert not [path for path in sent if path.startswith(("data/", "resolver/fixtures/"))]


def test_a_venv_or_baked_data_never_goes_up():
    for path in ("structure/venv/lib/python3.12/site-packages/torch/__init__.py",
                 "impact/venv/bin/activate_this.py",
                 "data/assets/mock/gene_ins.json", "data/notes.py",
                 "trafficking/fixtures/P01308.json", "resolver/__pycache__/resolve.cpython-312.pyc"):
        assert modal_app.not_code(Path(path)), path
    assert not modal_app.not_code(Path("resolver/worker.py"))
