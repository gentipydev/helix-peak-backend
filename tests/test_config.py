"""The app must not boot without a contact address for NCBI."""

import subprocess
import sys


def test_missing_email_blocks_startup(tmp_path):
    """NCBI_EMAIL unset -> importing the app raises, so uvicorn refuses to run.

    Run in a subprocess with a scratch cwd so neither the ambient environment
    nor a developer's local .env can satisfy the setting.
    """
    result = subprocess.run(
        [sys.executable, "-c", "import app.main"],
        cwd=tmp_path,
        env={"PATH": "/usr/bin:/bin", "PYTHONPATH": str(_project_root())},
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "ncbi_email" in result.stderr.lower()


def _project_root():
    import pathlib

    return pathlib.Path(__file__).resolve().parent.parent
