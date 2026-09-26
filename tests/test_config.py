"""The app must not boot without a contact address for NCBI."""

import os
import subprocess
import sys


def test_missing_email_blocks_startup(tmp_path):
    """NCBI_EMAIL unset -> importing the app raises, so uvicorn refuses to run.

    Run in a subprocess with a scratch cwd so neither the ambient environment
    nor a developer's local .env can satisfy the setting.

    ``SYSTEMROOT`` is the one variable passed through, where there is one:
    Windows Python cannot initialise its sockets without it, and the import
    then dies in asyncio before it reaches the setting this test is about.
    """
    env = {"PATH": "/usr/bin:/bin", "PYTHONPATH": str(_project_root())}
    if "SYSTEMROOT" in os.environ:
        env["SYSTEMROOT"] = os.environ["SYSTEMROOT"]
    result = subprocess.run(
        [sys.executable, "-c", "import app.main"],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "ncbi_email" in result.stderr.lower()


def _project_root():
    import pathlib

    return pathlib.Path(__file__).resolve().parent.parent
