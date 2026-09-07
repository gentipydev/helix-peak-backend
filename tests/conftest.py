"""Shared test fixtures.

NCBI_EMAIL is set before importing the app: app.config instantiates Settings at
import time, so without it the import itself fails. That is deliberate
production behaviour (see test_missing_email_blocks_startup), which the rest of
the suite has to work around.
"""

import os
import pathlib

os.environ.setdefault("NCBI_EMAIL", "tests@example.com")

import pytest
from fastapi.testclient import TestClient

from app.main import app

FIXTURES = pathlib.Path(__file__).parent / "fixtures"


@pytest.fixture
def client():
    return TestClient(app)


@pytest.fixture
def genbank_text():
    """The real NG_007114 GenBank record, saved from NCBI."""
    return (FIXTURES / "ng_007114.gb").read_text()
