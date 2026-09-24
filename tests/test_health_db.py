"""Tests for GET /health/db.

No test here opens a socket. What is worth pinning is the reporting contract:
an unconfigured service says so plainly, and a database that cannot be reached
degrades to a report instead of an exception, because gene reads do not need
the cache to work.
"""

import pytest

from app import db, main
from app.config import settings


@pytest.fixture
def database_url(monkeypatch):
    """Point the settings at a database without connecting to one."""

    def _set(url):
        monkeypatch.setattr(settings, "database_url", url)

    return _set


def test_unconfigured_is_not_an_error(client, database_url):
    database_url(None)
    assert client.get("/health/db").json() == {"database": "unconfigured"}


def test_reports_server_version_when_reachable(client, database_url, monkeypatch):
    database_url("postgresql://example/db")
    monkeypatch.setattr(db, "server_version", lambda: "PostgreSQL 15.8")

    assert client.get("/health/db").json() == {
        "database": "ok",
        "server": "PostgreSQL 15.8",
    }


def test_unreachable_database_reports_instead_of_raising(client, database_url, monkeypatch):
    database_url("postgresql://example/db")

    def boom():
        raise OSError("connection to server at 10.0.0.1, port 5432 failed")

    monkeypatch.setattr(db, "server_version", boom)
    response = client.get("/health/db")

    assert response.status_code == 200
    assert response.json() == {"database": "error", "error": "OSError"}


def test_error_response_does_not_leak_the_connection_details(
    client, database_url, monkeypatch
):
    """The endpoint is public; the host and user belong in the log, not the body."""
    database_url("postgresql://postgres.abcdefgh:hunter2@aws-0-eu-central-1.pooler:5432/postgres")

    def boom():
        raise OSError(
            "connection to server at aws-0-eu-central-1.pooler failed: "
            "user postgres.abcdefgh"
        )

    monkeypatch.setattr(db, "server_version", boom)
    body = client.get("/health/db").text

    assert "pooler" not in body
    assert "postgres.abcdefgh" not in body
    assert "hunter2" not in body


def test_liveness_stays_free_of_the_database(client, database_url, monkeypatch):
    """/health must not touch the pool: it is the probe that runs every minute."""
    database_url("postgresql://example/db")

    def boom():
        raise AssertionError("/health queried the database")

    monkeypatch.setattr(db, "server_version", boom)

    assert client.get("/health").json() == {"status": "ok"}


def test_pool_is_not_opened_without_a_database_url(database_url):
    database_url(None)
    db.close_pool()
    db.open_pool()

    assert db.pool is None


def test_server_version_without_a_pool_is_an_error(database_url):
    database_url(None)
    db.close_pool()

    with pytest.raises(RuntimeError):
        db.server_version()
