"""Tests for GET /health."""


def test_health(client):
    assert client.get("/health").json() == {"status": "ok"}
