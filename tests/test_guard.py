"""What every request is held to before a route sees it (`app/guard.py`).

The suite runs with the per-process limits off (the root `conftest.py`); each
test here sets the limits it is about, and starts from no counts.
"""

import pytest

from app import db, guard
from app.config import settings


@pytest.fixture(autouse=True)
def fresh_counts():
    guard.reset()
    yield
    guard.reset()


@pytest.fixture
def no_database(monkeypatch):
    monkeypatch.setattr(db, "pool", None)


def test_every_answer_says_not_to_sniff_or_frame_it(client):
    response = client.get("/health")
    assert response.headers["x-content-type-options"] == "nosniff"
    assert response.headers["x-frame-options"] == "DENY"
    assert response.headers["referrer-policy"] == "no-referrer"


def test_a_clients_requests_are_limited_a_minute(client, monkeypatch, no_database):
    monkeypatch.setattr(settings, "requests_per_minute", 3)
    said = [client.get("/catalog").status_code for _ in range(4)]
    # The first three reach the route (503: no database); the fourth does not.
    assert said == [503, 503, 503, 429]
    refused = client.get("/catalog")
    assert refused.status_code == 429
    assert refused.headers["retry-after"] == "60"
    assert refused.json() == {"detail": "Too many requests from one place. Ask again in a minute."}


def test_another_address_has_a_minute_of_its_own(client, monkeypatch, no_database):
    monkeypatch.setattr(settings, "requests_per_minute", 1)
    one = {"x-forwarded-for": "203.0.113.7, 10.0.0.1"}
    two = {"x-forwarded-for": "198.51.100.2"}
    assert client.get("/catalog", headers=one).status_code == 503
    assert client.get("/catalog", headers=one).status_code == 429
    assert client.get("/catalog", headers=two).status_code == 503


def test_the_service_as_a_whole_is_limited_whatever_address_is_claimed(
        client, monkeypatch, no_database):
    monkeypatch.setattr(settings, "service_requests_per_minute", 2)
    said = [client.get("/catalog", headers={"x-forwarded-for": f"203.0.113.{n}"}).status_code
            for n in range(3)]
    assert said == [503, 503, 429]


def test_posts_have_a_lower_limit_than_reads(client, monkeypatch, no_database):
    monkeypatch.setattr(settings, "posts_per_minute", 2)
    said = [client.post("/proteins/resolve", json={"gene": "BRCA1"}).status_code
            for _ in range(3)]
    assert said == [503, 503, 429]
    # Reads from the same place are still answered.
    assert client.get("/catalog").status_code == 503


def test_the_health_checks_are_never_counted(client, monkeypatch):
    monkeypatch.setattr(settings, "requests_per_minute", 1)
    monkeypatch.setattr(settings, "service_requests_per_minute", 1)
    assert [client.get("/health").status_code for _ in range(5)] == [200] * 5


@pytest.mark.parametrize("headers, peer, said", [
    ({"cf-connecting-ip": "198.51.100.9", "x-forwarded-for": "203.0.113.7"}, None,
     "198.51.100.9"),
    ({"x-forwarded-for": "203.0.113.7, 10.0.0.1, 10.0.0.2"}, ("10.0.0.3", 1), "203.0.113.7"),
    ({}, ("192.0.2.4", 5000), "192.0.2.4"),
    ({}, None, "unknown"),
])
def test_a_client_is_the_address_it_asks_from(headers, peer, said):
    scope = {"headers": [(k.encode(), v.encode()) for k, v in headers.items()], "client": peer}
    assert guard.client_of(scope) == said


def test_a_body_larger_than_the_service_reads_is_refused_unread(client, monkeypatch, no_database):
    monkeypatch.setattr(settings, "max_body_bytes", 64)
    response = client.post("/proteins/resolve", json={"gene": "B" * 100})
    assert response.status_code == 413
    assert response.json() == {"detail": "The request body is larger than this service reads."}


def test_a_body_with_no_length_is_still_held_to_the_limit(client, monkeypatch, no_database):
    monkeypatch.setattr(settings, "max_body_bytes", 64)

    def chunks():
        yield b'{"gene": "'
        yield b"B" * 100
        yield b'"}'

    response = client.post("/proteins/resolve", content=chunks(),
                           headers={"content-type": "application/json"})
    assert response.status_code == 413


def test_a_body_within_the_limit_reaches_the_route_whole(client, monkeypatch, no_database):
    monkeypatch.setattr(settings, "max_body_bytes", 64)
    # No database: the route read the body, validated it and answered 503.
    assert client.post("/proteins/resolve", json={"gene": "BRCA1"}).status_code == 503
    # And one the route reads as invalid is the route's 422, not the guard's.
    assert client.post("/proteins/resolve", json={"gene": ""}).status_code == 422


def test_a_window_forgets_its_counts_each_minute(monkeypatch):
    window = guard.Minute()
    now = [120.0]
    monkeypatch.setattr(guard.time, "monotonic", lambda: now[0])
    assert window.take("a", 1) and not window.take("a", 1)
    now[0] = 180.0
    assert window.take("a", 1)


def test_a_window_with_too_many_keys_lets_new_ones_through_uncounted():
    window = guard.Minute(keys_at_most=2)
    assert window.take("a", 1) and window.take("b", 1)
    # A third key is not kept: the service's own limit is what holds then.
    assert window.take("c", 1) and window.take("c", 1)
    assert not window.take("a", 1)
