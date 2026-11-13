"""What every request is held to before a route sees it.

- Its body: none the service reads is more than a small JSON object (a gene
  and an install's id), so a body past `MAX_BODY_BYTES` is a 413, read no
  further. Without this a client could make the service hold any size of
  body in memory while it parses it.
- How often: each client `REQUESTS_PER_MINUTE` requests a minute, of which
  `POSTS_PER_MINUTE` may be POSTs, and the whole service
  `SERVICE_REQUESTS_PER_MINUTE`. Past one, a 429 with `Retry-After`. The
  health checks are never counted: the keep-alive pings them.

There are no accounts, so a client is the address it asks from: Cloudflare's
`CF-Connecting-IP` where it is there, else the first address in
`X-Forwarded-For` (which Render's proxy sets), else the socket's peer. A client
can name another address in those headers, so a client's own limit is a
courtesy that a script can step round; the service's limit holds whatever is
claimed, and what costs money is capped in the database (`resolves.py`,
`budget.py`), not here. The counts are this process's, kept for a minute.

And every answer says, as an API's should, not to sniff its type or frame it.
"""

import json
import threading
import time
from typing import Dict, Optional

from .config import settings

# Never counted: the keep-alive pings both (`.github/workflows/keepalive.yml`).
EXEMPT = frozenset(("/health", "/health/db"))

_HEADERS = (
    (b"x-content-type-options", b"nosniff"),
    (b"x-frame-options", b"DENY"),
    (b"referrer-policy", b"no-referrer"),
)

_WITH_BODY = frozenset(("POST", "PUT", "PATCH"))


class Minute:
    """Requests counted by key in the current minute (a fixed window)."""

    def __init__(self, keys_at_most: int = 50_000):
        self._lock = threading.Lock()
        self._minute: Optional[int] = None
        self._counts: Dict[str, int] = {}
        self._keys_at_most = keys_at_most

    def take(self, key: str, limit: int) -> bool:
        """Count one request for `key`, where its minute has room. 0 is no limit.

        Past `keys_at_most` keys in one minute, a new key is let through
        uncounted rather than kept: many addresses at once are what the
        service's own limit is for, and memory stays bounded."""
        if not limit:
            return True
        minute = int(time.monotonic() // 60)
        with self._lock:
            if minute != self._minute:
                self._minute = minute
                self._counts.clear()
            count = self._counts.get(key)
            if count is None:
                if len(self._counts) >= self._keys_at_most:
                    return True
                count = 0
            if count >= limit:
                return False
            self._counts[key] = count + 1
            return True

    def clear(self) -> None:
        with self._lock:
            self._counts.clear()
            self._minute = None


clients = Minute()
posts = Minute()
service = Minute()


def reset() -> None:
    """Forget every count (for tests)."""
    for window in (clients, posts, service):
        window.clear()


def _header(scope, name: bytes) -> Optional[str]:
    for key, value in scope.get("headers") or ():
        if key == name:
            return value.decode("latin-1")
    return None


def client_of(scope) -> str:
    """The address a request is counted under (see the module's docstring)."""
    said = _header(scope, b"cf-connecting-ip")
    if said and said.strip():
        return said.strip()
    forwarded = _header(scope, b"x-forwarded-for")
    if forwarded:
        first = forwarded.split(",")[0].strip()
        if first:
            return first
    peer = scope.get("client")
    return peer[0] if peer else "unknown"


async def _refuse(send, status: int, detail: str, retry_after: Optional[int] = None) -> None:
    body = json.dumps({"detail": detail}).encode()
    headers = [(b"content-type", b"application/json"),
               (b"content-length", str(len(body)).encode()), *_HEADERS]
    if retry_after is not None:
        headers.append((b"retry-after", str(retry_after).encode()))
    await send({"type": "http.response.start", "status": status, "headers": headers})
    await send({"type": "http.response.body", "body": body})


class Guard:
    """The ASGI middleware that holds every request to the above."""

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return

        if scope.get("path") not in EXEMPT:
            who = client_of(scope)
            busy = None
            if not service.take("*", settings.service_requests_per_minute):
                busy = "The service is answering as many requests as it can. Ask again in a minute."
            elif not clients.take(who, settings.requests_per_minute):
                busy = "Too many requests from one place. Ask again in a minute."
            elif (scope.get("method") in _WITH_BODY
                  and not posts.take(who, settings.posts_per_minute)):
                busy = "Too many requests from one place. Ask again in a minute."
            if busy is not None:
                await _refuse(send, 429, busy, retry_after=60)
                return

        if scope.get("method") in _WITH_BODY:
            limit = settings.max_body_bytes
            declared = _header(scope, b"content-length")
            if declared is not None:
                try:
                    too_long = int(declared) > limit
                except ValueError:
                    await _refuse(send, 400, "Content-Length is not a number.")
                    return
                if too_long:
                    await _refuse(send, 413, "The request body is larger than this service reads.")
                    return
            # Read whole, and never past the limit, whatever the header said
            # (or where none did); the route is handed it as it came.
            body = b""
            more = True
            while more:
                message = await receive()
                if message["type"] == "http.disconnect":
                    return
                body += message.get("body", b"")
                if len(body) > limit:
                    await _refuse(send, 413, "The request body is larger than this service reads.")
                    return
                more = message.get("more_body", False)
            replay = [{"type": "http.request", "body": body, "more_body": False}]

            async def receive_again():
                if replay:
                    return replay.pop()
                return await receive()

            downstream = receive_again
        else:
            downstream = receive

        async def send_with_headers(message):
            if message["type"] == "http.response.start":
                held = {key.lower() for key, _ in message.get("headers") or ()}
                message = {**message, "headers": [
                    *(message.get("headers") or ()),
                    *(pair for pair in _HEADERS if pair[0] not in held)]}
            await send(message)

        await self.app(scope, downstream, send_with_headers)
