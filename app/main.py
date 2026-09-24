"""FastAPI application wiring."""

import logging
import socket
from contextlib import asynccontextmanager

from Bio import Entrez
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from . import db, record_cache
from .config import settings
from .router import router

Entrez.email = settings.ncbi_email
Entrez.tool = "helixpeek-backend"

socket.setdefaulttimeout(settings.ncbi_timeout_seconds)

logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Hold the database pool open for the life of the worker.

    Opening here rather than at import time means a wrong or unreachable
    DATABASE_URL surfaces on /health/db instead of killing the worker during
    startup, when nothing can report why.
    """
    db.open_pool()
    record_cache.create_schema()
    try:
        yield
    finally:
        db.close_pool()


app = FastAPI(
    title="Helix Peek Backend",
    description="Minimal NCBI GenBank fetch service.",
    version="0.1.0",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(router)


@app.get("/health")
def health() -> dict:
    """Liveness check that does not spend an NCBI request."""
    return {"status": "ok"}


@app.get("/health/db")
def health_db() -> dict:
    """Readiness check for the cache database: one real round trip.

    Separate from /health so the liveness probe stays free. A failure here is
    reported, not raised: the service can still serve genes from NCBI without a
    database.

    The endpoint is public, and psycopg spells the host, port and user out in
    its connection errors -- so the message goes to the log, where operating
    the service can read it, and only the exception class comes back over HTTP.
    """
    if not settings.database_url:
        return {"database": "unconfigured"}
    try:
        return {"database": "ok", "server": db.server_version()}
    except Exception as exc:
        logger.exception("Cache database unreachable")
        return {"database": "error", "error": type(exc).__name__}
