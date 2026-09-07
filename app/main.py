"""FastAPI application wiring."""

import socket

from Bio import Entrez
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .config import settings
from .router import router

# Set once at startup rather than per call. NCBI uses this to contact you
# before throttling or blocking, so it must be a real address.
Entrez.email = settings.ncbi_email
Entrez.tool = "helixpeak-backend"

# Entrez.efetch exposes no timeout argument, so bound it at the socket layer.
# The resulting socket.timeout is an OSError and lands in the router's 502 path.
socket.setdefaulttimeout(settings.ncbi_timeout_seconds)

app = FastAPI(
    title="HelixPeak Backend",
    description="Minimal NCBI GenBank fetch service.",
    version="0.1.0",
)

# Permissive by design: this is a local, single-user dev service, and Flutter
# web (flutter run -d chrome) is otherwise blocked by the browser. Tighten this
# before the service is exposed anywhere real.
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
