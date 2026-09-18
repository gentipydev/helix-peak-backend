"""FastAPI application wiring."""

import socket

from Bio import Entrez
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .config import settings
from .router import router

Entrez.email = settings.ncbi_email
Entrez.tool = "helixpeak-backend"

socket.setdefaulttimeout(settings.ncbi_timeout_seconds)

app = FastAPI(
    title="HelixPeak Backend",
    description="Minimal NCBI GenBank fetch service.",
    version="0.1.0",
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
