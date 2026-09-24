"""Application settings, loaded from the environment / .env file."""

from pydantic_settings import BaseSettings, SettingsConfigDict
from pathlib import Path
from typing import Optional


class Settings(BaseSettings):
    """Settings for the Helix Peek backend.

    ``ncbi_email`` has no default on purpose. NCBI rejects (and eventually
    blocks) traffic without a contact address, so an unset ``NCBI_EMAIL``
    should stop the app from booting rather than produce failures at request
    time.
    """

    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )

    ncbi_email: str

    ncbi_timeout_seconds: float = 20.0

    # Unset means "no cache": the service reads through to NCBI for every
    # request. Supabase's direct host (db.<ref>.supabase.co) resolves to IPv6
    # only, which Render cannot reach, so this must be the pooler host and
    # must carry sslmode=require.
    database_url: Optional[str] = None

    database_pool_max_size: int = 4

    database_timeout_seconds: float = 10.0

    # A GenBank record for a given accession is effectively immutable, so this
    # is about picking up NCBI's own corrections rather than about staleness.
    cache_ttl_days: int = 30

    # Generated artifacts can be mounted here in a deployed service. The local
    # default shares the exact payloads shipped by the companion mobile app.
    # A directory that is not there is reported as a fault rather than read as
    # an answer about the gene; see ``impact_explanations``.
    impact_explanations_dir: Path = (
        Path(__file__).resolve().parents[2] / "helix-peek" / "assets" / "impact_explanations"
    )

    # The Supabase project, as the REST origin rather than the database host:
    # https://<ref>.supabase.co. Track blobs are served to clients from its
    # public storage endpoint, which they fetch directly -- this service builds
    # the URL and never proxies the bytes, so neither its memory nor its egress
    # is in the path of a 9 MB ClinVar snapshot.
    supabase_url: Optional[str] = None

    # There is deliberately no service_role key here. This service only ever
    # builds public storage URLs, and it reaches its rows through
    # DATABASE_URL, so it needs no Supabase credential at all -- `supabase_url`
    # above is a public project URL. The key that bypasses row level security
    # lives only in the environment of `tool/upload_tracks.py`, which runs from
    # a developer's machine. A deployed process that never holds it cannot leak
    # it.

    tracks_bucket: str = "tracks"

    models_bucket: str = "models"


settings = Settings()
