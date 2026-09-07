"""Application settings, loaded from the environment / .env file."""

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Settings for the HelixPeak backend.

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

    # Entrez.efetch offers no timeout argument and urllib has no default, so a
    # stalled NCBI connection would otherwise hang a request forever.
    ncbi_timeout_seconds: float = 20.0


# Instantiated at import time so a missing NCBI_EMAIL raises a ValidationError
# during startup instead of on the first request.
settings = Settings()
