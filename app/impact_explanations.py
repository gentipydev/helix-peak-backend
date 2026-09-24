"""Serve generated evidence by exact accession/gene identity."""
import json

from fastapi import HTTPException

from .config import settings


def read_impact_explanations(accession: str, gene: str) -> dict:
    # Request values never become paths. The header identifies each generated
    # bundle, including the transcript, assembly and original impact digest.
    for path in sorted(settings.impact_explanations_dir.glob("*.json")):
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError):
            raise HTTPException(503, "AVI explanations are temporarily unavailable.")
        if not isinstance(data, dict):
            raise HTTPException(503, "AVI explanations are temporarily unavailable.")
        if data.get("accession") != accession or data.get("gene") != gene:
            continue
        if data.get("schema_version") != 1 or data.get("scorer") != "AVI_SCORE_FEATURE_IMPORTANCE":
            raise HTTPException(503, "AVI explanations are temporarily unavailable.")
        return data
    raise HTTPException(404, "AVI explanations are not included for this gene.")
