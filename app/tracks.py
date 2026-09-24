"""Where each track for a protein is, and what state it is in.

This names bytes; it never carries them. A ready track resolves to a public
Supabase storage URL that the client fetches directly, so a 9 MB ClinVar
snapshot never passes through this service's 512 MB of memory or its egress.

The four booleans the app shipped with -- ``scored``, ``impactScored``,
``clinvarAvailable``, ``impactExplanationsAvailable`` -- become four states
here. A boolean can say a gene has no snapshot, but it cannot tell a gene
nobody has baked from one whose bake lands in four minutes, and saying "not yet
included" about the second is the one claim R9.3 forbids.
"""

import logging
from typing import Dict, List, Optional

from . import db
from .config import settings

logger = logging.getLogger(__name__)

# Every track a protein can have, in the order the walk meets them. A protein
# with no row for a kind reports ``absent``, so a missing row and an explicit
# "nothing is coming" say the same thing and neither is an error.
KINDS = (
    "record",
    "constraint",
    "impact",
    "clinvar",
    "structure",
    "impact_explanations",
)

_BUCKET_OF_KIND = {"structure": "models"}

_SELECT = """
select kind, state, reason, bucket, object_path, bytes, sha256,
       content_encoding, format, provenance
from protein_track
where slug = %s
"""

_SELECT_MANY = """
select slug, kind, state
from protein_track
where slug = any(%s)
"""


def bucket_for(kind: str) -> str:
    """Which bucket holds a track of this kind.

    Models live apart from the JSON tracks because they are large, immutable
    and deserve their own cache policy -- not because anything reads them
    differently.
    """
    return _BUCKET_OF_KIND.get(kind, settings.tracks_bucket)


def public_url(bucket: str, object_path: str) -> Optional[str]:
    """The URL a client fetches this object from, or None if unconfigured.

    Buckets are public-read, so this is a plain HTTPS GET and the client ships
    no credentials. Signing them would buy nothing: the data is public science,
    and a signed URL expires while a cached one should not.
    """
    if not settings.supabase_url:
        return None
    origin = settings.supabase_url.rstrip("/")
    return "{}/storage/v1/object/public/{}/{}".format(origin, bucket, object_path)


def _absent(reason: Optional[str] = None) -> dict:
    return {"state": "absent", "reason": reason, "provenance": {}}


def _row_to_track(row) -> dict:
    (kind, state, reason, bucket, object_path, size,
     sha256, content_encoding, fmt, provenance) = row
    url = None
    if state == "ready" and bucket and object_path:
        url = public_url(bucket, object_path)
        if url is None:
            # The row says the object is there and this service cannot say
            # where. Reporting it ready without a URL would hand the client a
            # track it cannot fetch, so it is a fault, named as one.
            logger.error(
                "SUPABASE_URL is unset; %s track for a protein is unreachable", kind
            )
            return {
                "state": "refused",
                "reason": "This track is stored but the service cannot address it.",
                "provenance": provenance or {},
            }
    return {
        "state": state,
        "reason": reason,
        "url": url,
        "format": fmt,
        "bytes": size,
        "sha256": sha256,
        "content_encoding": content_encoding,
        "provenance": provenance or {},
    }


def for_slug(slug: str) -> Dict[str, dict]:
    """Every track for one protein, keyed by kind, with every kind present.

    A database that is missing or broken reports every track absent rather than
    raising. The walk then draws a protein with no tracks, which is a state it
    already knows how to draw -- the same trade `record_cache` makes, where a
    dead database costs what it can show and never availability.
    """
    tracks = {kind: _absent() for kind in KINDS}
    if db.pool is None:
        return tracks
    try:
        with db.pool.connection() as conn:
            rows = conn.execute(_SELECT, (slug,)).fetchall()
    except Exception:
        logger.exception("Track read failed for %s", slug)
        return tracks
    for row in rows:
        kind = row[0]
        if kind in tracks:
            tracks[kind] = _row_to_track(row)
    return tracks


def states_for(slugs: List[str]) -> Dict[str, Dict[str, str]]:
    """Just the state of each track for many proteins, for search cards.

    A card shows what is ready; it does not need a URL for a track nobody has
    opened. One query for a page of results rather than one per protein.
    """
    states = {slug: {kind: "absent" for kind in KINDS} for slug in slugs}
    if db.pool is None or not slugs:
        return states
    try:
        with db.pool.connection() as conn:
            rows = conn.execute(_SELECT_MANY, (list(slugs),)).fetchall()
    except Exception:
        logger.exception("Track state read failed for %d proteins", len(slugs))
        return states
    for slug, kind, state in rows:
        if slug in states and kind in states[slug]:
            states[slug][kind] = state
    return states
