"""Where the resolver reads and writes: the rows, and the storage bucket.

Every statement here runs on a psycopg connection opened with autocommit, so a
claim is visible to every other worker the moment it is made, and the writes
that have to land together run inside one `conn.transaction()`. Claims take
`FOR UPDATE SKIP LOCKED`, so two sweeps running at once never work the same
request or the same job.

The statements that already exist are reused rather than copied: a resolved
protein is upserted with `seed_catalog`'s statement, guarded so it can never
overwrite one of the twenty, and its tracks with `upload_tracks`'s.
"""

from __future__ import annotations

import urllib.request
from typing import Optional

from pipeline import seed_catalog, upload_tracks
from pipeline.resolver.resolve import IndexRow, Resolution

TRACKS_BUCKET = "tracks"

# A worker that has held a claim this long has died: its container was
# reclaimed, or it timed out. The claim goes back on the queue, or fails once
# it has been tried as often as anything is.
STALE_RESOLVE = "30 minutes"
STALE_BAKE = "90 minutes"
# A model is baked in seconds and its baker is stopped at fifteen minutes, so
# one that has run this long is a dead worker's too.
STALE_STRUCTURE = "20 minutes"
# AlphaGenome's scores and ClinVar's records are fetched in minutes, and an AVI
# bake is stopped at fifty (`local_worker.EVIDENCE_TIMEOUT`).
STALE_EVIDENCE = "60 minutes"
MAX_ATTEMPTS = 3

# What a track that was never finished says, by what was baking it.
_NEVER_SCORED = "The scorer stopped before it finished, every time it was tried."
_NEVER_MODELLED = "The model's bake stopped before it finished, every time it was tried."
_NEVER_PLACED = "AlphaGenome's bake stopped before it finished, every time it was tried."
_NEVER_FETCHED = "ClinVar's bake stopped before it finished, every time it was tried."

# How long a claim on each kind is held before it is a dead worker's, and what
# its track then says. Any other kind is held as long as a scoring is.
_STALE = {"structure": (STALE_STRUCTURE, _NEVER_MODELLED),
          "impact": (STALE_EVIDENCE, _NEVER_PLACED),
          "clinvar": (STALE_EVIDENCE, _NEVER_FETCHED)}

# The variant evidence a protein built on demand is given beside its scores:
# AlphaGenome's per-base scores, and ClinVar's records, which are placed on the
# protein's bases by the coordinate map the first carries.
EVIDENCE = ("impact", "clinvar")

# What a bake waits for, protein by protein: none of these kinds queued or
# running for it. ClinVar needs AVI's map, and ESM-2, the long step, goes last,
# so a build that opens once its scores are in has its evidence too.
WAITS_FOR = {"clinvar": ("impact",), "constraint": EVIDENCE}

# `seed_catalog`'s upsert, which ends by setting `resolved_at`, made unable to
# touch a curated row: those have a reading order, and nothing resolved does.
assert seed_catalog._UPSERT_PROTEIN.rstrip().endswith("resolved_at = now()")
UPSERT_RESOLVED = (seed_catalog._UPSERT_PROTEIN.rstrip()
                   + "\nwhere protein.catalog_order is null\n")

_JSONB_COLUMNS = ("regions", "disulfides", "structure", "provenance")


def jsonb(value):
    from psycopg.types.json import Jsonb
    return None if value is None else Jsonb(value)


def connect(url: str):
    import psycopg
    # The pooler Render and Modal reach Supabase through holds no prepared
    # statements, as `seed_catalog.py --check` found.
    return psycopg.connect(url, autocommit=True, prepare_threshold=None)


class TrackStorage(upload_tracks.Storage):
    """The uploader's storage client, and the public read every client makes."""

    def __init__(self, url: str, key: str):
        super().__init__(url, key)
        self.public = url.rstrip("/") + "/storage/v1/object/public"

    def get(self, bucket: str, path: str) -> bytes:
        with urllib.request.urlopen(f"{self.public}/{bucket}/{path}", timeout=120) as response:
            return response.read()


# ------------------------------------------------------------ requests


def reap(conn) -> None:
    """Put back what a dead worker held, or fail it once it has had its tries."""
    conn.execute(
        """
        update resolve_request
        set state = case when attempts >= %(max)s then 'failed' else 'queued' end,
            reason = case when attempts >= %(max)s
                          then 'The worker stopped before it finished, every time it was tried.'
                          else reason end,
            finished_at = case when attempts >= %(max)s then now() else null end,
            started_at = case when attempts >= %(max)s then started_at else null end
        where state = 'running' and started_at < now() - %(stale)s::interval
        """,
        {"max": MAX_ATTEMPTS, "stale": STALE_RESOLVE},
    )
    # Each kind with a patience of its own, then every other at a scoring's.
    groups = [("kind = %(kind)s", {"kind": kind}, stale, never)
              for kind, (stale, never) in _STALE.items()]
    groups.append(("kind <> all(%(kinds)s::text[])", {"kinds": list(_STALE)},
                   STALE_BAKE, _NEVER_SCORED))
    with conn.transaction():
        for which, said, stale, never in groups:
            dead = conn.execute(
                f"""
                update bake_job
                set state = 'failed', error = 'The worker stopped before it finished.',
                    finished_at = now()
                where state = 'running' and attempts >= %(max)s and {which}
                  and started_at < now() - %(stale)s::interval
                returning slug, kind
                """,
                {"max": MAX_ATTEMPTS, "stale": stale, **said},
            ).fetchall()
            for slug, kind in dead:
                refuse_track(conn, slug, kind, never)
            conn.execute(
                f"""
                update bake_job set state = 'queued', started_at = null
                where state = 'running' and {which}
                  and started_at < now() - %(stale)s::interval
                """,
                {"stale": stale, **said},
            )


def claim_request(conn) -> Optional[dict]:
    row = conn.execute(
        """
        update resolve_request
        set state = 'running', attempts = attempts + 1, started_at = now()
        where id = (
            select id from resolve_request where state = 'queued'
            order by requested_at
            for update skip locked
            limit 1
        )
        returning id, gene, uniprot, slug, attempts
        """
    ).fetchone()
    if row is None:
        return None
    return dict(zip(("id", "gene", "uniprot", "slug", "attempts"), row))


def queued_requests(conn) -> int:
    return conn.execute(
        "select count(*) from resolve_request where state = 'queued'"
    ).fetchone()[0]


def protein_for_gene(conn, gene: str) -> Optional[str]:
    """The slug of the protein this gene already makes, listed or resolved."""
    row = conn.execute(
        "select slug from protein where gene = %s and taxon_id = 9606", (gene,)
    ).fetchone()
    return row[0] if row else None


def index_row(conn, uniprot: str, gene: str) -> Optional[IndexRow]:
    row = conn.execute(
        """
        select uniprot, gene, name, length, refseq_nuc, refseq_prot,
               chrom_acc, chrom_start, chrom_end, refseqgene, buildable, unavailable_reason
        from protein_index where uniprot = %s and gene = %s
        """,
        (uniprot, gene),
    ).fetchone()
    return IndexRow(*row) if row else None


def mane_release(conn) -> Optional[str]:
    row = conn.execute("select mane from protein_index_release").fetchone()
    return row[0] if row else None


def write_resolution(conn, request_id: int, resolution: Resolution, record: dict,
                     resolver_version: int, structures: bool = False,
                     evidence: bool = False) -> None:
    """Everything a resolved protein is, in one transaction: its row and aliases,
    its record track ready, its constraint track pending behind a queued bake,
    and the request done.

    With `structures`, its structure track is pending behind a bake of its own
    as well, and with `evidence` its AVI and ClinVar tracks are
    (`queue_evidence`). Only a worker that can make them asks for either (the
    Mac's): a track is never left pending where nothing will bake it."""
    slug = resolution.target.slug
    protein = {key: (jsonb(value) if key in _JSONB_COLUMNS else value)
               for key, value in resolution.protein.items()}
    with conn.transaction():
        written = conn.execute(UPSERT_RESOLVED + " returning slug", protein).fetchone()
        if written is None:
            raise RuntimeError(f"{slug} is a curated protein; the resolver does not write over it.")
        conn.execute("delete from protein_alias where slug = %s", (slug,))
        for alias, kind in resolution.aliases:
            conn.execute(seed_catalog._UPSERT_ALIAS, (slug, alias, kind))
        conn.execute(upload_tracks._UPSERT, {**record, "provenance": jsonb(record["provenance"])})
        conn.execute(
            """
            insert into protein_track (slug, kind, state, reason, format, provenance)
            values (%s, 'constraint', 'pending', null, 'json', '{}'::jsonb)
            on conflict (slug, kind) do update
            set state = 'pending', reason = null, updated_at = now()
            where protein_track.state <> 'ready'
            """,
            (slug,),
        )
        # The reader who asked may stop the scoring too (0011).
        conn.execute(
            """
            insert into bake_job (slug, kind, asker)
            select %s, 'constraint', asker from resolve_request where id = %s
            on conflict (slug, kind) where state in ('queued', 'running') do nothing
            """,
            (slug, request_id),
        )
        if structures:
            conn.execute(
                """
                insert into protein_track (slug, kind, state, reason, format, provenance)
                values (%s, 'structure', 'pending', null, 'fsceneb', '{}'::jsonb)
                on conflict (slug, kind) do update
                set state = 'pending', reason = null, updated_at = now()
                where protein_track.state <> 'ready'
                """,
                (slug,),
            )
            conn.execute(
                """
                insert into bake_job (slug, kind, asker)
                select %s, 'structure', asker from resolve_request where id = %s
                on conflict (slug, kind) where state in ('queued', 'running') do nothing
                """,
                (slug, request_id),
            )
        if evidence:
            asker = conn.execute(
                "select asker from resolve_request where id = %s", (request_id,)).fetchone()
            queue_evidence(conn, slug, asker[0] if asker else None)
        finish_request(conn, request_id, "done", slug=slug, resolver_version=resolver_version)


def queue_evidence(conn, slug: str, asker: Optional[str] = None) -> list:
    """A protein's AVI and ClinVar tracks pending, each behind a queued bake.
    Returns the kinds left pending: every one not ready already.

    A kind whose track is ready is left as it is, and one with a bake queued
    or running keeps that bake: asking twice queues nothing twice. Only a protein
    resolved on demand is ever queued; one of the twenty is baked by hand.
    Called inside `write_resolution`'s transaction, or in one of its own for a
    protein built before its evidence was (a backfill).
    """
    with conn.transaction():
        found = conn.execute(
            "select 1 from protein where slug = %s and catalog_order is null "
            "and resolver_version > 0", (slug,)).fetchone()
        if found is None:
            raise RuntimeError(f"{slug} is not a protein resolved on demand; its evidence "
                               f"is baked by hand.")
        ready = {kind for (kind,) in conn.execute(
            "select kind from protein_track where slug = %s and kind = any(%s) "
            "and state = 'ready'", (slug, list(EVIDENCE))).fetchall()}
        queued = []
        for kind in EVIDENCE:
            if kind in ready:
                continue
            conn.execute(
                """
                insert into protein_track (slug, kind, state, reason, format, provenance)
                values (%s, %s, 'pending', null, 'json', '{}'::jsonb)
                on conflict (slug, kind) do update
                set state = 'pending', reason = null, updated_at = now()
                where protein_track.state <> 'ready'
                """,
                (slug, kind),
            )
            conn.execute(
                """
                insert into bake_job (slug, kind, asker) values (%s, %s, %s)
                on conflict (slug, kind) where state in ('queued', 'running') do nothing
                """,
                (slug, kind, asker),
            )
            queued.append(kind)
        return queued


def finish_request(conn, request_id: int, state: str, *, slug: Optional[str] = None,
                   reason: Optional[str] = None, resolver_version: Optional[int] = None) -> None:
    conn.execute(
        """
        update resolve_request
        set state = %(state)s, reason = %(reason)s, slug = coalesce(%(slug)s, slug),
            resolver_version = coalesce(%(version)s, resolver_version), finished_at = now()
        where id = %(id)s
        """,
        {"state": state, "reason": reason, "slug": slug, "version": resolver_version,
         "id": request_id},
    )


def requeue_request(conn, request_id: int) -> None:
    conn.execute(
        "update resolve_request set state = 'queued', started_at = null where id = %s",
        (request_id,),
    )


# ------------------------------------------------------------ bakes


def claim_bake(conn, kind: str) -> Optional[dict]:
    """The oldest queued bake of this kind for a protein resolved on demand,
    once nothing it waits for (`WAITS_FOR`) is queued or running for that
    protein.

    The twenty are baked by hand, and their rows carry no table to rebuild a
    `Target` from, so a job for one of them is left where it is.
    """
    row = conn.execute(
        """
        update bake_job
        set state = 'running', attempts = attempts + 1, started_at = now()
        where id = (
            select j.id from bake_job j
            join protein p on p.slug = j.slug
            where j.state = 'queued' and j.kind = %(kind)s and p.resolver_version > 0
              and not exists (
                  select 1 from bake_job w
                  where w.slug = j.slug and w.kind = any(%(waits)s::text[])
                    and w.state in ('queued', 'running'))
            order by j.requested_at
            for update of j skip locked
            limit 1
        )
        returning id, slug, kind, attempts
        """,
        {"kind": kind, "waits": list(WAITS_FOR.get(kind, ()))},
    ).fetchone()
    if row is None:
        return None
    return dict(zip(("id", "slug", "kind", "attempts"), row))


def queued_bakes(conn, kind: str) -> int:
    return conn.execute(
        "select count(*) from bake_job where state = 'queued' and kind = %s", (kind,)
    ).fetchone()[0]


def note_progress(conn, job_id: int, done: int, total: int, left: float) -> bool:
    """How far a running bake has got: `done` of `total`, and the seconds its
    baker estimates are left. Read by `GET /proteins/resolve/{gene}`, never by
    the worker; a job no longer running is left as it ended.

    Returns whether the job is still running: false once the reader who asked
    has stopped it (0011)."""
    written = conn.execute(
        """
        update bake_job
        set progress_done = %s, progress_total = %s, progress_left = %s, progress_at = now()
        where id = %s and state = 'running'
        returning id
        """,
        (done, total, left, job_id),
    ).fetchone()
    return written is not None


def still_running(conn, job_id: int) -> bool:
    """Whether a bake is still running: false once it has been stopped."""
    row = conn.execute("select state from bake_job where id = %s", (job_id,)).fetchone()
    return row is not None and row[0] == "running"


def protein(conn, slug: str) -> Optional[dict]:
    columns = seed_catalog._PROTEIN_COLUMNS
    row = conn.execute(
        f"select {', '.join(columns)} from protein where slug = %s", (slug,)
    ).fetchone()
    return dict(zip(columns, row)) if row else None


def track_state(conn, slug: str, kind: str) -> tuple:
    """A track's state and the sentence it carries: ('absent', None) where the
    protein has no row of that kind."""
    row = conn.execute(
        "select state, reason from protein_track where slug = %s and kind = %s", (slug, kind),
    ).fetchone()
    return (row[0], row[1]) if row else ("absent", None)


def ready_track(conn, slug: str, kind: str) -> Optional[dict]:
    row = conn.execute(
        """
        select bucket, object_path, sha256 from protein_track
        where slug = %s and kind = %s and state = 'ready'
        """,
        (slug, kind),
    ).fetchone()
    return dict(zip(("bucket", "object_path", "sha256"), row)) if row else None


def finish_bake(conn, job_id: int, track: dict) -> None:
    """The track ready and its job done, together."""
    with conn.transaction():
        conn.execute(upload_tracks._UPSERT, {**track, "provenance": jsonb(track["provenance"])})
        conn.execute(
            "update bake_job set state = 'done', error = null, finished_at = now() where id = %s",
            (job_id,),
        )


def finish_structure(conn, job_id: int, track: dict, structure: dict) -> None:
    """A model's track ready, the fold page's words and chains on its protein's
    row, and its job done, together.

    `structure` is the row's `structure` column as the twenty carry it,
    `{"chrome": ..., "chains": ...}`. It is written only on a row resolved on
    demand: one of the twenty is never written over.
    """
    with conn.transaction():
        written = conn.execute(
            """
            update protein set structure = %s
            where slug = %s and catalog_order is null
            returning slug
            """,
            (jsonb(structure), track["slug"]),
        ).fetchone()
        if written is None:
            raise RuntimeError(
                f"{track['slug']} is a curated protein; the resolver does not write over it.")
        conn.execute(upload_tracks._UPSERT, {**track, "provenance": jsonb(track["provenance"])})
        conn.execute(
            "update bake_job set state = 'done', error = null, finished_at = now() where id = %s",
            (job_id,),
        )


def refuse_track(conn, slug: str, kind: str, reason: str) -> None:
    conn.execute(
        """
        update protein_track set state = 'refused', reason = %s, updated_at = now()
        where slug = %s and kind = %s and state <> 'ready'
        """,
        (reason, slug, kind),
    )


def refuse_bake(conn, job_id: int, slug: str, kind: str, reason: str,
                error: Optional[str] = None) -> None:
    """The bake declined: its track says why, and its job is over.

    `error` is what the job keeps where that is more than a reader is told: a
    reader sees the track's reason, and an operator the job's error."""
    with conn.transaction():
        refuse_track(conn, slug, kind, reason)
        conn.execute(
            "update bake_job set state = 'failed', error = %s, finished_at = now() where id = %s",
            (error or reason, job_id),
        )


def requeue_bake(conn, job_id: int, error: str) -> None:
    conn.execute(
        "update bake_job set state = 'queued', error = %s, started_at = null where id = %s",
        (error, job_id),
    )
