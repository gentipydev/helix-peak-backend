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
        returning id, slug, kind, attempts, refresh
        """,
        {"kind": kind, "waits": list(WAITS_FOR.get(kind, ()))},
    ).fetchone()
    if row is None:
        return None
    return dict(zip(("id", "slug", "kind", "attempts", "refresh"), row))


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


# ------------------------------------------------------------ folds (0012)

# A search that found nothing suitable is believed this long, and then made
# again: four of the PDB's weekly releases.
NONE_EXPIRES = "28 days"

# What a fold is drawn over beside the record's letters: the protein's table of
# regions and its disulfides, as one digest the database takes, so that the
# worker's reading and the publishing transaction's compare the same thing.
_REGIONS = ("md5(coalesce(p.provenance -> 'regions_from', 'null'::jsonb)::text"
            " || coalesce(p.disulfides, 'null'::jsonb)::text)")

_TRACK_COLUMNS = ("slug", "kind", "state", "reason", "bucket", "object_path", "bytes", "sha256",
                  "format", "provenance")

# A track row in any state, written whole: a fold's tracks are written as they
# were made, or put back as a generation recorded them.
_PUT_TRACK = """
insert into protein_track
    (slug, kind, state, reason, bucket, object_path, bytes, sha256,
     content_encoding, format, provenance, updated_at)
values (%(slug)s, %(kind)s, %(state)s, %(reason)s, %(bucket)s, %(object_path)s, %(bytes)s,
        %(sha256)s, null, %(format)s, %(provenance)s, now())
on conflict (slug, kind) do update set
    state = excluded.state, reason = excluded.reason, bucket = excluded.bucket,
    object_path = excluded.object_path, bytes = excluded.bytes,
    sha256 = excluded.sha256, content_encoding = null,
    format = excluded.format, provenance = excluded.provenance,
    updated_at = now()
"""


class Superseded(RuntimeError):
    """A bake that is no longer this worker's to finish, or whose protein
    changed while it was made. Nothing it made is written."""


def folds_ready(conn) -> bool:
    """Whether the database holds what a fold is published with (0012)."""
    return conn.execute("select to_regclass('fold_generation') is not null").fetchone()[0]


def regions_of(conn, slug: str) -> Optional[str]:
    """The digest of a protein's regions and disulfides (`_REGIONS`)."""
    row = conn.execute(f"select {_REGIONS} from protein p where p.slug = %s", (slug,)).fetchone()
    return row[0] if row else None


def fold_prior(conn, slug: str) -> dict:
    """What earlier searches decided for a protein, for the next one:

    - `none`: the latest whole search that found nothing suitable and has not
      expired, which a search for the same inputs under the same rule need not
      make again;
    - `progress`: where the latest search stopped, where it was cut short
      (`incomplete`) and nothing has been decided since: a search for the same
      inputs, finding the same entries, goes on from there;
    - `active`: the fold the protein draws now, where one was published since
      0012: its source, the record and regions it was made from, and its rule.
    """
    said: dict = {"none": None, "progress": None, "active": None}
    none = conn.execute(
        """
        select id, policy, inputs, record_sha256, decided_at, (search ->> 'entities')::int
        from fold_selection
        where slug = %s and outcome = 'none' and expires_at > now()
        order by decided_at desc, id desc limit 1
        """, (slug,)).fetchone()
    if none is not None:
        said["none"] = {"id": none[0], "policy": none[1], "inputs": none[2],
                        "record_sha256": none[3], "decided_at": none[4].isoformat(),
                        "entities": none[5]}
    latest = conn.execute(
        """
        select id, outcome, policy, inputs, search, refused from fold_selection
        where slug = %s order by decided_at desc, id desc limit 1
        """, (slug,)).fetchone()
    if latest is not None and latest[1] == "incomplete":
        said["progress"] = {"id": latest[0], "policy": latest[2], "inputs": latest[3],
                            "search": (latest[4] or {}).get("sha256"),
                            "refused": list(latest[5] or [])}
    active = conn.execute(
        """
        select id, source, record_sha256, regions, policy from fold_generation
        where slug = %s and active
        """, (slug,)).fetchone()
    if active is not None:
        said["active"] = dict(zip(("id", "source", "record_sha256", "regions", "policy"), active))
    return said


def pinned(prior: dict, record_sha256: str, regions: str, policy: int) -> bool:
    """Whether the experimental fold a protein draws stands for its inputs: an
    entry chosen for the same stored record, regions and disulfides, under the
    same rule. Such a fold is kept until a search is asked to start afresh."""
    active = prior.get("active")
    return (active is not None and active["source"] == "pdb"
            and (active["record_sha256"], active["regions"], active["policy"])
            == (record_sha256, regions, policy))


def _track(conn, slug: str, kind: str) -> Optional[dict]:
    row = conn.execute(
        f"select {', '.join(_TRACK_COLUMNS)} from protein_track "
        f"where slug = %s and kind = %s for update", (slug, kind)).fetchone()
    return dict(zip(_TRACK_COLUMNS, row)) if row else None


def _put_track(conn, row: dict) -> None:
    conn.execute(_PUT_TRACK, {**{key: row.get(key) for key in _TRACK_COLUMNS},
                              "provenance": jsonb(row.get("provenance") or {})})


def _clear_folding(conn, slug: str) -> None:
    """A protein's folding track says nothing is coming: a fold drawn still,
    or none at all, has no animation."""
    conn.execute(
        """
        update protein_track
        set state = 'absent', reason = null, bucket = null, object_path = null, bytes = null,
            sha256 = null, content_encoding = null, provenance = '{}'::jsonb, updated_at = now()
        where slug = %s and kind = 'folding' and state <> 'absent'
        """, (slug,))


def _owned(conn, job: dict) -> None:
    """The bake is still running at this claim's attempt, and held for this
    transaction; otherwise `Superseded`: the reaper put it back, or another
    worker took it."""
    owned = conn.execute(
        "select id from bake_job where id = %s and state = 'running' and attempts = %s "
        "for update", (job["id"], job["attempts"])).fetchone()
    if owned is None:
        raise Superseded("The bake is no longer this worker's: it was put back on the queue, "
                         "or taken again.")


def _record_selection(conn, slug: str, job_id: int, said: Optional[dict],
                      record_sha256: str) -> Optional[int]:
    """What a search decided, kept: a new row, or the earlier one it relied on."""
    if not said:
        return None
    if said.get("cached"):
        return said["cached"]
    chosen = said.get("chosen") or {}
    return conn.execute(
        f"""
        insert into fold_selection (slug, job_id, policy, inputs, record_sha256, outcome, search,
                                    candidate, entry, refused, detail, expires_at)
        values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
                case when %s = 'none' then now() + interval '{NONE_EXPIRES}' end)
        returning id
        """,
        (slug, job_id, said["policy"], said["inputs"], record_sha256, said["outcome"],
         jsonb(said.get("search")), said.get("candidate"), chosen.get("entry"),
         jsonb(said.get("refused") or []), jsonb(said), said["outcome"]),
    ).fetchone()[0]


def _snapshot(conn, slug: str, structure: Optional[dict]) -> None:
    """The fold a protein drew before any was published here, kept as a
    generation of its own, so that it can be put back: an AlphaFold model, or
    the sentence why there was none. Only once, and only a settled one."""
    if conn.execute("select 1 from fold_generation where slug = %s limit 1",
                    (slug,)).fetchone() is not None:
        return
    track = _track(conn, slug, "structure")
    if track is None or track["state"] not in ("ready", "refused"):
        return
    folding = _track(conn, slug, "folding")
    source = "none" if track["state"] == "refused" else (
        "pdb" if folding is not None and folding["state"] == "ready" else "alphafold")
    conn.execute(
        """
        insert into fold_generation (slug, source, state, structure, structure_track,
                                     folding_track, active)
        values (%s, %s, %s, %s, %s, %s, false)
        """,
        (slug, source, track["state"], jsonb(structure), jsonb(_plain(track)),
         jsonb(_plain(folding) if folding is not None and folding["state"] == "ready" else None)))


def _plain(row: Optional[dict]) -> Optional[dict]:
    """A track row as a generation keeps it: JSON, its numbers as numbers."""
    if row is None:
        return None
    return {key: row.get(key) for key in _TRACK_COLUMNS}


def publish_fold(conn, job: dict, fold: dict) -> Optional[int]:
    """Publish a protein's fold whole, or nothing: the structure track and the
    row's words, the folding track, the generation that records them, the
    search that chose them, and the bake done, in one transaction. Their
    objects are stored before this, under names of their own digests, and
    nothing stored is ever deleted, so a generation written over stays whole.

    `job` is the claim (`claim_bake`): the bake must still be running at that
    attempt (`_owned`). And the protein must be as it was read: one resolved on
    demand, with the same stored record (`fold["record_sha256"]`) and the same
    regions and disulfides (`fold["regions"]`, `regions_of`). Otherwise
    `Superseded` is raised and nothing is written.

    `fold` says what was made: its `source` ("pdb", "alphafold" or "none"),
    `structure_row` and `structure` (the row's `structure` column) for a
    model, `folding_row` for an experimental one, `reason` where there is no
    model (the fold page's sentence), `audit` (where the bake's audit is),
    `selection` (`fold_report.selection`, or `{"cached": id}`) and `policy`.

    A protein that draws a model keeps it where nothing new can be drawn
    ("none"): its bake ends, saying why on its job, and nothing a reader sees
    changes. Returns the generation published, or None where none was.
    """
    slug = job["slug"]
    source = fold["source"]
    with conn.transaction():
        _owned(conn, job)
        current = conn.execute(
            f"""
            select p.structure, r.sha256, {_REGIONS}
            from protein p
            left join protein_track r on r.slug = p.slug and r.kind = 'record'
                                     and r.state = 'ready'
            where p.slug = %s and p.catalog_order is null and p.resolver_version > 0
            for update of p
            """, (slug,)).fetchone()
        if current is None:
            raise RuntimeError(f"{slug} is a curated protein; the resolver does not write over it.")
        structure, record_sha256, regions = current
        if (record_sha256, regions) != (fold["record_sha256"], fold["regions"]):
            raise Superseded(f"{slug}'s record or regions changed while its fold was made.")

        _snapshot(conn, slug, structure)
        selection_id = _record_selection(conn, slug, job["id"], fold.get("selection"),
                                         record_sha256)
        if source == "none":
            reason = fold["reason"]
            track = _track(conn, slug, "structure")
            if track is not None and track["state"] == "ready":
                # It keeps the model it has.
                conn.execute("update bake_job set state = 'failed', error = %s, "
                             "finished_at = now() where id = %s", (reason, job["id"]))
                return None
            structure_row = {"slug": slug, "kind": "structure", "state": "refused",
                             "reason": reason, "format": "fsceneb", "provenance": {}}
            _put_track(conn, structure_row)
            _clear_folding(conn, slug)
            folding_row = None
            conn.execute("update bake_job set state = 'failed', error = %s, finished_at = now() "
                         "where id = %s", (reason, job["id"]))
        else:
            structure = fold["structure"]
            structure_row = {**fold["structure_row"], "state": "ready", "reason": None}
            folding_row = (None if source != "pdb"
                           else {**fold["folding_row"], "state": "ready", "reason": None})
            written = conn.execute(
                "update protein set structure = %s where slug = %s and catalog_order is null "
                "returning slug", (jsonb(structure), slug)).fetchone()
            if written is None:
                raise RuntimeError(
                    f"{slug} is a curated protein; the resolver does not write over it.")
            _put_track(conn, structure_row)
            if folding_row is not None:
                _put_track(conn, folding_row)
            else:
                _clear_folding(conn, slug)
            conn.execute("update bake_job set state = 'done', error = null, finished_at = now() "
                         "where id = %s", (job["id"],))
        conn.execute("update fold_generation set active = false where slug = %s and active",
                     (slug,))
        return conn.execute(
            """
            insert into fold_generation (slug, source, state, job_id, selection_id, structure,
                                         structure_track, folding_track, audit, record_sha256,
                                         regions, policy, active)
            values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, true)
            returning id
            """,
            (slug, source, structure_row["state"], job["id"], selection_id,
             jsonb(structure if source != "none" else None), jsonb(_plain(structure_row)),
             jsonb(_plain(folding_row)), jsonb(fold.get("audit")), record_sha256, regions,
             fold.get("policy")),
        ).fetchone()[0]


def save_progress(conn, job: dict, selection: dict, record_sha256: str) -> None:
    """A search cut short at its bounds: how far it got kept
    (`fold_selection`, `incomplete`), and its bake put back on the queue
    without spending one of its tries. The next attempt goes on from there."""
    with conn.transaction():
        _owned(conn, job)
        _record_selection(conn, job["slug"], job["id"], selection, record_sha256)
        conn.execute(
            """
            update bake_job
            set state = 'queued', started_at = null, attempts = greatest(attempts - 1, 0),
                error = 'The search for an experimental structure goes on next cycle.'
            where id = %s
            """, (job["id"],))


def end_fold_bake(conn, job: dict, *, error: str, requeue: bool, refund: bool = False,
                  reason: Optional[str] = None) -> bool:
    """End a structure bake that published nothing, if it is still this
    worker's: back on the queue (`requeue`, with its try given back where
    `refund`), or failed, its track refused with `reason` where it has no model
    yet. Returns whether it was this worker's to end."""
    with conn.transaction():
        try:
            _owned(conn, job)
        except Superseded:
            return False
        if requeue:
            conn.execute(
                "update bake_job set state = 'queued', error = %s, started_at = null, "
                "attempts = greatest(attempts - %s, 0) where id = %s",
                (error, 1 if refund else 0, job["id"]))
        else:
            if reason is not None:
                refuse_track(conn, job["slug"], "structure", reason)
            conn.execute("update bake_job set state = 'failed', error = %s, finished_at = now() "
                         "where id = %s", (error, job["id"]))
    return True


def keep_fold(conn, job: dict) -> bool:
    """A structure bake whose protein's experimental fold stands (`pinned`):
    done, with nothing written. Returns whether it was this worker's."""
    with conn.transaction():
        try:
            _owned(conn, job)
        except Superseded:
            return False
        conn.execute(
            "update bake_job set state = 'done', finished_at = now(), "
            "error = 'Its experimental fold stands: the same record, regions and rule.' "
            "where id = %s", (job["id"],))
    return True


def queue_structure(conn, slug: str, *, refresh: bool = False, policy: Optional[int] = None,
                    asker: Optional[str] = None) -> str:
    """A protein's fold searched for again, as the worker searches for a new
    one's: its structure bake queued. Returns what was done:

    - `queued`: a bake is queued, and a track with no model yet says one is
      coming (`pending`). A model it draws is kept until another is published.
    - `waiting`: one is queued or running already.
    - `kept`: its experimental fold stands for the same record, regions and
      rule (`pinned`); `refresh` searches all the same.

    Only a protein resolved on demand; one of the twenty is baked by hand.
    """
    with conn.transaction():
        found = conn.execute(
            f"""
            select r.sha256, {_REGIONS} from protein p
            left join protein_track r on r.slug = p.slug and r.kind = 'record'
                                     and r.state = 'ready'
            where p.slug = %s and p.catalog_order is null and p.resolver_version > 0
            for update of p
            """, (slug,)).fetchone()
        if found is None:
            raise RuntimeError(f"{slug} is not a protein resolved on demand; its fold is "
                               f"baked by hand.")
        if found[0] is None:
            raise RuntimeError(f"{slug} has no ready record to draw a fold for.")
        if conn.execute(
                "select 1 from bake_job where slug = %s and kind = 'structure' "
                "and state in ('queued', 'running')", (slug,)).fetchone() is not None:
            return "waiting"
        if not refresh and policy is not None and pinned(fold_prior(conn, slug), found[0],
                                                         found[1], policy):
            return "kept"
        conn.execute(
            """
            insert into protein_track (slug, kind, state, reason, format, provenance)
            values (%s, 'structure', 'pending', null, 'fsceneb', '{}'::jsonb)
            on conflict (slug, kind) do update
            set state = 'pending', reason = null, updated_at = now()
            where protein_track.state <> 'ready'
            """, (slug,))
        conn.execute(
            "insert into bake_job (slug, kind, asker, refresh) values (%s, 'structure', %s, %s)",
            (slug, asker, refresh))
        return "queued"


def generations(conn, slug: str) -> list:
    """A protein's published folds, newest first."""
    columns = ("id", "source", "state", "active", "created_at", "job_id", "selection_id",
               "structure", "structure_track", "folding_track", "policy")
    rows = conn.execute(
        f"select {', '.join(columns)} from fold_generation where slug = %s order by id desc",
        (slug,)).fetchall()
    return [dict(zip(columns, row)) for row in rows]


def roll_back_fold(conn, slug: str, to: Optional[int] = None) -> dict:
    """Put back a fold a protein drew before, whole, in one transaction: the
    row's words, the structure track and the folding track as a generation
    recorded them (no folding track: none, read as absent), and that
    generation the active one. `to` names it; by default it is the one
    before the fold drawn now.

    Refused while a structure bake for it is queued or running, which would
    publish over what is put back. Returns which generation was left, and
    which put back."""
    with conn.transaction():
        found = conn.execute(
            "select 1 from protein where slug = %s and catalog_order is null "
            "and resolver_version > 0 for update", (slug,)).fetchone()
        if found is None:
            raise RuntimeError(f"{slug} is not a protein resolved on demand.")
        if conn.execute(
                "select 1 from bake_job where slug = %s and kind = 'structure' "
                "and state in ('queued', 'running')", (slug,)).fetchone() is not None:
            raise RuntimeError(f"A structure bake for {slug} is queued or running: let it "
                               f"finish, or stop the worker, before a fold is put back.")
        every = generations(conn, slug)
        active = next((g for g in every if g["active"]), None)
        if to is None:
            before = [g for g in every if active is None or g["id"] < active["id"]]
            if not before:
                raise RuntimeError(f"{slug} has no earlier fold to put back.")
            chosen = before[0]
        else:
            chosen = next((g for g in every if g["id"] == to), None)
            if chosen is None:
                raise RuntimeError(f"{slug} has no generation {to}.")
        if active is not None and chosen["id"] == active["id"]:
            raise RuntimeError(f"Generation {to} is the fold {slug} draws now.")
        written = conn.execute(
            "update protein set structure = %s where slug = %s and catalog_order is null "
            "returning slug",
            (jsonb(chosen["structure"] if chosen["state"] == "ready" else None), slug)).fetchone()
        if written is None:
            raise RuntimeError(f"{slug} is a curated protein; the resolver does not write over it.")
        _put_track(conn, chosen["structure_track"])
        if chosen["folding_track"] is not None:
            _put_track(conn, chosen["folding_track"])
        else:
            _clear_folding(conn, slug)
        # Two statements: the index allows one active a protein at every row.
        conn.execute("update fold_generation set active = false where slug = %s and active",
                     (slug,))
        conn.execute("update fold_generation set active = true where id = %s", (chosen["id"],))
        return {"from": None if active is None else active["id"], "to": chosen["id"],
                "source": chosen["source"], "state": chosen["state"]}

