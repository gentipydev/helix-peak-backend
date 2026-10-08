"""Proteins the catalog does not list, asked for by a reader: Phase 6.

The service does not resolve anything itself. A request is a row in
`resolve_request`; the resolver on Modal (`pipeline/resolver/`) works it, and
writes the protein row and its record track when it is done. What this module
does is decide what an ask means, write the row when it is a new one, and wake
the resolver so the reader does not wait for its schedule.

What a protein can be, as `ResolveResponse.state` says it:

- ``ready``: there is a protein row for this gene, listed or resolved. ``slug``
  is where it opens, which for the twenty is the slug they shipped with.
- ``pending``: a request is queued or running.
- ``refused``: the resolver declined, and ``reason`` is why. Final, for the
  resolver version that said it; asking again changes nothing.
- ``failed``: it broke more often than it is retried. Asking again queues a new
  request.
- ``unavailable``: the index cannot build this protein, and ``reason`` is why
  (no MANE Select transcript, or one for another isoform).
- ``buildable``: nobody has asked yet, and nothing stands in the way. Only a
  read says this; an ask turns it into ``pending``.

A gene the index does not hold at all is not a state but a 404, and an
unreadable database is a 503, never an answer: the catalog's rule.

The one thing that costs money is a request, so a day's are capped
(``RESOLVES_PER_DAY``). Two readers asking for one gene make one request, which
the database holds to (`resolve_request_inflight_idx`).

Beside the state, ``build`` says how far a build has got (`BuildReport`), so a
reader waiting minutes on ESM-2 sees which step it is at, the residues scored
and the time left. It is read from the request and the constraint bake, whose
progress the worker writes as the scorer prints it (`store.note_progress`).

The reader who asked may stop a build (`stop`): there are no accounts, so who
asked is the asking install's own random id, kept on the request and on the
bake it queues (0011). Stopped while queued, nothing was built and asking
again starts over. Stopped while ESM-2 was scoring, the protein stays built
without its track (``absent``), and asking again queues the scoring afresh.
"""

import logging
import urllib.request
from typing import Optional, Tuple

from .catalog import CatalogUnavailable, _require_pool
from .config import settings

logger = logging.getLogger(__name__)

HUMAN = 9606


class DailyCapReached(Exception):
    """Today's requests are used up."""


class NotYours(Exception):
    """Only the install that asked for a build may stop it."""


class NothingToStop(Exception):
    """No build a stop could end is under way; the message says what there is."""


_PROTEIN = """
select slug from protein where lower(gene) = lower(%s) and taxon_id = 9606
"""

# A buildable entry first: a gene can name several entries, at most one of
# them buildable (`protein_index_buildable_gene`).
_INDEX = """
select uniprot, gene, buildable, unavailable_reason from protein_index
where lower(gene) = lower(%s) and gene <> ''
order by buildable desc, uniprot
limit 1
"""

_LATEST = """
select state, reason from resolve_request where lower(gene) = lower(%s)
order by requested_at desc, id desc
limit 1
"""

_TODAY = """
select count(*) from resolve_request where requested_at > now() - interval '1 day'
"""

_QUEUE = """
insert into resolve_request (gene, uniprot, slug, asker) values (%s, %s, %s, %s)
on conflict (lower(gene)) where state in ('queued', 'running') do nothing
returning id
"""

# The latest request for a gene, the seconds since it was asked for (to now, or
# to when it ended), and the proteins the worker takes before it: requests
# asked for earlier and not finished, and every constraint bake still to run,
# which belongs to a request that was resolved first.
_REQUEST_BUILD = """
select r.state,
       extract(epoch from coalesce(r.finished_at, now()) - r.requested_at)::float8,
       (select count(*) from resolve_request o
        where o.state in ('queued', 'running')
          and (o.requested_at, o.id) < (r.requested_at, r.id))
       + (select count(*) from bake_job j join protein p on p.slug = j.slug
          where j.kind = 'constraint' and j.state in ('queued', 'running')
            and p.resolver_version > 0)
from resolve_request r
where lower(r.gene) = lower(%s)
order by r.requested_at desc, r.id desc
limit 1
"""

# A protein resolved on demand and its latest constraint bake: where the bake
# is, how far the scorer has got in the run that holds it (progress written by
# an earlier run does not count, 0010), the bakes taken before a queued one,
# the seconds since the request that queued it was asked for, the transcript
# its record was built from, and what its track says. The request is the one
# finished in the transaction that queued the bake, so the two share `now()`;
# a bake queued by hand counts from itself.
_BAKE_BUILD = """
select j.state,
       case when j.progress_at >= j.started_at then j.progress_done end,
       case when j.progress_at >= j.started_at then j.progress_total end,
       case when j.progress_at >= j.started_at and j.state = 'running'
            then greatest(j.progress_left - extract(epoch from now() - j.progress_at), 0)::float8
       end,
       case when j.state = 'queued' then
           (select count(*) from bake_job o join protein q on q.slug = o.slug
            where o.kind = 'constraint' and q.resolver_version > 0 and o.id <> j.id
              and (o.state = 'running'
                   or (o.state = 'queued' and (o.requested_at, o.id) < (j.requested_at, j.id))))
       end,
       extract(epoch from coalesce(j.finished_at, now()) - coalesce(
           (select r.requested_at from resolve_request r
            where lower(r.gene) = lower(p.gene) and r.finished_at = j.requested_at
            order by r.id desc limit 1),
           j.requested_at))::float8,
       p.transcript_id,
       t.state,
       t.reason,
       case when j.state = 'running'
            then extract(epoch from now() - j.started_at)::float8
       end
from protein p
join bake_job j on j.slug = p.slug and j.kind = 'constraint'
left join protein_track t on t.slug = p.slug and t.kind = 'constraint'
where p.slug = %s and p.catalog_order is null
order by j.requested_at desc, j.id desc
limit 1
"""


def _said(state: str, slug: Optional[str], reason: Optional[str] = None) -> dict:
    return {"slug": slug, "state": state, "reason": reason}


# The latest request for a gene, and the install that asked for it.
_LATEST_ASK = """
select id, state, asker from resolve_request where lower(gene) = lower(%s)
order by requested_at desc, id desc
limit 1
"""

# A request not yet taken, withdrawn; one a worker holds is left to it.
_STOP_REQUEST = """
update resolve_request
set state = 'stopped', reason = 'Stopped before it was built.', finished_at = now()
where id = %s and state = 'queued'
returning extract(epoch from finished_at - requested_at)::float8
"""

# A protein's latest constraint bake, and the install that asked for it.
_LATEST_BAKE = """
select id, state, asker from bake_job where slug = %s and kind = 'constraint'
order by requested_at desc, id desc
limit 1
"""

# The bake ended where it is: the worker hears it at its next line and ends the
# scorer, and the track says that nothing is coming.
_STOP_BAKE = """
update bake_job
set state = 'stopped', error = 'Stopped by the reader who asked.', finished_at = now()
where id = %s and state in ('queued', 'running')
returning id
"""
_STOPPED_TRACK = """
update protein_track set state = 'absent', reason = null, updated_at = now()
where slug = %s and kind = 'constraint' and state = 'pending'
"""

# A protein resolved on demand whose scoring was stopped, and whose track still
# has nothing: its scoring queued afresh, for the install asking now.
_RESUME = """
insert into bake_job (slug, kind, asker)
select p.slug, 'constraint', %(asker)s
from protein p
join protein_track t on t.slug = p.slug and t.kind = 'constraint' and t.state = 'absent'
where p.slug = %(slug)s and p.catalog_order is null
  and (select j.state from bake_job j where j.slug = p.slug and j.kind = 'constraint'
       order by j.requested_at desc, j.id desc limit 1) = 'stopped'
on conflict (slug, kind) where state in ('queued', 'running') do nothing
returning id
"""
_RESUMED_TRACK = """
update protein_track set state = 'pending', reason = null, updated_at = now()
where slug = %s and kind = 'constraint' and state = 'absent'
"""


def _report(step: str, elapsed: float, **said) -> dict:
    return {"step": step, "ahead": None, "scored": None, "residues": None, "left": None,
            "elapsed": elapsed, "transcript": None, "reason": None, "stopped": False,
            "ran": None, **said}


def _baked(row: tuple) -> dict:
    """A resolved protein's build, from its latest constraint bake (`_BAKE_BUILD`)."""
    job, done, total, left, ahead, elapsed, transcript, track, reason, ran = row
    if job in ("queued", "running"):
        if job == "queued":
            return _report("scoring", elapsed, transcript=transcript, ahead=ahead)
        # Every residue scored, and the bake not over: the gate and the upload.
        step = "check" if done is not None and done == total else "scoring"
        return _report(step, elapsed, transcript=transcript, scored=done, residues=total,
                       left=left, ran=ran)
    if job == "stopped":
        said = ("Stopped at {:,} of {:,} residues.".format(done, total)
                if done is not None and total else "Stopped before a residue was scored.")
        return _report("scoring", elapsed, transcript=transcript, scored=done, residues=total,
                       reason=said, stopped=True)
    if track == "ready":
        return _report("done", elapsed, transcript=transcript, scored=done, residues=total)
    # Refused: by the alignment gate once every residue was scored, or by the
    # scorer before then. Where the run left no progress, it is the scoring.
    step = "check" if done is not None and done == total else "scoring"
    return _report(step, elapsed, transcript=transcript, scored=done, residues=total,
                   reason=reason or "The ESM-2 track was not made.")


def _build(conn, gene: str, said: dict) -> Optional[dict]:
    """How far the build of what `said` names has got. None where there is none:
    one of the twenty, or a protein nobody has asked for."""
    state = said["state"]
    if state == "ready":
        row = conn.execute(_BAKE_BUILD, (said["slug"],)).fetchone()
        return None if row is None else _baked(row)
    if state not in ("pending", "refused", "failed"):
        return None
    row = conn.execute(_REQUEST_BUILD, (gene,)).fetchone()
    if row is None:
        return None
    request_state, elapsed, ahead = row
    if request_state == "queued":
        return _report("queued", elapsed, ahead=ahead)
    # Running, or ended here: the response's own state and reason say how.
    return _report("record", elapsed)


def _current(conn, gene: str) -> Tuple[Optional[dict], Optional[tuple]]:
    """What this gene is now, and its index entry when a request could be made.

    The first is None for a gene the index does not hold; the second is set
    only where nothing stands in the way of a new request.
    """
    found = conn.execute(_PROTEIN, (gene,)).fetchone()
    if found is not None:
        return _said("ready", found[0]), None
    entry = conn.execute(_INDEX, (gene,)).fetchone()
    if entry is None:
        return None, None
    uniprot, symbol, buildable, unavailable = entry
    slug = symbol.lower()
    if not buildable:
        return _said("unavailable", None, unavailable), None
    latest = conn.execute(_LATEST, (gene,)).fetchone()
    if latest is not None:
        state, reason = latest
        if state in ("queued", "running"):
            return _said("pending", slug), None
        if state == "refused":
            return _said("refused", None, reason), None
        if state == "failed":
            return _said("failed", slug, reason), (uniprot, symbol, slug)
    return None, (uniprot, symbol, slug)


def current(gene: str, taxon: int = HUMAN) -> Optional[dict]:
    """What a gene's protein is now. None where the index does not hold it."""
    if taxon != HUMAN:
        return None
    pool = _require_pool()
    try:
        with pool.connection() as conn:
            said, askable = _current(conn, gene)
            if said is not None:
                said["build"] = _build(conn, gene, said)
    except Exception as exc:
        logger.exception("Resolve state read failed for %r", gene)
        raise CatalogUnavailable("The catalog could not be read.") from exc
    if said is None and askable is not None:
        # Buildable, and never asked for.
        return _said("buildable", askable[2])
    return said


def _resume(conn, slug: str, asker: Optional[str]) -> bool:
    """Queue afresh the scoring of a protein whose scoring was stopped.

    Whether it did: only a protein resolved on demand, whose latest bake was
    stopped and whose track still has nothing, is scored again. It is no new
    request, so it does not count towards the day's.
    """
    if conn.execute(_RESUME, {"slug": slug, "asker": asker}).fetchone() is None:
        return False
    conn.execute(_RESUMED_TRACK, (slug,))
    return True


def request(gene: str, taxon: int = HUMAN,
            asker: Optional[str] = None) -> Tuple[Optional[dict], bool]:
    """Ask for a gene's protein. Returns what it is now, and whether this ask
    queued new work, a request or a stopped protein's scoring (and so should
    wake the resolver). `asker` is the asking install's id: only it may stop
    what this ask queues.

    Raises DailyCapReached where a new request would be one too many today.
    """
    if taxon != HUMAN:
        return None, False
    pool = _require_pool()
    try:
        with pool.connection() as conn:
            said, askable = _current(conn, gene)
            if askable is None:
                resumed = False
                if said is not None:
                    if said["state"] == "ready":
                        resumed = _resume(conn, said["slug"], asker)
                    said["build"] = _build(conn, gene, said)
                return said, resumed
            if conn.execute(_TODAY).fetchone()[0] >= settings.resolves_per_day:
                raise DailyCapReached()
            uniprot, symbol, slug = askable
            queued = conn.execute(_QUEUE, (symbol, uniprot, slug, asker)).fetchone()
            # Not queued means another reader's request got there in the same
            # moment: it is pending either way, and only the one that queued
            # it wakes anyone.
            said = _said("pending", slug)
            said["build"] = _build(conn, gene, said)
    except DailyCapReached:
        raise
    except Exception as exc:
        logger.exception("Resolve request failed for %r", gene)
        raise CatalogUnavailable("The request could not be recorded.") from exc
    return said, queued is not None


def stop(gene: str, asker: str, taxon: int = HUMAN) -> Optional[dict]:
    """Stop the build of a gene's protein, for the install that asked for it.

    Returns what it is now, ``stopped``, with the build as it ended; None where
    the index does not hold the gene. Raises NotYours where another install
    asked (or none that could be named, before 0011), and NothingToStop where
    no build a stop could end is under way: none at all, its gene record being
    written (seconds), or its track made a moment ago.
    """
    if taxon != HUMAN:
        return None
    pool = _require_pool()
    try:
        with pool.connection() as conn:
            said, askable = _current(conn, gene)
            if said is None and askable is None:
                return None
            latest = conn.execute(_LATEST_ASK, (gene,)).fetchone()
            if latest is not None and latest[1] in ("queued", "running"):
                request_id, state, owner = latest
                if owner is None or owner != asker:
                    raise NotYours()
                ended = None
                if state == "queued":
                    ended = conn.execute(_STOP_REQUEST, (request_id,)).fetchone()
                if ended is None:
                    raise NothingToStop("Its gene record is being written. Stop it in a moment.")
                build = _report("queued", ended[0], stopped=True,
                                reason="Stopped before it was built.")
                # Nothing was built, so nothing opens.
                return {**_said("stopped", None), "build": build}
            if said is not None and said["state"] == "ready":
                slug = said["slug"]
                bake = conn.execute(_LATEST_BAKE, (slug,)).fetchone()
                if bake is not None and bake[1] in ("queued", "running"):
                    job_id, _, owner = bake
                    if owner is None or owner != asker:
                        raise NotYours()
                    if conn.execute(_STOP_BAKE, (job_id,)).fetchone() is None:
                        raise NothingToStop("Its ESM-2 track has just been made.")
                    conn.execute(_STOPPED_TRACK, (slug,))
                    row = conn.execute(_BAKE_BUILD, (slug,)).fetchone()
                    return {**_said("stopped", slug), "build": _baked(row)}
            raise NothingToStop("Nothing is building it.")
    except (NotYours, NothingToStop):
        raise
    except Exception as exc:
        logger.exception("Stopping the build failed for %r", gene)
        raise CatalogUnavailable("The stop could not be recorded.") from exc


def wake() -> None:
    """Tell the resolver on Modal a request is waiting. Best effort.

    Unconfigured, or unanswered, this does nothing a reader would notice except
    wait: the resolver's own schedule takes every queued request within five
    minutes. So a failure is logged and never raised.
    """
    url = settings.modal_wake_url
    if not url:
        return
    headers = {}
    if settings.modal_key and settings.modal_secret:
        headers = {"Modal-Key": settings.modal_key, "Modal-Secret": settings.modal_secret}
    call = urllib.request.Request(url, data=b"", method="POST", headers=headers)
    try:
        with urllib.request.urlopen(call, timeout=10):
            return
    except Exception:
        logger.warning("Could not wake the resolver at %s; its schedule will take the request",
                       url, exc_info=True)
