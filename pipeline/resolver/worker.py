"""The resolver's work loop: take a request, resolve it, bake what it queued.

Two steps, on two kinds of machine. `resolve_next` takes one `resolve_request`,
resolves its protein (`resolve.py`), stores the record track and writes the rows
it needs; that is CPU and network, and takes about as long as NCBI does.
`score_next` takes one queued constraint bake and runs ESM-2 over the protein
the record holds; that wants a GPU. `modal_app.py` runs each where it belongs,
and either can be run by hand on any machine with the same credentials.

What ends a request or a bake, and how the reader hears of it:

- A refusal -- the resolver's `Refused`, the uploader's gate declining the
  record it built (`Declined`), or the scorer's own `ValueError` (its
  alignment gate, its budget, a sequence that is not the record's) -- is final
  for this resolver version. The request or the track says `refused`, with
  the sentence why.
- Anything else (NCBI or UniProt not answering, a storage hiccup) is tried
  again, up to `store.MAX_ATTEMPTS` times, and then fails, saying so.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import re
import sys
from typing import Callable, Optional, Tuple

from pipeline import uniprot, upload_tracks
from pipeline.resolver import store
from pipeline.resolver.resolve import RESOLVER_VERSION, Refused, resolve, target_of
from pipeline.resolver.scoring import score_with_esm
from pipeline.targets import Target

# How far a scoring has got: residues scored, of how many, and the seconds left.
Report = Callable[[int, int, float], None]

# The constraint track's bytes, from a target and its record. It tells `Report`
# how far it has got as it goes, for a reader watching the protein build.
Score = Callable[[Target, bytes, Report], bytes]

# What the scorer prints every 25 residues (`score_protein.py`).
_SCORED = re.compile(r"Scored (\d+)/(\d+) \(([\d.]+)s, [\d.]+ min left\)")


def progress_of(text: str) -> Optional[Tuple[int, int, float]]:
    """The last progress line the scorer printed in `text`, as `Report` takes it.

    The seconds left are the scorer's own rate over what is still to score,
    worked from its elapsed seconds rather than read from its minutes, which
    it rounds to a tenth.
    """
    found = _SCORED.findall(text)
    if not found:
        return None
    done, total, elapsed = int(found[-1][0]), int(found[-1][1]), float(found[-1][2])
    if done <= 0 or total <= 0 or done > total:
        return None
    return done, total, elapsed / done * (total - done)


class _ProgressTee(io.TextIOBase):
    """Standard output as it was, with each progress line also told to `report`."""

    def __init__(self, out, report: Report):
        self._out = out
        self._report = report
        self._line = ""

    def write(self, text: str) -> int:
        self._out.write(text)
        self._line += text
        *whole, self._line = self._line.split("\n")
        for line in whole:
            found = progress_of(line)
            if found is not None:
                self._report(*found)
        return len(text)

    def flush(self) -> None:
        self._out.flush()


def score_in_process(target: Target, record: bytes, report: Report) -> bytes:
    """`scoring.score_with_esm` in this process, as Modal runs it, with the
    progress the scorer prints passed to `report` as it prints it."""
    with contextlib.redirect_stdout(_ProgressTee(sys.stdout, report)):
        return score_with_esm(target, record)


class Declined(ValueError):
    """The uploader's gate turned a payload down. It would turn the same bytes
    down again, so this is never worth another try."""


def _said(error: BaseException) -> str:
    """An exception as one line a request or a track can carry."""
    text = " ".join(str(error).split()) or type(error).__name__
    return text if len(text) <= 300 else text[:297] + "..."


def stage(storage, kind: str, target: Target, payload: bytes) -> dict:
    """Validate a payload, store its bytes, and return the track row to write.

    The uploader's own gate (`upload_tracks.validate`) decides what is
    storable and what its provenance says, and the object is named as the
    uploader names it, `<kind>/<slug>.<sha12>.json`, so a track resolved on
    demand cannot be told from one uploaded by hand. A payload the gate turns
    down raises `Declined`, and nothing is stored.
    """
    try:
        provenance = upload_tracks.validate(kind, target, payload)
    except ValueError as exc:
        raise Declined(str(exc)) from exc
    digest = hashlib.sha256(payload).hexdigest()
    path = f"{kind}/{target.slug}.{digest[:12]}.json"
    storage.put(store.TRACKS_BUCKET, path, payload, "application/json")
    return {
        "slug": target.slug, "kind": kind, "bucket": store.TRACKS_BUCKET,
        "object_path": path, "bytes": len(payload), "sha256": digest,
        "format": "json", "provenance": provenance,
    }


def resolve_next(conn, storage, *, fetch_entry=uniprot.fetch_entry) -> Optional[dict]:
    """Resolve the oldest queued request. None when there is none."""
    request = store.claim_request(conn)
    if request is None:
        return None
    outcome = {"id": request["id"], "gene": request["gene"]}
    try:
        existing = store.protein_for_gene(conn, request["gene"])
        if existing is not None:
            # Listed, or resolved by a request that finished first.
            store.finish_request(conn, request["id"], "done", slug=existing)
            return {**outcome, "state": "done", "slug": existing}
        row = store.index_row(conn, request["uniprot"], request["gene"])
        if row is None:
            raise Refused(f"The protein index no longer holds {request['uniprot']} "
                          f"made by {request['gene']}.")
        resolution = resolve(row, fetch_entry(row.uniprot), mane_release=store.mane_release(conn))
        try:
            record = stage(storage, "record", resolution.target, resolution.record)
        except Declined as exc:
            # This resolver builds the same record for this gene every time,
            # so the gate's verdict on it is a refusal, not a failure to retry.
            raise Refused(str(exc)) from exc
        store.write_resolution(conn, request["id"], resolution, record, RESOLVER_VERSION)
        return {**outcome, "state": "done", "slug": resolution.target.slug}
    except Refused as exc:
        store.finish_request(conn, request["id"], "refused", reason=_said(exc),
                             resolver_version=RESOLVER_VERSION)
        return {**outcome, "state": "refused", "reason": _said(exc)}
    except Exception as exc:  # noqa: BLE001 -- anything else is worth another try
        if request["attempts"] < store.MAX_ATTEMPTS:
            store.requeue_request(conn, request["id"])
            return {**outcome, "state": "queued", "reason": _said(exc)}
        store.finish_request(conn, request["id"], "failed", reason=_said(exc))
        return {**outcome, "state": "failed", "reason": _said(exc)}


def _refusal(error: ValueError) -> str:
    """The scorer's refusals, as the sentence a reader is shown.

    A failed alignment gate is told as what was measured, never as scores
    misfiled. The gate cannot tell a track shifted by a position from a
    protein ESM-2 knows too little about to prefer its own residues, and among
    proteins resolved on demand it is the second that turns up: micropeptides
    and orphan proteins fall short, while a familiar one of 25 residues clears
    it (`README.md` has the measurement).
    """
    found = re.search(r"Alignment gate FAILED: ([\d.]+%) before, ([\d.]+%) after", str(error))
    if found:
        return (f"ESM-2 prefers the residue that is there to the one before it "
                f"{found.group(1)} of the time, and to the one after it {found.group(2)}. "
                f"More than half is needed to check that its scores sit on their own "
                f"residues, and a protein the model knows little about falls short. "
                f"The scores are not drawn.")
    return _said(error)


def _reporter(conn, job_id: int) -> Report:
    """`Report` for one job: its progress written where a reader asking after
    the protein reads it. A write that fails costs the bake nothing; the next
    line tries again, and the track is what the bake is for."""
    def report(done: int, total: int, left: float) -> None:
        try:
            store.note_progress(conn, job_id, done, total, left)
        except Exception:  # noqa: BLE001 -- see the docstring
            pass
    return report


def score_next(conn, storage, *, score: Score = score_in_process) -> Optional[dict]:
    """Bake the oldest queued constraint track. None when there is none."""
    job = store.claim_bake(conn, "constraint")
    if job is None:
        return None
    outcome = {"id": job["id"], "slug": job["slug"]}
    try:
        protein = store.protein(conn, job["slug"])
        record_row = store.ready_track(conn, job["slug"], "record")
        if protein is None or record_row is None:
            raise Refused("The protein has no ready record to score.")
        record = storage.get(record_row["bucket"], record_row["object_path"])
        if hashlib.sha256(record).hexdigest() != record_row["sha256"]:
            raise RuntimeError("The stored record does not match the sha256 its row carries.")
        target = target_of(protein)
        report = _reporter(conn, job["id"])
        track = stage(storage, "constraint", target, score(target, record, report))
        store.finish_bake(conn, job["id"], track)
        return {**outcome, "state": "ready"}
    except (Refused, ValueError) as exc:
        reason = _refusal(exc) if isinstance(exc, ValueError) else _said(exc)
        store.refuse_bake(conn, job["id"], job["slug"], "constraint", reason)
        return {**outcome, "state": "refused", "reason": reason}
    except Exception as exc:  # noqa: BLE001
        if job["attempts"] < store.MAX_ATTEMPTS:
            store.requeue_bake(conn, job["id"], _said(exc))
            return {**outcome, "state": "queued", "reason": _said(exc)}
        reason = f"Scoring failed {store.MAX_ATTEMPTS} times: {_said(exc)}"
        store.refuse_bake(conn, job["id"], job["slug"], "constraint", reason)
        return {**outcome, "state": "refused", "reason": reason}


def sweep(conn, storage, *, limit: int = 10, fetch_entry=uniprot.fetch_entry) -> dict:
    """Reap what dead workers held, then resolve up to `limit` requests.

    Returns what it did, and how many constraint bakes are waiting, so the
    caller knows whether to start a GPU.
    """
    store.reap(conn)
    resolved = []
    for _ in range(limit):
        outcome = resolve_next(conn, storage, fetch_entry=fetch_entry)
        if outcome is None:
            break
        resolved.append(outcome)
    return {"resolved": resolved, "constraint_queued": store.queued_bakes(conn, "constraint")}


def score_all(conn, storage, *, limit: int = 20, score: Score = score_in_process) -> list[dict]:
    """Bake queued constraint tracks until none are left, or `limit` are done."""
    done = []
    for _ in range(limit):
        outcome = score_next(conn, storage, score=score)
        if outcome is None:
            break
        done.append(outcome)
    return done
