"""The resolver's work loop: take a request, resolve it, bake what it queued.

Two steps, on two kinds of machine. `resolve_next` takes one `resolve_request`,
resolves its protein (`resolve.py`), stores the record track and writes the rows
it needs; that is CPU and network, and takes about as long as NCBI does.
`score_next` takes one queued constraint bake and runs ESM-2 over the protein
the record holds; that wants a GPU. `modal_app.py` runs each where it belongs,
and either can be run by hand on any machine with the same credentials.

A third, where a worker asks for it (`structures`): `structure_next` takes one
queued structure bake and makes the fold page's model from AlphaFold DB's
(`structure/alphafold.py`). That wants PyMOL and the app's scene importer,
which the Mac's worker has and Modal's image does not, so only the Mac asks.

And the variant evidence, where a worker asks for it (`evidence`), with the
bakers that made the twenty's, unchanged: `impact_next` takes one queued AVI
bake (`impact/bake_impact.py`, AlphaGenome's scores and the map from the
record's bases to GRCh38), and `clinvar_next` one ClinVar bake
(`clinvar/bake_clinvar.py`), whose records are placed by that map, so it waits
for it (`store.WAITS_FOR`). The first wants the AlphaGenome client and its
skill's GENCODE lookup, which only the Mac has. ESM-2 waits for both, so a
build that opens once its scores are in has its evidence too.

What ends a request or a bake, and how the reader hears of it:

- A refusal -- the resolver's `Refused`, the uploader's gate declining the
  record it built (`Declined`), or the scorer's own `ValueError` (its
  alignment gate, its budget, a sequence that is not the record's) -- is final
  for this resolver version. The request or the track says `refused`, with
  the sentence why.
- Anything else (NCBI or UniProt not answering, a storage hiccup) is tried
  again, up to `store.MAX_ATTEMPTS` times, and then fails, saying so.
- A stop by the reader who asked (`Stopped`) ends a scoring where it is. The
  service has already written it on the rows, so nothing more is. A model is
  made in seconds, as a record is, and like a record is not stopped.
- A protein AlphaFold DB has no model of, or none sure enough to draw
  (`Unmodelled`), has its structure track `refused`, with the sentence why.
- One of the AVI bake's own gates saying no (`Unplaced`: the exons do not
  pair, or the record's bases are not GRCh38's where GENCODE puts the gene)
  refuses its AVI track, and the ClinVar track with it: there is no map to
  place ClinVar's records by. NCBI or the Atlas not answering (`Unfetched`,
  or any other error) is tried again.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import re
import shutil
import struct
import sys
from pathlib import Path
from typing import Callable, NamedTuple, Optional, Tuple

from pipeline import paths, uniprot, upload_tracks
from pipeline.resolver import store
from pipeline.resolver.resolve import (
    RESOLVER_VERSION, VARIANT_FLOOR, VARIANT_SHARE, Refused, resolve, target_of)
from pipeline.resolver.scoring import score_with_esm
from pipeline.targets import Target

# How far a scoring has got: residues scored, of how many, and the seconds left.
Report = Callable[[int, int, float], None]

# The constraint track's bytes, from a target and its record. It tells `Report`
# how far it has got as it goes, for a reader watching the protein build.
Score = Callable[[Target, bytes, Report], bytes]



class Modelled(NamedTuple):
    """One protein's model as its baker hands it over (`alphafold.Built`): the
    `.glb`, the scene compiled from it, the fold page's words and chains, and
    the bake's provenance."""

    glb: bytes
    scene: bytes
    chrome: dict
    chains: list
    provenance: dict


# The model of a protein, from its target and its record.
Model = Callable[[Target, bytes], Modelled]

# The AVI track's bytes, from a target and its record (`bake_impact.bake`).
Avi = Callable[[Target, bytes], bytes]

# The ClinVar track's bytes, for a target whose record and AVI track are laid
# out where its baker reads them (`lay_out`).
ClinVar = Callable[[Target], bytes]

STOPPED_FETCHING = "The worker was stopped while ClinVar's records were being fetched."

# The fold page's words as every row carries them (`StructureChrome` in the
# app), and the tints a model made on demand may name: pLDDT's four bands
# (`alphafold.BANDS`) and its bridges.
_CHROME = ("pdb", "modelled", "label", "count", "unit", "sentence", "semantics")
_TINTS = frozenset({"plddtVeryHigh", "plddtConfident", "plddtLow", "plddtVeryLow", "cysteine"})

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


class Stopped(Exception):
    """The reader who asked stopped the bake (`POST /proteins/resolve/{gene}/stop`).

    The service has already said so on its rows: the job is `stopped` and its
    track `absent`. It is raised out through the scorer, which ends, and
    nothing is written after it.
    """


class Unmodelled(Exception):
    """No model is drawn for this protein (`alphafold.Refused`): AlphaFold DB
    has none, or none of this sequence, or none sure enough. The message is
    the sentence the fold page shows."""


class Unplaced(Exception):
    """AlphaGenome's scores cannot be placed on this protein's gene: one of the
    AVI bake's own gates said no (`bake_impact.BakeError`). The message is the
    gate's own words, which the job keeps; the track is told `_unplaced`."""


class Unfetched(RuntimeError):
    """ClinVar did not answer in full: its search, or a batch of its records,
    or the worker was stopped between two. Never a verdict on the protein."""


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


def model_job(target: Target, record: bytes) -> dict:
    """What a model is made for (`alphafold.Job`), from a resolved protein's
    target and its record: the protein the walk draws, its kept regions and
    disulfides, and how far the model may differ from it, which is as far as
    the resolver let the record differ from UniProt."""
    protein = json.loads(record)["protein"]["translation"]
    return {
        "slug": target.slug, "accession": target.uniprot, "display": target.display,
        "protein": protein,
        "kept": [[region.start, region.end] for region in target.regions if region.kept],
        "disulfides": [list(pair) for pair in target.disulfides],
        "allowed": max(VARIANT_FLOOR, int(len(protein) * VARIANT_SHARE)),
    }


def glb_nodes(payload: bytes) -> list:
    """The nodes of a `.glb` that carry a mesh, by name, in file order: read
    from its JSON chunk, which the format puts first."""
    if payload[:4] != b"glTF" or payload[16:20] != b"JSON":
        raise ValueError("not a .glb")
    length = struct.unpack_from("<I", payload, 12)[0]
    document = json.loads(payload[20:20 + length])
    return [node.get("name") for node in document.get("nodes", []) if "mesh" in node]


def stage_structure(storage, target: Target, built: Modelled) -> Tuple[dict, dict]:
    """Check a model, store its two objects, and return the track row and the
    protein row's `structure` column.

    The objects are named as the uploader names a model's
    (`upload_tracks.structure_row`): the scene the app loads, and the `.glb`
    it was compiled from beside it. What is checked is what the app would
    trip on: the seven words it reads, a tint for every node, and each node
    the row names in the model. A model that fails raises `Declined`, and
    nothing is stored.
    """
    chrome, chains = built.chrome, built.chains
    try:
        nodes = glb_nodes(built.glb)
    except ValueError as exc:
        raise Declined(str(exc)) from exc
    span = chrome.get("modelled") if isinstance(chrome, dict) else None
    if (not isinstance(chrome, dict) or tuple(chrome) != _CHROME
            or not all(isinstance(chrome[key], str) and chrome[key]
                       for key in ("pdb", "label", "unit", "sentence", "semantics"))
            or not isinstance(chrome["count"], int)
            or not (span is None or (isinstance(span, list) and len(span) == 2
                                     and all(isinstance(end, int) for end in span)))):
        raise Declined("the fold page's words are not the seven the app reads")
    if (not chains or any(not isinstance(chain, dict) or set(chain) != {"node", "tint"}
                          or chain["tint"] not in _TINTS for chain in chains)):
        raise Declined("a chain names no tint the app has")
    if [chain["node"] for chain in chains] != nodes:
        raise Declined(f"the model's nodes are {nodes}, and its row names "
                       f"{[chain['node'] for chain in chains]}")
    if not built.scene:
        raise Declined("the model has no compiled scene")
    if built.provenance.get("pdb") != chrome["pdb"]:
        raise Declined("the provenance names another entry than the fold page does")

    row, uploads = upload_tracks.structure_row(target.slug, built.glb, built.scene,
                                               built.provenance)
    for path, payload, content_type in uploads:
        storage.put(row["bucket"], path, payload, content_type)
    return row, {"chrome": chrome, "chains": chains}


def resolve_next(conn, storage, *, fetch_entry=uniprot.fetch_entry,
                 structures: bool = False, evidence: bool = False) -> Optional[dict]:
    """Resolve the oldest queued request. None when there is none.

    `structures` also queues the protein's structure bake, and `evidence` its
    AVI and ClinVar bakes (`store.write_resolution`): for a worker that goes
    on to make them."""
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
        store.write_resolution(conn, request["id"], resolution, record, RESOLVER_VERSION,
                               structures=structures, evidence=evidence)
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


class _Reporter:
    """`Report` for one job: its progress written where a reader asking after
    the protein reads it, and a stop noticed there.

    A write that fails costs the bake nothing; the next line tries again, and
    the track is what the bake is for. A job no longer running has been
    stopped, and `Stopped` is raised.
    """

    def __init__(self, conn, job_id: int):
        self._conn = conn
        self._job = job_id

    def __call__(self, done: int, total: int, left: float) -> None:
        try:
            running = store.note_progress(self._conn, self._job, done, total, left)
        except Exception:  # noqa: BLE001 -- see the docstring
            return
        if not running:
            raise Stopped()

    def check(self) -> None:
        """Between two progress lines: `Stopped` once the job has been stopped.
        A scorer that can ask (the Mac's) asks this while the model loads."""
        try:
            running = store.still_running(self._conn, self._job)
        except Exception:  # noqa: BLE001 -- as a write that fails
            return
        if not running:
            raise Stopped()


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
        report = _Reporter(conn, job["id"])
        track = stage(storage, "constraint", target, score(target, record, report))
        store.finish_bake(conn, job["id"], track)
        return {**outcome, "state": "ready"}
    except Stopped:
        return {**outcome, "state": "stopped"}
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


def structure_next(conn, storage, *, model: Model) -> Optional[dict]:
    """Bake the oldest queued structure track. None when there is none."""
    job = store.claim_bake(conn, "structure")
    if job is None:
        return None
    outcome = {"id": job["id"], "slug": job["slug"]}
    try:
        protein = store.protein(conn, job["slug"])
        record_row = store.ready_track(conn, job["slug"], "record")
        if protein is None or record_row is None:
            raise Refused("The protein has no ready record to hold a model to.")
        record = storage.get(record_row["bucket"], record_row["object_path"])
        if hashlib.sha256(record).hexdigest() != record_row["sha256"]:
            raise RuntimeError("The stored record does not match the sha256 its row carries.")
        target = target_of(protein)
        track, structure = stage_structure(storage, target, model(target, record))
        store.finish_structure(conn, job["id"], track, structure)
        return {**outcome, "state": "ready"}
    except (Refused, Unmodelled) as exc:
        store.refuse_bake(conn, job["id"], job["slug"], "structure", _said(exc))
        return {**outcome, "state": "refused", "reason": _said(exc)}
    except Exception as exc:  # noqa: BLE001
        if job["attempts"] < store.MAX_ATTEMPTS:
            store.requeue_bake(conn, job["id"], _said(exc))
            return {**outcome, "state": "queued", "reason": _said(exc)}
        # The fold page shows a refused track's reason, so it is a reader's
        # sentence; what broke stays on the job, for whoever mends it.
        reason = (f"AlphaFold's model could not be made: the bake broke each of the "
                  f"{store.MAX_ATTEMPTS} times it was tried.")
        store.refuse_bake(conn, job["id"], job["slug"], "structure", reason, error=_said(exc))
        return {**outcome, "state": "refused", "reason": reason}


# ------------------------------------------------------------ the evidence


def _stored(storage, row: dict, what: str) -> bytes:
    """A ready track's bytes, held to the digest its row carries."""
    payload = storage.get(row["bucket"], row["object_path"])
    if hashlib.sha256(payload).hexdigest() != row["sha256"]:
        raise RuntimeError(f"The stored {what} does not match the sha256 its row carries.")
    return payload


def lay_out(target: Target, record: bytes, impact: Optional[bytes] = None) -> None:
    """Put a protein's record, and its AVI track where one is given, where the
    bakers and the upload gate read them: `paths.DATA`, under the names the
    twenty's have there. That is the worker's own directory, never
    `pipeline/data/`."""
    for name, payload in ((target.mock_asset, record), (target.impact_asset, impact)):
        if payload is None:
            continue
        path = paths.DATA / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)


def _told(what: str, target: Target, error: BaseException) -> str:
    """A bake's verdict as a reader is told it: what was not done, and why in
    the gate's own words, without the name it starts with or the examples the
    job keeps for an operator."""
    lines = str(error).strip().splitlines()
    said = lines[0].strip() if lines else ""
    for name in (f"{target.slug} ClinVar", target.gene, target.slug):
        if said.startswith(f"{name}: "):
            said = said[len(name) + 2:]
            break
    said = said.rstrip(" :.")
    return _said(f"{what}: {said}." if said else f"{what}.")


def _unmapped(gene: str, why: Optional[str]) -> str:
    """Why ClinVar is refused where AVI was: its records are placed by AVI's
    map. AVI's own reason follows, where it gave one."""
    said = (f"ClinVar's records are placed by AlphaGenome's coordinate map, which "
            f"{gene} has none of.")
    return _said(f"{said} {why}" if why else said)


def impact_next(conn, storage, *, avi: Avi) -> Optional[dict]:
    """Bake the oldest queued AVI track. None when there is none."""
    job = store.claim_bake(conn, "impact")
    if job is None:
        return None
    outcome = {"id": job["id"], "slug": job["slug"]}
    try:
        protein = store.protein(conn, job["slug"])
        record_row = store.ready_track(conn, job["slug"], "record")
        if protein is None or record_row is None:
            raise Refused("The protein has no ready record to place AlphaGenome's scores on.")
        record = _stored(storage, record_row, "record")
        target = target_of(protein)
        # A verdict is the bake's own gate (`Unplaced`) or the uploader's
        # (`Declined`); the Atlas or storage failing is tried again.
        verdict = None
        try:
            payload = avi(target, record)
        except Unplaced as exc:
            verdict = exc
        else:
            # The gate holds the track to the record it is filed under.
            lay_out(target, record)
            try:
                track = stage(storage, "impact", target, payload)
            except Declined as exc:
                verdict = exc
        if verdict is not None:
            # The track carries a reader's sentence, the job the gate's words.
            reason = _told(f"AlphaGenome's scores were not placed on {target.gene}'s bases",
                           target, verdict)
            store.refuse_bake(conn, job["id"], job["slug"], "impact", reason,
                              error=_said(verdict))
            return {**outcome, "state": "refused", "reason": reason}
        store.finish_bake(conn, job["id"], track)
        return {**outcome, "state": "ready"}
    except Refused as exc:
        store.refuse_bake(conn, job["id"], job["slug"], "impact", _said(exc))
        return {**outcome, "state": "refused", "reason": _said(exc)}
    except Exception as exc:  # noqa: BLE001 -- the Atlas, GENCODE's lookup, storage
        if job["attempts"] < store.MAX_ATTEMPTS:
            store.requeue_bake(conn, job["id"], _said(exc))
            return {**outcome, "state": "queued", "reason": _said(exc)}
        reason = (f"AlphaGenome's scores could not be fetched: the bake broke each of the "
                  f"{store.MAX_ATTEMPTS} times it was tried.")
        store.refuse_bake(conn, job["id"], job["slug"], "impact", reason, error=_said(exc))
        return {**outcome, "state": "refused", "reason": reason}


def clinvar_next(conn, storage, *, clinvar: ClinVar) -> Optional[dict]:
    """Bake the oldest queued ClinVar track whose protein's AVI bake is over.
    None when there is none.

    Its records are placed on the protein's bases by AVI's coordinate map
    (`runs`, `chromosome`, `complemented`), so where AVI was refused, so is it,
    with AVI's reason after its own.
    """
    job = store.claim_bake(conn, "clinvar")
    if job is None:
        return None
    outcome = {"id": job["id"], "slug": job["slug"]}
    try:
        protein = store.protein(conn, job["slug"])
        record_row = store.ready_track(conn, job["slug"], "record")
        if protein is None or record_row is None:
            raise Refused("The protein has no ready record to place ClinVar's records on.")
        target = target_of(protein)
        map_row = store.ready_track(conn, job["slug"], "impact")
        if map_row is None:
            state, why = store.track_state(conn, job["slug"], "impact")
            raise Refused(_unmapped(target.gene, why if state == "refused" else None))
        record = _stored(storage, record_row, "record")
        lay_out(target, record, _stored(storage, map_row, "AVI track"))
        # A verdict is the baker's own `ValueError` on the records, or the
        # uploader's (`Declined`). NCBI not answering comes back `Unfetched`,
        # and storage failing as itself: both are tried again.
        verdict = None
        try:
            payload = clinvar(target)
        except ValueError as exc:
            verdict = exc
        else:
            try:
                track = stage(storage, "clinvar", target, payload)
            except Declined as exc:
                verdict = exc
        if verdict is not None:
            reason = _told(f"ClinVar's records were not placed on {target.gene}'s bases",
                           target, verdict)
            store.refuse_bake(conn, job["id"], job["slug"], "clinvar", reason,
                              error=_said(verdict))
            return {**outcome, "state": "refused", "reason": reason}
        store.finish_bake(conn, job["id"], track)
        return {**outcome, "state": "ready"}
    except Refused as exc:
        store.refuse_bake(conn, job["id"], job["slug"], "clinvar", _said(exc))
        return {**outcome, "state": "refused", "reason": _said(exc)}
    except Exception as exc:  # noqa: BLE001 -- NCBI or storage: worth another try
        if job["attempts"] < store.MAX_ATTEMPTS:
            store.requeue_bake(conn, job["id"], _said(exc))
            return {**outcome, "state": "queued", "reason": _said(exc)}
        reason = (f"ClinVar's records could not be fetched: the bake broke each of the "
                  f"{store.MAX_ATTEMPTS} times it was tried.")
        store.refuse_bake(conn, job["id"], job["slug"], "clinvar", reason, error=_said(exc))
        return {**outcome, "state": "refused", "reason": reason}


def clinvar_in_process(target: Target, cache: Path,
                       stopped: Callable[[], bool] = lambda: False) -> bytes:
    """`bake_clinvar.bake` in this process, unchanged, for a target laid out in
    `paths.DATA` (`lay_out`): the ClinVar track's bytes.

    Two of the baker's functions are wrapped while it runs, never edited. Its
    raw responses go to `cache`, not `pipeline/clinvar/cache/` beside the
    twenty's, and are deleted once the bake is over (Phase 0's fifth
    decision): the track names the query, the day and each batch's sha256,
    and a bake made again fetches again. NCBI not answering in full, or
    `stopped` saying so between two batches, is `Unfetched`. What the baker
    prints is kept out of the worker's own output.
    """
    from pipeline.clinvar import bake_clinvar

    fetch, fetch_batch = bake_clinvar.fetch, bake_clinvar.fetch_batch
    folder = cache / target.gene

    def fetching(gene, _cache, replay):
        try:
            return fetch(gene, folder, replay)
        except Unfetched:
            raise
        except Exception as exc:  # noqa: BLE001 -- NCBI's, whatever it was
            raise Unfetched(f"ClinVar did not answer in full: {_said(exc)}") from exc

    def batch(efetch, ids, *args, **kwargs):
        if stopped():
            raise Unfetched(STOPPED_FETCHING)
        return fetch_batch(efetch, ids, *args, **kwargs)

    bake_clinvar.fetch, bake_clinvar.fetch_batch = fetching, batch
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            bake_clinvar.bake(target, replay=False)
    finally:
        bake_clinvar.fetch, bake_clinvar.fetch_batch = fetch, fetch_batch
        shutil.rmtree(folder, ignore_errors=True)
    return (paths.DATA / f"assets/clinvar/{target.slug}_clinvar.json").read_bytes()


def impact_all(conn, storage, *, limit: int = 20, avi: Avi) -> list[dict]:
    """Bake queued AVI tracks until none are left, or `limit` are done."""
    done = []
    for _ in range(limit):
        outcome = impact_next(conn, storage, avi=avi)
        if outcome is None:
            break
        done.append(outcome)
    return done


def clinvar_all(conn, storage, *, limit: int = 20, clinvar: ClinVar) -> list[dict]:
    """Bake queued ClinVar tracks until none are left, or `limit` are done."""
    done = []
    for _ in range(limit):
        outcome = clinvar_next(conn, storage, clinvar=clinvar)
        if outcome is None:
            break
        done.append(outcome)
    return done


def sweep(conn, storage, *, limit: int = 10, fetch_entry=uniprot.fetch_entry,
          structures: bool = False, evidence: bool = False) -> dict:
    """Reap what dead workers held, then resolve up to `limit` requests.

    Returns what it did, and how many constraint bakes are waiting, so the
    caller knows whether to start a GPU, and how many structure bakes.
    `structures` queues one for each protein resolved, and `evidence` its AVI
    and ClinVar bakes (`resolve_next`).
    """
    store.reap(conn)
    resolved = []
    for _ in range(limit):
        outcome = resolve_next(conn, storage, fetch_entry=fetch_entry, structures=structures,
                               evidence=evidence)
        if outcome is None:
            break
        resolved.append(outcome)
    return {"resolved": resolved, "constraint_queued": store.queued_bakes(conn, "constraint"),
            "structure_queued": store.queued_bakes(conn, "structure")}


def structure_all(conn, storage, *, limit: int = 20, model: Model) -> list[dict]:
    """Bake queued structure tracks until none are left, or `limit` are done."""
    done = []
    for _ in range(limit):
        outcome = structure_next(conn, storage, model=model)
        if outcome is None:
            break
        done.append(outcome)
    return done


def score_all(conn, storage, *, limit: int = 20, score: Score = score_in_process) -> list[dict]:
    """Bake queued constraint tracks until none are left, or `limit` are done."""
    done = []
    for _ in range(limit):
        outcome = score_next(conn, storage, score=score)
        if outcome is None:
            break
        done.append(outcome)
    return done
