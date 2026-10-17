"""Upload baked track assets to Supabase storage and record where they landed.

    DATABASE_URL=... SUPABASE_URL=... SUPABASE_SERVICE_KEY=... \
        python3 pipeline/upload_tracks.py --kind impact_explanations
    python3 pipeline/upload_tracks.py --kind impact_explanations --dry-run

The files come from `pipeline/data/` (`paths.DATA`), where the bakes write them.

One family at a time, because each is validated differently and each is worth
watching land on its own. The bytes go to a public bucket the client fetches
directly; this writes the `protein_track` row that names them.

**Stored plain, never pre-gzipped.** Measured 2026-09-24: Supabase's CDN
compresses in transit whenever the client sends `Accept-Encoding: gzip`, on
files of any size -- dystrophin's 9,563,683-byte ClinVar snapshot goes over the
wire in 572,718 bytes, within 0.7% of our own `gzip -9`. Storing a `.gz` makes
the CDN gzip an already gzipped payload and forces the client to decompress
twice for nothing.

**Paths carry the digest**: `<kind>/<slug>.<sha12>.json`. The client is handed
the URL out of the track row and never builds it, so the path is free to change
with the content -- which is what makes `cache-control: immutable` true rather
than merely convenient. A rebake writes a new object and the superseded one is
deleted once its row no longer points at it.

**Validated before it is written**, the same order `record_cache` uses: a
payload that would fail on the way out never becomes a row that says `ready`.
A target refused for that is said and makes the run exit 1, after the others.

**A structure is compiled on the way**: the `.glb` the bake made goes up with
the `.fsceneb` the app loads, compiled by the app's own importer
(`alphafold.compile_scene`, the flutter_scene the app's checkout pins), as the
worker compiles a model made on demand. `--app` names the checkout and
`--dart` the Dart that runs the importer.

`--against <service>` (with `--dry-run`) says, for each track planned, whether
the service already serves exactly those bytes: the proof that a re-bake, or
the compiler, makes what is stored.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import check_assets  # noqa: E402
from pipeline.impact.check_explanations import validate as validate_explanations  # noqa: E402
from pipeline.paths import CLIENT, DATA  # noqa: E402
from pipeline.targets import TARGETS  # noqa: E402

MODELS_BUCKET = "models"

# A year, and immutable: the digest is in the path, so these bytes never change.
CACHE_CONTROL = "public, max-age=31536000, immutable"


def compiled_scene(glb: Path, app: Path, dart: Path) -> bytes:
    """The `.fsceneb` the app loads, compiled from [glb] by the app's own
    importer, as the worker compiles a model made on demand. Raises
    RuntimeError where no scene can be compiled from that checkout."""
    import tempfile

    from pipeline.structure import alphafold

    unready = alphafold.importer_unready(app)
    if unready:
        raise RuntimeError(unready)
    with tempfile.TemporaryDirectory() as work:
        return alphafold.compile_scene(glb, Path(work) / "model.fsceneb", app, str(dart))


def default_dart() -> Path:
    """The Dart the worker runs the importer with (`local_worker.default_dart`),
    unless `RESOLVER_DART` names another, as it does for the worker."""
    from pipeline.resolver.local_worker import default_dart as worker_dart

    named = os.environ.get("RESOLVER_DART")
    return Path(named) if named else worker_dart(Path.home())


def served_as(service: str, row: dict) -> str:
    """Whether [service] already serves exactly this row's bytes: `same`,
    `differs`, or why it cannot say."""
    url = f"{service.rstrip('/')}/protein/{row['slug']}/tracks"
    try:
        with urllib.request.urlopen(url, timeout=120) as response:
            live = json.loads(response.read()).get(row["kind"]) or {}
    except urllib.error.HTTPError as exc:
        return f"not served ({exc.code})"
    if live.get("state") != "ready":
        return f"not ready ({live.get('state')})"
    same = live.get("sha256") == row["sha256"]
    if row["kind"] == "structure":
        glb = (row["provenance"].get("glb") or {}).get("sha256")
        same = same and ((live.get("provenance") or {}).get("glb") or {}).get("sha256") == glb
    return "same" if same else "differs"


def asset_of(kind: str, target) -> Path | None:
    """The local file for this kind, or None where this protein has none."""
    if kind == "record":
        return DATA / target.mock_asset
    if kind == "impact_explanations":
        return DATA / f"assets/impact_explanations/{target.slug}.json"
    if kind == "constraint":
        return DATA / target.constraint_asset if target.scored else None
    if kind == "impact":
        return DATA / target.impact_asset if target.impact_scored else None
    if kind == "clinvar":
        path = DATA / f"assets/clinvar/{target.slug}_clinvar.json"
        return path if target.clinvar_available else None
    if kind == "structure":
        return DATA / target.structure_asset
    if kind == "folding":
        # Where `folding/bake_folding.py` writes it (`folding_asset`).
        path = DATA / f"assets/folding/{target.slug}_folding.json"
        return path if target.structure else None
    if kind == "locus":
        # Where `locus/bake_locus.py` writes it (`locus_asset`).
        return DATA / f"assets/locus/{target.slug}_locus.json"
    raise SystemExit(f"unknown kind {kind!r}")


def _checked(check, target, data: dict) -> None:
    """One of `check_assets`' checks on a payload, against the record it is
    filed under; a ValueError naming every problem it found."""
    record = json.loads((DATA / target.mock_asset).read_text())
    before = len(check_assets.problems)
    try:
        check(target, record, data)
    except (KeyError, TypeError, IndexError, AttributeError) as exc:
        raise ValueError(f"the payload is not shaped as the checks read it: {exc!r}") from exc
    found = check_assets.problems[before:]
    if found:
        raise ValueError("; ".join(found))


def validate(kind: str, target, payload: bytes) -> dict:
    """Check the payload and return what its provenance should say.

    Raises ValueError with a reason. Nothing that raises here is uploaded, so a
    malformed asset cannot become a row that claims to be ready.
    """
    if kind == "structure":
        if payload[:4] != b"glTF":
            raise ValueError("not a .glb")
        return {"pdb": target.structure.pdb,
                "nodes": [c.node for c in target.structure.chains]}

    data = json.loads(payload)
    if not isinstance(data, dict):
        raise ValueError("payload is not an object")

    if kind == "record":
        if data.get("gene") != target.gene:
            raise ValueError(f"gene is {data.get('gene')!r}, expected {target.gene!r}")
        # The offline checker's own gate, so a record that would fail
        # `check_assets.py` never becomes a row that says ready: the CDS
        # translates to the protein, every peptide to its slice, and every
        # intron is a known splice class.
        before = len(check_assets.problems)
        check_assets.check_record(data, target.slug)
        found = check_assets.problems[before:]
        if found:
            raise ValueError("; ".join(found))
        return {"source": "NCBI Entrez", "accession": target.source.accession,
                "protein_id": target.source.protein_id,
                "transcript_id": target.source.transcript_id,
                "built_by": "pipeline/mock/build_gene_record.py"}

    if kind == "impact_explanations":
        impact = (DATA / target.impact_asset).read_bytes()
        # The same check the endpoint used to make on every read, made once on
        # the way in: scorer, schema version, and the digest binding this
        # payload to the exact AVI track it explains.
        validate_explanations(impact, data)
        if data.get("accession") != target.source.accession or data.get("gene") != target.gene:
            raise ValueError("accession/gene do not match the target")
        return {key: data.get(key) for key in
                ("scorer", "units", "scope", "selection", "schema_version",
                 "generated_at", "client_version", "impact_sha256")}

    if kind == "constraint":
        if data.get("gene") != target.gene:
            raise ValueError(f"gene is {data.get('gene')!r}, expected {target.gene!r}")
        return {key: data.get(key) for key in
                ("model", "revision", "method", "normalization", "score_units",
                 "entropy_vocabulary", "vocabulary_size", "context")}

    if kind == "impact":
        if data.get("gene") != target.gene:
            raise ValueError(f"gene is {data.get('gene')!r}, expected {target.gene!r}")
        # The offline checker's own gate, held to the record the track is filed
        # under (`DATA / target.mock_asset`), as the record's is: its letters
        # and coordinates, a map covering every drawn base, three scores a
        # base, and the biology the sheet claims. A bake the worker made with
        # nobody watching meets it before it becomes a row that says ready.
        _checked(check_assets.check_impact, target, data)
        return {key: data.get(key) for key in
                ("scorer", "score_units", "assembly", "annotation", "chromosome",
                 "transcript", "orientation", "complemented")}

    if kind == "clinvar":
        if data.get("gene") != target.gene:
            raise ValueError(f"gene is {data.get('gene')!r}, expected {target.gene!r}")
        # The same, for the AVI map its records were placed by
        # (`DATA / target.impact_asset`): each allele, strand and codon
        # derived again from the record.
        _checked(check_assets.check_clinvar, target, data)
        return {key: data.get(key) for key in
                ("source", "schema_version", "assembly", "scope", "retrieved_at",
                 "searched_records", "excluded")}

    if kind == "folding":
        # The whole of `check_folding.py` on this protein: the mature chain
        # residue for residue, the record's letters, the entry's own CA atoms
        # in a chain, and each on its ribbon in the stored model. The
        # provenance is the entry, its chains and the frame they were put in.
        from pipeline.folding.check_folding import problems_of
        found = problems_of(target, data)
        if found:
            raise ValueError("; ".join(found))
        return {
            "pdb": data["pdb"],
            "chains": [{key: chain[key] for key in
                        ("node", "pdb_chain", "offset", "first", "last")}
                       for chain in data["chains"]],
            "frame": data["frame"],
            "secondary_structure": data["secondary_structure"],
            "schema_version": data["schema_version"],
            "built_by": data["built_by"],
        }

    if kind == "locus":
        # The whole of `check_locus.py` on this protein: its own names, the
        # chromosome whole and the gene at its bands, the record's slice where
        # it has one, and a version for every source. The provenance is where
        # the gene lies and what that was read from.
        from pipeline.locus.check_locus import problems_of
        found = problems_of(target, data)
        if found:
            raise ValueError("; ".join(found))
        return {key: data.get(key) for key in
                ("assembly", "genome", "chromosome", "sequence", "locus", "sources",
                 "retrieved", "schema_version", "built_by")}

    raise SystemExit(f"unknown kind {kind!r}")


def structure_row(slug: str, glb: bytes, compiled: bytes, provenance: dict) -> tuple[dict, list]:
    """A structure track's row and its two objects, from a model and the scene
    compiled from it.

    The row points at what the phone loads. `.fsceneb` is a versioned
    container tied to the flutter_scene version, so the `.glb` it was compiled
    from goes up beside it and is named in the provenance: an upgrade that
    invalidates every scene can then re-compile from storage rather than
    needing the repo and another PyMOL run. Returns the row, and the objects
    as (path, bytes, content type), the `.glb` first.
    """
    digest = hashlib.sha256(glb).hexdigest()
    scene_digest = hashlib.sha256(compiled).hexdigest()
    glb_path = f"structure/{slug}.{digest[:12]}.glb"
    scene_path = f"structure/{slug}.{scene_digest[:12]}.fsceneb"
    row = {
        "slug": slug, "kind": "structure", "bucket": MODELS_BUCKET,
        "object_path": scene_path,
        "bytes": len(compiled), "sha256": scene_digest, "format": "fsceneb",
        "provenance": {**provenance, "glb": {
            "path": glb_path,
            "sha256": digest,
            "bytes": len(glb),
        }},
    }
    return row, [(glb_path, glb, "model/gltf-binary"),
                 (scene_path, compiled, "application/octet-stream")]


class Storage:
    def __init__(self, url: str, key: str):
        self.base = url.rstrip("/") + "/storage/v1"
        self.headers = {"Authorization": f"Bearer {key}", "apikey": key}

    def put(self, bucket: str, path: str, payload: bytes, content_type: str) -> None:
        request = urllib.request.Request(
            f"{self.base}/object/{bucket}/{path}", method="POST", data=payload,
            headers={**self.headers, "Content-Type": content_type,
                     "Cache-Control": CACHE_CONTROL, "x-upsert": "true"})
        with urllib.request.urlopen(request, timeout=300):
            return

    def delete(self, bucket: str, paths: list[str]) -> None:
        if not paths:
            return
        request = urllib.request.Request(
            f"{self.base}/object/{bucket}", method="DELETE",
            data=json.dumps({"prefixes": paths}).encode(),
            headers={**self.headers, "Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(request, timeout=120):
                return
        except urllib.error.HTTPError as exc:
            print(f"  could not delete {paths}: {exc.code}", file=sys.stderr)


_UPSERT = """
insert into protein_track
    (slug, kind, state, reason, bucket, object_path, bytes, sha256,
     content_encoding, format, provenance, updated_at)
values (%(slug)s, %(kind)s, 'ready', null, %(bucket)s, %(object_path)s, %(bytes)s,
        %(sha256)s, null, %(format)s, %(provenance)s, now())
on conflict (slug, kind) do update set
    state = 'ready', reason = null, bucket = excluded.bucket,
    object_path = excluded.object_path, bytes = excluded.bytes,
    sha256 = excluded.sha256, content_encoding = null,
    format = excluded.format, provenance = excluded.provenance,
    updated_at = now()
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Upload baked tracks to Supabase.")
    parser.add_argument("--kind", required=True, choices=[
        "record", "impact_explanations", "constraint", "impact", "clinvar", "structure",
        "folding", "locus"])
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--app", type=Path, default=CLIENT,
                        help=f"the app's checkout, whose importer compiles a scene (default {CLIENT})")
    parser.add_argument("--dart", type=Path, default=None,
                        help="the Dart that runs it (default: the worker's, or $RESOLVER_DART)")
    parser.add_argument("--against", default=None, metavar="SERVICE",
                        help="with --dry-run: whether the service already serves these bytes")
    arguments = parser.parse_args()
    if arguments.against and not arguments.dry_run:
        parser.error("--against compares a dry run with what is served; add --dry-run")

    kind = arguments.kind
    wanted = set(arguments.target) if arguments.target else None
    bucket, suffix, content_type = {
        "structure": (MODELS_BUCKET, "glb", "model/gltf-binary"),
    }.get(kind, ("tracks", "json", "application/json"))

    unknown = sorted((wanted or set()) - {t.slug for t in TARGETS})
    if unknown:
        raise SystemExit(f"no such target: {unknown}")
    planned, refused = [], []
    for target in TARGETS:
        if wanted and target.slug not in wanted:
            continue
        path = asset_of(kind, target)
        if path is None or not path.exists():
            continue
        payload = path.read_bytes()
        try:
            provenance = validate(kind, target, payload)
        except (ValueError, KeyError) as exc:
            print(f"  {target.slug}: REFUSED -- {exc}", file=sys.stderr)
            refused.append(target.slug)
            continue
        digest = hashlib.sha256(payload).hexdigest()
        uploads = [(f"{kind}/{target.slug}.{digest[:12]}.{suffix}",
                    payload, content_type)]
        row = {
            "slug": target.slug, "kind": kind, "bucket": bucket,
            "object_path": uploads[0][0],
            "bytes": len(payload), "sha256": digest, "format": suffix,
            "provenance": provenance,
        }
        if kind == "structure":
            try:
                scene = compiled_scene(path, arguments.app, arguments.dart or default_dart())
            except (RuntimeError, OSError) as exc:
                print(f"  {target.slug}: REFUSED -- no scene could be compiled: {exc}",
                      file=sys.stderr)
                refused.append(target.slug)
                continue
            row, uploads = structure_row(target.slug, payload, scene, provenance)
        row["_uploads"] = uploads
        planned.append(row)

    total = sum(len(payload) for row in planned for _, payload, _ in row["_uploads"])
    for row in planned:
        for path, payload, _ in row["_uploads"]:
            print(f"  {row['slug']:<16} {len(payload):>10,} B  {path}")
        if arguments.against:
            print(f"  {row['slug']:<16} {served_as(arguments.against, row)} "
                  f"as {arguments.against}")
    objects = sum(len(row["_uploads"]) for row in planned)
    print(f"\n{len(planned)} track(s), {objects} object(s), {total / 1e6:.2f} MB"
          + (f"; {len(refused)} refused: {', '.join(refused)}" if refused else ""),
          file=sys.stderr)
    if arguments.dry_run:
        return 1 if refused else 0

    url, key = os.environ.get("SUPABASE_URL"), os.environ.get("SUPABASE_SERVICE_KEY")
    database = os.environ.get("DATABASE_URL")
    if not (url and key and database):
        raise SystemExit("SUPABASE_URL, SUPABASE_SERVICE_KEY and DATABASE_URL must be set.")

    import psycopg
    from psycopg.types.json import Jsonb

    storage = Storage(url, key)
    superseded = []
    with psycopg.connect(database) as conn:
        with conn.cursor() as cur:
            for row in planned:
                old = cur.execute(
                    "select object_path, provenance from protein_track "
                    "where slug = %s and kind = %s and state = 'ready'",
                    (row["slug"], row["kind"])).fetchone()
                uploads = row.pop("_uploads")
                for path, payload, mime in uploads:
                    storage.put(bucket, path, payload, mime)
                cur.execute(_UPSERT, {**row, "provenance": Jsonb(row["provenance"])})
                if old:
                    written = {path for path, _, _ in uploads}
                    # Both of a structure row's objects, where it had two.
                    stale = [old[0], ((old[1] or {}).get("glb") or {}).get("path")]
                    superseded.extend(p for p in stale if p and p not in written)
                print(f"  uploaded {row['slug']}", file=sys.stderr)
        conn.commit()

    # Only after the rows point elsewhere: nothing is deleted while it is still
    # what a client would be told to fetch.
    storage.delete(bucket, superseded)
    if superseded:
        print(f"  removed {len(superseded)} superseded object(s)", file=sys.stderr)
    print(f"{len(planned)} {kind} track(s) ready"
          + (f"; {len(refused)} refused: {', '.join(refused)}" if refused else ""),
          file=sys.stderr)
    return 1 if refused else 0


if __name__ == "__main__":
    raise SystemExit(main())
