"""Read-only check of one protein the resolver built, against the live service,
the stored bytes, the database rows and UniProt's own features.

    .venv/bin/python scripts/check_built.py ANTXR1        # any device, printed
    .venv/bin/python scripts/check_built.py ANTXR1 mps    # scored by the Mac's worker
    .venv/bin/python scripts/check_built.py ANTXR1 cuda   # scored on Modal

From the repository root: the rows are read through `.env`'s DATABASE_URL.
Nothing here writes. It prints what it found and exits 1 if anything is wrong.

It checks a protein whose ESM-2 track is ready. One whose scores the scorer
refused has no track to check here; the app's test/live_resolved_check.dart
reads it as the app does (LIVE_EXPECT=refused).

Its model is checked as it stands: ready (the stored scene and the `.glb` it
was compiled from, each to its digest, and the fold page's words and chains to
the model), refused (the sentence why, and nothing on the row), or never asked
for. A model still pending is not a finished build.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import urllib.request

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import psycopg  # noqa: E402

from app.config import settings  # noqa: E402
from pipeline import uniprot  # noqa: E402
from pipeline.constraint.score_protein import MODEL, REVISION  # noqa: E402
from pipeline.resolver import worker  # noqa: E402

SERVICE = "https://helix-peak-backend.onrender.com"
wrong: list[str] = []


def get(path: str):
    with urllib.request.urlopen(SERVICE + path, timeout=90) as response:
        return json.load(response)


def fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def hold(ok: bool, what: str) -> None:
    print(("  ok    " if ok else "  WRONG ") + what)
    if not ok:
        wrong.append(what)


def structure(row: dict, track: dict) -> None:
    """The protein's model, as far as it has got (`pipeline/structure/alphafold.py`)."""
    state, said = track["state"], track.get("provenance") or {}
    if state != "ready":
        print(f"  structure: {state}" + (f": {track['reason']}" if track.get("reason") else ""))
        hold(state in ("refused", "absent"), f"its model is not still on its way ({state})")
        hold(row["structure"] is None and row["chains"] == [],
             "no model, so no words or chains for one on the row")
        return

    scene = fetch(track["url"])
    hold(track["format"] == "fsceneb" and hashlib.sha256(scene).hexdigest() == track["sha256"]
         and len(scene) == track["bytes"], "scene bytes match their sha256 and size")
    # The `.glb` it was compiled from is beside it in the bucket.
    glb = fetch(track["url"].rsplit("/structure/", 1)[0] + "/" + said["glb"]["path"])
    hold(hashlib.sha256(glb).hexdigest() == said["glb"]["sha256"] and len(glb) == said["glb"]["bytes"],
         "the .glb beside it matches the sha256 and size its provenance names")
    chrome, chains = row["structure"], row["chains"]
    print(f"  structure: {len(scene):,} bytes at {track['url'].rsplit('/', 2)[-2:]}, {said.get('entry')}"
          f" v{said.get('model_version')}, span {said.get('span')}, mean pLDDT {said.get('mean_plddt')},"
          f" shares {said.get('plddt_shares')}, sampling {said.get('sampling')},"
          f" {said.get('representation')}, frame {said.get('frame')}")
    print(f"  structure: bridges {said.get('bridges')}, dropped {said.get('bridges_dropped')},"
          f" differences {said.get('sequence_differences')}, importer {said.get('importer')},"
          f" PyMOL {said.get('pymol')}")
    print(f"  fold page: {chrome}")
    print(f"  chains: {chains}")
    hold(said.get("source") == "AlphaFold DB" and said.get("licence") == "CC BY 4.0"
         and said.get("entry") == f"AF-{row['uniprot']}-F1" == (chrome or {}).get("pdb"),
         "AlphaFold DB's canonical entry for this protein, credited CC BY 4.0")
    hold(tuple(chrome or {}) == worker._CHROME, "the fold page's seven words")
    hold(worker.glb_nodes(glb) == [chain["node"] for chain in chains] == said.get("nodes")
         and all(chain["tint"] in worker._TINTS for chain in chains),
         "every node the row names is in the model, with a tint the app has")
    kept = [(r["start"], r["end"]) for r in row["regions"] if r["kept"]]
    span = [min(a for a, _ in kept), max(b for _, b in kept)] if kept else [1, row["facts"]["residues"]]
    hold(said.get("span") == span and chrome["count"] == span[1] - span[0] + 1
         and chrome["modelled"] == (None if span == [1, row["facts"]["residues"]] else span),
         f"drawn over the mature span {span}, and counted so")
    hold(said.get("mean_plddt", 0) >= said.get("gate", 50), "mean pLDDT over the span clears the gate")
    hold(all(pair in row["disulfides"] for pair in said.get("bridges", []))
         and len(said.get("bridges", [])) + len(said.get("bridges_dropped", [])) == len(row["disulfides"])
         and (("bonds" in said.get("nodes", [])) == bool(said.get("bridges"))),
         "its bridges are the row's, each drawn or accounted for")


def main(gene: str, device: str | None) -> int:
    slug = gene.lower()
    print(f"== {gene} ==")

    said = get(f"/proteins/resolve/{gene}")
    hold(said == {"slug": slug, "state": "ready", "reason": None}, f"resolve state: {said}")

    row = get(f"/protein/{slug}")
    tracks = get(f"/protein/{slug}/tracks")
    facts, provenance = row["facts"], row["provenance"]
    print(f"  row: {row['display']!r}, {row['gene']} / {row['uniprot']}, accession {row['accession']},"
          f" transcript {row.get('transcript_id')}, protein {row.get('protein_id')}")
    print(f"  facts: {facts}; mature_peptides {row['mature_peptides']}; chain {row.get('chain')!r};"
          f" structure {row.get('structure')}")
    print(f"  summary: {row['summary']}")
    print(f"  regions: {[(r['short'], r['label'], r['start'], r['end'], r['kept']) for r in row['regions']]}")
    print(f"  disulfides: {row['disulfides']}")
    print(f"  provenance: rule {provenance.get('regions_rule')!r}, cleaved {provenance.get('cleaved')},"
          f" chain_label {provenance.get('chain_label')!r}, variants {provenance.get('uniprot_variants')},"
          f" skipped {provenance.get('skipped_features')}, dropped {provenance.get('disulfides_dropped')}")
    print(f"  provenance: uniprot {provenance.get('uniprot')}, mane {provenance.get('mane')}")
    hold(row["catalog_order"] is None, "catalog_order is null (never joins the list)")
    hold(row["resolver_version"] == 1 and provenance.get("prose") == "templated"
         and provenance.get("source") == "pipeline/resolver", "resolver version 1, templated prose")
    states = {kind: track["state"] for kind, track in tracks.items()}
    hold(states.get("record") == "ready" and states.get("constraint") == "ready",
         f"record and constraint ready: {states}")
    hold(all(state == "absent" for kind, state in states.items()
             if kind not in ("record", "constraint", "structure")),
         "every other track absent")
    hold(row["tracks"] == states, "the row's bare states agree with its track rows")
    structure(row, tracks["structure"])

    # The stored bytes, held to the digests their rows carry.
    record_bytes = fetch(tracks["record"]["url"])
    constraint_bytes = fetch(tracks["constraint"]["url"])
    hold(hashlib.sha256(record_bytes).hexdigest() == tracks["record"]["sha256"]
         and len(record_bytes) == tracks["record"]["bytes"], "record bytes match their sha256 and size")
    hold(hashlib.sha256(constraint_bytes).hexdigest() == tracks["constraint"]["sha256"]
         and len(constraint_bytes) == tracks["constraint"]["bytes"], "constraint bytes match their sha256 and size")
    record = json.loads(record_bytes)
    constraint = json.loads(constraint_bytes)
    translation = record["protein"]["translation"]
    print(f"  record: {len(record_bytes):,} bytes at {tracks['record']['url'].rsplit('/', 2)[-2:]},"
          f" gene {record['gene']}, {len(record['exons'])} exons, protein {len(translation)} aa")
    print(f"  record provenance: {tracks['record'].get('provenance')}")
    hold(record["gene"] == gene and len(translation) == facts["residues"]
         and len(record["exons"]) == facts["exons"], "record is this gene; protein and exons as the row counts them")

    generation = constraint.get("generation", {})
    print(f"  constraint: {len(constraint_bytes):,} bytes, model {constraint.get('model')}@{str(constraint.get('revision'))[:12]},"
          f" generation {generation}, context {constraint.get('context')}")
    print(f"  constraint provenance: {tracks['constraint'].get('provenance')}")
    hold(constraint.get("model") == MODEL and constraint.get("revision") == REVISION,
         "the twenty's model at its pinned revision")
    if device:
        hold(generation.get("device") == device, f"scored on {device}")
    hold(constraint.get("gene") == gene and constraint.get("sequence") == translation
         and len(constraint.get("positions", [])) == len(translation),
         "one position a residue, over the record's own protein")
    hold([list(pair) for pair in constraint.get("disulfides", [])] == row["disulfides"],
         "the track carries the row's disulfides")
    print(f"  constraint regions: {constraint.get('regions')}")

    # UniProt's own word, read independently of the resolver.
    entry = uniprot.fetch_entry(row["uniprot"])
    canonical = entry.body["sequence"]["value"]
    features = [(f["type"], f["location"]["start"].get("value"), f["location"]["end"].get("value"),
                 f.get("description", "")) for f in entry.body.get("features", [])
                if f["type"] in ("Initiator methionine", "Signal", "Transit peptide", "Propeptide",
                                 "Chain", "Peptide", "Disulfide bond")]
    print(f"  UniProt {entry.accession} ({entry.release}): {len(canonical)} aa; features {features}")
    differences = [[i, a, b] for i, (a, b) in enumerate(zip(translation, canonical), 1) if a != b]
    hold(len(translation) == len(canonical) and differences == [list(v) for v in provenance.get("uniprot_variants", [])],
         f"record protein vs UniProt: {len(differences)} difference(s), as the row declares")
    bonds = sorted([a, b] for kind, a, b, _ in features if kind == "Disulfide bond")
    hold(sorted(row["disulfides"]) == bonds, f"disulfides are UniProt's: {bonds}")
    hold(all(translation[a - 1] == "C" and translation[b - 1] == "C" for a, b in row["disulfides"]),
         "every bridge joins two cysteines of the record's protein")
    regions = row["regions"]
    hold(regions[0]["start"] == 1 and regions[-1]["end"] == len(translation)
         and all(a["end"] + 1 == b["start"] for a, b in zip(regions, regions[1:])),
         "regions tile the protein, end to end")
    leaders = [(a, b) for kind, a, b, _ in features if kind in ("Initiator methionine", "Signal", "Transit peptide")]
    removed = [(r["start"], r["end"]) for r in regions if not r["kept"]]
    hold(all(leader in removed for leader in leaders), f"UniProt's leader(s) {leaders} are removed regions {removed}")

    # The rows behind it, and the list beside it.
    with psycopg.connect(settings.database_url, prepare_threshold=None) as conn:
        with conn.transaction():
            conn.execute("set transaction read only")
            request = conn.execute(
                "select state, attempts, resolver_version, reason, finished_at - requested_at "
                "from resolve_request where lower(gene) = lower(%s) order by id desc limit 1", (gene,)).fetchone()
            bake = conn.execute(
                "select state, attempts, error, finished_at - started_at from bake_job where slug = %s "
                "and kind = 'constraint' order by id desc limit 1", (slug,)).fetchone()
            model = conn.execute(
                "select state, attempts, error, finished_at - started_at from bake_job where slug = %s "
                "and kind = 'structure' order by id desc limit 1", (slug,)).fetchone()
            where = conn.execute(
                "select accession, slice_start, slice_end, slice_strand, transcript_id, protein_id, chain_name "
                "from protein where slug = %s", (slug,)).fetchone()
            curated = conn.execute("select count(*) from protein where catalog_order is not null").fetchone()[0]
            aliases = conn.execute(
                "select alias, kind from protein_alias where slug = %s order by kind, alias", (slug,)).fetchall()
    print(f"  protein row source: {where}")
    print(f"  request: {request}")
    print(f"  bake: {bake}")
    print(f"  structure bake: {model}")
    print(f"  aliases: {aliases}")
    hold(request is not None and request[0] == "done" and request[2] == 1, "request done by resolver version 1")
    hold(bake is not None and bake[0] == "done", "constraint bake done")
    catalog = get("/catalog")
    hold(len(catalog["proteins"]) == 20 and curated == 20
         and slug not in {p["slug"] for p in catalog["proteins"]}, "/catalog is still the twenty")
    suggestions = get(f"/proteins/suggest?q={gene}")["suggestions"]
    mine = [s for s in suggestions if s["gene"] == gene]
    hold(bool(mine) and mine[0]["status"] == "ready" and mine[0]["slug"] == slug,
         f"suggest calls it ready: {[(s['gene'], s['status'], s['slug']) for s in mine[:1]]}")

    print(f"== {gene}: {'ALL OK' if not wrong else str(len(wrong)) + ' WRONG'} ==")
    return 1 if wrong else 0


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        raise SystemExit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None))
