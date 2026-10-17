"""Draft a curated row from a gene: everything the pipeline derives, for a person to finish.

    NCBI_EMAIL=you@example.com .venv/bin/python scripts/draft_curated.py --out /tmp/draft \\
        PLN=phospholamban@1ZLL S100B=s100b

For each `GENE=slug`, and the PDB entry it must draw where one is named after `@`:

1. Resolves it as the resolver does (`resolve.resolve`, live NCBI and UniProt;
   its index row and MANE release read from `.env`'s database in a read-only
   transaction), under the curated slug.
2. Chooses its fold by the on-demand path's rule (`fold_job.choose`, policy 1),
   bakes the choice and holds it to the app's own check, by running
   `structure/experimental.py --job` in the structure bake's environment, as
   `dry_fold.py --bake` does. A gene whose search ends anything but `chosen`,
   or on another entry than the one named, is not drafted.
3. Stores the working copy that bake was made from as
   `pipeline/structure/structures/working/<slug>.pdb` (the file the curated
   bakers then cut the row's model and fold from, `pdb.source_of`), and beside
   it `<slug>.json`: where it came from (entry, model, policy, the entry's own
   chains and numbering) and the digests of the record, model, scene and fold
   the app's check passed. A stored working copy is never overwritten unless
   `--replace` says so.
4. Writes into `--out/<slug>/` the `targets.py` row (`row.py`) and the
   `curated/catalog.json` entry (`entry.json`), beside the bake's own files,
   and prints both.

What is left to a person: every word, which starts `DRAFT:` until it is
written by hand (R4.5; `check_assets.py` and `seed_catalog.py` refuse a row
still saying it), and the decisions the rules leave to the record, marked
`# HAND:` in the row (R3.3, R3.4). Nothing is queued, stored in the database
or uploaded. The GenBank cache and the bakers' data directory are temporary.
"""

from __future__ import annotations

import argparse
import atexit
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

DRAFT = "DRAFT:"
SPEC = re.compile(r"^([A-Za-z0-9-]+)=([a-z0-9][a-z0-9-]*)(?:@([0-9A-Za-z]{4}))?$")
# The families every curated row so far is baked with (`curated/catalog.json`'s
# `tracks`); AVI explanations are a pilot of three, not drafted.
TRACKS = ["clinvar", "constraint", "impact", "structure"]


@dataclass(frozen=True)
class Spec:
    gene: str
    slug: str
    entry: str | None = None


def parse_spec(text: str) -> Spec:
    """`GENE=slug` or `GENE=slug@ENTRY`."""
    found = SPEC.match(text)
    if not found:
        raise ValueError(f"{text!r}: write GENE=slug, or GENE=slug@ENTRY")
    gene, slug, entry = found.groups()
    return Spec(gene.upper(), slug, entry.upper() if entry else None)


# ------------------------------------------------------------ what is drafted


def structure_of(provenance: dict):
    """The row's `Structure`, as the bake drew it: the working copy's own chain
    letters, a span only where the entry holds a domain, and the cartoon's
    sampling and representation the bake settled on."""
    from pipeline.structure.normalize import LETTERS
    from pipeline.targets import Chain, Structure

    chains = provenance["chains"]
    span = (chains[0]["first"], chains[0]["last"]) if provenance["scope"] == "domain" else None
    return Structure(
        pdb=provenance["pdb"],
        chains=tuple(Chain(c["node"], LETTERS[c["node"]]) for c in chains),
        bonds=bool(provenance["bridges"]),
        residues=span,
        sampling=int(provenance["sampling"]),
        representation=provenance["representation"],
    )


def _number(n: int) -> str:
    return f"{n:_}"


def _region(region) -> str:
    kept = "" if region.kept else ", kept=False"
    return f"Region({region.label!r}, {region.short!r}, {region.start}, {region.end}{kept})"


def row_source(target) -> str:
    """The `targets.py` row for [target], as the twenty's are written, with the
    decisions the rules leave to a person marked `# HAND:`."""
    source, structure = target.source, target.structure
    where = [repr(source.accession)]
    if source.seq_start is not None:
        where += [_number(source.seq_start), _number(source.seq_stop)]
    if source.protein_id:
        where.append(f"protein_id={source.protein_id!r}")
    if source.transcript_id:
        where.append(f"transcript_id={source.transcript_id!r}")
    lines = [
        "    Target(",
        f"        slug={target.slug!r}, gene={target.gene!r}, uniprot={target.uniprot!r}, "
        f"aa={target.aa},",
        f"        display={target.display!r},  # HAND: as curated/catalog.json says it",
        f"        source=Source({', '.join(where)}),",
        f"        cleaved={target.cleaved}, chain_label={target.chain_label!r},",
    ]
    if target.cleaved:
        lines.append("        # HAND: proprotein='...' where the record names none (R3.3)")
    lines.append(f"        mature_peptides={target.mature_peptides},  # HAND: True where the "
                 f"record's peptides are one cleavage series (R3.4)")
    if target.uniprot_variants:
        lines.append(f"        uniprot_variants={tuple(target.uniprot_variants)!r},")
    lines.append("        regions=(  # HAND: UniProt's processing; read against the record")
    lines += [f"            {_region(r)}," for r in target.regions]
    lines.append("        ),")
    if target.disulfides:
        lines.append(f"        disulfides={tuple(tuple(p) for p in target.disulfides)!r},")
    drawn = [f"pdb={structure.pdb!r}"]
    if structure.residues is not None:
        drawn.append(f"residues={tuple(structure.residues)!r}")
    if structure.sampling != 8:
        drawn.append(f"sampling={structure.sampling}")
    if structure.representation != "cartoon":
        drawn.append(f"representation={structure.representation!r}")
    chains = ", ".join(f"Chain({c.node!r}, {c.pdb_chain!r})" for c in structure.chains)
    lines += [
        "        structure=Structure(",
        f"            # Cut from structures/working/{target.slug}.pdb, a working copy of "
        f"{structure.pdb}.",
        f"            {', '.join(drawn)},",
        f"            chains=({chains},),",
        f"            bonds={structure.bonds},",
        "        ),",
        "    ),",
    ]
    return "\n".join(lines) + "\n"


def catalog_entry(target, record: dict, order: int, outcome: dict, summary: str) -> dict:
    """The `curated/catalog.json` entry: data as derived, words drafted.

    The facts are the record's as the walk derives them (`check_assets.check_facts`):
    they follow the record the row is baked from, so a change of
    `mature_peptides` changes `chains`, and `chain`, which is named exactly
    where the record names no peptide."""
    peptides = record.get("peptides") or []
    chrome = dict(outcome["chrome"])
    for word in ("label", "sentence", "semantics"):
        chrome[word] = f"{DRAFT} {chrome[word]}"
    return {
        "slug": target.slug,
        "catalog_order": order,
        "display": f"{DRAFT} {target.display}",
        "gene": target.gene,
        "uniprot": target.uniprot,
        "accession": target.source.accession,
        "summary": f"{DRAFT} {summary}",
        "chain": None if peptides else f"{DRAFT} {target.chain_label.lower()}",
        "facts": {
            "residues": len(record["protein"]["translation"]),
            "exons": len(record["exons"]),
            "chains": max(1, len(peptides)),
            "bridges": len(target.disulfides),
        },
        "chains": outcome["chains"],
        "structure": chrome,
        "tracks": TRACKS,
    }


def sidecar(outcome: dict, record: bytes, working_copy: str) -> dict:
    """Where a stored working copy came from, and what the app's check passed."""
    provenance, folding = outcome["provenance"], outcome["folding_provenance"]
    return {
        "entry": provenance["entry"],
        "model": provenance["model"],
        "method": provenance["method"],
        "resolution": provenance["resolution"],
        "revision": provenance["revision"],
        "file": provenance["file"],
        "policy": outcome["policy"],
        "candidate": outcome["candidate"],
        "scope": provenance["scope"],
        "domains": provenance["domains"],
        # The entry's own chains, and its authors' numbering where one offset
        # does for every residue: what the working copy's letters and record
        # numbering stand for.
        "chains": [{**chain, "offset": said.get("offset")}
                   for chain, said in zip(provenance["chains"], folding["chains"])],
        "bridges": provenance["bridges"],
        "bridges_dropped": provenance["bridges_dropped"],
        "secondary_structure": provenance["secondary_structure"],
        "representation": provenance["representation"],
        "sampling": provenance["sampling"],
        "residues": provenance["residues"],
        "sequence_differences": provenance["sequence_differences"],
        "sha256": {
            "record": hashlib.sha256(record).hexdigest(),
            "translation": hashlib.sha256(
                json.loads(record)["protein"]["translation"].encode()).hexdigest(),
            "working_copy": hashlib.sha256(working_copy.encode()).hexdigest(),
            **outcome["sha256"],
        },
        "measured": outcome["measured"],
        "drafted_by": "scripts/draft_curated.py",
    }


# ------------------------------------------------------------ the run


def _env_file(given: str | None) -> Path:
    for candidate in ([Path(given)] if given else [BACKEND / ".env"]):
        if candidate.exists():
            return candidate
    raise SystemExit("No .env with DATABASE_URL; pass --env.")


def _index(env: dict, specs: list[Spec]) -> tuple[dict, str | None]:
    """Each gene's buildable index row, and the MANE release, read only."""
    import psycopg

    from pipeline.resolver import store

    rows = {}
    with psycopg.connect(env["DATABASE_URL"], prepare_threshold=None) as conn:
        with conn.transaction():
            conn.execute("set transaction read only")
            mane = store.mane_release(conn)
            for spec in specs:
                hit = conn.execute(
                    "select uniprot from protein_index where gene = %s and buildable",
                    (spec.gene,)).fetchone()
                rows[spec.gene] = store.index_row(conn, hit[0], spec.gene) if hit else None
    return rows, mane


def _bake(args, target, record: bytes, protein: dict, folder: Path) -> dict:
    """The choice baked and held to the app's check, by `experimental.py`."""
    from pipeline.structure import experimental

    python = args.structure_python or BACKEND / "pipeline" / "structure" / "venv" / "bin" / "python"
    if not Path(python).exists():
        raise SystemExit(f"No structure bake environment at {python}; pass --structure-python.")
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    (folder / "record.json").write_bytes(record)
    job = {"target": experimental.target_json(target), "record": "record.json",
           "row": experimental.served_row(protein)}
    (folder / "job.json").write_text(json.dumps(job, ensure_ascii=False), encoding="utf-8")
    command = [str(python), "-u", str(BACKEND / "pipeline" / "structure" / "experimental.py"),
               "--job", str(folder / "job.json"), "--out", str(folder), "--app", str(args.app),
               "--describe", target.gene, "drafted"]
    for flag, value in (("--dart", args.dart), ("--flutter", args.flutter),
                        ("--limit", args.limit), ("--draws", args.draws)):
        if value:
            command += [flag, str(value)]
    ran = subprocess.run(command, cwd=str(BACKEND), stdin=subprocess.DEVNULL,
                         capture_output=True, text=True)
    (folder / "stdout.txt").write_text(ran.stdout)
    (folder / "stderr.txt").write_text(ran.stderr)
    said = folder / "outcome.json"
    if ran.returncode != 0 or not said.exists():
        raise RuntimeError(f"the bake exited {ran.returncode}: {' '.join(ran.stderr.split())[-600:]}")
    return json.loads(said.read_text(encoding="utf-8"))


def draft(args, spec: Spec, row, mane: str | None, order: int) -> tuple[bool, str]:
    """Draft one gene: whether it was drafted, and what was made of it, in a line."""
    from pipeline import uniprot
    from pipeline.paths import CURATED
    from pipeline.resolver.resolve import Refused, resolve
    from pipeline.structure import pdb

    if row is None:
        return False, f"{spec.gene}: no buildable index entry"
    curated = {p["slug"] for p in json.loads(CURATED.read_text())["proteins"]}
    if spec.slug in curated:
        return False, f"{spec.gene}: {spec.slug} is already a curated row"
    copy = pdb.WORKING / f"{spec.slug}.pdb"
    if copy.exists() and not args.replace:
        return False, f"{spec.gene}: {copy} is stored already (--replace to draft it again)"
    try:
        resolution = resolve(row, uniprot.fetch_entry(row.uniprot), mane_release=mane)
    except Refused as exc:
        return False, f"{spec.gene}: the resolver refuses it: {exc}"
    target = replace(resolution.target, slug=spec.slug)
    protein = {**resolution.protein, "slug": spec.slug}
    folder = args.out / spec.slug
    outcome = _bake(args, target, resolution.record, protein, folder / "bake")
    if outcome["outcome"] != "chosen":
        return False, f"{spec.gene}: the search ended {outcome['outcome']!r}; nothing drafted"
    entry = outcome["provenance"]["pdb"]
    if spec.entry and entry != spec.entry:
        return False, f"{spec.gene}: chose {entry}, not {spec.entry} as asked; nothing drafted"

    working = (folder / "bake" / "working_copy.pdb").read_text(encoding="utf-8")
    drawn = replace(target, structure=structure_of(outcome["provenance"]))
    record = json.loads(resolution.record)
    entry_json = catalog_entry(drawn, record, order, outcome, protein["summary"])
    row_text = row_source(drawn)
    pdb.WORKING.mkdir(parents=True, exist_ok=True)
    copy.write_text(working, encoding="utf-8")
    (pdb.WORKING / f"{spec.slug}.json").write_text(
        json.dumps(sidecar(outcome, resolution.record, working), ensure_ascii=False, indent=1)
        + "\n", encoding="utf-8")
    (folder / "row.py").write_text(row_text, encoding="utf-8")
    (folder / "entry.json").write_text(
        json.dumps(entry_json, ensure_ascii=False, indent=2, separators=(",", ": ")) + "\n",
        encoding="utf-8")
    measured = outcome["measured"]
    return True, (f"{spec.gene}: {spec.slug} at order {order}, {entry} "
            f"({outcome['chrome']['sentence']}) scene {measured.get('scene_bytes', 0):,} B; "
            f"working copy {copy.relative_to(BACKEND)}")


def main(arguments: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("genes", nargs="+", metavar="GENE=slug[@ENTRY]")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--env", default=None, help="the .env with DATABASE_URL (this checkout's)")
    parser.add_argument("--order", type=int, default=None,
                        help="the first new catalog_order (default: after the curated rows)")
    parser.add_argument("--structure-python", type=Path, default=None)
    parser.add_argument("--app", type=Path, default=BACKEND.parent / "helix-peek")
    parser.add_argument("--dart", default=None)
    parser.add_argument("--flutter", default=None)
    parser.add_argument("--limit", type=int, default=None, help="entry files examined at most")
    parser.add_argument("--draws", type=int, default=None, help="candidates baked at most")
    parser.add_argument("--replace", action="store_true",
                        help="draft again over a stored working copy")
    args = parser.parse_args(arguments)
    try:
        specs = [parse_spec(text) for text in args.genes]
    except ValueError as exc:
        parser.error(str(exc))
    if len({s.slug for s in specs}) != len(specs):
        parser.error("two genes asked for one slug")

    from dotenv import dotenv_values

    env = dotenv_values(_env_file(args.env))
    email = os.environ.get("NCBI_EMAIL") or env.get("NCBI_EMAIL")
    if not email:
        raise SystemExit("Set NCBI_EMAIL: the record is fetched from NCBI.")
    # Named before any baker is imported: each reads its directory once, then.
    scratch = Path(tempfile.mkdtemp(prefix="helixpeek-draft-"))
    atexit.register(shutil.rmtree, scratch, True)
    os.environ["NCBI_EMAIL"] = email
    os.environ["HELIXPEEK_GB_CACHE"] = str(scratch / "genbank")
    os.environ["HELIXPEEK_DATA"] = str(scratch / "data")

    from pipeline.paths import CURATED

    orders = [p["catalog_order"] for p in json.loads(CURATED.read_text())["proteins"]]
    first = args.order if args.order is not None else max(orders) + 1
    rows, mane = _index(env, specs)
    args.out.mkdir(parents=True, exist_ok=True)
    made = 0
    for order, spec in enumerate(specs, first):
        print(f"{spec.gene} -> {spec.slug}", flush=True)
        try:
            drafted, said = draft(args, spec, rows[spec.gene], mane, order)
        except Exception as exc:  # noqa: BLE001 -- said, and the others go on
            drafted, said = False, f"{spec.gene}: {type(exc).__name__}: {exc}"
        print(f"  {said}", flush=True)
        folder = args.out / spec.slug
        if drafted:
            made += 1
            print((folder / "row.py").read_text(), flush=True)
            print((folder / "entry.json").read_text(), flush=True)
    print(f"{made} of {len(specs)} drafted into {args.out}")
    return 0 if made == len(specs) else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
