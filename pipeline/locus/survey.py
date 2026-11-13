"""The census: every protein the search offers to build, its locus baked from
the pinned snapshot and checked. Writes only under `--out`.

    .venv/bin/python -I pipeline/locus/survey.py --live --out OUT [--save-dump DIR]
    .venv/bin/python -I pipeline/locus/survey.py --dump DIR --out OUT
        [--gene-info Homo_sapiens.gene_info.gz]

`--live` reads the protein index through `.env`'s DATABASE_URL in one
transaction set read only; `--dump` reads a directory a live run saved
(`--save-dump`): `protein_index.json`, `protein_index_release.json` and
`protein.json`. Every buildable row is baked exactly as the resolver bakes a
protein it builds (`bake_locus.bake`: the row and the snapshot, nothing
fetched), held to `check_locus` with the record the resolver would cut for it
(`resolve.source_for`), and, with `--gene-info`, to NCBI Gene's place
(`verify_locus`). A curated protein keeps its own slug; any other is named as
the resolver names it.

Under OUT: `locus/<slug>_locus.json` for each track, `refused.json`,
`failed.json`, and `summary.json`, whose numbers it also prints. With
`--gff` (MANE v1.5's genomic GFF, `MANE_GFF`), also `records/<slug>_record.json`:
a record per protein in the record track's own shape, its exons MANE's and
its bases made up, for the app's check of every protein's zoom
(`helix-peek/test/zoom_catalogue_check.dart`).
"""

from __future__ import annotations

import argparse
import collections
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.locus import bake_locus, check_locus, verify_locus  # noqa: E402
from pipeline.locus import inputs as locus_inputs  # noqa: E402
from pipeline.locus.tissues import pair_tissue as tissue_of_pair  # noqa: E402
from pipeline.resolver.resolve import IndexRow, Refused, slug_for, source_for  # noqa: E402
from pipeline.targets import Target  # noqa: E402

INDEX_COLUMNS = ("uniprot", "gene", "name", "length", "hgnc_id", "refseq_nuc", "refseq_prot",
                 "ensembl_nuc", "mane_isoform", "chrom_acc", "chrom_start", "chrom_end",
                 "chrom_strand", "refseqgene", "buildable", "unavailable_reason")


def read_live(env: Path) -> tuple[list, dict, list]:
    """The index's rows, its release row and the protein table's names, read in
    one transaction set read only."""
    import psycopg
    from dotenv import dotenv_values

    url = dotenv_values(env)["DATABASE_URL"]
    with psycopg.connect(url, prepare_threshold=None) as conn:
        with conn.transaction():
            conn.execute("set transaction read only")
            rows = conn.execute(
                f"select {', '.join(INDEX_COLUMNS)} from protein_index").fetchall()
            cursor = conn.execute("select * from protein_index_release")
            names = [column.name for column in cursor.description]
            release = dict(zip(names, cursor.fetchone()))
            proteins = conn.execute(
                "select slug, gene, uniprot, catalog_order from protein").fetchall()
    return ([dict(zip(INDEX_COLUMNS, row)) for row in rows],
            {key: str(value) for key, value in release.items()},
            [dict(zip(("slug", "gene", "uniprot", "catalog_order"), row)) for row in proteins])


def read_dump(folder: Path) -> tuple[list, dict, list]:
    release = json.loads((folder / "protein_index_release.json").read_text())
    return (json.loads((folder / "protein_index.json").read_text()),
            release[0] if isinstance(release, list) else release,
            json.loads((folder / "protein.json").read_text()))


def save_dump(folder: Path, rows: list, release: dict, proteins: list) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "protein_index.json").write_text(json.dumps(rows))
    (folder / "protein_index_release.json").write_text(json.dumps([release]))
    (folder / "protein.json").write_text(json.dumps(proteins))


MANE_GFF = ("https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/"
            "MANE.GRCh38.v1.5.refseq_genomic.gff.gz")

# A RefSeq transcript as MANE's GFF files a copy of it: `NM_000937.5-2`.
_COPY = re.compile(r"^(N[MR]_\d+\.\d+)-\d+$")


def _accession(named: str) -> str:
    """The transcript's accession, without the number of the copy."""
    found = _COPY.match(named)
    return found.group(1) if found else named


def read_gff(path: Path) -> dict:
    """MANE's transcripts, by RefSeq accession: the mRNA's span and strand,
    and its exons and its coding stretches, as the genome numbers them."""
    import gzip

    found: dict = {}
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 9 or columns[2] not in ("mRNA", "exon", "CDS"):
                continue
            attributes = dict(item.split("=", 1) for item in columns[8].split(";")
                              if "=" in item)
            named = attributes.get("ID" if columns[2] == "mRNA" else "Parent", "")
            if not named.startswith("rna-"):
                continue
            # A transcript MANE annotates on a patch scaffold is filed as a
            # copy of its accession (`rna-NM_000937.5-2`).
            entry = found.setdefault(_accession(named[4:]), {"exons": [], "cds": []})
            start, end = int(columns[3]), int(columns[4])
            if columns[2] == "mRNA":
                entry.update(start=start, end=end, strand=1 if columns[6] == "+" else -1)
            else:
                entry["exons" if columns[2] == "exon" else "cds"].append((start, end))
    return found


# A run of bases, made up: the census's records are drawn at the DNA's stop,
# whose letters are only ever the record's own.
_BASES = "ACGTTGCAAGCTGATCCTAGGATCCAGTCAGTTCAAGCTTGGCCATGCATCGATCGTACGTAAGGCTTACGG"


def record_of(row: dict, transcript: dict) -> dict:
    """A record for the app's check of the census, in the record track's own
    shape (`mock/build_gene_record.py`'s): MANE Select's exons and coding
    stretches on the gene's own strand, numbered from its 5′ end, its bases
    made up, and its introns shortened past the gene page's budget exactly as
    the builder shortens a real record's (`intron_scale`, `compress`). Raises
    ValueError where its exons alone are over the budget, as the builder does."""
    from pipeline.mock import build_gene_record as builder
    from pipeline.targets import GENE_PAGE_BUDGET_BP

    start, end, strand = transcript["start"], transcript["end"], transcript["strand"]

    def own(a: int, b: int) -> tuple[int, int]:
        return (a - start + 1, b - start + 1) if strand > 0 else (end - b + 1, end - a + 1)

    exons = sorted(own(a, b) for a, b in transcript["exons"])
    coding = sorted(own(a, b) for a, b in transcript["cds"])
    span = end - start + 1
    payload = {
        "gene": row["gene"],
        "location": {"start": 1, "end": span, "strand": 1},
        "sequence": (_BASES * (span // len(_BASES) + 1))[:span],
        "exons": [{"number": n, "start": a, "end": b} for n, (a, b) in enumerate(exons, 1)],
        "transcript": {"segments": [{"start": a, "end": b} for a, b in exons]},
        "protein": {"segments": [{"start": a, "end": b} for a, b in coding],
                    "translation": "M" * max(1, row.get("length") or 1),
                    "product": row.get("name")},
        "peptides": [],
    }
    scale = builder.intron_scale(payload, GENE_PAGE_BUDGET_BP)
    return builder.compress(payload, scale) if scale < 1.0 else payload


def target_for(row: dict, slug: str) -> Target:
    """The `Target` the resolver makes for [row], as far as the check reads it:
    its names and the record it would cut (`source_for`)."""
    index = IndexRow(**{key: row.get(key) for key in IndexRow.__dataclass_fields__})
    return Target(slug=slug, gene=row["gene"], uniprot=row["uniprot"], display=row["name"],
                  source=source_for(index), structure=None, cleaved=False, chain_label="")


def survey(rows: list, release: dict, proteins: list, snapshot: locus_inputs.Inputs,
           out: Path, gff: dict | None = None) -> dict:
    """Bake and check every buildable row; write each track under [out], and
    with [gff] (`read_gff`) a record for each, for the app's check."""
    curated = {p["gene"]: p["slug"] for p in proteins if p.get("catalog_order") is not None}
    offered = [row for row in rows if row["buildable"]]
    folder = out / "locus"
    folder.mkdir(parents=True, exist_ok=True)
    uniprot_source = {"release": release["uniprot"]}
    tracks, refused, failed, problems = {}, [], [], []
    records, unrecorded = 0, []
    if gff is not None:
        (out / "records").mkdir(parents=True, exist_ok=True)
    for row in offered:
        slug = curated.get(row["gene"]) or slug_for(row["gene"])
        try:
            target = target_for(row, slug)
            track = bake_locus.bake(row, snapshot, slug=slug, uniprot_source=uniprot_source,
                                    mane_release=release["mane"])
        except bake_locus.LocusRefused as exc:
            refused.append({"slug": slug, "gene": row["gene"], "uniprot": row["uniprot"],
                            "reason": str(exc)})
            continue
        except (Refused, LookupError, ValueError, KeyError) as exc:
            failed.append({"slug": slug, "gene": row["gene"],
                           "error": f"{type(exc).__name__}: {exc}"})
            continue
        problems += check_locus.problems_of(target, track)
        (folder / f"{slug}_locus.json").write_bytes(bake_locus.encode(track))
        tracks[slug] = track
        if gff is not None:
            transcript = gff.get(row["refseq_nuc"] or "")
            if transcript is None:
                unrecorded.append((slug, "MANE's GFF has no transcript " + str(row["refseq_nuc"])))
                continue
            try:
                record = record_of(row, transcript)
            except ValueError as exc:
                unrecorded.append((slug, str(exc)))
                continue
            (out / "records" / f"{slug}_record.json").write_text(
                json.dumps(record, separators=(",", ":")))
            records += 1
    (out / "refused.json").write_text(json.dumps(refused, indent=1))
    (out / "failed.json").write_text(json.dumps(failed, indent=1))
    return {"offered": len(offered), "curated": sorted(curated.values()), "tracks": tracks,
            "refused": refused, "failed": failed, "problems": problems,
            "records": records, "unrecorded": unrecorded}


def _tied(levels: dict) -> bool:
    """Whether the two highest of [levels] are level, and above zero."""
    top = sorted(levels.values(), reverse=True)[:2]
    return len(top) == 2 and top[0] == top[1] and top[0] > 0


def ties_of(track: dict) -> list[str]:
    """Where rule v3's choice for [track] was settled by its tie-break."""
    expression, path = track["expression"], track["path"]
    found = []
    if path["tissue_from"] == "tissue_level" and _tied(expression["tissue_levels"]):
        found.append("tissue")
    if path["cell_from"] == "single_cell_level" and _tied(expression.get("cell_levels") or {}):
        found.append("cell")
    if path["tissue_from"] == "tissue_cell_type" and expression["tissue_cell_type"]:
        levels = expression["tissue_levels"]
        by_tissue = {}
        for pair in expression["tissue_cell_type"]:
            tissue = tissue_of_pair(pair["tissue"])
            by_tissue[tissue] = levels.get(tissue, 0)
        if len(by_tissue) > 1 and _tied(by_tissue):
            found.append("pair")
    return found


def summary_of(found: dict) -> dict:
    """The census's numbers, in the shape of the handoff's findings."""
    tracks = found["tracks"]
    count = collections.Counter
    paths = [t["path"] for t in tracks.values()]
    pairs = {(p["tissue"], p["cell_type"]) for p in paths if p["tissue"]}
    v2 = {slug: bake_locus.path_of_v2(t["expression"]) for slug, t in tracks.items()}
    changed = sorted(slug for slug, t in tracks.items() if t["path"] != v2[slug])
    curated = set(found["curated"])
    ties = {slug: ties_of(t) for slug, t in tracks.items()}
    return {
        "records": found.get("records", 0),
        "unrecorded": found.get("unrecorded", []),
        "v3_changed": len(changed),
        "v3_changed_from": dict(count(str(v2[slug]["tissue_from"]) for slug in changed)),
        "v3_changed_matched_by_mane": sum(
            tracks[slug]["expression"]["atlas_match"] == "mane" for slug in changed),
        "v3_changed_organ": sum(tracks[slug]["path"]["tissue"] != v2[slug]["tissue"]
                                for slug in changed),
        "v3_changed_curated": sorted(slug for slug in changed if slug in curated),
        "v3_ties": dict(count(kind for found_ties in ties.values() for kind in found_ties)),
        "v3_tied_genes": sorted(tracks[slug]["gene"] for slug, t in ties.items() if t),
        "offered": found["offered"],
        "tracks": len(tracks),
        "refused": [(r["gene"], r["reason"]) for r in found["refused"]],
        "failed": found["failed"],
        "check_problems": found["problems"],
        "chromosomes": dict(sorted(count(t["chromosome"] for t in tracks.values()).items())),
        "bands_per_gene": dict(count(len(t["band"]["names"]) for t in tracks.values())),
        "placed": sorted((t["gene"], t["locus"], t["placed_from"]["ucsc"])
                         for t in tracks.values() if "placed_from" in t),
        "atlas_match": dict(count(t["expression"]["atlas_match"] for t in tracks.values())),
        "atlas_named_otherwise": sorted(
            (t["gene"], t["expression"]["atlas_gene"]) for t in tracks.values()
            if t["expression"]["atlas_gene"] != t["gene"]),
        "tissue_from": dict(count(str(p["tissue_from"]) for p in paths)),
        "cell_from": dict(count(str(p["cell_from"]) for p in paths)),
        "no_tissue": sorted(t["gene"] for t in tracks.values() if t["path"]["tissue"] is None),
        "no_cell": sum(p["cell_type"] is None for p in paths),
        "lands_in": dict(count((p["lands_in"] or {}).get("cell") for p in paths if p["lands_in"])),
        "organs": len({p["tissue"] for p in paths if p["tissue"]}),
        "organ_cell_pairs": len(pairs),
        "top_organs": count(p["tissue"] for p in paths).most_common(8),
        "spans": sorted(t["span"]["end"] - t["span"]["start"] + 1 for t in tracks.values())[::max(
            1, len(tracks) // 4)] if tracks else [],
    }


def main(arguments: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--live", action="store_true")
    source.add_argument("--dump", type=Path)
    parser.add_argument("--save-dump", type=Path, default=None)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--gene-info", type=Path, default=None)
    parser.add_argument("--gff", type=Path, default=None,
                        help="MANE v1.5's genomic GFF: writes a record for the app's check")
    args = parser.parse_args(arguments)

    rows, release, proteins = (read_live(BACKEND / ".env") if args.live
                               else read_dump(args.dump))
    if args.save_dump:
        save_dump(args.save_dump, rows, release, proteins)
    snapshot = locus_inputs.load()
    args.out.mkdir(parents=True, exist_ok=True)
    found = survey(rows, release, proteins, snapshot, args.out,
                   read_gff(args.gff) if args.gff else None)
    summary = {"index": release, "locus_inputs": snapshot.sha256, **summary_of(found)}
    if args.gene_info:
        verdicts = verify_locus.verify(list(found["tracks"].values()), snapshot,
                                       verify_locus.read_gene_info(args.gene_info))
        summary["ncbi"] = dict(verdicts["classes"])
        summary["ncbi_differences"] = verdicts["rows"]
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
    for key, value in summary.items():
        if key in ("ncbi_differences", "atlas_named_otherwise", "placed", "check_problems",
                   "v3_tied_genes", "no_tissue", "unrecorded"):
            print(f"{key}: {len(value)}")
            for item in value[:40]:
                print(f"    {item}")
        else:
            print(f"{key}: {value}")
    bad = summary["failed"] or summary["check_problems"] or (summary.get("ncbi") or {}).get("wrong")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
