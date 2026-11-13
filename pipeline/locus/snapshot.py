"""Build the locus bake's pinned inputs: one snapshot of every source it reads.

    .venv/bin/python -I pipeline/locus/snapshot.py
    .venv/bin/python -I pipeline/locus/snapshot.py --work /tmp/x --keep

Downloads every source from a pinned URL into a directory of its own
(`--work`, else a new temporary one, removed after), reads them, and writes
`pipeline/data/locus_inputs/locus-inputs.<sha12>.jsonl.gz` (`--out`), where
the bake, the census and the resolver read it (`inputs.py`). It prints the
sha256 of the snapshot's text and says whether it is the one `inputs.py`
pins. Nothing is uploaded.

The sources, each pinned to a release that no longer changes:

- the Human Protein Atlas 25.1, from its own archive site
  (`v25.proteinatlas.org`): every gene's summary, its consensus RNA levels by
  tissue and its single cell type levels, and the release line its releases
  page names (`bake_locus.atlas_of`). The live site serves 26.0;
- UCSC's hg38 cytoBand, chromAlias, altLocations and fixLocations, through its
  REST API as `bake_locus.fetch_table` reads them, with the time UCSC last
  updated each (cytoBand 2022-10-28);
- MANE v1.5's summary, from its release directory (not `current/`);
- HGNC's complete set as its monthly archive kept it on 2026-10-06 (the live
  file changes every week; this one gives every offered gene the Ensembl gene
  the live one did on 2026-11-13);
- NCBI's placements of GRCh38.p14's patch and alternate locus scaffolds on
  their chromosomes (`all_alt_scaffold_placement.txt`), p14 being the last
  patch release GRCh38 will have.

The text is made the same way every time: genes in Ensembl order, keys
sorted, levels as numbers, and the date it is made on given (`--date`, the
pinned `INPUTS_CREATED` unless a new snapshot is meant), so the same sources
give the same digest. Each download's own sha256 is in the header.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from app.protein_index import parse_mane  # noqa: E402
from pipeline.locus import inputs  # noqa: E402
from pipeline.locus.bake_locus import (  # noqa: E402
    HPA_LICENCE,
    Table,
    _get,
    atlas_of,
    fetch_table,
)

HPA_VERSION = "25.1"
HPA_SITE = "https://v25.proteinatlas.org"
HPA_RELEASES = HPA_SITE + "/about/releases"
HPA_FILES = {
    "proteinatlas.json.gz": HPA_SITE + "/download/proteinatlas.json.gz",
    "rna_tissue_consensus.tsv.zip": HPA_SITE + "/download/tsv/rna_tissue_consensus.tsv.zip",
    "rna_single_cell_type.tsv.zip": HPA_SITE + "/download/tsv/rna_single_cell_type.tsv.zip",
}

MANE_RELEASE = "v1.5"
MANE_SUMMARY = ("https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/"
                "MANE.GRCh38.v1.5.summary.txt.gz")

HGNC_DATE = "2026-10-06"
HGNC = ("https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/"
        f"tsv/hgnc_complete_set_{HGNC_DATE}.txt")

PLACEMENT_ASSEMBLY = "GRCh38.p14"
PLACEMENT_ACCESSION = "GCF_000001405.40"
PLACEMENT = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/"
             "GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/"
             "all_alt_scaffold_placement.txt")

UCSC_TABLES = {"cytoBand": "cytoband", "chromAlias": "chrom_alias",
               "altLocations": "alt_locations", "fixLocations": "fix_locations"}

# The fields of an Atlas gene the snapshot keeps: its names, and every field
# `bake_locus.expression_of` reads.
ATLAS_FIELDS = (
    "Gene", "Ensembl", "Uniprot",
    "RNA tissue specificity", "RNA tissue distribution", "RNA tissue specific nTPM",
    "RNA single cell type specificity", "RNA single cell type distribution",
    "RNA single cell type specific nCPM",
    "RNA tissue cell type enrichment",
    "Subcellular main location", "Subcellular additional location",
    "Secretome location",
)


def _number(text: str):
    """A level as a number: a whole one without its point."""
    value = float(text)
    return int(value) if value.is_integer() else value


def _levels(tsv_zip: Path, key: str, unit: str) -> tuple[list[str], dict[str, dict]]:
    """An Atlas level table (`Gene`, `Gene name`, [key], [unit]): its names,
    sorted, and each gene's level by name."""
    by_gene: dict[str, dict] = {}
    with zipfile.ZipFile(tsv_zip) as archive:
        names = archive.namelist()
        if len(names) != 1:
            raise ValueError(f"{tsv_zip.name} holds {names}, not one table")
        with archive.open(names[0]) as raw:
            for row in csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8"), delimiter="\t"):
                by_gene.setdefault(row["Gene"], {})[row[key]] = _number(row[unit])
    names = sorted({name for levels in by_gene.values() for name in levels})
    return names, by_gene


def _hgnc(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8") as stream:
        return {row["hgnc_id"]: row["ensembl_gene_id"]
                for row in csv.DictReader(stream, delimiter="\t")
                if row.get("ensembl_gene_id")}


def _placements(path: Path) -> dict[str, list[dict]]:
    """NCBI's placement rows, by the scaffold's RefSeq accession."""
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].lstrip("#").split("\t")
    found: dict[str, list[dict]] = {}
    for line in lines[1:]:
        values = line.split("\t")
        if len(values) != len(header):
            continue
        row = dict(zip(header, values))
        found.setdefault(row["alt_scaf_acc"], []).append(row)
    return found


def _file(path: Path, url: str) -> dict:
    payload = path.read_bytes()
    return {"url": url, "sha256": hashlib.sha256(payload).hexdigest(), "bytes": len(payload)}


def snapshot_text(*, folder: Path, release_page: str, tables: dict[str, Table],
                  created: str) -> bytes:
    """The snapshot's text, from the downloads in [folder] (named as
    `HPA_FILES`, `mane_summary.txt.gz`, `hgnc_complete_set.txt` and
    `alt_scaffold_placement.txt`), the Atlas's releases page and the UCSC
    tables as read."""
    atlas = atlas_of(release_page)
    if atlas.version != HPA_VERSION:
        raise ValueError(f"the Atlas site names release {atlas.version}, not {HPA_VERSION}")
    with gzip.open(folder / "proteinatlas.json.gz", "rt", encoding="utf-8") as stream:
        genes = json.load(stream)
    tissue_names, tissues = _levels(folder / "rna_tissue_consensus.tsv.zip", "Tissue", "nTPM")
    cell_names, cells = _levels(folder / "rna_single_cell_type.tsv.zip", "Cell type", "nCPM")
    mane = parse_mane(io.StringIO(
        gzip.decompress((folder / "mane_summary.txt.gz").read_bytes()).decode("utf-8")))
    hgnc = _hgnc(folder / "hgnc_complete_set.txt")
    placements = _placements(folder / "alt_scaffold_placement.txt")

    lines: list[str] = []
    uniprot: dict[str, list[str]] = {}
    for gene in sorted(genes, key=lambda found: found["Ensembl"]):
        ensembl = gene["Ensembl"]
        for accession in gene.get("Uniprot") or []:
            uniprot.setdefault(accession, []).append(ensembl)
        entry = {
            "atlas": {name: gene.get(name) for name in ATLAS_FIELDS},
            "tissues": [tissues.get(ensembl, {}).get(name, 0) for name in tissue_names],
            "cells": [cells.get(ensembl, {}).get(name, 0) for name in cell_names],
        }
        lines.append(json.dumps([ensembl, entry], sort_keys=True, separators=(",", ":"),
                                ensure_ascii=False))

    files = {name: _file(folder / name, url) for name, url in HPA_FILES.items()}
    header = {
        "kind": inputs.KIND,
        "format": inputs.FORMAT,
        "created": created,
        "sources": {
            "hpa": {
                "version": atlas.version,
                "release_date": atlas.release_date,
                "ensembl": atlas.ensembl,
                "licence": HPA_LICENCE,
                "site": HPA_SITE,
                "releases": HPA_RELEASES,
                "genes": len(genes),
                "files": files,
            },
            **{key: tables[track].provenance() for track, key in UCSC_TABLES.items()},
            "mane": {"release": MANE_RELEASE, "select_rows": len(mane),
                     **_file(folder / "mane_summary.txt.gz", MANE_SUMMARY)},
            "hgnc": {"date": HGNC_DATE, "ids": len(hgnc),
                     **_file(folder / "hgnc_complete_set.txt", HGNC)},
            "placement": {"assembly": PLACEMENT_ASSEMBLY, "accession": PLACEMENT_ACCESSION,
                          "scaffolds": len(placements),
                          **_file(folder / "alt_scaffold_placement.txt", PLACEMENT)},
        },
        "tables": {track: {"genome": table.genome, "updated": table.updated, "url": table.url,
                           "rows": list(table.rows)}
                   for track, table in tables.items()},
        "placements": placements,
        "mane": mane,
        "hgnc": hgnc,
        "uniprot": uniprot,
        "consensus_tissues": tissue_names,
        "cell_types": cell_names,
    }
    head = json.dumps(header, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return ("\n".join([head, *lines]) + "\n").encode("utf-8")


def write(text: bytes, folder: Path) -> tuple[Path, str]:
    """The snapshot written, gzipped with no name or time in its header, and
    the sha256 of its text."""
    digest = hashlib.sha256(text).hexdigest()
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / inputs.file_name(digest)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0,
                           compresslevel=9) as stream:
            stream.write(text)
    return path, digest


def download(folder: Path) -> str:
    """Every file source into [folder], under the names `snapshot_text`
    reads; the Atlas's releases page as text."""
    for name, url in {**HPA_FILES,
                      "mane_summary.txt.gz": MANE_SUMMARY,
                      "hgnc_complete_set.txt": HGNC,
                      "alt_scaffold_placement.txt": PLACEMENT}.items():
        body, _ = _get(url)
        (folder / name).write_bytes(body)
        print(f"  {name:<32} {len(body):>12,} B  {hashlib.sha256(body).hexdigest()[:12]}")
    return _get(HPA_RELEASES)[0].decode("utf-8", "replace")


def main(arguments: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--work", type=Path, default=None,
                        help="where the downloads go (default: a new temporary directory)")
    parser.add_argument("--keep", action="store_true", help="keep the downloads")
    parser.add_argument("--out", type=Path, default=inputs.FOLDER)
    parser.add_argument("--date", default=inputs.INPUTS_CREATED,
                        help=f"the snapshot's date (default the pinned {inputs.INPUTS_CREATED})")
    args = parser.parse_args(arguments)

    work = args.work or Path(tempfile.mkdtemp(prefix="helixpeek-locus-inputs-"))
    work.mkdir(parents=True, exist_ok=True)
    try:
        page = download(work)
        tables = {track: fetch_table(track) for track in UCSC_TABLES}
        for track, table in tables.items():
            print(f"  {track:<32} {len(table.rows):>12,} rows, updated {table.updated}")
        text = snapshot_text(folder=work, release_page=page, tables=tables, created=args.date)
    finally:
        if not args.keep and args.work is None:
            shutil.rmtree(work, ignore_errors=True)
    path, digest = write(text, args.out)
    pinned = digest == inputs.INPUTS_SHA256
    print(f"\n{path}\n  text {len(text):,} B, sha256 {digest}\n  "
          + ("the pinned snapshot" if pinned else
             f"NOT the pinned one ({inputs.INPUTS_SHA256[:12]}): pin it in inputs.py "
             f"only if a new snapshot is meant"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
