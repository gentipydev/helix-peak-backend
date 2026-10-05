"""Bake `locus`: where on its chromosome each protein's gene lies, by band.

    .venv/Scripts/python pipeline/locus/bake_locus.py --all
    .venv/Scripts/python pipeline/locus/bake_locus.py --target insulin

The app's zoom ends on a chromosome with the gene marked at its band. This
track carries what that needs, per protein:

- `span`: where the gene lies on GRCh38, 1-based and inclusive, and its
  strand. It is the MANE Select transcript's span, found the way the protein
  index finds it: the protein's UniProt entry joined to MANE by
  `rows_for_entry` (`app/protein_index.py`), the row for the table's gene.
- `bands`: every cytogenetic band of that chromosome, in order from the end
  of the short arm, each with its Giemsa stain, from the UCSC Genome Browser's
  cytoBand table for hg38. UCSC's 0-based starts are made 1-based, so every
  position in the track is counted the way MANE and GenBank count.
- `band`: the band or bands the span lies in, and where they begin and end;
  and `locus`, the place as cytogenetics writes it: 11p15.5, or Xp21.2-p21.1
  for a gene that runs across two bands.
- `expression`: where in the body the gene is read, as the Human Protein
  Atlas measures it: how specific its RNA is to a tissue and to a cell type,
  in the Atlas's own categories, and the tissues and cell types it names,
  highest first. Since schema 2, also the Atlas's tissue cell type pairs
  (the cell types it finds the gene enriched in within a tissue), where in a
  cell it finds the protein, and where it is secreted to.
- `path` (schema 2): the one way down the zoom takes for this gene, chosen
  by one rule from those readings (`path_of`): the organ, the kind of cell
  in it, and where that cell has no nucleus, the cell the zoom lands in
  instead. Every cell is one that lives in the organ, so no gene's zoom
  lands in a cell of another organ.

A band is a stain pattern seen down a microscope at low resolution, millions
of bases long. It says where a gene lies, not that the gene can be seen there.

Each payload names its sources and their versions: the cytoBand and
chromAlias tables by genome and by the time UCSC's API says each was last
updated, with a digest of the rows read; the MANE release; the UniProt release
and entry version; and the Protein Atlas version, its release date and the
Ensembl version it is built on, with its licence, CC BY 4.0. UCSC's download server does not answer from every network,
so the tables are read through its REST API, which serves the same tables.

Which proteins get the track: all twenty (`LOCUS_TARGETS`), recorded here,
never in targets.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from app.protein_index import rows_for_entry, versionless  # noqa: E402
from pipeline.locus.cell_types import anucleate, cell_class, homes, single  # noqa: E402
from pipeline.locus.tissues import pair_tissue, tissue_name  # noqa: E402
from pipeline.mane import current_summary  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import BY_SLUG, TARGETS, Target  # noqa: E402
from pipeline.uniprot import Entry, fetch_entry  # noqa: E402

# Every protein gets the track: each has a chromosome to zoom to.
LOCUS_TARGETS = TARGETS

# 2 adds keys only: the Atlas's tissue cell type pairs, subcellular location
# and secretome location under `expression`, and the zoom's `path`. A schema 1
# reader ignores them.
SCHEMA_VERSION = 2

ASSEMBLY = "GRCh38"
GENOME = "hg38"
UCSC_API = "https://api.genome.ucsc.edu/getData/track?genome={};track={}"

_AGENT = "helixpeek-locus-bake"

# The Human Protein Atlas: one gene's summary, by Ensembl gene, and the page
# that names the release it is from.
HPA_GENE = "https://www.proteinatlas.org/{}.json"
HPA_RELEASES = "https://www.proteinatlas.org/about/releases"
HPA_LICENCE = "CC BY 4.0"

# How specific a gene's RNA is, in the Atlas's categories: to a tissue, and to
# a single cell type. The last two name nothing specific.
TISSUE_SPECIFICITY = ("Tissue enriched", "Group enriched", "Tissue enhanced",
                      "Low tissue specificity", "Not detected")
CELL_SPECIFICITY = ("Cell type enriched", "Group enriched", "Cell type enhanced",
                    "Low cell type specificity", "Not detected")
UNSPECIFIC = ("Low tissue specificity", "Low cell type specificity", "Not detected")

# The stains cytoBand gives a band: G-negative, four depths of G-positive,
# the centromere, the variable heterochromatin, and the stalks of the
# acrocentric chromosomes.
STAINS = ("gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar", "stalk")


@dataclass(frozen=True)
class Table:
    """One UCSC table as its API served it."""

    track: str
    genome: str
    # When UCSC last updated the table, as its API reports it: the table's
    # version, since UCSC tables carry no other.
    updated: str
    rows: tuple[dict, ...]
    url: str

    @property
    def sha256(self) -> str:
        """A digest of the rows, in a canonical order: the same table gives the
        same digest however the API happens to order its answer."""
        canonical = sorted(json.dumps(row, sort_keys=True) for row in self.rows)
        return hashlib.sha256("\n".join(canonical).encode()).hexdigest()

    def provenance(self) -> dict:
        return {
            "table": self.track,
            "genome": self.genome,
            "updated": self.updated,
            "rows": len(self.rows),
            "sha256": self.sha256,
            "url": self.url,
        }


@dataclass(frozen=True)
class Atlas:
    """The Human Protein Atlas release the genes were read from."""

    version: str
    release_date: str
    ensembl: str

    def provenance(self, gene: str) -> dict:
        return {
            "version": self.version,
            "release_date": self.release_date,
            "ensembl": self.ensembl,
            "licence": HPA_LICENCE,
            "url": HPA_GENE.format(gene),
        }


_RELEASE = re.compile(
    r"Protein Atlas version (\d+\.\d+) Release date: (\d{4})\.(\d{2})\.(\d{2}) "
    r"Ensembl version: (\d+)")


def atlas_of(page: str) -> Atlas:
    """The latest release its releases page lists: the first it names."""
    text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", page))
    found = _RELEASE.search(text)
    if not found:
        raise LookupError("the Protein Atlas's releases page names no release")
    version, year, month, day, ensembl = found.groups()
    return Atlas(version, f"{year}-{month}-{day}", ensembl)


def ensembl_gene_of(row: dict, mane: dict) -> str:
    """The Ensembl gene MANE puts the row's transcript on, without version."""
    return versionless(mane[versionless(row["ensembl_nuc"])]["Ensembl_Gene"])


def _ranked(values: dict | None, unit: str) -> list[dict]:
    return sorted(
        ({"name": name, unit: float(value)} for name, value in (values or {}).items()),
        key=lambda found: (-found[unit], found["name"]))


def _pairs(values: list | None) -> list[dict]:
    """The Atlas's `Tissue - Cell type` pairs, in the order it lists them."""
    out = []
    for value in values or []:
        tissue, sep, cell = str(value).partition(" - ")
        if not sep or not tissue.strip() or not cell.strip():
            raise ValueError(f"the Atlas's tissue cell type {value!r} is not 'Tissue - Cell type'")
        out.append({"tissue": tissue.strip(), "cell_type": cell.strip()})
    return out


def expression_of(gene: dict) -> dict:
    """What the Atlas says of where a gene is read: its specificity to a
    tissue and to a cell type, and the ones it names, highest first; the cell
    types it finds the gene enriched in within a tissue; where in a cell it
    finds the protein; and where it is secreted to."""
    return {
        "tissue": {
            "specificity": gene.get("RNA tissue specificity"),
            "distribution": gene.get("RNA tissue distribution"),
            "specific": _ranked(gene.get("RNA tissue specific nTPM"), "ntpm"),
        },
        "cell_type": {
            "specificity": gene.get("RNA single cell type specificity"),
            "distribution": gene.get("RNA single cell type distribution"),
            "specific": _ranked(gene.get("RNA single cell type specific nCPM"), "ncpm"),
        },
        "tissue_cell_type": _pairs(gene.get("RNA tissue cell type enrichment")),
        "subcellular": {
            "main": list(gene.get("Subcellular main location") or []),
            "additional": list(gene.get("Subcellular additional location") or []),
        },
        "secretome": gene.get("Secretome location"),
    }


def path_of(expression: dict) -> dict:
    """The one way down the zoom takes, from the Atlas's readings alone.

    The organ is the tissue the gene's RNA is highest in. The cell is one
    that lives there: the first cell type the Atlas finds the gene enriched
    in within that tissue, else the single cell type it is highest in among
    those whose home is that tissue, else none, and the app draws the
    tissue's own cells and says so. A gene specific to no tissue goes to the
    home of the single cell type it is highest in, else to its first tissue
    cell type pair, else nowhere in particular. A cell with no nucleus lands
    in the cell that still has one.

    Refuses a reading whose names the Atlas's vocabularies do not have, so a
    new Atlas release that renames a tissue or a cell type stops the bake
    rather than lead a zoom astray.
    """
    tissues = []
    for found in expression["tissue"]["specific"]:
        name = tissue_name(found["name"])
        if name is None:
            raise LookupError(f"the Atlas's tissue {found['name']!r} is not one it lists")
        tissues.append(name)
    cells = []
    for found in expression["cell_type"]["specific"]:
        if single(found["name"]) is None:
            raise LookupError(f"the Atlas's cell type {found['name']!r} is not one it lists")
        cells.append(found["name"])
    pairs = []
    for pair in expression["tissue_cell_type"]:
        tissue = pair_tissue(pair["tissue"])
        if tissue is None:
            raise LookupError(f"the Atlas's tissue cell type tissue {pair['tissue']!r} "
                              "is not one it lists")
        pairs.append((tissue, pair["cell_type"]))

    tissue = cell = tissue_from = cell_from = None
    if tissues:
        tissue, tissue_from = tissues[0], "tissue"
        cell = next((c for t, c in pairs if t == tissue), None)
        if cell is not None:
            cell_from = "tissue_cell_type"
        else:
            cell = next((c for c in cells if tissue in homes(c)), None)
            cell_from = "single_cell_type" if cell is not None else None
    elif any(homes(c) for c in cells):
        cell = next(c for c in cells if homes(c))
        tissue, tissue_from, cell_from = homes(cell)[0], "cell_type", "single_cell_type"
    elif pairs:
        (tissue, cell), tissue_from, cell_from = pairs[0], "tissue_cell_type", "tissue_cell_type"

    kind = None
    if cell is not None:
        kind = cell_class(cell)
        if kind is None:
            raise LookupError(f"the path's cell type {cell!r} has no class the Atlas gives")
    lands = anucleate(cell) if cell is not None else None
    return {
        "tissue": tissue,
        "tissue_from": tissue_from,
        "cell_type": cell,
        "cell_class": kind,
        "cell_from": cell_from,
        "lands_in": None if lands is None else {**lands, "why": "no nucleus"},
    }


def locus_asset(target: Target) -> str:
    return f"assets/locus/{target.slug}_locus.json"


def _get(url: str) -> tuple[bytes, dict[str, str]]:
    request = urllib.request.Request(url, headers={"User-Agent": _AGENT})
    with urllib.request.urlopen(request, timeout=300) as response:
        return response.read(), {k.lower(): v for k, v in response.headers.items()}


def table_of(answer: dict, url: str) -> Table:
    """A UCSC API answer as a table: its rows, flattened out of the per-
    chromosome lists the API groups a whole-genome answer into."""
    track = answer["track"]
    found = answer[track]
    rows = found if isinstance(found, list) else [
        row for chrom in found.values() for row in chrom]
    return Table(
        track=track,
        genome=answer["genome"],
        updated=answer["dataTime"],
        rows=tuple(rows),
        url=url,
    )


def fetch_table(track: str) -> Table:
    url = UCSC_API.format(GENOME, track)
    body, _ = _get(url)
    return table_of(json.loads(body), url)


def chromosome_of(accession: str, aliases: Table) -> str:
    """UCSC's name for a RefSeq sequence ("NC_000011.10" -> "chr11"), as its
    own chromAlias table gives it."""
    for row in aliases.rows:
        if row["alias"] == accession and row.get("source") == "refseq":
            return row["chrom"]
    raise LookupError(f"UCSC's chromAlias names no chromosome {accession}")


def mane_row(target: Target, entry: Entry, mane: dict, release: str) -> dict:
    """The protein index's row for [target]: its entry joined to MANE, the
    buildable row for the table's gene."""
    rows = [row for row in rows_for_entry(entry.body, mane, {}, release)
            if row["buildable"] and row["gene"] == target.gene]
    if len(rows) != 1:
        raise LookupError(
            f"{target.slug}: {len(rows)} buildable MANE Select rows for {target.gene} "
            f"in {entry.accession}")
    return rows[0]


def bands_of(chrom: str, cytoband: Table) -> list[dict]:
    """Every band of [chrom], from the end of the short arm, 1-based."""
    found = sorted(
        (row for row in cytoband.rows if row["chrom"] == chrom),
        key=lambda row: row["chromStart"])
    return [{
        "name": row["name"],
        "start": row["chromStart"] + 1,
        "end": row["chromEnd"],
        "stain": row["gieStain"],
    } for row in found]


def locus_of(chromosome: str, names: list[str]) -> str:
    """How cytogenetics writes the place: 11p15.5; Xp21.2-p21.1 across two."""
    return chromosome + names[0] + ("" if len(names) == 1 else "-" + names[-1])


def payload(target: Target, entry: Entry, mane: dict, mane_release: str,
            cytoband: Table, aliases: Table, atlas: Atlas, atlas_gene: dict,
            retrieved: str) -> dict:
    """The track for one protein, from the sources as they were read."""
    body = entry.body
    if body.get("primaryAccession") != target.uniprot:
        raise ValueError(f"{target.slug}: asked for {target.uniprot}, "
                         f"UniProt answered {body.get('primaryAccession')}")
    row = mane_row(target, entry, mane, mane_release)
    ensembl = ensembl_gene_of(row, mane)
    if atlas_gene.get("Ensembl") != ensembl or atlas_gene.get("Gene") != target.gene:
        raise ValueError(f"{target.slug}: asked the Atlas for {ensembl}, it answered "
                         f"{atlas_gene.get('Ensembl')} ({atlas_gene.get('Gene')})")
    chrom = chromosome_of(row["chrom_acc"], aliases)
    bands = bands_of(chrom, cytoband)
    if not bands:
        raise LookupError(f"{target.slug}: cytoBand has no bands for {chrom}")
    start, end = row["chrom_start"], row["chrom_end"]
    inside = [b for b in bands if b["start"] <= end and b["end"] >= start]
    if not inside:
        raise LookupError(f"{target.slug}: {chrom}:{start}-{end} falls in no band")
    chromosome = chrom[3:] if chrom.startswith("chr") else chrom
    names = [b["name"] for b in inside]
    audit = body.get("entryAudit", {})
    expression = {"ensembl_gene": ensembl, **expression_of(atlas_gene)}
    return {
        "slug": target.slug,
        "gene": target.gene,
        "uniprot": target.uniprot,
        "assembly": ASSEMBLY,
        "genome": cytoband.genome,
        "chromosome": chromosome,
        "sequence": row["chrom_acc"],
        "length": bands[-1]["end"],
        "span": {"start": start, "end": end, "strand": row["chrom_strand"]},
        "transcript": {"refseq": row["refseq_nuc"], "ensembl": row["ensembl_nuc"]},
        "band": {"names": names, "start": inside[0]["start"], "end": inside[-1]["end"]},
        "locus": locus_of(chromosome, names),
        "bands": bands,
        "expression": expression,
        "path": path_of(expression),
        "sources": {
            "cytoband": cytoband.provenance(),
            "chrom_alias": aliases.provenance(),
            "mane": {"release": mane_release},
            "uniprot": {
                "release": entry.release,
                "release_date": entry.release_date,
                "entry_version": audit.get("entryVersion"),
            },
            "hpa": atlas.provenance(ensembl),
        },
        "retrieved": retrieved,
        "schema_version": SCHEMA_VERSION,
        "built_by": "pipeline/locus/bake_locus.py",
    }


def encode(track: dict) -> bytes:
    return (json.dumps(track, indent=2, ensure_ascii=False) + "\n").encode()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if not args.all and not args.target:
        parser.error("name --target or --all")
    unknown = sorted(set(args.target or []) - set(BY_SLUG))
    if unknown:
        parser.error(f"no such target: {unknown}")
    chosen = LOCUS_TARGETS if args.all else tuple(BY_SLUG[s] for s in args.target)

    cytoband = fetch_table("cytoBand")
    aliases = fetch_table("chromAlias")
    mane, mane_release = current_summary(_get)
    atlas = atlas_of(_get(HPA_RELEASES)[0].decode("utf-8", "replace"))
    retrieved = datetime.now(timezone.utc).date().isoformat()
    print(f"cytoBand {cytoband.genome}, updated {cytoband.updated}, {len(cytoband.rows)} rows; "
          f"chromAlias updated {aliases.updated}; MANE {mane_release}; "
          f"Protein Atlas {atlas.version} ({atlas.release_date}, Ensembl {atlas.ensembl})")

    failed = []
    for target in chosen:
        try:
            entry = fetch_entry(target.uniprot)
            ensembl = ensembl_gene_of(mane_row(target, entry, mane, mane_release), mane)
            atlas_gene = json.loads(_get(HPA_GENE.format(ensembl))[0])
            track = payload(target, entry, mane, mane_release, cytoband, aliases,
                            atlas, atlas_gene, retrieved)
        except (LookupError, ValueError) as exc:
            failed.append(f"{target.slug}: {exc}")
            print(f"{target.slug:<16} FAILED: {exc}")
            continue
        destination = DATA / locus_asset(target)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(encode(track))
        span = track["span"]
        path = track["path"]
        lands = path["lands_in"]
        print(f"{target.slug:<16} {track['gene']:<6} {track['locus']:<14} "
              f"{track['sequence']}:{span['start']}-{span['end']}  "
              f"{path['tissue'] or '-'} ({path['tissue_from'] or '-'}) / "
              f"{path['cell_type'] or 'drawn'} ({path['cell_from'] or '-'})"
              f"{' -> ' + lands['cell'] if lands else ''}")
    if failed:
        print(f"\n{len(failed)} protein(s) did not resolve to a band:")
        for line in failed:
            print("  " + line)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
