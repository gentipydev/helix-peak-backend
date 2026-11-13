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

Everything but UniProt is read from one pinned snapshot (`inputs.py`, built by
`snapshot.py`): the Atlas 25.1, the UCSC tables and MANE v1.5, as they were
read once, and nothing fetched gene by gene. The Atlas's site serves its
newest release, and a bake that read it would move a curated zoom to another
release without anyone deciding to. `retrieved` is the snapshot's date.

Which proteins get the track: every curated one (`LOCUS_TARGETS`), recorded
here, never in targets.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.request
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from app.protein_index import rows_for_entry, versionless  # noqa: E402
from pipeline.locus import inputs as locus_inputs  # noqa: E402
from pipeline.locus.cell_types import anucleate, cell_class, homes, single  # noqa: E402
from pipeline.locus.tissues import pair_tissue, tissue_name  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import BY_SLUG, TARGETS, Target  # noqa: E402
from pipeline.uniprot import Entry, fetch_entry  # noqa: E402

# Every protein gets the track: each has a chromosome to zoom to.
LOCUS_TARGETS = TARGETS

# 2 adds keys only: the Atlas's tissue cell type pairs, subcellular location
# and secretome location under `expression`, and the zoom's `path`. A schema 1
# reader ignores them.
# 3 adds keys only, for every protein the search offers: how the Atlas entry
# was found and the Atlas's own name for the gene (`expression.atlas_match`,
# `atlas_gene`), where a gene MANE puts on a patch scaffold was carried from
# (`placed_from`), and the snapshot it was made from (`sources.locus_inputs`).
# A schema 2 reader ignores them. The curated keep their stored schema 2.
SCHEMA_VERSION = 3
SCHEMAS = (2, 3)

ASSEMBLY = "GRCh38"
GENOME = "hg38"
UCSC_API = "https://api.genome.ucsc.edu/getData/track?genome={};track={}"

_AGENT = "helixpeek-locus-bake"

# The Human Protein Atlas's site, where a gene's summary is `<site>/<ENSG>.json`.
# The snapshot names the site of its own release (`v25.proteinatlas.org`).
HPA_SITE = "https://www.proteinatlas.org"
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

    @cached_property
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
    site: str = HPA_SITE

    def provenance(self, gene: str) -> dict:
        return {
            "version": self.version,
            "release_date": self.release_date,
            "ensembl": self.ensembl,
            "licence": HPA_LICENCE,
            "url": f"{self.site}/{gene}.json",
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


def table_from(snapshot: locus_inputs.Inputs, track: str) -> Table:
    """One UCSC table as the snapshot keeps it: the rows its API served, so
    its digest is the one a live read of the same table gives."""
    found = snapshot.table(track)
    return Table(track=track, genome=found["genome"], updated=found["updated"],
                 rows=tuple(found["rows"]), url=found["url"])


def atlas_from(snapshot: locus_inputs.Inputs) -> Atlas:
    """The Atlas release the snapshot was read from, at that release's site."""
    hpa = snapshot.atlas
    return Atlas(hpa["version"], hpa["release_date"], hpa["ensembl"], site=hpa["site"])


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


class LocusRefused(Exception):
    """No zoom can be drawn for this gene. The message is the reader's
    sentence, in the app's data-first voice; the track is `refused` with it."""


# What a gene the Atlas has no entry for is refused with. GSTT1, KIR2DS2,
# HLA-DRB3 and LILRA3 live only on alternate haplotypes, which the Atlas's
# Ensembl does not annotate.
NO_ATLAS_ENTRY = "The Human Protein Atlas has no reading for this gene."

# The accessions MANE puts a transcript on: the primary chromosomes 1 to 22,
# X (23) and Y (24), and otherwise a GRCh38 patch or alternate locus scaffold
# (an `NW_` or `NT_` accession), which `placed` carries to its chromosome.
_PRIMARY = re.compile(r"^NC_0000(\d\d)\.\d+$")

# How a scaffold's span is carried to its chromosome: each position of the
# part of the scaffold NCBI aligns, at the same share of the way along the
# region of the chromosome it is aligned to, the other way along where the
# scaffold lies reversed. Exact at the ends of the alignment, and in between
# off by at most the difference of the two lengths (none for a fix patch,
# 1.5% for the longest alternate locus): far finer than a band.
PLACEMENT_METHOD = "proportional"
PLACEMENT_ASSEMBLY = "GRCh38.p14"

# How the Atlas entry that reads a gene was found (`atlas_match`).
MATCHES = ("mane", "hgnc", "uniprot")


def on_primary(accession: str | None) -> bool:
    found = _PRIMARY.match(accession or "")
    return found is not None and 1 <= int(found.group(1)) <= 24


def carried(position: int, scaffold: tuple[int, int], chromosome: tuple[int, int],
            orientation: str) -> int:
    """A scaffold position on its chromosome (`PLACEMENT_METHOD`)."""
    (a0, a1), (c0, c1) = scaffold, chromosome
    share = (position - a0) / (a1 - a0) if a1 > a0 else 0.0
    if orientation == "-":
        share = 1 - share
    return c0 + round(share * (c1 - c0))


def placed(row: dict, placements: list, aliases: Table) -> tuple[dict, dict]:
    """[row]'s span carried from the GRCh38 patch or alternate locus scaffold
    MANE puts its transcript on to the chromosome NCBI places that scaffold on
    (GRCh38.p14's `alt_scaffold_placement`): where the gene lies, and the
    `placed_from` block that says where it was read and how it was carried.
    The strand turns where the scaffold lies reversed."""
    accession = row["chrom_acc"]
    if len(placements) != 1:
        raise LookupError(f"{row['gene']}: NCBI places scaffold {accession} "
                          f"{len(placements)} times, not once")
    found = placements[0]
    parent = found["parent_acc"]
    if found["parent_type"] != "CHROMOSOME" or not on_primary(parent):
        raise LookupError(f"{row['gene']}: scaffold {accession} is placed on "
                          f"{found['parent_type']} {parent}, not a chromosome")
    orientation = found["ori"]
    if orientation not in ("+", "-"):
        raise LookupError(f"{row['gene']}: scaffold {accession} has orientation {orientation!r}")
    scaffold = (int(found["alt_scaf_start"]), int(found["alt_scaf_stop"]))
    chromosome = (int(found["parent_start"]), int(found["parent_stop"]))
    start, end, strand = row["chrom_start"], row["chrom_end"], row["chrom_strand"]
    if not scaffold[0] <= start <= end <= scaffold[1]:
        raise LookupError(f"{row['gene']}: {accession}:{start}-{end} is outside the part "
                          f"of the scaffold NCBI aligns ({scaffold[0]}-{scaffold[1]})")
    ends = sorted(carried(p, scaffold, chromosome, orientation) for p in (start, end))
    place = {"sequence": parent, "start": ends[0], "end": ends[1],
             "strand": strand if orientation == "+" else -strand}
    placed_from = {
        "sequence": accession,
        "ucsc": chromosome_of(accession, aliases),
        "scaffold": found["alt_scaf_name"],
        "span": {"start": start, "end": end, "strand": strand},
        "placement": {
            "assembly": PLACEMENT_ASSEMBLY,
            "orientation": orientation,
            "scaffold": {"start": scaffold[0], "end": scaffold[1]},
            "chromosome": {"start": chromosome[0], "end": chromosome[1]},
        },
        "method": PLACEMENT_METHOD,
    }
    return place, placed_from


def atlas_match(row: dict, snapshot: locus_inputs.Inputs) -> tuple[str, str] | None:
    """The Atlas entry that reads an index row's gene, and how it was found:
    the Ensembl gene MANE puts its transcript on; else the one HGNC gives its
    HGNC id (a gene MANE annotates on a patch scaffold has a gene id of the
    patch's own, which the Atlas, built on the primary assembly, never saw);
    else the one entry whose UniProt accessions name its protein, if exactly
    one does. None where none does."""
    mane = snapshot.mane.get(versionless(row.get("ensembl_nuc") or ""))
    if mane is not None:
        gene = versionless(mane["Ensembl_Gene"])
        if snapshot.atlas_gene(gene) is not None:
            return gene, "mane"
    gene = snapshot.hgnc_ensembl(row.get("hgnc_id"))
    if gene and snapshot.atlas_gene(gene) is not None:
        return gene, "hgnc"
    found = snapshot.by_uniprot(row.get("uniprot") or "")
    if len(found) == 1:
        return found[0], "uniprot"
    return None


def track_of(*, slug: str, gene: str, uniprot: str, row: dict, atlas_gene: dict,
             match: tuple[str, str], cytoband: Table, aliases: Table, atlas: Atlas,
             mane_release: str, uniprot_source: dict, retrieved: str,
             placements: list = (), inputs_source: dict | None = None) -> dict:
    """The track for one protein: its index row's span, on its chromosome
    (carried there from a patch scaffold where MANE puts it on one), and the
    reading of the Atlas entry [match] chose. Refuses an Atlas answer for
    another entry than the one chosen."""
    ensembl, how = match
    if how not in MATCHES:
        raise ValueError(f"{slug}: no such way to match the Atlas: {how!r}")
    if atlas_gene.get("Ensembl") != ensembl:
        raise ValueError(f"{slug}: asked the Atlas for {ensembl}, it answered "
                         f"{atlas_gene.get('Ensembl')} ({atlas_gene.get('Gene')})")
    if on_primary(row["chrom_acc"]):
        place = {"sequence": row["chrom_acc"], "start": row["chrom_start"],
                 "end": row["chrom_end"], "strand": row["chrom_strand"]}
        placed_from = None
    else:
        place, placed_from = placed(row, list(placements), aliases)
    chrom = chromosome_of(place["sequence"], aliases)
    bands = bands_of(chrom, cytoband)
    if not bands:
        raise LookupError(f"{slug}: cytoBand has no bands for {chrom}")
    start, end = place["start"], place["end"]
    inside = [b for b in bands if b["start"] <= end and b["end"] >= start]
    if not inside:
        raise LookupError(f"{slug}: {chrom}:{start}-{end} falls in no band")
    chromosome = chrom[3:] if chrom.startswith("chr") else chrom
    names = [b["name"] for b in inside]
    expression = {"ensembl_gene": ensembl, **expression_of(atlas_gene),
                  "atlas_match": how, "atlas_gene": atlas_gene.get("Gene")}
    track = {
        "slug": slug,
        "gene": gene,
        "uniprot": uniprot,
        "assembly": ASSEMBLY,
        "genome": cytoband.genome,
        "chromosome": chromosome,
        "sequence": place["sequence"],
        "length": bands[-1]["end"],
        "span": {"start": start, "end": end, "strand": place["strand"]},
        "transcript": {"refseq": row["refseq_nuc"], "ensembl": row["ensembl_nuc"]},
    }
    if placed_from is not None:
        track["placed_from"] = placed_from
    sources = {
        "cytoband": cytoband.provenance(),
        "chrom_alias": aliases.provenance(),
        "mane": {"release": mane_release},
        "uniprot": uniprot_source,
        "hpa": atlas.provenance(ensembl),
    }
    if inputs_source is not None:
        sources["locus_inputs"] = inputs_source
    track.update({
        "band": {"names": names, "start": inside[0]["start"], "end": inside[-1]["end"]},
        "locus": locus_of(chromosome, names),
        "bands": bands,
        "expression": expression,
        "path": path_of(expression),
        "sources": sources,
        "retrieved": retrieved,
        "schema_version": SCHEMA_VERSION,
        "built_by": "pipeline/locus/bake_locus.py",
    })
    return track


def payload(target: Target, entry: Entry, mane: dict, mane_release: str,
            cytoband: Table, aliases: Table, atlas: Atlas, atlas_gene: dict,
            retrieved: str, *, match: tuple[str, str] | None = None,
            inputs_source: dict | None = None) -> dict:
    """A curated protein's track, from its UniProt entry joined to MANE, the
    way the protein index joins them, and the Atlas entry [match] chose
    (MANE's Ensembl gene, unless said)."""
    body = entry.body
    if body.get("primaryAccession") != target.uniprot:
        raise ValueError(f"{target.slug}: asked for {target.uniprot}, "
                         f"UniProt answered {body.get('primaryAccession')}")
    row = mane_row(target, entry, mane, mane_release)
    audit = body.get("entryAudit", {})
    return track_of(
        slug=target.slug, gene=target.gene, uniprot=target.uniprot, row=row,
        atlas_gene=atlas_gene, match=match or (ensembl_gene_of(row, mane), "mane"),
        cytoband=cytoband, aliases=aliases, atlas=atlas, mane_release=mane_release,
        uniprot_source={"release": entry.release, "release_date": entry.release_date,
                        "entry_version": audit.get("entryVersion")},
        retrieved=retrieved, inputs_source=inputs_source)


def inputs_source(snapshot: locus_inputs.Inputs) -> dict:
    """What a track names of the snapshot it was made from."""
    return {"sha256": snapshot.sha256, "created": snapshot.created}


_READ: dict = {}


def _read(snapshot: locus_inputs.Inputs) -> tuple[Table, Table, Atlas]:
    """The snapshot's cytoBand and chromAlias tables and its Atlas, made once."""
    if snapshot.sha256 not in _READ:
        _READ[snapshot.sha256] = (table_from(snapshot, "cytoBand"),
                                  table_from(snapshot, "chromAlias"), atlas_from(snapshot))
    return _READ[snapshot.sha256]


def bake(row: dict, snapshot: locus_inputs.Inputs, *, slug: str, uniprot_source: dict,
         mane_release: str | None = None) -> dict:
    """The track for one protein the index offers, from its index row and the
    pinned snapshot alone, nothing fetched: the same bytes wherever and
    whenever it is made from the same row. [uniprot_source] names the UniProt
    release the index row was read from, and [mane_release] its MANE release
    (the snapshot's unless said). Raises `LocusRefused` for a gene no zoom can
    be drawn for, with the reader's sentence."""
    match = atlas_match(row, snapshot)
    if match is None:
        raise LocusRefused(NO_ATLAS_ENTRY)
    cytoband, aliases, atlas = _read(snapshot)
    return track_of(
        slug=slug, gene=row["gene"], uniprot=row["uniprot"], row=row,
        atlas_gene=snapshot.atlas_gene(match[0]), match=match, cytoband=cytoband,
        aliases=aliases, atlas=atlas, mane_release=mane_release or snapshot.mane_release,
        uniprot_source=uniprot_source, retrieved=snapshot.created,
        placements=snapshot.placements(row["chrom_acc"]),
        inputs_source=inputs_source(snapshot))


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

    snapshot = locus_inputs.load()
    cytoband = table_from(snapshot, "cytoBand")
    aliases = table_from(snapshot, "chromAlias")
    mane, mane_release = snapshot.mane, snapshot.mane_release
    atlas = atlas_from(snapshot)
    retrieved = snapshot.created
    print(f"locus inputs {snapshot.sha256[:12]} ({snapshot.created}): "
          f"cytoBand {cytoband.genome}, updated {cytoband.updated}, {len(cytoband.rows)} rows; "
          f"chromAlias updated {aliases.updated}; MANE {mane_release}; "
          f"Protein Atlas {atlas.version} ({atlas.release_date}, Ensembl {atlas.ensembl})")

    failed = []
    for target in chosen:
        try:
            entry = fetch_entry(target.uniprot)
            match = atlas_match(mane_row(target, entry, mane, mane_release), snapshot)
            if match is None:
                raise LookupError(f"the Protein Atlas {atlas.version} has no entry for "
                                  f"{target.gene}")
            track = payload(target, entry, mane, mane_release, cytoband, aliases,
                            atlas, snapshot.atlas_gene(match[0]), retrieved, match=match,
                            inputs_source=inputs_source(snapshot))
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
