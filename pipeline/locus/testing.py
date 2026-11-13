"""Offline snapshots for the locus tests: the pinned inputs, made small.

A snapshot is built from the locus fixtures (six Atlas genes, their MANE rows,
five chromosomes' UCSC tables, the Atlas's releases page), with what the
fixtures lack made up small for it: level tables, HGNC ids, NCBI's placement
of two scaffolds (a fix patch and an alternate locus placed reversed), the
UCSC rows naming those scaffolds and chromosome 19, and an Atlas entry for each
placed gene, made from p53's reading. Nothing fetches. Not a test module: the
tests import it.
"""

from __future__ import annotations

import gzip
import json
import zipfile
from pathlib import Path

from pipeline.locus import bake_locus, inputs, snapshot

FIXTURES = Path(__file__).resolve().parent / "fixtures"

# The six fixture genes' UniProt accessions, which the cut-down Atlas
# fixtures leave out.
UNIPROT = {
    "ENSG00000001626": ["P13569"], "ENSG00000141510": ["P04637"],
    "ENSG00000142168": ["P00441"], "ENSG00000198947": ["P11532"],
    "ENSG00000244734": ["P68871"], "ENSG00000254647": ["P01308"],
}

# Two genes MANE puts on scaffolds, as the protein index rows them: POLR2A on
# a fix patch of chromosome 17, KIR2DS4 on an alternate locus of chromosome 19
# that NCBI places reversed. Their Atlas entries (made here) are the
# primary-assembly genes HGNC gives their HGNC ids.
POLR2A = {"uniprot": "P24928", "gene": "POLR2A", "name": "DNA-directed RNA polymerase II subunit RPB1",
          "length": 1970, "hgnc_id": "HGNC:9187", "refseq_nuc": "NM_000937.5",
          "refseq_prot": "NP_000928.1", "ensembl_nuc": "ENST00000674977.2",
          "mane_isoform": None, "chrom_acc": "NW_016107299.1", "chrom_start": 104823,
          "chrom_end": 135071, "chrom_strand": 1, "refseqgene": None, "buildable": True,
          "unavailable_reason": None}
KIR2DS4 = {"uniprot": "P43632", "gene": "KIR2DS4",
           "name": "Killer cell immunoglobulin-like receptor 2DS4", "length": 304,
           "hgnc_id": "HGNC:6336", "refseq_nuc": "NM_012314.6", "refseq_prot": "NP_036446.3",
           "ensembl_nuc": "ENST00000612454.1", "mane_isoform": None,
           "chrom_acc": "NT_187669.1", "chrom_start": 20467, "chrom_end": 36357,
           "chrom_strand": -1, "refseqgene": None, "buildable": True,
           "unavailable_reason": None}
PLACED_ATLAS = {"ENSG00000181222": ("POLR2A", "P24928"),
                "ENSG00000221957": ("KIR2DS4", "P43632")}

PLACEMENT = (
    "#alt_asm_name\tprim_asm_name\talt_scaf_name\talt_scaf_acc\tparent_type\tparent_name\t"
    "parent_acc\tregion_name\tori\talt_scaf_start\talt_scaf_stop\tparent_start\tparent_stop\t"
    "alt_start_tail\talt_stop_tail\n"
    "PATCHES\tPrimary Assembly\tHG2046_PATCH\tNW_016107299.1\tCHROMOSOME\t17\tNC_000017.11\t"
    "REGION108\t+\t1\t154723\t7379544\t7534268\t0\t0\n"
    "ALT_REF_LOCI_21\tPrimary Assembly\tHSCHR19KIR_T7526_A_HAP_CTG3_1\tNT_187669.1\tCHROMOSOME\t"
    "19\tNC_000019.10\tLRC_KIR\t-\t1\t170665\t54698377\t54869056\t0\t0\n"
)

# Chromosome 19 cut down to three bands, whole from base 1 to its last, and the
# UCSC names of it and of the two scaffolds.
CHR19_BANDS = [
    {"chrom": "chr19", "chromStart": 0, "chromEnd": 26200000, "name": "p13", "gieStain": "gneg"},
    {"chrom": "chr19", "chromStart": 26200000, "chromEnd": 28100000, "name": "q11",
     "gieStain": "acen"},
    {"chrom": "chr19", "chromStart": 28100000, "chromEnd": 58617616, "name": "q13.42",
     "gieStain": "gpos25"},
]
EXTRA_ALIASES = [
    {"alias": "NC_000019.10", "chrom": "chr19", "source": "refseq"},
    {"alias": "NW_016107299.1", "chrom": "chr17_KV575245v1_fix", "source": "refseq"},
    {"alias": "NT_187669.1", "chrom": "chr19_KI270915v1_alt", "source": "refseq"},
]

RELEASES = (FIXTURES / "hpa_releases.html").read_text(encoding="utf-8")


def atlas_genes() -> list:
    """The six fixture genes with their UniProt accessions, and an entry for
    each placed gene made from p53's reading."""
    genes = []
    for path in sorted((FIXTURES / "hpa").glob("ENSG*.json")):
        gene = json.loads(path.read_text(encoding="utf-8"))
        genes.append({**gene, "Uniprot": UNIPROT[gene["Ensembl"]]})
    p53 = next(gene for gene in genes if gene["Gene"] == "TP53")
    for ensembl, (name, accession) in PLACED_ATLAS.items():
        genes.append({**p53, "Ensembl": ensembl, "Gene": name, "Uniprot": [accession]})
    return genes


def _answer(track: str, rows: list) -> bake_locus.Table:
    return bake_locus.table_of({"genome": "hg38", "dataTime": "2022-11-03T14:54:06",
                                "track": track, track: rows},
                               f"https://api.genome.ucsc.edu/{track}")


def tables() -> dict:
    def fixture(name: str, extra: list) -> bake_locus.Table:
        answer = json.loads((FIXTURES / name).read_text(encoding="utf-8"))
        table = bake_locus.table_of(answer, "https://api.genome.ucsc.edu/test")
        return bake_locus.Table(track=table.track, genome=table.genome, updated=table.updated,
                                rows=table.rows + tuple(extra), url=table.url)

    return {
        "cytoBand": fixture("cytoband_hg38_five.json", CHR19_BANDS),
        "chromAlias": fixture("chromalias_hg38_five.json", EXTRA_ALIASES),
        "altLocations": _answer("altLocations", [
            {"chrom": "chr19", "chromStart": 54698376, "chromEnd": 54869056,
             "name": "chr19_KI270915v1_alt:1-170665"}]),
        "fixLocations": _answer("fixLocations", [
            {"chrom": "chr17", "chromStart": 7379543, "chromEnd": 7534268,
             "name": "chr17_KV575245v1_fix:1-154723"}]),
    }


def _zip(path: Path, name: str, text: str) -> None:
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(name, text)


def downloads(folder: Path) -> Path:
    """The files `snapshot.snapshot_text` reads, made from the fixtures."""
    folder.mkdir(parents=True, exist_ok=True)
    genes = atlas_genes()
    with gzip.open(folder / "proteinatlas.json.gz", "wt", encoding="utf-8") as stream:
        json.dump(genes, stream)
    tissues = ["Gene\tGene name\tTissue\tnTPM"]
    cells = ["Gene\tGene name\tCell type\tnCPM"]
    for n, gene in enumerate(genes):
        tissues += [f"{gene['Ensembl']}\t{gene['Gene']}\tpancreas\t{10.5 * n}",
                    f"{gene['Ensembl']}\t{gene['Gene']}\tcolon\t{3.0 + n}",
                    f"{gene['Ensembl']}\t{gene['Gene']}\tsmall intestine\t0.0"]
        cells += [f"{gene['Ensembl']}\t{gene['Gene']}\tpancreatic islet cells\t{100 + n}",
                  f"{gene['Ensembl']}\t{gene['Gene']}\tenterocytes\t2.5"]
    _zip(folder / "rna_tissue_consensus.tsv.zip", "rna_tissue_consensus.tsv",
         "\n".join(tissues) + "\n")
    _zip(folder / "rna_single_cell_type.tsv.zip", "rna_single_cell_type.tsv",
         "\n".join(cells) + "\n")
    (folder / "mane_summary.txt.gz").write_bytes(
        gzip.compress((FIXTURES / "mane_v1.5_six.tsv").read_bytes(), mtime=0))
    (folder / "hgnc_complete_set.txt").write_text(
        "hgnc_id\tsymbol\tensembl_gene_id\n"
        "HGNC:6081\tINS\tENSG00000254647\n"
        "HGNC:9187\tPOLR2A\tENSG00000181222\n"
        "HGNC:6336\tKIR2DS4\tENSG00000221957\n"
        "HGNC:999999\tNOENSEMBL\t\n", encoding="utf-8")
    (folder / "alt_scaffold_placement.txt").write_text(PLACEMENT, encoding="utf-8")
    return folder


def text(folder: Path, created: str = "2026-11-13", release_page: str = RELEASES) -> bytes:
    return snapshot.snapshot_text(folder=folder, release_page=release_page, tables=tables(),
                                  created=created)


def built(folder: Path) -> tuple[Path, str, inputs.Inputs]:
    """A snapshot written under [folder] and read back: its path, its digest,
    and it read."""
    path, digest = snapshot.write(text(downloads(folder / "downloads")), folder / "snapshot")
    return path, digest, inputs.read(path, digest)
