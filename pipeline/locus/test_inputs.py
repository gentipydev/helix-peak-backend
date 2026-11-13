"""The locus bake's pinned inputs: the snapshot, made and read, offline.

A snapshot is built here from the locus fixtures (six Atlas genes, their MANE
rows, five chromosomes' UCSC tables, the Atlas's releases page), with level
tables, HGNC ids and NCBI placements made up small for it. Nothing fetches.
One test reads the real pinned snapshot, and skips where it has not been built
(`snapshot.py`), as the tests that need stored tracks skip until
`fetch_tracks.py` has filled `pipeline/data/`.
"""

from __future__ import annotations

import gzip
import io
import json
import sys
import zipfile
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.locus import bake_locus, inputs, snapshot  # noqa: E402
from pipeline.paths import DATA  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures"

# The six fixture genes' UniProt accessions, which the cut-down Atlas
# fixtures leave out.
UNIPROT = {
    "ENSG00000001626": ["P13569"], "ENSG00000141510": ["P04637"],
    "ENSG00000142168": ["P00441"], "ENSG00000198947": ["P11532"],
    "ENSG00000244734": ["P68871"], "ENSG00000254647": ["P01308"],
}

PLACEMENT = (
    "#alt_asm_name\tprim_asm_name\talt_scaf_name\talt_scaf_acc\tparent_type\tparent_name\t"
    "parent_acc\tregion_name\tori\talt_scaf_start\talt_scaf_stop\tparent_start\tparent_stop\t"
    "alt_start_tail\talt_stop_tail\n"
    "PATCHES\tPrimary Assembly\tHG2046_PATCH\tNW_016107299.1\tCHROMOSOME\t17\tNC_000017.11\t"
    "REGION108\t+\t1\t154723\t7379544\t7534268\t0\t0\n"
    "ALT_REF_LOCI_21\tPrimary Assembly\tHSCHR19KIR_T7526_A_HAP_CTG3_1\tNT_187669.1\tCHROMOSOME\t"
    "19\tNC_000019.10\tLRC_KIR\t-\t1\t170665\t54698377\t54869056\t0\t0\n"
)


def _atlas_genes() -> list:
    genes = []
    for path in sorted((FIXTURES / "hpa").glob("ENSG*.json")):
        gene = json.loads(path.read_text(encoding="utf-8"))
        genes.append({**gene, "Uniprot": UNIPROT[gene["Ensembl"]]})
    return genes


def _zip(path: Path, name: str, text: str) -> None:
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(name, text)


def _answer(track: str, rows: list) -> bake_locus.Table:
    return bake_locus.table_of({"genome": "hg38", "dataTime": "2022-11-03T14:54:06",
                                "track": track, track: rows},
                               f"https://api.genome.ucsc.edu/{track}")


def _tables() -> dict:
    def fixture(name: str) -> bake_locus.Table:
        answer = json.loads((FIXTURES / name).read_text(encoding="utf-8"))
        return bake_locus.table_of(answer, "https://api.genome.ucsc.edu/test")

    return {
        "cytoBand": fixture("cytoband_hg38_five.json"),
        "chromAlias": fixture("chromalias_hg38_five.json"),
        "altLocations": _answer("altLocations", [
            {"chrom": "chr19", "chromStart": 54698376, "chromEnd": 54869056,
             "name": "chr19_KI270915v1_alt:1-170665"}]),
        "fixLocations": _answer("fixLocations", [
            {"chrom": "chr17", "chromStart": 7379543, "chromEnd": 7534268,
             "name": "chr17_KV575245v1_fix:1-154723"}]),
    }


def _downloads(folder: Path) -> Path:
    """The files `snapshot_text` reads, made from the fixtures."""
    folder.mkdir(parents=True, exist_ok=True)
    genes = _atlas_genes()
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
        "HGNC:999999\tNOENSEMBL\t\n", encoding="utf-8")
    (folder / "alt_scaffold_placement.txt").write_text(PLACEMENT, encoding="utf-8")
    return folder


RELEASES = (FIXTURES / "hpa_releases.html").read_text(encoding="utf-8")


def _text(folder: Path, created: str = "2026-11-13") -> bytes:
    return snapshot.snapshot_text(folder=folder, release_page=RELEASES, tables=_tables(),
                                  created=created)


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """A snapshot written and read back: its path, its digest, and it read."""
    folder = _downloads(tmp_path_factory.mktemp("downloads"))
    path, digest = snapshot.write(_text(folder), tmp_path_factory.mktemp("snapshot"))
    return path, digest, inputs.read(path, digest)


def test_the_same_sources_make_the_same_snapshot(tmp_path):
    folder = _downloads(tmp_path / "downloads")
    first, digest = snapshot.write(_text(folder), tmp_path / "a")
    second, again = snapshot.write(_text(folder), tmp_path / "b")
    assert digest == again
    assert first.read_bytes() == second.read_bytes()
    assert first.name == f"locus-inputs.{digest[:12]}.jsonl.gz"
    # The digest is the text's, so another compressor still matches it.
    assert inputs.read(first, digest).sha256 == digest


def test_another_date_is_another_snapshot(tmp_path):
    folder = _downloads(tmp_path / "downloads")
    assert _text(folder, "2026-11-13") != _text(folder, "2026-11-14")


def test_any_other_snapshot_is_refused(built, tmp_path):
    path, digest, _ = built
    with pytest.raises(inputs.InputsUnavailable, match="not the pinned"):
        inputs.read(path, "f" * 64)
    with pytest.raises(inputs.InputsUnavailable, match="snapshot.py"):
        inputs.read(tmp_path / "nowhere.jsonl.gz", digest)
    # A file that is not one at all.
    junk = tmp_path / "junk.jsonl.gz"
    junk.write_bytes(gzip.compress(b'{"kind": "something else", "format": 1}\n'))
    with pytest.raises(inputs.InputsUnavailable):
        inputs.read(junk, None)


def test_a_gene_reads_as_the_atlas_answered_it(built):
    _, _, read = built
    for gene in _atlas_genes():
        found = read.atlas_gene(gene["Ensembl"])
        assert found == {name: gene.get(name) for name in snapshot.ATLAS_FIELDS}
        assert bake_locus.expression_of(found) == bake_locus.expression_of(gene)
    assert read.atlas_gene("ENSG00000000000") is None


def test_its_levels_are_the_tables_own(built):
    _, _, read = built
    tissues, cells = read.levels("ENSG00000254647")
    assert set(tissues) == {"colon", "pancreas", "small intestine"}
    assert tissues["small intestine"] == 0
    assert cells == {"enterocytes": 2.5, "pancreatic islet cells": 105}
    assert read.levels("ENSG00000000000") is None


def test_it_finds_a_gene_by_uniprot_and_by_hgnc(built):
    _, _, read = built
    assert read.by_uniprot("P01308") == ("ENSG00000254647",)
    assert read.by_uniprot("Q00000") == ()
    assert read.hgnc_ensembl("HGNC:9187") == "ENSG00000181222"
    assert read.hgnc_ensembl("HGNC:999999") is None
    assert read.hgnc_ensembl(None) is None


def test_it_places_a_scaffold_by_its_refseq_accession(built):
    _, _, read = built
    (kir,) = read.placements("NT_187669.1")
    assert (kir["parent_acc"], kir["ori"], kir["parent_start"]) == ("NC_000019.10", "-", "54698377")
    assert read.placements("NW_000000000.1") == []


def test_a_table_from_it_is_the_table_as_read(built):
    _, _, read = built
    for track, table in _tables().items():
        kept = bake_locus.table_from(read, track)
        assert kept.sha256 == table.sha256
        assert kept.provenance() == table.provenance()


def test_it_names_its_sources(built):
    _, _, read = built
    assert (read.atlas["version"], read.atlas["site"]) == ("25.1", snapshot.HPA_SITE)
    assert read.mane_release == "v1.5"
    assert read.sources["hgnc"]["date"] == snapshot.HGNC_DATE
    assert read.sources["placement"]["assembly"] == "GRCh38.p14"
    atlas = bake_locus.atlas_from(read)
    assert atlas.provenance("ENSG00000254647")["url"] == \
        "https://v25.proteinatlas.org/ENSG00000254647.json"


def test_it_refuses_another_atlas_release(tmp_path):
    folder = _downloads(tmp_path / "downloads")
    page = RELEASES.replace("version 25.1", "version 26.0", 1)
    with pytest.raises(ValueError, match="26.0"):
        snapshot.snapshot_text(folder=folder, release_page=page, tables=_tables(),
                               created="2026-11-13")


def test_it_is_looked_for_where_it_is_named(tmp_path, monkeypatch):
    monkeypatch.setenv(inputs.ENV, str(tmp_path))
    assert inputs.default_path("ab" * 32) == tmp_path / "locus-inputs.abababababab.jsonl.gz"
    named = tmp_path / "elsewhere.jsonl.gz"
    monkeypatch.setenv(inputs.ENV, str(named))
    assert inputs.default_path() == named
    monkeypatch.delenv(inputs.ENV)
    assert inputs.default_path().parent == inputs.FOLDER


_PINNED = inputs.default_path()


@pytest.mark.skipif(not _PINNED.is_file(), reason="the pinned snapshot is not built here")
def test_the_pinned_snapshot_gives_every_stored_curated_reading():
    """Every curated protein's stored `expression`, as the snapshot gives it:
    the Atlas 25.1 the twenty and those after them were baked from."""
    pinned = inputs.load()
    assert pinned.sha256 == inputs.INPUTS_SHA256
    assert pinned.created == inputs.INPUTS_CREATED
    assert pinned.atlas["version"] == "25.1"
    stored = sorted((DATA / "assets" / "locus").glob("*_locus.json"))
    if not stored:
        pytest.skip("no stored locus tracks: run pipeline/fetch_tracks.py")
    for path in stored:
        track = json.loads(path.read_bytes())
        expression = dict(track["expression"])
        ensembl = expression.pop("ensembl_gene")
        assert bake_locus.expression_of(pinned.atlas_gene(ensembl)) == expression, path.name


def test_the_snapshot_folder_is_ignored_by_git():
    ignore = (BACKEND / ".gitignore").read_text(encoding="utf-8").splitlines()
    assert "pipeline/data/" in ignore
    assert inputs.FOLDER.is_relative_to(BACKEND / "pipeline" / "data")


def test_the_bake_fetches_nothing_but_uniprot():
    """The Atlas, UCSC and MANE come from the snapshot alone: the bake names
    no fetch of them."""
    source = Path(bake_locus.__file__).read_text(encoding="utf-8")
    main = source[source.index("def main("):]
    assert "fetch_table(" not in main and "current_summary" not in source
    assert "proteinatlas.org/about" not in source
