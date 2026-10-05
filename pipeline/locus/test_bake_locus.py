"""The locus bake and its check, offline.

The fixtures are the sources' own answers, cut down: UniProt's entries for
six of the twenty (release 2026_03), keeping what the protein index's join
reads; their six MANE Select rows from the v1.5 summary; UCSC's cytoBand
and chromAlias tables for hg38, as its API served them, for the five
chromosomes those genes are on; and the Human Protein Atlas's summaries of the
six genes (version 25.1), keeping the fields the bake reads (since schema 2,
its tissue cell type pairs, subcellular and secretome locations too), with the
text of its releases page that the parser reads. Nothing here fetches.
"""

from __future__ import annotations

import copy
import io
import json
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from app.protein_index import parse_mane  # noqa: E402
from pipeline import fetch_tracks, upload_tracks  # noqa: E402
from pipeline.locus import bake_locus, cell_types, check_locus, tissues, verify_locus  # noqa: E402
from pipeline.targets import BY_SLUG  # noqa: E402
from pipeline.uniprot import Entry  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures"

# The places the plan names to hold the bake to, as the literature gives
# them. They are expectations here, never inputs to the bake.
KNOWN = {
    "insulin": "11p15.5",
    "hemoglobin": "11p15.4",
    "p53": "17p13.1",
    "cftr": "7q31.2",
    "dystrophin": "Xp21",
    "sod1": "21q22",
}


def _table(name: str) -> bake_locus.Table:
    answer = json.loads((FIXTURES / name).read_text(encoding="utf-8"))
    return bake_locus.table_of(answer, "https://api.genome.ucsc.edu/test")


@pytest.fixture(scope="module")
def sources():
    mane = parse_mane(io.StringIO((FIXTURES / "mane_v1.5_six.tsv").read_text()))
    return (mane, _table("cytoband_hg38_five.json"), _table("chromalias_hg38_five.json"))


ATLAS = bake_locus.atlas_of((FIXTURES / "hpa_releases.html").read_text(encoding="utf-8"))


def _atlas_gene(ensembl: str) -> dict:
    return json.loads((FIXTURES / "hpa" / f"{ensembl}.json").read_text(encoding="utf-8"))


def _entry(accession: str) -> Entry:
    body = json.loads((FIXTURES / "uniprot" / f"{accession}.json").read_text(encoding="utf-8"))
    return Entry(accession, body, "2026_03", "2026-09-02", "2026-09-27")


def _track(slug: str, sources) -> dict:
    mane, cytoband, aliases = sources
    target = BY_SLUG[slug]
    entry = _entry(target.uniprot)
    ensembl = bake_locus.ensembl_gene_of(
        bake_locus.mane_row(target, entry, mane, "v1.5"), mane)
    return bake_locus.payload(target, entry, mane, "v1.5", cytoband, aliases,
                              ATLAS, _atlas_gene(ensembl), "2026-09-27")


# -- the places -----------------------------------------------------------------


@pytest.mark.parametrize("slug", sorted(KNOWN))
def test_each_gene_lies_in_the_band_the_literature_gives_it(slug, sources):
    track = _track(slug, sources)
    # The known place, read the way HGNC's is: the track's own place lies in
    # it, whether it is written as finely (11p15.5) or more coarsely (21q22).
    assert verify_locus.within(track, KNOWN[slug]) is True
    assert track["locus"].startswith(KNOWN[slug])


def test_a_gene_across_two_bands_is_written_across_them(sources):
    track = _track("dystrophin", sources)
    assert track["band"]["names"] == ["p21.2", "p21.1"]
    assert track["locus"] == "Xp21.2-p21.1"
    bands = {b["name"]: b for b in track["bands"]}
    assert track["band"]["start"] == bands["p21.2"]["start"]
    assert track["band"]["end"] == bands["p21.1"]["end"]
    # Its 2.1 million bases, on the minus strand, inside those two bands.
    span = track["span"]
    assert span["end"] - span["start"] + 1 > 2_000_000
    assert span["strand"] == -1
    assert track["band"]["start"] <= span["start"] <= span["end"] <= track["band"]["end"]


def test_the_span_is_the_protein_index_rows_own(sources):
    # The join the index makes: the entry's MANE Select cross-reference, the
    # MANE row for the table's gene, its GRCh38 span.
    track = _track("insulin", sources)
    assert track["sequence"] == "NC_000011.10"
    assert track["span"] == {"start": 2159779, "end": 2161209, "strand": -1}
    assert track["transcript"] == {"refseq": "NM_000207.3", "ensembl": "ENST00000381330.5"}


def test_positions_are_made_one_based(sources):
    mane, cytoband, aliases = sources
    track = _track("insulin", sources)
    first = min((row for row in cytoband.rows if row["chrom"] == "chr11"),
                key=lambda row: row["chromStart"])
    assert first["chromStart"] == 0
    assert track["bands"][0] == {"name": first["name"], "start": 1,
                                 "end": first["chromEnd"], "stain": first["gieStain"]}


def test_the_chromosome_is_named_by_ucsc_own_alias_table(sources):
    _, _, aliases = sources
    assert bake_locus.chromosome_of("NC_000023.11", aliases) == "chrX"
    assert bake_locus.chromosome_of("NC_000011.10", aliases) == "chr11"
    with pytest.raises(LookupError):
        bake_locus.chromosome_of("NC_000011.9", aliases)


# -- what it refuses ---------------------------------------------------------------


def test_refuses_an_entry_with_no_row_for_the_gene(sources):
    mane, cytoband, aliases = sources
    target = BY_SLUG["insulin"]
    body = copy.deepcopy(_entry(target.uniprot).body)
    body["uniProtKBCrossReferences"] = [
        x for x in body["uniProtKBCrossReferences"] if x["database"] != "MANE-Select"]
    with pytest.raises(LookupError, match="0 buildable MANE Select rows"):
        bake_locus.payload(target, Entry(target.uniprot, body, "2026_03", "2026-09-02",
                                         "2026-09-27"),
                           mane, "v1.5", cytoband, aliases, ATLAS,
                           _atlas_gene("ENSG00000254647"), "2026-09-27")


def test_refuses_an_entry_uniprot_did_not_answer_for(sources):
    mane, cytoband, aliases = sources
    with pytest.raises(ValueError, match="UniProt answered"):
        bake_locus.payload(BY_SLUG["hemoglobin"], _entry("P01308"), mane, "v1.5",
                           cytoband, aliases, ATLAS, _atlas_gene("ENSG00000244734"),
                           "2026-09-27")


def test_refuses_a_chromosome_with_no_bands(sources):
    mane, cytoband, aliases = sources
    bare = bake_locus.Table("cytoBand", "hg38", cytoband.updated,
                            tuple(r for r in cytoband.rows if r["chrom"] != "chr11"),
                            cytoband.url)
    with pytest.raises(LookupError, match="no bands"):
        bake_locus.payload(BY_SLUG["insulin"], _entry("P01308"), mane, "v1.5",
                           bare, aliases, ATLAS, _atlas_gene("ENSG00000254647"),
                           "2026-09-27")


# -- where in the body, as the Atlas reads it ------------------------------------------


def test_names_where_the_atlas_finds_the_gene_read(sources):
    hemoglobin = _track("hemoglobin", sources)["expression"]
    assert hemoglobin["ensembl_gene"] == "ENSG00000244734"
    assert hemoglobin["tissue"]["specificity"] == "Tissue enriched"
    assert hemoglobin["tissue"]["specific"][0]["name"] == "bone marrow"
    # The Atlas's own reading: the cells it finds the RNA in are red cells,
    # which is what makes the zoom land in a precursor.
    assert hemoglobin["cell_type"]["specific"][0]["name"] == "Erythrocytes"
    insulin = _track("insulin", sources)["expression"]
    assert insulin["tissue"]["specific"][0]["name"] == "pancreas"
    assert insulin["cell_type"]["specific"][0]["name"] == "Pancreatic islet cells"


def test_a_gene_read_everywhere_names_no_tissue(sources):
    p53 = _track("p53", sources)["expression"]
    assert p53["tissue"] == {"specificity": "Low tissue specificity",
                             "distribution": "Detected in all", "specific": []}
    assert p53["cell_type"]["specific"] == []
    dystrophin = _track("dystrophin", sources)["expression"]
    assert dystrophin["tissue"]["specific"] == []
    assert dystrophin["cell_type"]["specific"], "specific to cell types, not tissues"


def test_names_them_highest_first(sources):
    for slug in KNOWN:
        expression = _track(slug, sources)["expression"]
        for key, unit in (("tissue", "ntpm"), ("cell_type", "ncpm")):
            levels = [found[unit] for found in expression[key]["specific"]]
            assert levels == sorted(levels, reverse=True)


def test_refuses_an_atlas_answer_for_another_gene(sources):
    mane, cytoband, aliases = sources
    with pytest.raises(ValueError, match="asked the Atlas"):
        bake_locus.payload(BY_SLUG["insulin"], _entry("P01308"), mane, "v1.5",
                           cytoband, aliases, ATLAS, _atlas_gene("ENSG00000244734"),
                           "2026-09-27")


def test_reads_the_atlas_release_from_its_releases_page():
    assert ATLAS == bake_locus.Atlas("25.1", "2026-05-25", "109")
    provenance = ATLAS.provenance("ENSG00000254647")
    assert provenance["licence"] == "CC BY 4.0"
    assert provenance["url"] == "https://www.proteinatlas.org/ENSG00000254647.json"
    with pytest.raises(LookupError):
        bake_locus.atlas_of("<html>no release here</html>")


# -- schema 2: the readings it adds, and the zoom's path --------------------------------


def test_schema_two_adds_the_atlas_pairs_locations_and_secretome(sources):
    assert bake_locus.SCHEMA_VERSION == 2
    insulin = _track("insulin", sources)
    assert insulin["schema_version"] == 2
    expression = insulin["expression"]
    assert expression["tissue_cell_type"] == [
        {"tissue": "Adrenal gland", "cell_type": "Adrenal medulla cells"},
        {"tissue": "Pancreas", "cell_type": "Beta cells"},
    ]
    assert expression["subcellular"] == {"main": [], "additional": []}
    assert expression["secretome"] == "Secreted to blood"
    p53 = _track("p53", sources)["expression"]
    assert p53["subcellular"]["main"] == ["Nucleoplasm"]
    assert p53["secretome"] is None
    # A schema 1 reader finds every key it read where it was.
    for key in ("ensembl_gene", "tissue", "cell_type"):
        assert key in expression


PATHS = {
    # slug: (tissue, tissue_from, cell_type, cell_from, lands_in)
    "insulin": ("pancreas", "tissue", "Beta cells", "tissue_cell_type", None),
    "hemoglobin": ("bone marrow", "tissue", "Erythrocytes", "single_cell_type",
                   {"cell": "erythroblasts", "place": "bone marrow", "why": "no nucleus"}),
    "cftr": ("pancreas", "tissue", "Pancreatic duct cells", "single_cell_type", None),
    "sod1": ("liver", "tissue", "Hepatocytes", "tissue_cell_type", None),
    "dystrophin": ("skeletal muscle", "cell_type", "Myonuclei", "single_cell_type", None),
    "p53": ("stomach", "tissue_cell_type", "Mitotic cells (Stomach)", "tissue_cell_type",
            None),
}


@pytest.mark.parametrize("slug", sorted(PATHS))
def test_each_gene_takes_one_path_through_an_organ_and_a_cell_of_it(slug, sources):
    path = _track(slug, sources)["path"]
    tissue, tissue_from, cell, cell_from, lands = PATHS[slug]
    assert (path["tissue"], path["tissue_from"]) == (tissue, tissue_from)
    assert (path["cell_type"], path["cell_from"]) == (cell, cell_from)
    assert path["lands_in"] == lands
    assert path["cell_class"] in cell_types.CLASSES
    # The cell lives in the organ: by its home, or by the Atlas's own pair.
    if cell_from == "single_cell_type":
        assert tissue in cell_types.homes(cell)


def _reading(tissues=(), cells=(), pairs=()):
    return {
        "tissue": {"specific": [{"name": n, "ntpm": 10.0} for n in tissues]},
        "cell_type": {"specific": [{"name": n, "ncpm": 10.0} for n in cells]},
        "tissue_cell_type": [{"tissue": t, "cell_type": c} for t, c in pairs],
    }


def test_the_rule_takes_a_pair_in_the_organ_before_a_single_cell_type():
    path = bake_locus.path_of(_reading(
        tissues=["pancreas"], cells=["Pancreatic islet cells"],
        pairs=[("Liver", "Hepatocytes"), ("Pancreas", "Beta cells")]))
    assert (path["tissue"], path["cell_type"], path["cell_from"]) == (
        "pancreas", "Beta cells", "tissue_cell_type")
    assert path["cell_class"] == "Endocrine cells"


def test_the_rule_never_takes_a_cell_of_another_organ():
    # Bone marrow first, microglia the only cell type: microglia live in the
    # brain, so no cell is named and the app draws the marrow's own.
    path = bake_locus.path_of(_reading(
        tissues=["bone marrow", "lymphoid tissue"], cells=["Microglia"],
        pairs=[("Breast", "T-cells")]))
    assert path == {"tissue": "bone marrow", "tissue_from": "tissue", "cell_type": None,
                    "cell_class": None, "cell_from": None, "lands_in": None}


def test_a_sampled_tissue_is_the_tissue_without_its_sample_number():
    path = bake_locus.path_of(_reading(tissues=["stomach 1"], cells=["Foveolar cells"]))
    assert (path["tissue"], path["cell_type"]) == ("stomach", "Foveolar cells")


def test_a_gene_specific_to_no_tissue_goes_to_its_cell_type_home():
    # The first cell type with a home: lacrimal acinar cells live where the
    # consensus reading does not sample, so the next one leads.
    path = bake_locus.path_of(_reading(
        cells=["Lacrimal acinar cells", "Late primary spermatocytes"]))
    assert (path["tissue"], path["tissue_from"], path["cell_type"]) == (
        "testis", "cell_type", "Late primary spermatocytes")


def test_a_gene_specific_to_nothing_follows_its_first_pair_or_nowhere():
    path = bake_locus.path_of(_reading(pairs=[("Stomach", "Mitotic cells (Stomach)")]))
    assert (path["tissue"], path["cell_type"], path["cell_class"]) == (
        "stomach", "Mitotic cells (Stomach)", "Stem and proliferating cells")
    assert bake_locus.path_of(_reading()) == {
        "tissue": None, "tissue_from": None, "cell_type": None, "cell_class": None,
        "cell_from": None, "lands_in": None}


def test_a_cell_with_no_nucleus_lands_in_the_one_that_has_one():
    path = bake_locus.path_of(_reading(tissues=["bone marrow"], cells=["Platelets"]))
    assert path["lands_in"] == {"cell": "megakaryocytes", "place": "bone marrow",
                                "why": "no nucleus"}


@pytest.mark.parametrize("reading, word", [
    (dict(tissues=["thymus gland"]), "tissue 'thymus gland'"),
    (dict(cells=["Glow cells"]), "cell type 'Glow cells'"),
    (dict(pairs=[("Pancreas islets", "Beta cells")]), "tissue cell type tissue"),
    (dict(pairs=[("Pancreas", "Moon cells")]), "has no class"),
])
def test_the_rule_refuses_a_name_the_atlas_does_not_list(reading, word):
    with pytest.raises(LookupError, match=word):
        bake_locus.path_of(_reading(**reading))


def test_the_pairs_are_read_as_the_atlas_writes_them():
    assert bake_locus._pairs(["Adipose subcutaneous - Adipocytes (Subcutaneous)"]) == [
        {"tissue": "Adipose subcutaneous", "cell_type": "Adipocytes (Subcutaneous)"}]
    with pytest.raises(ValueError, match="Tissue - Cell type"):
        bake_locus._pairs(["Pancreas: Beta cells"])


def test_the_vocabularies_are_the_atlas_own():
    assert len(cell_types.CELL_TYPES) == 154
    assert len(cell_types.CLASSES) == 15
    assert {kind for kind, _ in cell_types.CELL_TYPES.values()} == set(cell_types.CLASSES)
    assert set(cell_types.PAIR_CELL_CLASSES.values()) <= set(cell_types.CLASSES)
    assert len(tissues.TISSUES) == 37
    for kind, homes in cell_types.CELL_TYPES.values():
        assert set(homes) <= set(tissues.TISSUES)
    assert set(tissues.PAIR_TISSUES.values()) <= set(tissues.TISSUES)
    # Case is the Atlas's, not the reader's.
    assert cell_types.single("monocytes") == "Monocytes"
    assert cell_types.cell_class("cDC") == "Blood and immune cells"


@pytest.mark.parametrize("damage, finding", [
    (lambda t: t["path"].update(cell_type="Pancreatic islet cells"), "the rule takes"),
    (lambda t: t.pop("path"), "no path"),
    (lambda t: t["expression"]["subcellular"]["main"].append("Cytoplasm"), "Atlas's word"),
    (lambda t: t["expression"].update(secretome="Secreted to lymph"), "Atlas's word"),
    (lambda t: t["expression"]["tissue_cell_type"].append(
        {"tissue": "Pancreas islets", "cell_type": "Beta cells"}), "one the Atlas lists"),
    (lambda t: t["expression"]["cell_type"]["specific"][0].update(name="Glow cells"),
     "not a name the Atlas lists"),
])
def test_the_check_holds_schema_two_to_the_atlas_and_the_rule(damage, finding, sources):
    track = _track("insulin", sources)
    damage(track)
    found = check_locus.problems_of(BY_SLUG["insulin"], track)
    assert any(finding in problem for problem in found), found


# -- the table and its version -------------------------------------------------------


def test_the_table_is_read_flat_and_named_by_its_version(sources):
    _, cytoband, _ = sources
    assert cytoband.track == "cytoBand"
    assert cytoband.genome == "hg38"
    assert cytoband.updated == "2022-10-28T15:06:46"
    assert {row["chrom"] for row in cytoband.rows} == {"chr7", "chr11", "chr17", "chr21", "chrX"}
    provenance = cytoband.provenance()
    assert provenance["rows"] == len(cytoband.rows)
    assert len(provenance["sha256"]) == 64


def test_the_digest_is_the_rows_not_their_order(sources):
    _, cytoband, _ = sources
    shuffled = bake_locus.Table(cytoband.track, cytoband.genome, cytoband.updated,
                                tuple(reversed(cytoband.rows)), cytoband.url)
    assert shuffled.sha256 == cytoband.sha256
    changed = [dict(row) for row in cytoband.rows]
    changed[0]["gieStain"] = "gpos100" if changed[0]["gieStain"] != "gpos100" else "gneg"
    edited = bake_locus.Table(cytoband.track, cytoband.genome, cytoband.updated,
                              tuple(changed), cytoband.url)
    assert edited.sha256 != cytoband.sha256


# -- the check -----------------------------------------------------------------------


@pytest.mark.parametrize("slug", sorted(KNOWN))
def test_a_true_track_passes_the_check(slug, sources):
    assert check_locus.problems_of(BY_SLUG[slug], _track(slug, sources)) == []


@pytest.mark.parametrize("damage, finding", [
    (lambda t: t.update(locus="11p15.4"), "locus"),
    (lambda t: t["band"].update(names=["p15.4"]), "band names"),
    (lambda t: t["bands"].pop(3), "expected to start at"),
    (lambda t: t["bands"][2].update(stain="blue"), "stain"),
    (lambda t: t["span"].update(end=t["length"] + 5), "not on the chromosome"),
    (lambda t: t.update(gene="HBB"), "gene is"),
    (lambda t: t["sources"]["cytoband"].update(updated="yesterday"), "updated"),
    (lambda t: t["sources"]["mane"].update(release="latest"), "MANE release"),
    (lambda t: t["expression"]["tissue"].update(specificity="Everywhere"), "specificity"),
    (lambda t: t["expression"]["cell_type"].update(
        specific=[{"name": "A", "ncpm": 1.0}, {"name": "B", "ncpm": 9.0}]), "highest first"),
    (lambda t: t["sources"]["hpa"].update(licence="CC BY-SA 4.0"), "licence"),
])
def test_the_check_catches_a_damaged_track(damage, finding, sources):
    track = _track("insulin", sources)
    damage(track)
    found = check_locus.problems_of(BY_SLUG["insulin"], track)
    assert any(finding in problem for problem in found), found


def test_the_check_holds_the_atlas_categories_to_their_lists(sources):
    track = _track("p53", sources)
    track["expression"]["tissue"]["specific"] = [{"name": "liver", "ntpm": 3.0}]
    found = check_locus.problems_of(BY_SLUG["p53"], track)
    assert any("Low tissue specificity, yet names 1" in problem for problem in found)
    track = _track("insulin", sources)
    track["expression"]["tissue"]["specific"] = []
    found = check_locus.problems_of(BY_SLUG["insulin"], track)
    assert any("Tissue enriched, yet names none" in problem for problem in found)


def test_the_check_holds_a_sliced_record_to_its_slice(sources):
    # Glucagon's record is cut from chromosome 2 around its gene; a track
    # placing it elsewhere, or on another sequence, disagrees with it.
    target = BY_SLUG["glucagon"]
    assert target.source.accession.startswith("NC_")
    track = _track("insulin", sources)
    track.update(slug="glucagon", gene=target.gene, uniprot=target.uniprot)
    found = check_locus.problems_of(target, track)
    assert any("the record is cut from" in problem for problem in found)
    assert any("outside the record's slice" in problem for problem in found)


# -- HGNC's places, read -----------------------------------------------------------------


def test_reads_hgnc_places_coarse_fine_and_across(sources):
    sod1 = _track("sod1", sources)
    assert verify_locus.within(sod1, "21q22") is True
    assert verify_locus.within(sod1, "21q22.11") is True
    assert verify_locus.within(sod1, "21q22.3") is False
    assert verify_locus.within(sod1, "20q22") is False
    dmd = _track("dystrophin", sources)
    assert verify_locus.within(dmd, "Xp21.2-p21.1") is True
    assert verify_locus.within(dmd, "Xp21.2") is False
    assert verify_locus.within(dmd, "mitochondria") is None


# -- the plumbing ----------------------------------------------------------------------


def test_the_uploader_and_the_fetcher_know_the_kind(sources, tmp_path, monkeypatch):
    target = BY_SLUG["insulin"]
    assert upload_tracks.asset_of("locus", target).name == "insulin_locus.json"
    assert fetch_tracks.asset_path("locus", target) == bake_locus.locus_asset(target)
    assert "locus" in fetch_tracks.KINDS
    track = _track("insulin", sources)
    provenance = upload_tracks.validate("locus", target, bake_locus.encode(track))
    assert provenance["locus"] == "11p15.5"
    assert provenance["sources"]["cytoband"]["updated"] == "2022-10-28T15:06:46"
    assert provenance["sources"]["hpa"]["version"] == "25.1"
    damaged = dict(track, locus="11p15.4")
    with pytest.raises(ValueError, match="locus"):
        upload_tracks.validate("locus", target, bake_locus.encode(damaged))
