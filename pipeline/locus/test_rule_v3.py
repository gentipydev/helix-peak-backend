"""Rule v3, the levels it reads, and the genes no zoom can be drawn for.

The user's decision of 2026-11-13: where a gene has only tissue cell type
pairs, the pair in the tissue its RNA is highest in, not the Atlas's
alphabetical first; where it has no reading at all, the tissue its RNA is
highest in and the cell type living there it is highest in; and a gene read
in no tissue is refused. Rule v2's first two ways stay as they are, so every
curated protein keeps its path.
"""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from app.protein_index import parse_mane  # noqa: E402
from pipeline.locus import bake_locus, check_locus  # noqa: E402
from pipeline.locus import inputs as locus_inputs  # noqa: E402
from pipeline.locus.testing import FIXTURES, tables  # noqa: E402
from pipeline.locus.tissues import TISSUES  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import Source, Target  # noqa: E402

LOW = "Low tissue specificity"


def reading(*, tissues=(), cells=(), pairs=(), levels=None, cell_levels=None) -> dict:
    """An Atlas reading as `expression_of` makes one, with the levels v3 reads."""
    expression = {
        "tissue": {"specificity": "Tissue enhanced" if tissues else LOW,
                   "distribution": "Detected in many",
                   "specific": [{"name": name, "ntpm": level} for name, level in tissues]},
        "cell_type": {"specificity": "Cell type enhanced" if cells else "Low cell type specificity",
                      "distribution": "Detected in many",
                      "specific": [{"name": name, "ncpm": level} for name, level in cells]},
        "tissue_cell_type": [{"tissue": tissue, "cell_type": cell} for tissue, cell in pairs],
        "subcellular": {"main": [], "additional": []},
        "secretome": None,
    }
    expression["tissue_levels"] = {**{tissue: 0 for tissue in TISSUES}, **(levels or {})}
    if cell_levels is not None:
        expression["cell_levels"] = cell_levels
    return expression


def _cells_of(tissue: str, **found) -> dict:
    """Every cell type of [tissue] at zero but those given."""
    return {cell: found.get(cell, 0) for cell in bake_locus.cell_candidates(tissue)}


# -- the first two ways, as v2 takes them -------------------------------------------------


def test_a_tissue_reading_goes_as_v2_takes_it():
    expression = reading(tissues=[("pancreas", 900.0)], pairs=[("Adrenal gland", "Adrenal medulla cells"),
                                                               ("Pancreas", "Beta cells")],
                         levels={"adrenal gland": 5000})
    assert bake_locus.path_of_v3(expression) == bake_locus.path_of_v2(expression)
    assert bake_locus.path_of_v3(expression)["tissue"] == "pancreas"


def test_a_cell_type_home_goes_as_v2_takes_it():
    expression = reading(cells=[("Late primary spermatocytes", 300.0)], levels={"brain": 900})
    path = bake_locus.path_of_v3(expression)
    assert path == bake_locus.path_of_v2(expression)
    assert (path["tissue"], path["tissue_from"]) == ("testis", "cell_type")


# -- pairs only: the pair in the tissue the RNA is highest in -------------------------------


def test_pairs_only_go_to_the_pair_in_the_highest_tissue():
    # v2 took the first, the Atlas's alphabet: adipose tissue.
    pairs = [("Adipose subcutaneous", "Adipocytes"), ("Adrenal gland", "Adrenal cortex cells"),
             ("Skin", "Melanocytes")]
    expression = reading(pairs=pairs, levels={"adipose tissue": 20, "adrenal gland": 5,
                                              "skin": 31.4})
    assert bake_locus.path_of_v2(expression)["tissue"] == "adipose tissue"
    path = bake_locus.path_of_v3(expression)
    assert (path["tissue"], path["cell_type"], path["tissue_from"], path["cell_from"]) == \
        ("skin", "Melanocytes", "tissue_cell_type", "tissue_cell_type")
    assert path["cell_class"] == "Pigment cells"


def test_a_single_pair_is_kept_whatever_the_levels():
    # p53's case: one pair, the stomach's mitotic cells, kept.
    expression = reading(pairs=[("Stomach", "Mitotic cells (Stomach)")],
                         levels={"lymphoid tissue": 80, "stomach": 12})
    path = bake_locus.path_of_v3(expression)
    assert (path["tissue"], path["cell_type"]) == ("stomach", "Mitotic cells (Stomach)")
    assert path == bake_locus.path_of_v2(expression)


def test_level_pairs_go_to_the_earlier_in_the_atlas_order():
    pairs = [("Kidney", "Proximal tubular cells"), ("Liver", "Hepatocytes")]
    expression = reading(pairs=pairs, levels={"kidney": 40, "liver": 40})
    assert bake_locus.path_of_v3(expression)["tissue"] == "kidney"


def test_two_pairs_in_one_tissue_take_the_first():
    pairs = [("Adipose subcutaneous", "Adipocytes"), ("Adipose visceral", "Adipose progenitor cells")]
    expression = reading(pairs=pairs, levels={"adipose tissue": 3})
    assert bake_locus.path_of_v3(expression)["cell_type"] == "Adipocytes"


# -- no reading: the highest tissue and its highest cell type ------------------------------


def test_no_reading_goes_to_the_highest_tissue_and_its_highest_cell():
    cells = _cells_of("bone marrow", **{"Neutrophils": 40.2, "Monocytes": 12.0})
    expression = reading(levels={"bone marrow": 75.3, "lymphoid tissue": 70.1},
                         cell_levels=cells)
    assert bake_locus.path_of_v2(expression)["tissue"] is None
    path = bake_locus.path_of_v3(expression)
    assert (path["tissue"], path["tissue_from"]) == ("bone marrow", "tissue_level")
    assert (path["cell_type"], path["cell_from"]) == ("Neutrophils", "single_cell_level")
    assert path["cell_class"] == "Blood and immune cells"


def test_no_reading_with_no_cell_read_there_draws_the_tissue_own_cells():
    expression = reading(levels={"skeletal muscle": 2.0},
                         cell_levels=_cells_of("skeletal muscle"))
    path = bake_locus.path_of_v3(expression)
    assert (path["tissue"], path["tissue_from"]) == ("skeletal muscle", "tissue_level")
    assert (path["cell_type"], path["cell_from"], path["cell_class"]) == (None, None, None)


def test_a_cell_with_no_nucleus_still_lands_in_its_precursor():
    expression = reading(levels={"bone marrow": 5},
                         cell_levels=_cells_of("bone marrow", Platelets=50.0, Neutrophils=3.0))
    path = bake_locus.path_of_v3(expression)
    assert path["cell_type"] == "Platelets"
    assert path["lands_in"] == {"cell": "megakaryocytes", "place": "bone marrow", "why": "no nucleus"}


def test_level_tissues_and_cells_go_to_the_later_name():
    # The planning census's rule, kept so its 396 changes are this rule's.
    expression = reading(levels={"brain": 10, "liver": 10},
                         cell_levels=_cells_of("liver", Hepatocytes=4, Cholangiocytes=4))
    path = bake_locus.path_of_v3(expression)
    assert (path["tissue"], path["cell_type"]) == ("liver", "Hepatocytes")


def test_a_gene_read_nowhere_has_no_path():
    path = bake_locus.path_of_v3(reading())
    assert path["tissue"] is None


# -- the levels, read ---------------------------------------------------------------------


def test_a_group_level_is_its_members_highest():
    levels = bake_locus.tissue_levels({"cerebellum": 3.0, "hypothalamus": 9.5, "colon": 2,
                                       "small intestine": 4, "spleen": 1, "liver": 7,
                                       "stomach 1": 2.5})
    assert list(levels) == list(TISSUES)
    assert (levels["brain"], levels["intestine"], levels["lymphoid tissue"]) == (9.5, 4, 1)
    assert (levels["liver"], levels["stomach"], levels["kidney"]) == (7, 2.5, 0)
    with pytest.raises(LookupError):
        bake_locus.tissue_levels({"tail fin": 1})


def test_the_cell_levels_are_the_cell_types_of_the_tissue():
    levels = bake_locus.cell_levels({"neutrophils": 4.0, "hepatocytes": 9.0}, "bone marrow")
    assert list(levels) == bake_locus.cell_candidates("bone marrow")
    assert levels["Neutrophils"] == 4.0 and "Hepatocytes" not in levels
    with pytest.raises(LookupError):
        bake_locus.cell_levels({"dragon cells": 1}, "bone marrow")


def test_cell_levels_are_kept_only_where_the_rule_goes_by_level():
    by_level = reading(levels={})
    kept = bake_locus.readings_v3(by_level, ({"liver": 3.0}, {"hepatocytes": 1.0}))
    assert list(kept["cell_levels"]) == bake_locus.cell_candidates("liver")
    with_pairs = reading(pairs=[("Liver", "Hepatocytes")])
    assert "cell_levels" not in bake_locus.readings_v3(with_pairs, ({"liver": 3.0}, {}))
    nowhere = reading()
    assert "cell_levels" not in bake_locus.readings_v3(nowhere, ({}, {}))


# -- a track, refused or checked --------------------------------------------------------


def _mane():
    return parse_mane(open(FIXTURES / "mane_v1.5_six.tsv"))


def _ins_row() -> dict:
    found = _mane()["ENST00000381330"]
    return {"uniprot": "P01308", "gene": "INS", "name": "Insulin", "hgnc_id": found["HGNC_ID"],
            "refseq_nuc": found["RefSeq_nuc"], "refseq_prot": found["RefSeq_prot"],
            "ensembl_nuc": found["Ensembl_nuc"], "chrom_acc": found["GRCh38_chr"],
            "chrom_start": int(found["chr_start"]), "chrom_end": int(found["chr_end"]),
            "chrom_strand": -1}


def _gene(**fields) -> dict:
    """INS's Atlas entry, its reading replaced by [fields]."""
    return {"Ensembl": "ENSG00000254647", "Gene": "INS", "Uniprot": ["P01308"],
            "RNA tissue specificity": "Not detected", "RNA tissue distribution": "Not detected",
            "RNA tissue specific nTPM": None, "RNA single cell type specificity": "Not detected",
            "RNA single cell type distribution": "Not detected",
            "RNA single cell type specific nCPM": None, "RNA tissue cell type enrichment": None,
            "Subcellular main location": None, "Subcellular additional location": None,
            "Secretome location": None, **fields}


def _track_of(gene: dict, levels: tuple) -> dict:
    t = tables()
    return bake_locus.track_of(
        slug="ins", gene="INS", uniprot="P01308", row=_ins_row(), atlas_gene=gene,
        match=("ENSG00000254647", "mane"), cytoband=bake_locus.table_of(
            json.loads((FIXTURES / "cytoband_hg38_five.json").read_text()), "u"),
        aliases=t["chromAlias"], atlas=bake_locus.Atlas("25.1", "2026-05-25", "109"),
        mane_release="v1.5", uniprot_source={"release": "2026_04"}, retrieved="2026-11-13",
        levels=levels, inputs_source={"sha256": "b" * 64, "created": "2026-11-13"})


INS = Target(slug="ins", gene="INS", uniprot="P01308", display="Insulin",
             source=Source("NC_000011.10", 2159779, 2161209, protein_id="NP_000198.1",
                           transcript_id="NM_000207.3"),
             structure=None, cleaved=False, chain_label="")


def test_a_gene_read_in_no_tissue_is_refused():
    with pytest.raises(bake_locus.LocusRefused) as refused:
        _track_of(_gene(), ({}, {}))
    assert str(refused.value) == "The Human Protein Atlas detects this gene's RNA in no tissue."


def test_a_gene_with_no_reading_is_baked_by_level_and_checked():
    track = _track_of(_gene(**{"RNA tissue specificity": LOW,
                               "RNA tissue distribution": "Detected in all",
                               "RNA single cell type specificity": "Low cell type specificity",
                               "RNA single cell type distribution": "Detected in all"}),
                      ({"liver": 12.5, "colon": 30.0, "duodenum": 41.0},
                       {"enterocytes": 22.0, "paneth cells": 22.0, "goblet cells": 3.0}))
    path = track["path"]
    assert (path["tissue"], path["tissue_from"]) == ("intestine", "tissue_level")
    # Enterocytes and Paneth cells are level: the later name.
    assert (path["cell_type"], path["cell_from"]) == ("Paneth cells", "single_cell_level")
    assert track["expression"]["tissue_levels"]["intestine"] == 41.0
    assert list(track["expression"]["cell_levels"]) == bake_locus.cell_candidates("intestine")
    assert check_locus.problems_of(INS, track) == []


@pytest.mark.parametrize("damage, finding", [
    (lambda e: e["tissue_levels"].pop("liver"), "every tissue"),
    (lambda e: e.update(tissue_levels=dict(reversed(list(e["tissue_levels"].items())))),
     "every tissue"),
    (lambda e: e["tissue_levels"].update(liver=-1), "not a level"),
    (lambda e: e.pop("cell_levels"), "cell types that live in intestine"),
    (lambda e: e["cell_levels"].pop("Enterocytes"), "cell types that live in intestine"),
    (lambda e: e["cell_levels"].update({"Paneth cells": 99.0}), "the rule takes"),
    (lambda e: e["tissue_levels"].update(liver=99.0), "the rule takes"),
])
def test_the_check_holds_the_levels_to_the_rule(damage, finding):
    track = copy.deepcopy(_track_of(_gene(), ({"duodenum": 41.0}, {"enterocytes": 22.0})))
    damage(track["expression"])
    found = check_locus.problems_of(INS, track)
    assert any(finding in problem for problem in found), found


def test_cell_levels_where_the_rule_does_not_go_by_level_are_caught():
    track = copy.deepcopy(_track_of(
        _gene(**{"RNA tissue cell type enrichment": ["Pancreas - Beta cells"]}),
        ({"pancreas": 50.0}, {})))
    track["expression"]["cell_levels"] = {"Hepatocytes": 1.0}
    assert any("does not go by level" in p for p in check_locus.problems_of(INS, track))


# -- every curated protein keeps its path ---------------------------------------------------


_PINNED = locus_inputs.default_path()
_STORED = sorted((DATA / "assets" / "locus").glob("*_locus.json"))


@pytest.mark.skipif(not _PINNED.is_file() or not _STORED,
                    reason="needs the pinned snapshot and the stored curated tracks")
def test_rule_v3_gives_every_curated_protein_its_stored_path():
    snapshot = locus_inputs.load()
    for path in _STORED:
        track = json.loads(path.read_bytes())
        assert track["schema_version"] == 2
        expression = dict(track["expression"])
        expression.update(bake_locus.readings_v3(
            expression, snapshot.levels(expression["ensembl_gene"])))
        assert bake_locus.path_of_v3(expression) == track["path"], path.name
