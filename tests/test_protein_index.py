"""Tests for how the protein index is built and how a query meets it.

Pure functions only: the entries below are cut down from UniProt's JSON for
release 2026_03, and the MANE rows from its v1.5 summary, keeping just the
fields the index reads. No network and no database.
"""

import io

import pytest

from app.protein_index import (
    ACCESSION,
    NAME,
    OTHER_NAME,
    SYMBOL,
    SYNONYM,
    WORD,
    bounds,
    curated_terms,
    normalize,
    parse_mane,
    parse_refseqgene,
    rows_for_entry,
    terms_for_row,
)

_MANE_TEXT = "\t".join([
    "#NCBI_GeneID", "Ensembl_Gene", "HGNC_ID", "symbol", "name", "RefSeq_nuc",
    "RefSeq_prot", "Ensembl_nuc", "Ensembl_prot", "MANE_status", "GRCh38_chr",
    "chr_start", "chr_end", "chr_strand",
]) + "\n" + "\n".join("\t".join(row) for row in [
    ["GeneID:3630", "ENSG00000254647.8", "HGNC:6081", "INS", "insulin",
     "NM_000207.3", "NP_000198.1", "ENST00000381330.5", "ENSP00000370731.5",
     "MANE Select", "NC_000011.10", "2159779", "2161209", "-"],
    ["GeneID:3039", "ENSG00000206172.8", "HGNC:4823", "HBA1", "hemoglobin subunit alpha 1",
     "NM_000558.5", "NP_000549.1", "ENST00000320868.9", "ENSP00000322421.5",
     "MANE Select", "NC_000016.10", "176680", "177522", "+"],
    ["GeneID:3040", "ENSG00000188536.13", "HGNC:4824", "HBA2", "hemoglobin subunit alpha 2",
     "NM_000517.6", "NP_000508.1", "ENST00000251595.11", "ENSP00000251595.6",
     "MANE Select", "NC_000016.10", "172876", "173710", "+"],
    ["GeneID:7159", "ENSG00000143514.17", "HGNC:12000", "TP53BP2", "tumor protein p53 binding protein 2",
     "NM_001031685.3", "NP_001026855.2", "ENST00000343537.12", "ENSP00000341957.7",
     "MANE Select", "NC_000001.11", "223779893", "223845966", "-"],
    # A second transcript MANE keeps for a gene's clinical variants; not a gene.
    ["GeneID:1", "ENSG00000000001.1", "HGNC:1", "EXAMPLE1", "example gene 1",
     "NM_000001.1", "NP_000001.1", "ENST00000000001.1", "ENSP00000000001.1",
     "MANE Plus Clinical", "NC_000001.11", "1000", "2000", "+"],
]) + "\n"


@pytest.fixture
def mane():
    return parse_mane(io.StringIO(_MANE_TEXT))


def _entry(accession, name, genes, mane_xrefs=(), length=100, alternative=(),
           contains=(), score=5.0):
    return {
        "primaryAccession": accession,
        "annotationScore": score,
        "proteinExistence": "1: Evidence at protein level",
        "sequence": {"length": length},
        "genes": [
            {"geneName": {"value": gene}, "synonyms": [{"value": s} for s in synonyms]}
            for gene, synonyms in genes
        ],
        "proteinDescription": {
            "recommendedName": {"fullName": {"value": name}},
            "alternativeNames": [
                {"fullName": {"value": full}, "shortNames": [{"value": s} for s in shorts]}
                for full, shorts in alternative
            ],
            "contains": [{"recommendedName": {"fullName": {"value": part}}} for part in contains],
        },
        "uniProtKBCrossReferences": [
            {"database": "MANE-Select", "id": enst, "isoformId": isoform}
            for enst, isoform in mane_xrefs
        ],
    }


INSULIN = _entry("P01308", "Insulin", [("INS", [])],
                 mane_xrefs=[("ENST00000381330.5", None)], length=110,
                 contains=["Insulin B chain", "Insulin A chain"])
HEMOGLOBIN_ALPHA = _entry(
    "P69905", "Hemoglobin subunit alpha", [("HBA1", []), ("HBA2", [])],
    mane_xrefs=[("ENST00000251595.11", None), ("ENST00000320868.9", None)],
    length=142, alternative=[("Alpha-globin", []), ("Hemoglobin alpha chain", [])])
ASPP2 = _entry("Q13625", "Apoptosis-stimulating of p53 protein 2",
               [("TP53BP2", ["ASPP2", "BBP"])],
               mane_xrefs=[("ENST00000343537.12", "Q13625-3")], length=1128,
               alternative=[("Bcl2-binding protein", ["Bbp"])])


# --- normalising ------------------------------------------------------------


@pytest.mark.parametrize("typed, stored", [
    ("Insulin", "insulin"),
    ("  TP53 ", "tp53"),
    ("Insulin-degrading enzyme", "insulin degrading enzyme"),
    ("NG_007114", "ng 007114"),
    ("p53%", "p53"),
    ("α-actinin", "alpha actinin"),
    ("IL-1β", "il 1 beta"),
    ("Ca(2+)/calmodulin", "ca 2 calmodulin"),
    ("Protéine", "proteine"),
    ("%%%", ""),
])
def test_a_query_and_a_term_normalise_the_same_way(typed, stored):
    assert normalize(typed) == stored


def test_the_range_holds_every_term_the_needle_starts_and_nothing_else():
    lo, hi = bounds("ins")
    inside = ["ins", "ins 1", "insr", "insulin", "insz"]
    outside = ["inr", "int", "in", "i", "ionosphere"]
    assert all(lo <= term < hi for term in inside)
    assert not any(lo <= term < hi for term in outside)
    assert bounds("hba9") == ("hba9", "hba:")


def test_an_empty_needle_has_no_range():
    with pytest.raises(ValueError):
        bounds("")


# --- the sources ------------------------------------------------------------


def test_mane_keeps_select_rows_keyed_by_transcript_without_version(mane):
    assert set(mane) == {"ENST00000381330", "ENST00000320868", "ENST00000251595",
                         "ENST00000343537"}
    assert mane["ENST00000381330"]["RefSeq_prot"] == "NP_000198.1"


def test_a_refseqgene_is_found_in_any_category_and_the_standard_wins():
    text = "\n".join([
        "#tax_id\tGeneID\tSymbol\tRSG\tLRG\tRNA\tt\tProtein\tp\tCategory",
        "9606\t1756\tDMD\tNG_012232.1\tLRG_199\tNM_004006.2\tt1\tNP_003997.1\tp1\treference standard",
        "9606\t1756\tDMD\tNG_012232.1\t\tNM_004006.3\t\tNP_003997.2\t\taligned: Selected",
        "9606\t3630\tINS\tNG_OTHER.1\t\tNM_000207.3\t\tNP_000198.1\t\taligned: Selected",
        "9606\t3630\tINS\tNG_007114.1\t\tNM_000207.3\t\tNP_000198.1\t\treference standard",
    ])
    found = parse_refseqgene(io.StringIO(text))
    assert found["NP_003997.2"] == "NG_012232.1"
    assert found["NP_000198.1"] == "NG_007114.1"


# --- rows -------------------------------------------------------------------


def test_a_canonical_mane_protein_is_one_buildable_row(mane):
    (row,) = rows_for_entry(INSULIN, mane, {"NP_000198.1": "NG_007114.1"}, "v1.5")
    assert row["gene"] == "INS" and row["buildable"]
    assert row["unavailable_reason"] is None
    assert (row["refseq_nuc"], row["refseq_prot"]) == ("NM_000207.3", "NP_000198.1")
    assert (row["chrom_acc"], row["chrom_start"], row["chrom_end"], row["chrom_strand"]) == (
        "NC_000011.10", 2159779, 2161209, -1)
    assert row["refseqgene"] == "NG_007114.1"
    assert row["gene_id"] == 3630 and row["length"] == 110


def test_one_entry_made_by_two_genes_is_two_rows(mane):
    """P69905 is HBA1's protein and HBA2's; each gene is its own walk."""
    rows = rows_for_entry(HEMOGLOBIN_ALPHA, mane, {}, "v1.5")
    assert sorted(row["gene"] for row in rows) == ["HBA1", "HBA2"]
    assert all(row["buildable"] and row["uniprot"] == "P69905" for row in rows)


def test_a_non_canonical_mane_isoform_is_not_buildable_yet(mane):
    (row,) = rows_for_entry(ASPP2, mane, {}, "v1.5")
    assert row["gene"] == "TP53BP2" and not row["buildable"]
    assert row["mane_isoform"] == "Q13625-3"
    assert "Q13625-3" in row["unavailable_reason"]


def test_an_entry_with_no_mane_transcript_is_listed_under_its_own_gene(mane):
    entry = _entry("A6NMY6", "Putative annexin A2-like protein", [("ANXA2P2", [])])
    (row,) = rows_for_entry(entry, mane, {}, "v1.5")
    assert (row["gene"], row["buildable"]) == ("ANXA2P2", False)
    assert row["unavailable_reason"] == "No MANE Select transcript."


def test_an_entry_with_no_gene_is_listed_under_an_empty_gene(mane):
    entry = _entry("Q8N8Q1", "Uncharacterized protein", [])
    (row,) = rows_for_entry(entry, mane, {}, "v1.5")
    assert (row["gene"], row["buildable"]) == ("", False)
    assert "no gene" in row["unavailable_reason"]


def test_a_transcript_mane_has_dropped_is_named_in_the_reason(mane):
    entry = _entry("P99999", "Gone protein", [("GONE", [])],
                   mane_xrefs=[("ENST00000999999.1", None)])
    (row,) = rows_for_entry(entry, mane, {}, "v1.5")
    assert not row["buildable"]
    assert row["unavailable_reason"] == "MANE v1.5 has no transcript ENST00000999999.1."


def test_mane_plus_clinical_is_not_a_gene_row(mane):
    assert "ENST00000000001" not in mane
    entry = _entry("P00001", "Example protein 1", [("EXAMPLE1", [])],
                   mane_xrefs=[("ENST00000000001.1", None)])
    (row,) = rows_for_entry(entry, mane, {}, "v1.5")
    assert not row["buildable"]


# --- terms ------------------------------------------------------------------


def test_a_row_is_found_by_accession_symbol_synonyms_names_and_their_words(mane):
    (row,) = rows_for_entry(ASPP2, mane, {}, "v1.5")
    terms = terms_for_row(ASPP2, row, "tumor protein p53 binding protein 2")
    assert ("q13625", ACCESSION) in terms
    assert ("tp53bp2", SYMBOL) in terms
    assert {("aspp2", SYNONYM), ("bbp", SYNONYM)} <= terms
    # The protein's own name and MANE's name for its gene rank as its name;
    # everything else UniProt calls it ranks below them.
    assert ("apoptosis stimulating of p53 protein 2", NAME) in terms
    assert ("tumor protein p53 binding protein 2", NAME) in terms
    assert ("bcl2 binding protein", OTHER_NAME) in terms
    assert {("apoptosis", WORD), ("stimulating", WORD), ("p53", WORD)} <= terms
    # Stopwords and short words are only ever part of a whole name.
    assert ("protein", WORD) not in terms
    assert ("of", WORD) not in terms


def test_a_word_is_never_stored_as_a_whole_name(mane):
    """Only a whole term may match exactly, so no word may pass as one."""
    (row,) = rows_for_entry(ASPP2, mane, {}, "v1.5")
    terms = terms_for_row(ASPP2, row)
    assert ("apoptosis", NAME) not in terms
    assert ("apoptosis", OTHER_NAME) not in terms


def test_what_a_precursor_is_cut_into_finds_it(mane):
    (row,) = rows_for_entry(INSULIN, mane, {}, "v1.5")
    terms = terms_for_row(INSULIN, row)
    assert ("insulin b chain", OTHER_NAME) in terms


def test_the_other_gene_of_a_shared_protein_is_a_synonym(mane):
    rows = {row["gene"]: row for row in rows_for_entry(HEMOGLOBIN_ALPHA, mane, {}, "v1.5")}
    terms = terms_for_row(HEMOGLOBIN_ALPHA, rows["HBA1"])
    assert ("hba1", SYMBOL) in terms
    assert ("hba2", SYNONYM) in terms
    assert ("hba1", SYNONYM) not in terms


def test_a_listed_protein_lends_its_own_names():
    terms = curated_terms("Hemoglobin (beta chain)", "hemoglobin", "NG_059281")
    assert ("hemoglobin", NAME) in terms
    assert ("hemoglobin beta chain", NAME) in terms
    assert ("ng 059281", ACCESSION) in terms
    assert ("beta", WORD) in terms
