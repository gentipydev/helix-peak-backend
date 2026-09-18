"""Tests for the pure feature extractor.

These run straight off the saved record with no HTTP and no mocking: the parser
takes an already-parsed SeqRecord and does no I/O.
"""

import pytest

from app.genbank_parser import GeneNotFound, extract_gene

PREPROINSULIN = (
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQ"
    "VGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)


@pytest.fixture
def ins(genbank_record):
    return extract_gene(genbank_record, "INS")


def test_reports_genbank_coordinates_not_biopython_offsets(ins):
    """The flat file says 4986..6416; Biopython holds that as 4985..6416."""
    assert ins["gene"] == "INS"
    assert ins["location"] == {"start": 4986, "end": 6416, "strand": 1}
    assert len(ins["sequence"]) == 1431


def test_exons_are_ordered_by_number_not_by_file_position(ins):
    """Exon 3 trails the mat_peptides in the flat file."""
    assert [exon["number"] for exon in ins["exons"]] == [1, 2, 3]
    assert [(exon["start"], exon["end"]) for exon in ins["exons"]] == [
        (4986, 5027),
        (5207, 5410),
        (6198, 6416),
    ]


def test_compound_locations_keep_their_introns(ins):
    """join(5224..5410,6198..6343) must not collapse to 5224..6343."""
    assert ins["protein"]["segments"] == [
        {"start": 5224, "end": 5410},
        {"start": 6198, "end": 6343},
    ]
    assert ins["transcript"]["segments"] == [
        {"start": 4986, "end": 5027},
        {"start": 5207, "end": 5410},
        {"start": 6198, "end": 6416},
    ]


def test_protein_is_the_preproinsulin_translation(ins):
    assert ins["protein"]["product"] == "insulin preproprotein"
    assert ins["protein"]["translation"] == PREPROINSULIN


def test_mature_peptides_are_the_three_insulin_chains(ins):
    assert [
        (peptide["product"], peptide["translation"]) for peptide in ins["peptides"]
    ] == [
        ("insulin B chain", "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"),
        ("C-peptide", "EAEDLQVGQVELGGGPGAGSLQPLALEGSLQ"),
        ("insulin A chain", "GIVEQCCTSICSLYQLENYCN"),
    ]


def test_signal_peptide_and_proinsulin(ins):
    """The two peptides GenBank states with no /translation of their own."""
    assert len(ins["signal_peptide"]["translation"]) == 24
    assert ins["proprotein"]["product"] == "proinsulin"
    # B chain + RR + C-peptide + KR + A chain.
    assert len(ins["proprotein"]["translation"]) == 86
    assert ins["proprotein"]["segments"] == [
        {"start": 5296, "end": 5410},
        {"start": 6198, "end": 6340},
    ]


def test_readthrough_features_are_excluded(genbank_record, ins):
    """INS-IGF2 starts at the very same base as INS and shares its signal peptide."""
    text = repr(ins)
    assert "INS-IGF2" not in text
    assert "insulin, isoform 2 precursor" not in text
    assert "tyrosine" not in text

    readthrough = extract_gene(genbank_record, "INS-IGF2")
    assert readthrough["protein"]["product"] == "insulin, isoform 2 precursor"
    assert readthrough["location"]["end"] == 8416


def test_gene_match_is_exact_not_a_prefix(genbank_record):
    """A prefix test would fold INS-IGF2's features into the INS answer."""
    with pytest.raises(GeneNotFound):
        extract_gene(genbank_record, "IN")
    with pytest.raises(GeneNotFound):
        extract_gene(genbank_record, "ins")


def test_absent_gene_raises(genbank_record):
    with pytest.raises(GeneNotFound):
        extract_gene(genbank_record, "BRCA1")
