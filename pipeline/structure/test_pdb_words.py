"""The fold page's words for an experimental entry: templated from facts.

Every claim is a number the entry gives; nothing knows which protein it is.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.structure import pdb_words  # noqa: E402
from pipeline.structure.normalize import Atom, Bridge, Drawn, Examined, Residue  # noqa: E402
from pipeline.structure.pdb_policy import Candidate, Inputs, Pick, Product  # noqa: E402
from pipeline.structure.pdb_search import Domain, Entity, Entry  # noqa: E402

CA = (Atom("CA", "C", 0.0, 0.0, 0.0, 1.0, 20.0, False),)


def chain(node, product, first, last, unplaced=(), entity="1AAA_1", asym="A") -> Drawn:
    return Drawn(node, product, entity, asym, asym, (first, last), tuple(
        Residue(n, n, str(n), "ALA", "A", None, () if n in unplaced else CA)
        for n in range(first, last + 1)))


def examined(chains, bridges=(), methods=("X-RAY DIFFRACTION",)) -> Examined:
    return Examined(entry="1AAA", sha256="0" * 64, model="1", methods=tuple(methods),
                    chains=tuple(chains), bridges=tuple(bridges), secondary={},
                    secondary_method=None)


def candidate(scope="product", domains=(), resolution=1.7, picks=1) -> Candidate:
    made = tuple(Pick(i, f"1AAA_{i + 1}", "AB"[i], "AB"[i], (1, 1, 10), 1, 10, 10, ())
                 for i in range(picks))
    return Candidate(entry="1AAA", scope=scope, picks=made, domains=tuple(domains),
                     method="X-RAY DIFFRACTION", resolution=resolution, useful=10, coverage=1.0,
                     share=1.0, precursor=False, mismatches=0, allowed=0, tier=0, support=1.0)


def protein(products=(("Transthyretin", 21, 147),), length=147, domains=()) -> Inputs:
    return Inputs(slug="x", gene="X", accession="P00001", display="Protein X",
                  translation="A" * length,
                  products=tuple(Product(i, label, a, b) for i, (label, a, b) in enumerate(products)),
                  domains=tuple(domains))


def entry(polymers=1, ligands=0) -> Entry:
    return Entry(entry="1AAA", methods=("X-RAY DIFFRACTION",), resolution=1.7, rfree=0.2, models=1,
                 nmr_conformer=None, clashscore=1.0, ramachandran_outliers=0.0,
                 rotamer_outliers=0.0, released=None, revised=None, revision=None, status="REL",
                 polymer_entities=polymers, nonpolymer_entities=ligands, oligomeric_count=None)


def entity(copies=1, alone=True, symmetry=None) -> Entity:
    return Entity(id="1AAA_1", entry=entry(), entity="1", description=None, polymer="Protein",
                  sequence="A" * 10, mutations=0, sources=(9606,), hosts=(), alignments=(),
                  chains=(), copies=copies, alone=alone, symmetry=symmetry)


def words(protein_, candidate_, examined_, entry_=None, entity_=None):
    chrome, chains = pdb_words.describe(protein_, candidate_, examined_, entry_ or entry(), entity_)
    assert tuple(chrome) == ("pdb", "modelled", "label", "count", "unit", "sentence", "semantics")
    assert chrome["pdb"] == "1AAA" and chrome["unit"] == "residues"
    assert len(chrome["sentence"]) <= pdb_words.LINE
    assert chrome["semantics"].endswith("Drag to turn it.")
    return chrome, chains


def test_a_whole_chain_alone_is_the_fold():
    chrome, chains = words(protein(), candidate(), examined([chain("chainA", 0, 21, 147)]),
                           entity_=entity())
    assert chrome["label"] == "the fold" and chrome["count"] == 127 and chrome["modelled"] is None
    assert chrome["sentence"] == "X-ray, 1.70 Å. Every residue placed."
    assert chrome["semantics"] == (
        "The Protein X structure in PDB 1AAA, by X-ray crystallography, 1.70 Å: its mature "
        "chain, residues 21 to 147. Every residue is placed. Drag to turn it.")
    assert chains == [{"node": "chainA", "tint": "mature1"}]


def test_the_count_is_the_span_and_the_sentence_says_how_much_of_it_is_placed():
    unplaced = set(range(21, 30)) | {146, 147}
    chrome, _ = words(protein(), candidate(), examined([chain("chainA", 0, 21, 147, unplaced)]),
                      entity_=entity())
    assert chrome["count"] == 127
    assert chrome["sentence"] == "X-ray, 1.70 Å. 116 of these 127 residues placed."
    assert "116 of its 127 residues are placed" in chrome["semantics"]
    # Never a claim about why: an unplaced residue is not thereby disordered.
    assert "disorder" not in (chrome["sentence"] + chrome["semantics"]).lower()


def test_a_chain_of_a_homo_oligomer_is_the_subunit_and_a_partnered_one_is_not():
    drawn = examined([chain("chainA", 0, 21, 147)])
    chrome, _ = words(protein(), candidate(), drawn,
                      entity_=entity(copies=4, alone=True, symmetry="D2"))
    assert chrome["label"] == "the subunit"
    assert chrome["sentence"] == "X-ray, 1.70 Å. One chain of four; every residue placed."
    assert "one of its four identical chains" in chrome["semantics"]
    # Beside a receptor, or the antibody fragment it was solved with, it is
    # one chain of nothing in particular: the others are only said not to be drawn.
    chrome, _ = words(protein(), candidate(), drawn, entry(polymers=3, ligands=2),
                      entity(copies=2, alone=False, symmetry="C2"))
    assert chrome["label"] == "the fold" and "One chain of" not in chrome["sentence"]
    assert "The entry's other chains and its ligands are not drawn." in chrome["semantics"]
    # Two molecules that share a crystal, with no symmetry between them, are
    # not a dimer (trypsin's 1TRN).
    chrome, _ = words(protein(), candidate(), drawn, entity_=entity(copies=2, symmetry="C1"))
    assert chrome["label"] == "the fold" and "One chain of" not in chrome["sentence"]


def test_what_the_entry_holds_and_is_not_drawn_is_said_to_a_screen_reader():
    drawn = examined([chain("chainA", 0, 21, 147)])
    assert "Its ligands are not drawn." in words(protein(), candidate(), drawn,
                                                 entry(ligands=1), entity())[0]["semantics"]
    assert "The entry's other chains are not drawn." in words(
        protein(), candidate(), drawn, entry(polymers=2), entity(alone=False))[0]["semantics"]
    assert "not drawn" not in words(protein(), candidate(), drawn, entry(), entity())[0]["semantics"]


TWO = (("Insulin B chain", 25, 54), ("Insulin A chain", 90, 110))


def test_two_chains_are_the_mature_chains_and_their_bridges_are_counted_by_kind():
    drawn = examined([chain("chainA", 0, 25, 54), chain("chainB", 1, 90, 110, entity="1AAA_2")],
                     [Bridge(31, "chainA", 96, "chainB", 2.03),
                      Bridge(43, "chainA", 109, "chainB", 2.04),
                      Bridge(95, "chainB", 100, "chainB", 2.04)])
    chrome, chains = words(protein(TWO, 110), candidate(picks=2), drawn)
    assert chrome["label"] == "the mature chains" and chrome["count"] == 51
    assert chrome["modelled"] is None
    assert chrome["sentence"] == "X-ray, 1.70 Å. Two chains, joined by two disulfide bridges."
    assert ("two chains, Insulin B chain (residues 25 to 54) and Insulin A chain (residues 90 to "
            "110)") in chrome["semantics"]
    assert "Two disulfide bridges join the chains; one more folds a chain." in chrome["semantics"]
    assert chains == [{"node": "chainA", "tint": "mature1"}, {"node": "chainB", "tint": "mature2"},
                      {"node": "bonds", "tint": "cysteine"}]


def test_a_chains_tint_is_its_products_place_so_the_legend_can_name_it():
    three = TWO + (("Third", 111, 120),)
    drawn = examined([chain("chainA", 0, 25, 54), chain("chainB", 2, 111, 120, entity="1AAA_2")])
    chrome, chains = words(protein(three, 120), candidate(picks=2), drawn)
    # The first and the third kept regions: mature1 and mature3, as the legend lists them.
    assert [c["tint"] for c in chains] == ["mature1", "mature3"]
    assert chrome["label"] == "two of its chains"
    assert chrome["sentence"] == "X-ray, 1.70 Å. Two chains, as the entry holds them."
    assert pdb_words.tint_of(5, drawn=1) == "mature1"
    with pytest.raises(ValueError):
        pdb_words.tint_of(3, drawn=2)


def test_one_of_a_precursors_several_products_is_named():
    six = (("Glicentin-related polypeptide", 21, 50), ("Glucagon", 53, 81),
           ("Glucagon-like peptide 1", 92, 128))
    drawn = examined([chain("chainA", 1, 53, 81, unplaced={53})])
    chrome, chains = words(protein(six, 180), candidate(resolution=1.1), drawn, entity_=entity())
    assert chrome["label"] == "one of its products" and chrome["modelled"] is None
    assert chrome["sentence"] == ("X-ray, 1.10 Å. Glucagon, one of three products; "
                                  "28 of 29 placed.")
    assert "Glucagon, one of its three products, residues 53 to 81" in chrome["semantics"]
    assert chains == [{"node": "chainA", "tint": "mature1"}]
    # A long name gives way before the line does.
    long = examined([chain("chainA", 2, 92, 128, unplaced={92, 93})],
                    methods=("ELECTRON MICROSCOPY",))
    chrome, _ = words(protein(six, 180), candidate(resolution=2.9), long, entity_=entity())
    assert chrome["sentence"] == ("cryo-EM, 2.9 Å. Glucagon-like peptide 1, one of three "
                                  "products; 35 of 37 placed.")
    assert len(chrome["sentence"]) == 80
    # And a longer one still drops the count rather than the name.
    longer = (six[0], six[1], ("Glucagon-like peptide 1(7-36)", 92, 128))
    chrome, _ = words(protein(longer, 180), candidate(resolution=2.9), long, entity_=entity())
    assert chrome["sentence"] == ("cryo-EM, 2.9 Å. Glucagon-like peptide 1(7-36), one of "
                                  "three products.")


def test_a_domain_is_named_and_its_span_is_said():
    kinase = Domain("Domain", "Protein kinase", 100, 288)
    drawn = examined([chain("chainA", 0, 94, 293, unplaced=set(range(286, 294)))])
    chrome, _ = words(protein((("p53", 1, 393),), 393, [kinase]),
                      candidate("domain", ["Protein kinase"], resolution=1.32), drawn,
                      entry(polymers=2), entity(alone=False))
    assert chrome["label"] == "the protein kinase domain"
    assert chrome["modelled"] == [94, 293] and chrome["count"] == 200
    assert chrome["sentence"] == ("X-ray, 1.32 Å. Only this part is drawn; 192 of these 200 "
                                  "residues placed.")
    assert "the protein kinase domain, residues 94 to 293" in chrome["semantics"]


@pytest.mark.parametrize("name, kind, said", [
    ("VWFA", "Domain", "the VWFA domain"),
    ("Protein kinase", "Domain", "the protein kinase domain"),
    ("Ig-like C1-type", "Domain", "the Ig-like C1-type domain"),
    ("P53 DNA-binding domain", "Domain", "the P53 DNA-binding domain"),
    ("Prion/Doppel alpha-helical domain", "Domain", "the Prion/Doppel alpha-helical domain"),
    ("RCC1 3", "Repeat", "the RCC1 3 repeat"),
    ("Spectrin 12", "Repeat", "the spectrin 12 repeat"),
    ("C2H2-type 1", "Zinc finger", "the C2H2-type 1 zinc finger"),
])
def test_a_domains_name_is_written_as_uniprot_writes_it_but_for_a_plain_first_word(name, kind, said):
    assert pdb_words.domain_phrase(name, kind) == said


def test_two_domains_are_named_and_more_are_counted():
    kinds = {"Calponin-homology (CH) 1": "Domain", "Calponin-homology (CH) 2": "Domain",
             "WW": "Domain", "EF hand": "Domain", "RCC1 1": "Repeat", "RCC1 2": "Repeat",
             "RCC1 3": "Repeat"}
    assert pdb_words.domains_phrase(["WW", "EF hand"], kinds) == "the WW and EF hand domains"
    assert pdb_words.domains_phrase(["RCC1 1", "RCC1 2"], kinds) == "the RCC1 1 and RCC1 2 repeats"
    assert pdb_words.domains_phrase(["RCC1 1", "RCC1 2", "RCC1 3"], kinds) == "three of its repeats"
    # Repeats and domains mixed (mTOR's HEAT and TPR repeats, its FRB and FATC
    # domains) are not worth a count.
    assert pdb_words.domains_phrase(["WW", "RCC1 1", "RCC1 2"], kinds) == "part of the protein"
    assert pdb_words.domains_phrase(["WW", "RCC1 1"], kinds) == (
        "the WW domain and the RCC1 1 repeat")
    # Two long names are counted rather than run off the header.
    assert pdb_words.domains_phrase(["Calponin-homology (CH) 1", "Calponin-homology (CH) 2"],
                                    kinds) == "two of its domains"


@pytest.mark.parametrize("methods, resolution, short, long", [
    (("X-RAY DIFFRACTION",), 1.7, "X-ray, 1.70 Å", "X-ray crystallography, 1.70 Å"),
    (("ELECTRON MICROSCOPY",), 3.04, "cryo-EM, 3.0 Å", "cryo-electron microscopy, 3.0 Å"),
    (("SOLUTION NMR",), None, "NMR", "solution NMR"),
    (("SOLID-STATE NMR",), None, "solid-state NMR", "solid-state NMR"),
    (("NEUTRON DIFFRACTION", "X-RAY DIFFRACTION"), 1.1, "neutron and X-ray, 1.10 Å",
     "neutron crystallography and X-ray crystallography, 1.10 Å"),
    (("ELECTRON CRYSTALLOGRAPHY",), 2.5, "electron crystallography, 2.50 Å",
     "electron crystallography, 2.50 Å"),
])
def test_a_method_is_named_with_its_resolution_where_it_has_one(methods, resolution, short, long):
    assert pdb_words.method(methods, resolution) == short
    assert pdb_words.method(methods, resolution, long=True) == long


def test_the_sentence_keeps_to_two_lines_at_the_longest_numbers_there_are():
    # R4.5: about 80 characters on a phone.
    methods = (("NEUTRON DIFFRACTION", "X-RAY DIFFRACTION"), ("ELECTRON CRYSTALLOGRAPHY",),
               ("ELECTRON MICROSCOPY",), ("SOLID-STATE NMR",))
    long_name = "Gamma-aminobutyric acid receptor-associated protein-like 2"
    many = (("first", 1, 9_000), (long_name, 9_001, 9_999), ("third", 10_000, 10_500))
    for used in methods:
        whole = examined([chain("chainA", 0, 1, 9_000, unplaced={1})], methods=used)
        part = examined([chain("chainA", 1, 9_001, 9_999, unplaced={9_001})], methods=used)
        for protein_, candidate_, drawn, entity_ in (
                (protein((("x", 1, 9_000),), 9_000), candidate(resolution=2.55), whole, entity()),
                (protein((("x", 1, 9_000),), 9_000), candidate(resolution=2.55), whole,
                 entity(copies=12, symmetry="D6")),
                (protein(many, 10_500), candidate(resolution=2.55), part, entity()),
                (protein((("x", 1, 9_000),), 9_000,
                         [Domain("Domain", "D", 1, 9_000)]),
                 candidate("domain", ["D"], resolution=2.55), whole, entity())):
            chrome, _ = words(protein_, candidate_, drawn, entity_=entity_)
            assert len(chrome["sentence"]) <= 80, chrome["sentence"]
