"""Reading the twenty entries residue by residue, on their own quirks.

Offline: the entries are in `structures/`, in the repository.
"""

from __future__ import annotations

import sys
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.structure.pdb import (  # noqa: E402
    STRUCTURES,
    ca_atoms,
    chain_residues,
    inserted,
    secondary_structure,
    substituted,
)
from pipeline.targets import BY_SLUG, TARGETS  # noqa: E402


def _residues(slug: str, chain: str | None = None):
    structure = BY_SLUG[slug].structure
    return chain_residues(structure, chain or structure.chains[0].pdb_chain)


def test_a_residue_modelled_without_its_ca_is_listed_without_one():
    # 1HGU's Pro37 is a single nitrogen: in the entry, with no place for a CA.
    by_number = {r.number: r for r in _residues("somatotropin")}
    assert by_number[37].letter == "P" and by_number[37].ca is None
    assert by_number[36].ca is not None


def test_residues_the_experiment_never_located_are_listed_without_a_ca():
    by_number = {r.number: r for r in _residues("somatotropin")}
    assert [n for n in (1, 38, 39, 191) if by_number[n].ca is None] == [1, 38, 39, 191]


def test_the_missing_list_starts_at_its_own_header_not_the_prose_above_it():
    # 6LMK's chain is E, and the prose explaining the columns has an E in the
    # chain column; read as a residue it would not even parse.
    residues = _residues("glucagon")
    assert [r.number for r in residues] == list(range(1, 30))
    assert "".join(r.letter for r in residues) == "HSQGTFTSDYSKYLDSRRAQDFVQWLMNT"


def test_expression_tags_are_not_the_protein():
    # 4KML's chain A starts with a 33-residue tag, numbered -9 to 23.
    assert {("A", n) for n in range(-9, 24)} <= inserted(STRUCTURES / "4KML.pdb")
    residues = _residues("prion")
    assert residues[0].number == 24 and residues[-1].number == 231


def test_a_modified_residue_reads_as_the_one_it_modifies():
    first = _residues("relaxin", "A")[0]
    assert (first.number, first.letter) == (-3, "Q")      # PCA, a pyroglutamate
    assert first.ca is not None
    assert _residues("amylase")[0].letter == "Q"


def test_a_ligand_is_never_a_residue():
    # 1SMD holds a calcium ion in chain A, whose one atom is also named CA.
    numbers = [r.number for r in _residues("amylase")]
    assert numbers == list(range(1, 497))


def test_the_reading_agrees_with_the_bakes_own_ca_atoms():
    # Every CA `ca_atoms` exports, and nothing else but a modified residue.
    for target in TARGETS:
        structure = target.structure
        located = sum(1 for c in structure.chains
                      for r in chain_residues(structure, c.pdb_chain) if r.ca is not None)
        modified = {"relaxin": 1, "amylase": 1}.get(target.slug, 0)
        assert located == len(ca_atoms(structure)) + modified, target.slug


def test_helices_and_strands_are_the_entrys_own():
    sod1 = secondary_structure(STRUCTURES / "2C9V.pdb")
    assert sum(1 for (c, _), s in sod1.items() if c == "A" and s == "strand") > 50
    # SOD1's only helices are 3-10 helices, and they count.
    assert any(s == "helix" for (c, _), s in sod1.items() if c == "A")
    insulin = secondary_structure(STRUCTURES / "3I40.pdb")
    assert set(insulin.values()) == {"helix"}


def test_a_substitution_is_declared_and_a_tag_is_not_one():
    assert ("A", 100) in substituted(STRUCTURES / "1AX8.pdb")     # leptin's W121E
    assert ("B", 30) in substituted(STRUCTURES / "3I40.pdb")      # insulin's B30
    assert not substituted(STRUCTURES / "4KML.pdb") & inserted(STRUCTURES / "4KML.pdb")
