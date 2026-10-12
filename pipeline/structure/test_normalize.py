"""An entry held to its own file: what is read, chosen and refused.

Three real entries, PDBe's SIFTS-enriched files cut down under `fixtures/pdb/`:
relaxin's two chains (6RLX, with its bridges and a pyroglutamate), an NMR
ensemble (2K6T) and thrombin in chymotrypsin numbering (1PPB, insertion
codes). Each refusal is held on a small file written here, where one thing is
wrong at a time.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.structure import mmcif, normalize, pdb_policy, pdb_search  # noqa: E402
from pipeline.structure.normalize import Unsuitable  # noqa: E402
from pipeline.structure.pdb import ONE_LETTER  # noqa: E402
from pipeline.structure.pdb_policy import Candidate, Inputs, Pick, Product  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "pdb"
THREE = {one: three for three, one in ONE_LETTER.items()}


# ------------------------------------------------------------ real entries


def recorded(name: str) -> list:
    whole = json.loads((FIXTURES / name).read_bytes())
    return [pdb_search._entity(e) for e in whole["data"]["polymer_entities"]]


def protein_of(accession: str, products, removed) -> Inputs:
    body = json.loads((FIXTURES / f"uniprot_{accession}.json").read_bytes())
    made = tuple(Product(i, label, a, b) for i, (label, a, b) in enumerate(products))
    pairs = tuple((f["location"]["start"]["value"], f["location"]["end"]["value"])
                  for f in body["features"] if f["type"] == "Disulfide bond"
                  and f["location"]["start"]["value"] != f["location"]["end"]["value"])
    return Inputs(slug=accession.lower(), gene=accession, accession=accession, display=accession,
                  translation=body["sequence"]["value"], products=made, removed=removed,
                  domains=pdb_policy.domains_within(pdb_search.domains_of(body), made),
                  disulfides=pairs)


RELAXIN = protein_of("P04090", [("B chain", 25, 53), ("A chain", 162, 185)],
                     ((1, 24, "Signal"), (54, 161, "Connecting peptide")))
INSL3 = protein_of("P51460", [("Insulin-like 3 B chain", 21, 55), ("Insulin-like 3 A chain", 106, 131)],
                   ((1, 20, "Signal"), (58, 104, "C peptide like")))
PROTHROMBIN = protein_of("P00734", [("Prothrombin", 44, 622)],
                         ((1, 24, "Signal"), (25, 43, "Propeptide")))


def candidate_of(protein: Inputs, fixture: str, entry: str, chains: int = None,
                 scope: str = None) -> Candidate:
    ranking = pdb_policy.rank(protein, recorded(fixture))
    return next(c for c in ranking.candidates if c.entry == entry
                and (chains is None or len(c.picks) == chains)
                and (scope is None or c.scope == scope))


def test_relaxins_two_chains_are_read_as_the_model_draws_them():
    text = (FIXTURES / "6rlx_cut.cif").read_bytes()
    examined = normalize.examine(RELAXIN, candidate_of(RELAXIN, "graphql_relaxin.json", "6RLX", 2),
                                 text)
    assert examined.entry == "6RLX" and examined.model == "1"
    assert examined.methods == ("X-RAY DIFFRACTION",)
    assert examined.sha256 == hashlib.sha256(text).hexdigest()
    b_chain, a_chain = examined.chains
    # Nodes in the precursor's order: the B chain comes first in the precursor.
    assert (b_chain.node, b_chain.product, b_chain.asym, b_chain.auth) == ("chainA", 0, "B", "B")
    assert (a_chain.node, a_chain.product, a_chain.asym, a_chain.auth) == ("chainB", 1, "A", "A")
    # Each mature chain is spanned whole, as the fold page counts it. The
    # construct's residues are those the entry held: the B chain from Ser26
    # (6RLX never had Asp25), its last never placed.
    assert (b_chain.first, b_chain.last, a_chain.first, a_chain.last) == (25, 53, 162, 185)
    assert [r.number for r in b_chain.residues] == list(range(26, 54))
    assert [r.number for r in b_chain.residues if not r.placed] == [53]
    assert all(r.placed for r in a_chain.residues)
    assert (examined.placed, examined.spanned) == (51, 53)
    # Every residue is the record's, letter for letter, and numbered without a gap.
    for chain in examined.chains:
        numbers = [r.number for r in chain.residues]
        assert numbers == list(range(numbers[0], numbers[-1] + 1))
        assert all(r.letter == RELAXIN.translation[r.number - 1] for r in chain.residues)
    assert examined.audit["mismatches"] == []


def test_relaxins_bridges_are_the_entrys_where_its_sulfurs_meet():
    examined = normalize.examine(RELAXIN, candidate_of(RELAXIN, "graphql_relaxin.json", "6RLX", 2),
                                 (FIXTURES / "6rlx_cut.cif").read_bytes())
    assert [(b.a, b.a_node, b.b, b.b_node) for b in examined.bridges] == [
        (35, "chainA", 172, "chainB"), (47, "chainA", 185, "chainB"),
        (171, "chainB", 176, "chainB")]
    assert all(1.9 < b.distance < 2.2 for b in examined.bridges)
    assert examined.audit["bridges_dropped"] == []
    # UniProt's three are the entry's three (R5.6).
    assert sorted(RELAXIN.disulfides) == sorted((b.a, b.b) for b in examined.bridges)
    assert examined.audit["uniprot_bridges_not_drawn"] == []


def test_a_modified_residue_reads_as_the_residue_the_gene_makes():
    examined = normalize.examine(RELAXIN, candidate_of(RELAXIN, "graphql_relaxin.json", "6RLX", 2),
                                 (FIXTURES / "6rlx_cut.cif").read_bytes())
    first = examined.chains[1].residues[0]
    # The A chain starts on a pyroglutamate, a cyclised glutamine.
    assert (first.number, first.name, first.letter, first.parent) == (162, "PCA", "Q", "GLN")
    assert first.placed and any(a.hetero for a in first.atoms)
    assert examined.audit["modified"] == [{"node": "chainB", "n": 162, "residue": "PCA",
                                           "parent": "GLN"}]


def test_the_entrys_own_helices_and_strands_are_read_in_the_precursors_numbering():
    examined = normalize.examine(RELAXIN, candidate_of(RELAXIN, "graphql_relaxin.json", "6RLX", 2),
                                 (FIXTURES / "6rlx_cut.cif").read_bytes())
    assert examined.secondary_method == ("DSSP in the entry's file from PDBe (struct_conf, "
                                         "struct_sheet_range)")
    for node, chain in zip(("chainA", "chainB"), examined.chains):
        spans = examined.secondary[node]
        assert spans and all(kind in ("helix", "strand") for kind, _, _ in spans)
        assert all(chain.first <= first <= last <= chain.last for _, first, last in spans)
    said = normalize.summary(examined)
    assert said["helix"] >= 3 and said["chains"][0]["placed"] == 27
    assert said["chains"][0]["listed"] == 29 and said["chains"][0]["construct"] == [26, 53]
    assert "atoms" not in json.dumps(said)


def test_an_nmr_entry_is_drawn_from_its_representative_conformer_and_never_an_average():
    text = (FIXTURES / "2k6t_cut.cif").read_bytes()
    candidate = candidate_of(INSL3, "graphql_insl3.json", "2K6T", 2)
    first = normalize.examine(INSL3, candidate, text)
    assert first.model == "1" and first.methods == ("SOLUTION NMR",)
    second = normalize.examine(INSL3, candidate, text.replace(
        b"_pdbx_nmr_representative.conformer_id 1", b"_pdbx_nmr_representative.conformer_id 2"))
    assert second.model == "2"
    one, two = (e.chains[0].residues[5].atom("CA") for e in (first, second))
    assert (one.x, one.y, one.z) != (two.x, two.y, two.z)
    # The coordinates are a conformer's own, to the digit the file gives.
    block = mmcif.parse(text.decode())
    in_file = {(float(x), float(y), float(z)) for x, y, z, model in block["atom_site"].rows(
        "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num") if model == "2"}
    assert (two.x, two.y, two.z) in in_file
    # A conformer the file does not hold is not drawn: the first is.
    odd = normalize.examine(INSL3, candidate, text.replace(
        b"_pdbx_nmr_representative.conformer_id 1", b"_pdbx_nmr_representative.conformer_id 99"))
    assert odd.model == "1"
    # Hydrogens, which an NMR entry is full of, are dropped.
    assert all(a.element not in ("H", "D") for c in first.chains for r in c.residues
               for a in r.atoms)


def test_insertion_codes_are_kept_and_the_precursors_numbering_runs_through_them():
    candidate = candidate_of(PROTHROMBIN, "graphql_thrombin.json", "1PPB", scope="domain")
    examined = normalize.examine(PROTHROMBIN, candidate, (FIXTURES / "1ppb_cut.cif").read_bytes())
    (heavy,) = examined.chains
    numbers = [r.number for r in heavy.residues]
    assert numbers == list(range(numbers[0], numbers[-1] + 1))
    # Thrombin is numbered as chymotrypsinogen, with lettered insertions (60A to 60I).
    lettered = [r.auth for r in heavy.residues if r.auth[-1].isalpha()]
    assert {"60A", "60B", "60I"} <= set(lettered) and len(lettered) > 20
    assert all(r.letter == PROTHROMBIN.translation[r.number - 1] for r in heavy.residues)


# ------------------------------------------------------------ files written here


SEQUENCE = "MKTAYIAKQRCISFVKSHFSCQLEERLGLI"       # 30 residues; cysteines at 11 and 21


def protein(translation=SEQUENCE, products=((1, 30),), disulfides=()) -> Inputs:
    return Inputs(slug="x", gene="X", accession="P00001", display="X", translation=translation,
                  products=tuple(Product(i, f"chain {i + 1}", a, b)
                                 for i, (a, b) in enumerate(products)),
                  disulfides=tuple(disulfides))


def candidate(picks=None, scope="product", entry="1AAA") -> Candidate:
    picks = picks or [dict()]
    made = tuple(Pick(product=p.get("product", i), entity=f"{entry}_{p.get('entity', i + 1)}",
                      asym=p.get("asym", "AB"[i]), auth=p.get("asym", "AB"[i]),
                      segment=p.get("segment", (1, 1, 30)), first=p.get("first", 1),
                      last=p.get("last", 30), placed=30, mismatches=(),
                      mapped_by=p.get("mapped_by", "SIFTS")) for i, p in enumerate(picks))
    return Candidate(entry=entry, scope=scope, picks=made, domains=(), method="X-RAY DIFFRACTION",
                     resolution=2.0, useful=30, coverage=1.0, share=1.0, precursor=False,
                     mismatches=0, allowed=0, tier=0, support=1.0)


class File:
    """A small mmCIF file: chains of residues laid out 3.8 A apart, a chain a row."""

    def __init__(self, entry="1AAA", method="X-RAY DIFFRACTION"):
        self.entry, self.method = entry, method
        self.atoms, self.scheme, self.sifts, self.poly, self.conn = [], [], [], [], []
        self.extra = ""

    def chain(self, asym="A", entity="1", sequence=SEQUENCE, unp_start=1, accession="P00001",
              row=0.0, skip=(), names=None, model="1", shift=None, alt=None, sifts=True,
              auth_start=1, declare=True):
        """Lay a chain down. `skip` are label numbers the entry never placed;
        `names` renames residues ({label: "MSE"}); `shift` moves residues
        ({label: dx}); `alt` gives a residue alternate CA positions
        ({label: [(alt id, occupancy, dy), ...]})."""
        names, shift, alt = names or {}, shift or {}, alt or {}
        if declare:
            self.poly.append((entity, sequence))
        for i, letter in enumerate(sequence):
            label = i + 1
            name = names.get(label, THREE[letter])
            if declare:
                self.scheme.append((asym, entity, label, name, auth_start + i, None))
                if sifts:
                    self.sifts.append((entity, asym, label, accession, unp_start + i, letter))
            if label in skip:
                continue
            x = 3.8 * i + shift.get(label, 0.0)
            group = "HETATM" if name not in ONE_LETTER else "ATOM"
            for atom, element, dx, dy in (("N", "N", -1.2, 0.6), ("CA", "C", 0.0, 0.0),
                                          ("C", "C", 1.2, 0.6), ("O", "O", 1.4, 1.8),
                                          ("CB", "C", 0.0, -1.5), ("H", "H", -1.5, 1.4)):
                versions = alt.get(label, [(None, 1.0, 0.0)]) if atom == "CA" else [(None, 1.0, 0.0)]
                for alt_id, occupancy, lift in versions:
                    self.atoms.append((group, element, atom, alt_id, name, asym, entity, label,
                                       x + dx, row + dy + lift, 0.0, occupancy, model,
                                       auth_start + i))
            if letter == "C":
                self.atoms.append((group, "S", "SG", None, name, asym, entity, label,
                                   x, row - 3.3, 0.0, 1.0, model, auth_start + i))
        return self

    def sulfur(self, asym, label, x, y, model="1"):
        """Move a cysteine's sulfur."""
        self.atoms = [a if not (a[2] == "SG" and a[5] == asym and a[7] == label and a[12] == model)
                      else a[:8] + (x, y) + a[10:] for a in self.atoms]
        return self

    def bridge(self, asym1, label1, asym2, label2, sym1="1_555", sym2="1_555"):
        self.conn.append((asym1, label1, sym1, asym2, label2, sym2))
        return self

    def text(self) -> bytes:
        def value(v):
            return "?" if v is None else str(v)

        lines = [f"data_{self.entry}", f"_entry.id {self.entry}", f"_exptl.method '{self.method}'",
                 "loop_", "_entity_poly.entity_id", "_entity_poly.pdbx_seq_one_letter_code_can"]
        lines += [f"{entity} {sequence}" for entity, sequence in self.poly]
        lines += ["loop_"] + [f"_pdbx_poly_seq_scheme.{item}" for item in (
            "asym_id", "entity_id", "seq_id", "mon_id", "pdb_seq_num", "pdb_ins_code")]
        lines += [" ".join(value(v) for v in row) for row in self.scheme]
        if self.sifts:
            lines += ["loop_"] + [f"_pdbx_sifts_xref_db.{item}" for item in (
                "entity_id", "asym_id", "seq_id", "unp_acc", "unp_num", "unp_res")]
            lines += [" ".join(value(v) for v in row) for row in self.sifts]
        lines += ["loop_"] + [f"_atom_site.{item}" for item in (
            "group_PDB", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
            "label_asym_id", "label_entity_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z",
            "occupancy", "pdbx_PDB_model_num", "auth_seq_id")]
        lines += [" ".join(value(v) for v in row) for row in self.atoms]
        if self.conn:
            lines += ["loop_"] + [f"_struct_conn.{item}" for item in (
                "conn_type_id", "ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_symmetry",
                "ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_symmetry")]
            lines += ["disulf " + " ".join(value(v) for v in row) for row in self.conn]
        return ("\n".join(lines) + "\n" + self.extra).encode()


def test_a_plain_chain_is_read_whole():
    examined = normalize.examine(protein(), candidate(), File().chain().text())
    (chain,) = examined.chains
    assert (chain.node, chain.first, chain.last, examined.placed) == ("chainA", 1, 30, 30)
    assert [r.auth for r in chain.residues[:3]] == ["1", "2", "3"]
    assert examined.secondary_method is None and examined.secondary == {"chainA": ()}
    assert all(a.element != "H" for r in chain.residues for a in r.atoms)


def test_residues_the_entry_never_placed_are_listed_and_not_placed():
    examined = normalize.examine(protein(), candidate(), File().chain(skip={1, 2, 15, 30}).text())
    (chain,) = examined.chains
    assert len(chain.residues) == 30
    assert [r.number for r in chain.residues if not r.placed] == [1, 2, 15, 30]
    assert (examined.placed, examined.spanned) == (26, 30)
    # A construct that never held the product's ends leaves the span whole.
    short = File().chain(sequence=SEQUENCE[4:26], unp_start=5).text()
    cut = candidate([dict(segment=(1, 5, 22), first=5, last=26)])
    examined = normalize.examine(protein(), cut, short)
    (chain,) = examined.chains
    assert (chain.first, chain.last, examined.spanned, examined.placed) == (1, 30, 30, 22)
    assert [chain.residues[0].number, chain.residues[-1].number] == [5, 26]


def test_of_alternate_locations_the_likelier_is_kept_whole_and_ties_go_to_the_first():
    likelier = File().chain(alt={7: [("A", 0.4, 0.0), ("B", 0.6, 0.5)]}).text()
    examined = normalize.examine(protein(), candidate(), likelier)
    residue = examined.chains[0].residues[6]
    assert residue.alternate == "B" and residue.atom("CA").y == 0.5
    assert len([a for a in residue.atoms if a.name == "CA"]) == 1
    assert examined.audit["alternates"] == [{"node": "chainA", "n": 7, "kept": "B"}]
    tied = File().chain(alt={7: [("A", 0.5, 0.0), ("B", 0.5, 0.5)]}).text()
    assert normalize.examine(protein(), candidate(), tied).chains[0].residues[6].alternate == "A"


def test_a_selenomethionine_is_the_methionine_the_gene_makes():
    examined = normalize.examine(protein(), candidate(), File().chain(names={1: "MSE"}).text())
    first = examined.chains[0].residues[0]
    assert (first.name, first.letter, first.parent) == ("MSE", "M", "MET")
    assert first.placed and examined.audit["mismatches"] == []


def test_two_residue_types_at_one_position_are_refused():
    made = File().chain()
    made.scheme.append(("A", "1", 5, "SER", 5, None))       # label 5 is TYR, and SER
    with pytest.raises(Unsuitable) as refused:
        normalize.examine(protein(), candidate(), made.text())
    assert refused.value.reason == "microheterogeneity"


def test_the_files_sifts_map_must_be_rcsbs_segment_residue_for_residue():
    with pytest.raises(Unsuitable, match="no SIFTS map") as refused:
        normalize.examine(protein(), candidate(), File().chain(sifts=False).text())
    assert refused.value.reason == "mapping"
    shifted = File().chain(unp_start=2).text()
    with pytest.raises(Unsuitable, match="RCSB's segment says P00001 1"):
        normalize.examine(protein(), candidate(), shifted)
    with pytest.raises(Unsuitable, match="RCSB's segment says"):
        normalize.examine(protein(), candidate(), File().chain(accession="Q99999").text())


def test_a_chain_laid_on_the_protein_by_alignment_needs_no_sifts_and_no_contrary_one():
    aligned = candidate([dict(mapped_by="alignment")])
    examined = normalize.examine(protein(), aligned, File().chain(sifts=False).text())
    assert examined.placed == 30
    with pytest.raises(Unsuitable, match="SIFTS maps chain A to Q99999"):
        normalize.examine(protein(), aligned, File().chain(accession="Q99999").text())


def test_a_break_inside_the_stretch_is_refused_and_never_bridged():
    apart = File().chain(shift={n: 9.0 for n in range(16, 31)}).text()
    with pytest.raises(Unsuitable) as refused:
        normalize.examine(protein(), candidate(), apart)
    assert refused.value.reason == "chain break" and "CA 15 to CA 16" in refused.value.detail
    # Across residues with no place, the gap may be as long as they could reach.
    reachable = File().chain(skip={15, 16}).text()
    assert normalize.examine(protein(), candidate(), reachable).placed == 28
    too_far = File().chain(skip={15, 16}, shift={n: 9.0 for n in range(17, 31)}).text()
    with pytest.raises(Unsuitable, match="across 3 residue"):
        normalize.examine(protein(), candidate(), too_far)


def test_letters_are_counted_again_from_the_file():
    other = SEQUENCE[:4] + "W" + SEQUENCE[5:]
    with pytest.raises(Unsuitable) as refused:
        normalize.examine(protein(), candidate(), File().chain(sequence=other).text())
    assert refused.value.reason == "fidelity"
    # Where the entry never placed the residue, its letter is not drawn and not counted.
    hidden = normalize.examine(protein(), candidate(), File().chain(sequence=other, skip={5}).text())
    assert hidden.audit["mismatches"] == []


def test_a_domain_is_listed_from_its_first_placed_residue_to_its_last():
    domain = candidate([dict(first=5, last=24)], scope="domain")
    examined = normalize.examine(protein(), domain, File().chain(skip={5, 6, 24}).text())
    (chain,) = examined.chains
    assert (chain.first, chain.last) == (7, 23)


TWO = dict(translation=SEQUENCE + "GSGSGSGSGS" + SEQUENCE[::-1], products=((1, 30), (41, 70)))


def pair(row=6.0, **second):
    made = File().chain("A", "1", SEQUENCE, 1)
    return made.chain("B", "2", SEQUENCE[::-1], 41, row=row, **second)


def both() -> Candidate:
    return candidate([dict(asym="A", segment=(1, 1, 30)),
                      dict(asym="B", segment=(1, 41, 30), first=41, last=70)])


def test_two_chains_are_drawn_together_only_where_they_touch():
    touching = normalize.examine(protein(**TWO), both(), pair(row=6.0).text())
    assert [c.node for c in touching.chains] == ["chainA", "chainB"]
    assert [(c.first, c.last) for c in touching.chains] == [(1, 30), (41, 70)]
    with pytest.raises(Unsuitable) as refused:
        normalize.examine(protein(**TWO), both(), pair(row=40.0).text())
    assert refused.value.reason == "apart"


def test_a_disulfide_joins_chains_that_nothing_else_does():
    # Chain A's Cys11 and chain B's Cys10 (the protein's 50), forty angstroms
    # apart but for one sulfur each, brought together.
    made = pair(row=40.0).sulfur("A", 11, 38.0, 19.0).sulfur("B", 10, 38.0, 21.0)
    made.bridge("A", 11, "B", 10)
    examined = normalize.examine(protein(**TWO), both(), made.text())
    assert [(b.a, b.a_node, b.b, b.b_node, b.distance) for b in examined.bridges] == [
        (11, "chainA", 50, "chainB", 2.0)]


def test_a_bridge_is_drawn_only_where_the_sulfurs_meet_in_this_copy():
    made = File().chain().sulfur("A", 11, 38.0, -3.3).sulfur("A", 21, 38.0, -5.3)
    made.bridge("A", 11, "A", 21)
    examined = normalize.examine(protein(disulfides=((11, 21),)), candidate(), made.text())
    assert [(b.a, b.b, b.distance) for b in examined.bridges] == [(11, 21, 2.0)]
    assert examined.audit["uniprot_bridges_not_drawn"] == []

    apart = File().chain().bridge("A", 11, "A", 21)
    examined = normalize.examine(protein(disulfides=((11, 21),)), candidate(), apart.text())
    assert examined.bridges == ()
    (dropped,) = examined.audit["bridges_dropped"]
    assert dropped["reason"] == "sulfurs apart" and dropped["distance"] > 30
    assert examined.audit["uniprot_bridges_not_drawn"] == [[11, 21]]

    crystal = File().chain().sulfur("A", 21, 38.0, -5.3).bridge("A", 11, "A", 21, sym2="2_655")
    examined = normalize.examine(protein(), candidate(), crystal.text())
    assert examined.bridges == ()
    assert examined.audit["bridges_dropped"][0]["reason"] == "to a copy in the crystal"

    made = File().chain().bridge("A", 11, "A", 21)
    made.atoms = [a for a in made.atoms if not (a[2] == "SG" and a[7] == 21)]
    examined = normalize.examine(protein(), candidate(), made.text())
    assert examined.audit["bridges_dropped"][0]["reason"] == "no sulfur in the entry"

    # The rods run from CA through CB: a cysteine modelled without its CB has none.
    made = File().chain().sulfur("A", 11, 38.0, -3.3).sulfur("A", 21, 38.0, -5.3)
    made.bridge("A", 11, "A", 21)
    made.atoms = [a for a in made.atoms if not (a[2] == "CB" and a[7] == 21)]
    examined = normalize.examine(protein(), candidate(), made.text())
    assert examined.bridges == ()
    assert examined.audit["bridges_dropped"][0]["reason"] == "a cysteine without its CA or CB"


def test_another_entrys_file_and_a_file_with_no_atoms_are_no_verdict():
    with pytest.raises(mmcif.CifError, match="not 1AAA's"):
        normalize.examine(protein(), candidate(), File(entry="2BBB").chain().text())
    with pytest.raises(mmcif.CifError):
        normalize.examine(protein(), candidate(), b"data_1AAA\n_entry.id 1AAA\n")
    with pytest.raises(mmcif.CifError):
        normalize.examine(protein(), candidate(), b"not a file at all")


def test_a_residue_numbered_past_what_pdb_format_holds_is_refused():
    long = "A" * 9_990 + SEQUENCE
    far = candidate([dict(segment=(1, 9_991, 30), first=9_991, last=10_020)])
    made = File().chain(unp_start=9_991)
    with pytest.raises(Unsuitable) as refused:
        normalize.examine(protein(translation=long, products=((1, len(long)),)), far, made.text())
    assert refused.value.reason == "numbering"


# ------------------------------------------------------------ the working copy


def _written(examined, tmp_path, secondary=None) -> Path:
    path = tmp_path / "copy.pdb"
    path.write_text(normalize.working_copy(examined, secondary))
    return path


def test_the_working_copy_reads_as_an_entry_does_to_the_twentys_own_readers(tmp_path):
    from pipeline.structure import pdb
    from pipeline.targets import Chain, Structure

    examined = normalize.examine(RELAXIN, candidate_of(RELAXIN, "graphql_relaxin.json", "6RLX", 2),
                                 (FIXTURES / "6rlx_cut.cif").read_bytes())
    path = _written(examined, tmp_path)
    structure = Structure(pdb="6RLX", chains=(Chain("chainA", "A"), Chain("chainB", "B")),
                          bonds=True)

    # One chain a node, in the record's numbering, a residue the construct
    # held and the entry never placed listed without a CA.
    b_chain = pdb.chain_residues(structure, "A", pdb=path)
    assert [r.number for r in b_chain] == list(range(26, 54))
    assert [r.number for r in b_chain if r.ca is None] == [53]
    assert "".join(r.letter for r in b_chain) == RELAXIN.translation[25:53]
    a_chain = pdb.chain_residues(structure, "B", pdb=path)
    assert [r.number for r in a_chain] == list(range(162, 186))
    # The pyroglutamate is drawn as the glutamine the gene makes.
    assert a_chain[0].letter == "Q" and a_chain[0].ca is not None
    assert "PCA" not in path.read_text() and "HETATM" not in path.read_text()

    # Its bridges are the entry's, as `bake.build_bonds` and the folding bake read them.
    assert pdb.ssbonds(path, structure) == [("A", 35, "B", 172), ("A", 47, "B", 185),
                                            ("B", 171, "B", 176)]
    at = pdb.atoms(path)
    for chain, number in (("A", 35), ("B", 172), ("B", 171), ("B", 176)):
        assert {(chain, number, name) for name in ("CA", "CB", "SG")} <= set(at)
    # Coordinates are the entry's own, to the digit.
    placed = examined.chains[0].residues[0].atom("CA")
    assert tuple(at[("A", 26, "CA")]) == (placed.x, placed.y, placed.z)
    assert len(pdb.ca_atoms(structure, path)) == examined.placed

    # Its helices and strands are the ones examined.
    shape = pdb.secondary_structure(path)
    for letter, chain in zip("AB", examined.chains):
        for kind, first, last in examined.secondary[chain.node]:
            assert all(shape[(letter, n)] == kind for n in range(first, last + 1))
    assert len(shape) == sum(last - first + 1 for spans in examined.secondary.values()
                             for _, first, last in spans)
    assert not pdb.inserted(path) and not pdb.substituted(path)


def test_the_working_copy_takes_another_assignment_of_shape(tmp_path):
    from pipeline.structure import pdb

    examined = normalize.examine(protein(), candidate(), File().chain().text())
    assert pdb.secondary_structure(_written(examined, tmp_path)) == {}
    asked = {"chainA": (("helix", 3, 12), ("strand", 20, 24))}
    shape = pdb.secondary_structure(_written(examined, tmp_path, asked))
    assert {n for (_, n), kind in shape.items() if kind == "helix"} == set(range(3, 13))
    assert {n for (_, n), kind in shape.items() if kind == "strand"} == set(range(20, 25))


def test_the_working_copy_declares_a_residue_that_is_not_the_records(tmp_path):
    from pipeline.structure import pdb

    long = SEQUENCE * 4                                      # 120 residues: two may differ
    other = long[:4] + "W" + long[5:]
    drawn = candidate([dict(segment=(1, 1, 120), first=1, last=120)])
    examined = normalize.examine(protein(translation=long, products=((1, 120),)), drawn,
                                 File().chain(sequence=other).text())
    assert examined.audit["mismatches"] == [{"node": "chainA", "n": 5, "entry": "W", "record": "Y"}]
    path = _written(examined, tmp_path)
    assert pdb.substituted(path) == {("A", 5)} and not pdb.inserted(path)


def test_a_residue_with_no_standard_parent_cannot_be_drawn(tmp_path):
    made = File().chain(names={7: "XYZ"})
    made.poly = [(entity, sequence[:6] + "X" + sequence[7:]) for entity, sequence in made.poly]
    with pytest.raises(Unsuitable):
        examined = normalize.examine(protein(), candidate(), made.text())
        normalize.working_copy(examined)
