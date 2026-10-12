"""The suitability policy: what is eligible, why not, and what ranks first.

Each rule is held on entities made here, so that one thing differs at a time;
the last tests rank real entities, RCSB's answers recorded under
`fixtures/pdb/`. Nothing reaches the network.
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.structure import pdb_policy, pdb_search  # noqa: E402
from pipeline.structure.pdb_policy import Inputs, Product  # noqa: E402
from pipeline.structure.pdb_search import Alignment, Chain, Domain, Entity, Entry  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "pdb"

# A protein of 200 residues: a 20-residue signal peptide, then one chain.
LETTERS = "ACDEFGHIKLMNPQRSTVWY"
PROTEIN = "".join(LETTERS[(7 * i + i // 20) % 20] for i in range(200))


def inputs(products=((21, 200),), domains=(), removed=((1, 20, "Signal peptide"),),
           translation=PROTEIN) -> Inputs:
    return Inputs(slug="x", gene="X", accession="P00001", display="Protein X",
                  translation=translation,
                  products=tuple(Product(i, f"chain {i + 1}", a, b)
                                 for i, (a, b) in enumerate(products)),
                  removed=removed, domains=tuple(domains))


def entry(name="1AAA", methods=("X-RAY DIFFRACTION",), resolution=2.0, rfree=0.22, **said) -> Entry:
    base = dict(entry=name, methods=tuple(methods), resolution=resolution, rfree=rfree, models=1,
                nmr_conformer=None, clashscore=5.0, ramachandran_outliers=0.5,
                rotamer_outliers=1.0, released="2020-01-01", revised="2024-01-01",
                revision="1.2", status="REL", polymer_entities=1, nonpolymer_entities=0,
                oligomeric_count=1)
    base.update(said)
    return Entry(**base)


def chain(asym="A", unobserved=(), first=1, last=180, rscc=0.95, rsrz=0.0, q=None, clashes=(),
          measured=True) -> Chain:
    placed = [n for n in range(first, last + 1) if n not in set(unobserved)]
    return Chain(asym=asym, auth=asym, unobserved=frozenset(unobserved), helices=(), strands=(),
                 rscc={n: rscc for n in placed} if measured and rscc is not None else {},
                 rsrz={n: rsrz for n in placed} if measured and rsrz is not None else {},
                 q_score={n: q for n in placed} if q is not None else {},
                 clashes=frozenset(clashes))


def entity(name="1AAA", number=1, start=21, end=200, chains=None, sequence=None,
           sources=(9606,), isoform=None, accession="P00001", translation=PROTEIN,
           **entry_said) -> Entity:
    """An entity that is residues `start..end` of the protein, mapped by SIFTS."""
    sequence = translation[start - 1:end] if sequence is None else sequence
    length = len(sequence)
    return Entity(
        id=f"{name}_{number}", entry=entry(name, **entry_said), entity=str(number),
        description="protein x", polymer="Protein", sequence=sequence, mutations=0,
        sources=tuple(sources), hosts=(562,),
        alignments=(Alignment("UniProt", accession, isoform, ((1, start, length),), "SIFTS"),),
        chains=tuple(chains) if chains is not None else (chain(last=length),))


def ranked(protein: Inputs, *entities):
    return pdb_policy.rank(protein, list(entities))


def reasons(ranking) -> list:
    return [(r.reason, r.detail) for r in ranking.rejections]


# ------------------------------------------------------------ the protein


def region(label, start, end, kept=True):
    return SimpleNamespace(label=label, start=start, end=end, kept=kept)


def test_a_cleaved_precursors_products_are_its_kept_regions():
    products, removed = pdb_policy.products_of(
        [region("Signal peptide", 1, 24, False), region("Insulin B chain", 25, 54),
         region("C peptide", 57, 87, False), region("Insulin A chain", 90, 110)],
        cleaved=True, length=110, chain_label="Proinsulin")
    assert [(p.index, p.label, p.start, p.end) for p in products] == [
        (0, "Insulin B chain", 25, 54), (1, "Insulin A chain", 90, 110)]
    assert removed == ((1, 24, "Signal peptide"), (57, 87, "C peptide"))


def test_an_uncut_protein_is_one_product_less_what_its_ends_lose():
    products, _ = pdb_policy.products_of([region("Initiator methionine", 1, 1, False)],
                                         cleaved=False, length=154, chain_label="myoglobin")
    assert [(p.label, p.start, p.end) for p in products] == [("myoglobin", 2, 154)]
    # A table of domains, as the twenty's uncut rows carry, is not a table of products.
    products, _ = pdb_policy.products_of(
        [region("Transactivation", 1, 61), region("DNA-binding", 102, 292)],
        cleaved=False, length=393, chain_label="p53")
    assert [(p.start, p.end) for p in products] == [(1, 393)]


def test_a_precursor_with_nothing_named_as_kept_leaves_what_its_removed_regions_leave():
    # Apelin: its peptides overlap, so the resolver keeps only what is removed.
    products, removed = pdb_policy.products_of(
        [region("Signal peptide", 1, 22, False), region("Propeptide", 23, 41, False)],
        cleaved=True, length=77, chain_label="Apelin")
    assert [(p.index, p.label, p.start, p.end) for p in products] == [(0, "Apelin", 42, 77)]
    assert len(removed) == 2
    # And a removed stretch in the middle leaves a product either side of it.
    products, _ = pdb_policy.products_of(
        [region("Signal peptide", 1, 20, False), region("Connecting peptide", 51, 60, False)],
        cleaved=True, length=90, chain_label="X")
    assert [(p.start, p.end) for p in products] == [(21, 50), (61, 90)]


def test_a_domain_is_clipped_to_the_mature_product_it_lies_in():
    products = (Product(0, "Major prion protein", 23, 230),)
    pfam = Domain("Domain", "Prion/Doppel alpha-helical domain", 134, 251, "Pfam PF00377")
    (clipped,) = pdb_policy.domains_within([pfam], products)
    assert (clipped.start, clipped.end, clipped.source) == (134, 230, "Pfam PF00377")
    # Less than half of it inside a product, or none: not a domain of the mature protein.
    assert pdb_policy.domains_within([Domain("Domain", "Leader", 1, 40)], products) == ()
    assert pdb_policy.domains_within([Domain("Domain", "Tail", 240, 253)], products) == ()


def test_pfam_names_a_domain_only_where_uniprot_names_none():
    products = (Product(0, "x", 1, 400),)
    named = Domain("Domain", "Ig-like C1-type", 25, 114)
    same = Domain("Domain", "Immunoglobulin C1-set domain", 30, 110, "Pfam PF07654")
    other = Domain("Domain", "Kinase", 200, 380, "Pfam PF00069")
    kept = pdb_policy.domains_within([same, other, named], products)
    assert [d.name for d in kept] == ["Ig-like C1-type", "Kinase"]


# ------------------------------------------------------------ eligibility


def test_a_whole_mature_chain_is_eligible_and_says_what_it_draws():
    ranking = ranked(inputs(), entity())
    assert not ranking.rejections
    (found,) = ranking.candidates
    (pick,) = found.picks
    assert (found.entry, found.scope, found.useful) == ("1AAA", "product", 180)
    assert (pick.asym, pick.first, pick.last, pick.placed, pick.mismatches) == ("A", 21, 200, 180, ())
    assert found.share == 1.0 and found.coverage == 1.0 and found.shortfall == 0
    assert (found.tier, found.precursor, found.mismatches) == (0, False, 0)


def test_the_source_organism_decides_and_the_expression_host_does_not():
    human_in_coli = entity()                       # hosts (562,): Escherichia coli
    assert ranked(inputs(), human_in_coli).candidates
    for sources in ((10090,), (9606, 10665), ()):
        ranking = ranked(inputs(), entity(sources=sources))
        assert not ranking.candidates and reasons(ranking)[0][0] == "organism"


@pytest.mark.parametrize("methods, eligible", [
    (("X-RAY DIFFRACTION",), True), (("NEUTRON DIFFRACTION", "X-RAY DIFFRACTION"), True),
    (("ELECTRON MICROSCOPY",), True), (("SOLUTION NMR",), True), (("SOLID-STATE NMR",), True),
    (("ELECTRON CRYSTALLOGRAPHY",), True),
    (("SOLUTION SCATTERING",), False), (("SOLUTION NMR", "SOLUTION SCATTERING"), False),
    (("THEORETICAL MODEL",), False), (("INTEGRATIVE",), False), ((), False),
])
def test_only_experimental_methods_that_place_atoms_are_eligible(methods, eligible):
    ranking = ranked(inputs(), entity(methods=methods))
    assert bool(ranking.candidates) is eligible
    if not eligible:
        assert reasons(ranking)[0][0] == "method"


def test_an_entry_that_is_not_released_is_not_eligible():
    assert reasons(ranked(inputs(), entity(status="OBS")))[0] == ("status", "the entry is OBS")


def test_only_the_canonical_sequence_of_this_protein_is_mapped():
    isoform = ranked(inputs(), entity(isoform="P00001-2"))
    assert reasons(isoform) == [("isoform", "mapped to isoform P00001-2")]
    elsewhere = ranked(inputs(), entity(accession="Q99999"))
    assert reasons(elsewhere)[0][0] == "mapping" and "UniProt Q99999" in reasons(elsewhere)[0][1]
    beyond = ranked(inputs(), entity(start=150, end=200, sequence=PROTEIN[149:] + "AAAA"))
    assert reasons(beyond)[0][0] == "mapping" and "beyond the protein" in reasons(beyond)[0][1]


def test_coverage_counts_residues_placed_once_each_and_not_the_construct():
    # 180 residues in the construct, 100 of them never placed: 44% of the product.
    sparse = entity(chains=[chain(unobserved=range(41, 141))])
    ranking = ranked(inputs(), sparse)
    assert not ranking.candidates
    assert reasons(ranking) == [("coverage", "places 44% of a mature product and no whole "
                                             "annotated domain")]
    # Two copies of the chain do not add up: each is measured alone.
    twice = entity(chains=[chain("A", unobserved=range(41, 141)),
                           chain("B", unobserved=range(41, 141))])
    assert not ranked(inputs(), twice).candidates
    # At the threshold itself it is eligible.
    enough = entity(chains=[chain(unobserved=range(1, 55))])     # 126 of 180 = 70%
    assert ranked(inputs(), enough).candidates


def test_letters_are_held_over_what_is_drawn_and_a_short_peptide_may_differ_nowhere():
    def mutated(count, length):
        letters = list(PROTEIN[20:20 + length])
        for k in range(count):
            letters[5 + 3 * k] = "W" if letters[5 + 3 * k] != "W" else "A"
        return "".join(letters)

    short = inputs(products=((21, 60),))
    assert ranked(short, entity(end=60)).candidates
    one_off = ranked(short, entity(end=60, sequence=mutated(1, 40)))
    assert reasons(one_off) == [("fidelity", "1 residue(s) drawn are not the record's, "
                                             "0 allowed over 40")]
    long = inputs(products=((21, 170),))
    assert ranked(long, entity(end=170, sequence=mutated(3, 150))).candidates[0].mismatches == 3
    assert reasons(ranked(long, entity(end=170, sequence=mutated(4, 150))))[0][0] == "fidelity"


def test_a_mismatch_the_entry_never_places_is_not_counted():
    letters = list(PROTEIN[20:])
    letters[3] = "W" if letters[3] != "W" else "A"
    hidden = entity(sequence="".join(letters), chains=[chain(unobserved=[4])])
    assert ranked(inputs(), hidden).candidates[0].mismatches == 0


@pytest.mark.parametrize("said, why", [
    (dict(resolution=3.2), "resolution 3.20"),
    (dict(resolution=None), "no resolution"),
    (dict(rfree=0.34), "R-free 0.340"),
])
def test_a_crystal_is_held_to_its_resolution_and_r_free(said, why):
    ranking = ranked(inputs(), entity(**said))
    assert not ranking.candidates and why in reasons(ranking)[0][1]


def test_a_crystal_is_held_to_the_density_under_the_residues_drawn():
    poor_fit = entity(chains=[chain(rscc=0.6)])
    assert "fit their density poorly" in reasons(ranked(inputs(), poor_fit))[0][1]
    outliers = entity(chains=[chain(rsrz=3.5)])
    assert "RSRZ outliers" in reasons(ranked(inputs(), outliers))[0][1]


def test_an_em_map_is_held_to_its_resolution_and_its_fit():
    em = dict(methods=("ELECTRON MICROSCOPY",), rfree=None)
    good = entity(resolution=3.0, chains=[chain(q=0.6, rscc=None, rsrz=None)], **em)
    (found,) = ranked(inputs(), good).candidates
    assert found.tier == 0 and found.support == 1.0
    coarse = ranked(inputs(), entity(resolution=4.5, chains=[chain(q=0.6)], **em))
    assert "resolution 4.50" in reasons(coarse)[0][1]
    # At 3 A a residue is expected near Q 0.585; 0.3 is far under it.
    faint = ranked(inputs(), entity(resolution=3.0, chains=[chain(q=0.3)], **em))
    assert "fit their map poorly" in reasons(faint)[0][1]
    assert pdb_policy.expected_q(3.0) == pytest.approx(0.585)


def test_an_nmr_model_is_held_to_its_geometry():
    nmr = dict(methods=("SOLUTION NMR",), resolution=None, rfree=None)
    sound = entity(chains=[chain(measured=False)], **nmr)
    assert ranked(inputs(), sound).candidates[0].tier == 0
    strained = ranked(inputs(), entity(chains=[chain(measured=False)],
                                       ramachandran_outliers=7.5, **nmr))
    assert reasons(strained) == [("quality", "7.5% Ramachandran outliers")]


def test_too_many_clashing_residues_are_refused_whatever_the_method():
    crowded = entity(chains=[chain(clashes=range(1, 80))])
    assert "of the residues drawn clash" in reasons(ranked(inputs(), crowded))[0][1]


def test_a_region_nobody_measured_is_eligible_and_ranks_as_unmeasured():
    old = entity(rfree=None, chains=[chain(measured=False)])
    (found,) = ranked(inputs(), old).candidates
    assert found.tier == 1 and found.support is None


# ------------------------------------------------------------ scope


def domain(name="Kinase", start=40, end=139, kind="Domain") -> Domain:
    return Domain(kind, name, start, end)


def test_a_whole_annotated_domain_qualifies_where_the_product_does_not():
    piece = entity(start=35, end=145)               # 111 of the product's 180
    with_domain = ranked(inputs(domains=[domain()]), piece)
    (found,) = with_domain.candidates
    assert found.scope == "domain" and found.domains == ("Kinase",)
    assert reasons(ranked(inputs(), piece))[0][0] == "coverage"


def test_a_domain_is_drawn_with_no_more_than_its_flank():
    # The construct runs sixty residues before the domain; ten of them are drawn.
    long_tail = entity(start=21, end=139)
    protein = inputs(products=((21, 200),), domains=[domain(start=81, end=139)])
    # 119 placed of 180 is under the product's threshold, so the domain is what qualifies.
    (found,) = ranked(protein, long_tail).candidates
    (pick,) = found.picks
    assert found.scope == "domain"
    assert (pick.first, pick.last) == (81 - pdb_policy.DOMAIN_FLANK, 139)
    assert found.useful == pick.placed == 59 + pdb_policy.DOMAIN_FLANK


def test_a_domain_needs_nearly_all_of_itself_and_enough_residues():
    most = entity(start=55, end=145)                # 85 of the domain's 100
    assert not ranked(inputs(domains=[domain()]), most).candidates
    small = entity(start=40, end=65)
    assert not ranked(inputs(domains=[domain(start=40, end=65)]), small).candidates


def test_the_mature_product_ranks_before_any_domain():
    piece = entity("1AAA", start=35, end=145, resolution=1.0)
    whole = entity("2BBB", resolution=2.8)
    ranking = ranked(inputs(domains=[domain()]), piece, whole)
    assert [(c.entry, c.scope) for c in ranking.candidates] == [("2BBB", "product"),
                                                                ("1AAA", "domain")]
    assert pdb_policy.why(*ranking.candidates) == ("it places the mature product; the runner-up "
                                                   "only a domain")


# ------------------------------------------------------------ rank


def test_within_a_tenth_of_the_best_coverage_the_better_experiment_wins():
    whole_but_coarse = entity("1AAA", resolution=2.9)
    nearly_and_sharp = entity("2BBB", resolution=1.2, chains=[chain(unobserved=range(1, 10))])
    ranking = ranked(inputs(), whole_but_coarse, nearly_and_sharp)
    assert [c.entry for c in ranking.candidates] == ["2BBB", "1AAA"]
    assert [c.shortfall for c in ranking.candidates] == [0, 0]
    assert pdb_policy.why(*ranking.candidates) == "resolution: 1.20 Å against 2.90 Å"


def test_more_than_a_tenth_short_of_the_best_coverage_ranks_after_it():
    whole_but_coarse = entity("1AAA", resolution=2.9)
    short_and_sharp = entity("2BBB", resolution=1.2, chains=[chain(unobserved=range(1, 31))])
    ranking = ranked(inputs(), short_and_sharp, whole_but_coarse)
    assert [(c.entry, c.shortfall) for c in ranking.candidates] == [("1AAA", 0), ("2BBB", 1)]
    assert pdb_policy.why(*ranking.candidates) == ("coverage: 100% of what it draws is placed, "
                                                   "against 83%")


def test_a_mature_construct_ranks_before_a_precursors():
    mature = entity("2BBB", resolution=2.5)
    with_leader = entity("1AAA", start=1, end=200, resolution=1.5)
    ranking = ranked(inputs(), with_leader, mature)
    assert [(c.entry, c.precursor) for c in ranking.candidates] == [("2BBB", False),
                                                                    ("1AAA", True)]
    assert "nothing the precursor loses" in pdb_policy.why(*ranking.candidates)
    # A residue or two past a cut is not a precursor's structure; nor is the
    # initiator methionine a construct keeps.
    barely = entity("1AAA", start=18, end=200, resolution=1.5)
    assert [c.precursor for c in ranked(inputs(), barely).candidates] == [False]
    met = inputs(products=((2, 200),), removed=((1, 1, "Initiator methionine"),))
    assert [c.precursor for c in ranked(met, entity(start=1, end=200)).candidates] == [False]


def test_fewer_residues_that_are_not_the_genes_rank_first_then_quality():
    letters = list(PROTEIN[20:])
    letters[9] = "W" if letters[9] != "W" else "A"
    mutant_sharp = entity("1AAA", sequence="".join(letters), resolution=1.0)
    wild_coarse = entity("2BBB", resolution=2.9)
    ranking = ranked(inputs(), mutant_sharp, wild_coarse)
    assert [c.entry for c in ranking.candidates] == ["2BBB", "1AAA"]
    assert pdb_policy.why(*ranking.candidates) == ("fidelity: 0 residue(s) not the gene's "
                                                   "against 1")


def test_a_resolution_ranks_before_none_and_a_measured_region_before_an_unmeasured_one():
    unmeasured = entity("1AAA", resolution=1.5, rfree=None, chains=[chain(measured=False)])
    measured = entity("2BBB", resolution=2.5)
    nmr = entity("3CCC", methods=("SOLUTION NMR",), resolution=None, rfree=None,
                 chains=[chain(measured=False)])
    ranking = ranked(inputs(), unmeasured, nmr, measured)
    # A crystal nobody deposited the data of is still an experiment with a
    # resolution, and ranks before an NMR model; one measured residue by
    # residue ranks before it, however much coarser.
    assert [c.entry for c in ranking.candidates] == ["2BBB", "1AAA", "3CCC"]
    assert "measured residue by residue" in pdb_policy.why(ranking.candidates[0],
                                                          ranking.candidates[1])
    assert pdb_policy.why(ranking.candidates[1], ranking.candidates[2]) == (
        "its method has a resolution (1.50 Å); the runner-up's has none")


def test_ties_break_by_entry_then_chain_and_the_order_of_the_input_never_matters():
    both = entity("2BBB", chains=[chain("B"), chain("A")])
    first = entity("1AAA")
    entities = [both, first, entity("3CCC", resolution=1.9)]
    expected = [("3CCC", "A"), ("1AAA", "A"), ("2BBB", "A"), ("2BBB", "B")]
    for seed in range(6):
        shuffled = entities[:]
        random.Random(seed).shuffle(shuffled)
        ranking = ranked(inputs(), *shuffled)
        assert [(c.entry, c.picks[0].asym) for c in ranking.candidates] == expected
    assert pdb_policy.why(ranking.candidates[1], ranking.candidates[2]) == (
        "entry id: 1AAA before 2BBB")
    assert pdb_policy.why(ranking.candidates[0], None) == "the only eligible candidate"


def test_rejections_are_counted_by_reason():
    ranking = ranked(inputs(), entity("1AAA", sources=(10090,)), entity("2BBB", resolution=9.0),
                     entity("3CCC", sources=(10090,)), entity("4DDD"))
    assert ranking.reasons() == {"organism": 2, "quality": 1}


# ------------------------------------------------------------ several products


TWO = dict(products=((21, 60), (101, 130)), removed=((1, 20, "Signal peptide"),
                                                      (61, 100, "Connecting peptide")))


def test_two_products_from_two_entities_of_one_entry_are_drawn_together_first():
    b = entity("1AAA", 1, start=21, end=60, chains=[chain("A", last=40)])
    a = entity("1AAA", 2, start=101, end=130, chains=[chain("B", last=30)])
    alone = entity("2BBB", 1, start=21, end=60, resolution=1.0, chains=[chain("A", last=40)])
    ranking = ranked(inputs(**TWO), a, b, alone)
    together = ranking.candidates[0]
    assert together.entry == "1AAA" and together.useful == 70 and together.share == 1.0
    assert [(p.product, p.asym, p.first, p.last) for p in together.picks] == [
        (0, "A", 21, 60), (1, "B", 101, 130)]
    # The chains alone are candidates too, and rank after the two together,
    # however sharp.
    assert [len(c.picks) for c in ranking.candidates] == [2, 1, 1, 1]
    assert pdb_policy.why(together, ranking.candidates[1]) == (
        "it draws 2 of the products together; the runner-up 1")


def test_two_chains_that_tie_for_a_product_are_told_apart_by_their_ids():
    b_one = entity("1AAA", 1, start=21, end=60, chains=[chain("A", last=40), chain("C", last=40)])
    a = entity("1AAA", 2, start=101, end=130, chains=[chain("B", last=30)])
    together = ranked(inputs(**TWO), b_one, a).candidates[0]
    assert [(p.product, p.asym) for p in together.picks] == [(0, "A"), (1, "B")]


def test_one_chain_holding_both_products_is_a_precursor_and_not_two_chains():
    single_chain = entity("1AAA", 1, start=21, end=130, chains=[chain("A", last=110)])
    ranking = ranked(inputs(**TWO), single_chain)
    assert sorted(len(c.picks) for c in ranking.candidates) == [1, 1]
    assert all(c.precursor for c in ranking.candidates)


def test_copies_of_one_product_are_never_drawn_as_chains_that_belong_together():
    unit = PROTEIN[:40]
    repeats = inputs(products=((1, 40), (41, 80)), removed=(), translation=unit + unit)
    one = entity("1AAA", 1, start=1, end=40, translation=unit + unit, chains=[chain("A", last=40)])
    two = entity("1AAA", 2, start=41, end=80, translation=unit + unit,
                 chains=[chain("B", last=40)])
    ranking = ranked(repeats, one, two)
    assert all(len(c.picks) == 1 for c in ranking.candidates)
    assert ("copies", "the products drawn together are copies of one sequence") in reasons(ranking)


def test_of_several_hormones_the_longest_does_not_win_for_its_length():
    hormones = inputs(products=((21, 49), (61, 97)), removed=((1, 20, "Signal peptide"),))
    short_whole = entity("1AAA", start=21, end=49, resolution=1.1, chains=[chain(last=29)])
    long_partial = entity("2BBB", start=61, end=97, resolution=2.9,
                          chains=[chain(last=37, unobserved=range(31, 38))])
    ranking = ranked(hormones, long_partial, short_whole)
    # 29 of 29 against 30 of 37: more residues, less of its own product.
    assert [(c.entry, c.useful, round(c.share, 2)) for c in ranking.candidates] == [
        ("1AAA", 29, 1.0), ("2BBB", 30, 0.81)]


def test_no_more_than_three_products_are_drawn_and_only_the_first_three_have_colours():
    spans = ((1, 40), (51, 90), (101, 140), (151, 190))
    four = inputs(products=spans, removed=())
    entities = [entity("1AAA", i + 1, start=a, end=b, chains=[chain("ABCD"[i], last=40)])
                for i, (a, b) in enumerate(spans)]
    together = ranked(four, *entities).candidates[0]
    assert [p.product for p in together.picks] == [0, 1, 2]


# ------------------------------------------------------------ recovery


def unmapped(sequence: str) -> Entity:
    return Entity(id="9XYZ_1", entry=entry("9XYZ"), entity="1", description=None, polymer="Protein",
                  sequence=sequence, mutations=None, sources=(9606,), hosts=(), alignments=(),
                  chains=(chain(last=len(sequence)),))


def test_an_unmapped_entity_is_laid_on_the_protein_where_exactly_one_place_takes_it():
    tagged = "MHHHHHHSSG" + PROTEIN[40:140]
    alignment, why = pdb_policy.recover(inputs(), unmapped(tagged))
    assert why == "" and alignment.provenance == "alignment"
    # The tag is trimmed: entity residue 11 is the protein's 41.
    assert alignment.regions == ((11, 41, 100),)
    recovered = Entity(**{**unmapped(tagged).__dict__, "alignments": (alignment,)})
    (found,) = pdb_policy.rank(inputs(domains=[domain(start=41, end=140)]), [recovered]).candidates
    assert found.picks[0].mapped_by == "alignment" and found.picks[0].first == 41


def test_a_repeat_a_stranger_and_a_mapped_entity_are_not_recovered():
    unit = PROTEIN[:40]
    repeats = inputs(products=((1, 120),), removed=(), translation=unit * 3)
    assert pdb_policy.recover(repeats, unmapped(unit)) == (None, "it lies on the protein in 3 places")
    stranger = "".join(reversed(PROTEIN[:60]))
    assert pdb_policy.recover(inputs(), unmapped(stranger))[1] == (
        "no stretch of it lies on the protein")
    assert pdb_policy.recover(inputs(), entity())[1] == "mapped by SIFTS"
    # Too short to be told from chance.
    assert pdb_policy.recover(inputs(), unmapped(PROTEIN[50:60]))[0] is None


# ------------------------------------------------------------ real entities


def recorded(name: str) -> list:
    whole = json.loads((FIXTURES / name).read_bytes())
    return [pdb_search._entity(e) for e in whole["data"]["polymer_entities"]]


def protein_of(accession: str, products, removed, pfam=()) -> Inputs:
    body = json.loads((FIXTURES / f"uniprot_{accession}.json").read_bytes())
    made = tuple(Product(i, label, a, b) for i, (label, a, b) in enumerate(products))
    return Inputs(slug=accession.lower(), gene=accession, accession=accession, display=accession,
                  translation=body["sequence"]["value"], products=made, removed=removed,
                  domains=pdb_policy.domains_within(
                      list(pdb_search.domains_of(body)) + list(pfam), made))


def test_transthyretins_recorded_entries_rank_by_the_rule():
    ttr = protein_of("P02766", [("Transthyretin", 21, 147)], ((1, 20, "Signal"),))
    ranking = pdb_policy.rank(ttr, recorded("graphql_ttr.json"))
    assert not ranking.rejections
    best = ranking.candidates[0]
    # 8C86: 1.10 A, 116 of the 127 residues placed, which is within a tenth of
    # the full-length 1TTA (1.70 A, no R-free deposited, so unmeasured). 1BM7
    # (2.00 A) places 114, more than a tenth short, and ranks after both.
    assert (best.entry, best.resolution, best.useful, best.shortfall) == ("8C86", 1.1, 116, 0)
    assert [(c.entry, c.shortfall, c.tier) for c in ranking.candidates
            if c.picks[0].asym == "A"] == [("8C86", 0, 0), ("1TTA", 0, 1), ("1BM7", 1, 0)]


def test_relaxins_two_chains_are_drawn_together_from_one_entry():
    relaxin = protein_of("P04090", [("B chain", 25, 53), ("A chain", 162, 185)],
                         ((1, 24, "Signal"), (54, 161, "Connecting peptide")))
    ranking = pdb_policy.rank(relaxin, recorded("graphql_relaxin.json"))
    together = [c for c in ranking.candidates if len(c.picks) == 2]
    # 6RLX, a crystal at 1.50 A whose data nobody deposited, before the NMR
    # model: a resolution before none.
    assert [c.entry for c in together] == ["6RLX", "2MV1"]
    assert ranking.candidates[:2] == together
    assert pdb_policy.why(*together) == ("its method has a resolution (1.50 Å); the "
                                         "runner-up's has none")
    x_ray = together[0]
    assert [(p.product, p.asym, p.first, p.last, p.placed) for p in x_ray.picks] == [
        (0, "B", 26, 52, 27), (1, "A", 162, 185, 24)]


def test_thrombins_protease_domain_qualifies_where_prothrombin_does_not():
    prothrombin = protein_of("P00734", [("Prothrombin", 44, 622)],
                             ((1, 24, "Signal"), (25, 43, "Propeptide")))
    ranking = pdb_policy.rank(prothrombin, recorded("graphql_thrombin.json"))
    (found,) = [c for c in ranking.candidates if c.scope == "domain"]
    assert found.entry == "1PPB" and found.domains == ("Peptidase S1",)
    (pick,) = found.picks
    assert pick.first >= 364 - pdb_policy.DOMAIN_FLANK and pick.last <= 622
    # The light chain alone is neither most of prothrombin nor a domain.
    assert ("coverage" in ranking.reasons())
