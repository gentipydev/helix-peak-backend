"""The resolver, offline, on insulin.

`fixtures/P01308.json` is UniProt's entry cut down to what the resolver reads:
the signal and C-peptide from `trafficking/fixtures/` and the names from
`locus/fixtures/` (both release 2026_03), with the two chains and three
disulfides `targets.py` pins from the same entry, and the sequence, which
`targets.py` declares identical to NG_007114's. The record is the real
NG_007114, served from the builder's own cache by `conftest.py`'s `offline`,
so nothing here fetches.

Insulin is one of the twenty, which is the point: its hand-written row is what
an automatic one has to agree with.
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import seed_catalog, upload_tracks  # noqa: E402
from pipeline.constraint import score_protein  # noqa: E402
from pipeline.resolver import resolve as resolver  # noqa: E402
from pipeline.resolver.resolve import (  # noqa: E402
    Feature, IndexRow, Refused, chromosome_of, regions_of, resolve, source_for, target_of,
)
from pipeline.targets import BY_SLUG, Source, partition  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures"

# INS as `protein_index` holds it: MANE v1.5's row (`locus/fixtures/`).
INS = IndexRow(
    uniprot="P01308", gene="INS", name="Insulin", length=110,
    refseq_nuc="NM_000207.3", refseq_prot="NP_000198.1",
    chrom_acc="NC_000011.10", chrom_start=2_159_779, chrom_end=2_161_209,
    refseqgene="NG_007114.1", buildable=True,
)


def _body() -> dict:
    return json.loads((FIXTURES / "P01308.json").read_text(encoding="utf-8"))


def _mutated(body: dict, positions: list[int]) -> dict:
    sequence = list(body["sequence"]["value"])
    for position in positions:
        sequence[position - 1] = "W" if sequence[position - 1] != "W" else "Y"
    return {**body, "sequence": {**body["sequence"], "value": "".join(sequence)}}


# ---------------------------------------------------------------- insulin


def test_insulin_resolves_to_the_table_it_was_curated_as(offline):
    resolution = resolve(INS, offline(_body()), mane_release="v1.5")
    target, curated = resolution.target, BY_SLUG["insulin"]

    # The same cut, piece for piece: where, whether kept, and what each is called by.
    assert [(r.start, r.end, r.kept, r.short) for r in target.regions] == \
        [(r.start, r.end, r.kept, r.short) for r in curated.regions]
    assert [r.label for r in target.regions] == \
        ["Signal peptide", "Insulin B chain", "C peptide", "Insulin A chain"]
    assert target.disulfides == curated.disulfides
    assert target.cleaved is curated.cleaved is True
    assert target.uniprot_variants == curated.uniprot_variants == ()
    assert target.mature_peptides is False
    assert target.slug == "ins"

    # So the partition reads the same between the pieces: the cleavage sites
    # and every residue's numbering.
    ours, theirs = partition(target), partition(curated)
    assert [(p["start"], p["end"], p["origin"], p["kept"]) for p in ours] == \
        [(p["start"], p["end"], p["origin"], p["kept"]) for p in theirs]
    assert [p["label"] for p in ours if "cleavage" in p["label"]] == \
        ["B / C cleavage site", "C / A cleavage site"]


def test_the_record_is_the_builders_and_passes_the_upload_gate(offline):
    resolution = resolve(INS, offline(_body()))
    record = json.loads(resolution.record)

    assert record["gene"] == "INS"
    assert record["protein"]["translation"] == _body()["sequence"]["value"]
    assert len(record["exons"]) == 3
    assert record["peptides"] == []
    # Byte for byte what `build_gene_record.write` puts on disk.
    assert resolution.record == (json.dumps(record, indent=4) + "\n").encode()

    provenance = upload_tracks.validate("record", resolution.target, resolution.record)
    assert provenance["accession"] == "NG_007114.1"
    assert provenance["protein_id"] == "NP_000198.1"
    assert provenance["transcript_id"] == "NM_000207.3"


def test_the_row_has_the_seeds_shape_and_says_how_it_was_made(offline):
    resolution = resolve(INS, offline(_body()), mane_release="v1.5")
    row = resolution.protein

    assert set(row) == set(seed_catalog._PROTEIN_COLUMNS)
    assert row["slug"] == "ins" and row["gene"] == "INS" and row["uniprot"] == "P01308"
    assert row["accession"] == "NG_007114.1" and row["slice_start"] is None
    assert row["display"] == "Insulin"
    assert row["summary"] == ("Made by INS on chromosome 11. Built on demand from UniProt "
                              "P01308 and MANE Select NM_000207.3, not curated by hand.")
    assert row["chain_name"] is None  # cut, so it names its pieces instead
    assert (row["residues"], row["exons"], row["chains"], row["bridges"]) == (110, 3, 2, 3)
    assert row["structure"] is None
    assert row["catalog_order"] is None
    assert row["resolver_version"] == resolver.RESOLVER_VERSION == 1
    assert row["regions"] == partition(resolution.target)
    assert row["disulfides"] == [[31, 96], [43, 109], [95, 100]]

    provenance = row["provenance"]
    assert provenance["prose"] == "templated"
    assert provenance["regions_rule"] == "processing"
    assert provenance["uniprot"]["release"] == "2026_03"
    assert provenance["mane"] == {"release": "v1.5", "transcript": "NM_000207.3",
                                  "protein": "NP_000198.1"}

    assert set(resolution.aliases) == {
        ("Insulin", "display"), ("INS", "gene"), ("ins", "slug"),
        ("P01308", "uniprot"), ("NG_007114.1", "accession"),
    }


def test_a_resolved_row_rebuilds_the_target_it_was_made_as(offline):
    resolution = resolve(INS, offline(_body()))
    # As the database hands it back: JSON arrays, not tuples.
    stored = json.loads(json.dumps(resolution.protein))
    assert target_of(stored) == resolution.target


def test_the_scorer_reads_the_protein_it_was_resolved_with(offline, monkeypatch, tmp_path):
    resolution = resolve(INS, offline(_body()))
    monkeypatch.setattr(score_protein, "DATA", tmp_path)
    destination = tmp_path / resolution.target.mock_asset
    destination.parent.mkdir(parents=True)
    destination.write_bytes(resolution.record)

    assert score_protein.sequence_of(target_of(resolution.protein)) == _body()["sequence"]["value"]


# ---------------------------------------------------------------- refusals


def test_a_few_differences_from_uniprot_are_declared(offline):
    body = _mutated(_body(), [30, 70])
    resolution = resolve(INS, offline(body))

    record = _body()["sequence"]["value"]
    assert resolution.target.uniprot_variants == (
        (30, record[29], body["sequence"]["value"][29]),
        (70, record[69], body["sequence"]["value"][69]),
    )
    assert resolution.protein["provenance"]["uniprot_variants"] == \
        [list(v) for v in resolution.target.uniprot_variants]


def test_more_differences_than_alleles_are_refused(offline):
    with pytest.raises(Refused, match="at 4 residues, more than the 3"):
        resolve(INS, offline(_mutated(_body(), [30, 60, 70, 80])))


def test_a_different_length_is_refused(offline):
    body = _body()
    longer = {**body, "sequence": {**body["sequence"], "value": body["sequence"]["value"] + "K"}}
    with pytest.raises(Refused, match="110 residues long; UniProt P01308 is 111"):
        resolve(INS, offline(longer))


def test_an_unbuildable_entry_is_refused_with_the_index_reason(offline):
    row = replace(INS, buildable=False, unavailable_reason="No MANE Select transcript.")
    with pytest.raises(Refused, match="No MANE Select transcript."):
        resolve(row, offline(_body()))


def test_an_entry_without_a_sequence_is_refused(offline):
    with pytest.raises(Refused, match="carries no sequence"):
        resolve(INS, offline({**_body(), "sequence": {"length": 110}}))


def test_a_bond_to_a_residue_that_is_not_a_cysteine_is_left_out(offline):
    body = _body()
    body["features"].append({
        "type": "Disulfide bond",
        "location": {"start": {"value": 2, "modifier": "EXACT"},
                     "end": {"value": 96, "modifier": "EXACT"}},
        "description": "",
    })
    resolution = resolve(INS, offline(body))
    assert resolution.target.disulfides == ((31, 96), (43, 109), (95, 100))
    assert resolution.protein["provenance"]["disulfides_dropped"] == [[2, 96]]


def test_a_feature_with_an_uncertain_end_is_skipped_and_named(offline):
    body = _body()
    body["features"][1]["location"]["start"]["modifier"] = "UNKNOWN"
    resolution = resolve(INS, offline(body))
    # The A chain is now the one kept piece, so the cut before it is named.
    assert [r.label for r in resolution.target.regions] == \
        ["Signal peptide", "C peptide", "Cleavage site", "Insulin A chain"]
    assert resolution.protein["provenance"]["skipped_features"] == [
        {"type": "Chain", "start": 25, "end": 54, "description": "Insulin B chain"}
    ]


# ---------------------------------------------------------------- the rules


def test_a_record_comes_from_the_refseqgene_or_else_the_mane_span():
    assert source_for(INS) == Source("NG_007114.1", protein_id="NP_000198.1",
                                     transcript_id="NM_000207.3")
    sliced = replace(INS, refseqgene=None)
    assert source_for(sliced) == Source("NC_000011.10", 2_159_779, 2_161_209,
                                        protein_id="NP_000198.1", transcript_id="NM_000207.3")
    with pytest.raises(Refused, match="no transcript of INS"):
        source_for(replace(sliced, chrom_acc=None))


def test_an_initiator_methionine_alone_is_not_a_cut():
    regions, rule, cleaved, chain_label = regions_of(
        [Feature("Initiator methionine", 1, 1, "Removed"),
         Feature("Chain", 2, 147, "Hemoglobin subunit beta")],
        147, "Hemoglobin subunit beta")
    assert [(r.label, r.short, r.kept) for r in regions] == \
        [("Initiator methionine", "Met", False), ("Hemoglobin subunit beta", "", True)]
    assert (rule, cleaved, chain_label) == ("processing", False, "Hemoglobin subunit beta")


def test_no_processing_features_is_one_chain_called_by_its_name():
    assert regions_of([], 393, "Cellular tumor antigen p53") == \
        ((), "processing", False, "Cellular tumor antigen p53")


def test_peptides_inside_a_chain_are_left_out_first():
    regions, rule, cleaved, _ = regions_of(
        [Feature("Signal", 1, 20, ""), Feature("Chain", 21, 180, "Proglucagon"),
         Feature("Peptide", 53, 81, "Glucagon")],
        180, "Pro-glucagon")
    assert [(r.label, r.start, r.end) for r in regions] == \
        [("Signal peptide", 1, 20), ("Proglucagon", 21, 180)]
    assert rule == "processing, without the peptides inside its chains" and cleaved


def test_overlapping_chains_leave_the_leader_and_one_chain():
    regions, rule, cleaved, chain_label = regions_of(
        [Feature("Signal", 1, 17, ""), Feature("Chain", 18, 770, "Amyloid-beta precursor protein"),
         Feature("Chain", 672, 713, "Amyloid-beta protein 42"),
         Feature("Chain", 688, 770, "Gamma-secretase C-terminal fragment 83")],
        770, "Amyloid-beta precursor protein")
    assert [(r.label, r.start, r.end, r.kept) for r in regions] == [
        ("Signal peptide", 1, 17, False),
        ("Amyloid-beta precursor protein", 18, 770, True),
    ]
    assert rule == "leader and one chain: UniProt's processing features overlap"
    assert cleaved and chain_label == "Amyloid-beta precursor protein"


def test_pieces_without_letters_are_numbered():
    regions, *_ = regions_of(
        [Feature("Signal", 1, 19, ""), Feature("Peptide", 20, 28, "Oxytocin"),
         Feature("Chain", 32, 125, "Neurophysin 1"), Feature("Propeptide", 126, 140, "")],
        140, "Oxytocin-neurophysin 1")
    assert [(r.label, r.short) for r in regions] == [
        ("Signal peptide", "S"), ("Oxytocin", "1"), ("Neurophysin 1", "2"), ("Propeptide", "Pro1"),
    ]


def test_a_letter_two_pieces_share_is_numbered_instead():
    regions, *_ = regions_of(
        [Feature("Chain", 1, 30, "Toxin A chain"), Feature("Propeptide", 33, 40, "A peptide"),
         Feature("Chain", 43, 90, "Toxin B chain")],
        90, "Toxin")
    assert [r.short for r in regions] == ["1", "Pro1", "2"]


def test_a_gap_beside_a_lone_chain_is_a_cleavage_site_of_its_own():
    regions, _, cleaved, _ = regions_of(
        [Feature("Signal", 1, 20, ""), Feature("Chain", 23, 100, "Hormone")],
        100, "Prohormone")
    assert [(r.label, r.short, r.start, r.end, r.kept) for r in regions] == [
        ("Signal peptide", "S", 1, 20, False),
        ("Cleavage site", "", 21, 22, False),
        ("Hormone", "", 23, 100, True),
    ]
    target = replace(BY_SLUG["lysozyme"], slug="hormone", aa=100, regions=regions,
                     cleaved=cleaved, disulfides=())
    tiled = partition(target)
    assert [p["label"] for p in tiled] == ["Signal peptide", "Cleavage site", "Hormone"]
    # The chain is numbered from its own first residue, not from the gap.
    assert tiled[-1]["origin"] == 23


def test_a_transit_peptide_says_where_it_takes_the_protein():
    regions, _, cleaved, _ = regions_of(
        [Feature("Transit peptide", 1, 24, "Mitochondrion"), Feature("Chain", 25, 300, "Enzyme")],
        300, "Enzyme")
    assert [(r.label, r.short) for r in regions] == \
        [("Transit peptide (mitochondrion)", "T"), ("Enzyme", "")]
    assert cleaved


def test_chromosomes_are_named_from_their_accessions():
    assert chromosome_of("NC_000011.10") == "11"
    assert chromosome_of("NC_000023.11") == "X"
    assert chromosome_of("NC_000024.10") == "Y"
    assert chromosome_of("NC_012920.1") is None
    assert chromosome_of(None) is None
