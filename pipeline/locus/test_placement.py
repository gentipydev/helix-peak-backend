"""Where a gene MANE puts on a scaffold lies, and which Atlas entry reads a
gene: the census's two fixes, offline (`testing.py` makes the snapshot).

Sixteen of the 7,311 proteins the search offered on 2026-11-13 have their
MANE Select transcript on a GRCh38 patch or alternate locus scaffold, and the
bake once wrote them on the scaffold (`17_KV575245v1_fix`, with no band).
Eighteen had no Atlas entry by MANE's Ensembl gene, fourteen of which another
route finds.
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.locus import bake_locus, check_locus, verify_locus  # noqa: E402
from pipeline.locus.testing import KIR2DS4, POLR2A, built, tables  # noqa: E402
from pipeline.targets import Source, Target  # noqa: E402

RELEASE = {"release": "2026_04"}


@pytest.fixture(scope="module")
def snapshot(tmp_path_factory):
    return built(tmp_path_factory.mktemp("locus"))[2]


def _target(row: dict, accession: str | None = None, start=None, stop=None) -> Target:
    """The target the resolver makes for [row]: its record cut from the
    scaffold MANE puts it on, around its span."""
    source = Source(accession or row["chrom_acc"],
                    row["chrom_start"] if start is None else start,
                    row["chrom_end"] if stop is None else stop,
                    protein_id=row["refseq_prot"], transcript_id=row["refseq_nuc"])
    return Target(slug=row["gene"].lower(), gene=row["gene"], uniprot=row["uniprot"],
                  display=row["name"], source=source, structure=None, cleaved=False,
                  chain_label="")


def _bake(row: dict, snapshot) -> dict:
    return bake_locus.bake(row, snapshot, slug=row["gene"].lower(), uniprot_source=RELEASE)


# -- carried to the chromosome -----------------------------------------------------


def test_the_ends_of_the_alignment_are_carried_exactly():
    scaffold, chromosome = (1, 154723), (7379544, 7534268)
    assert bake_locus.carried(1, scaffold, chromosome, "+") == 7379544
    assert bake_locus.carried(154723, scaffold, chromosome, "+") == 7534268
    assert bake_locus.carried(1, scaffold, chromosome, "-") == 7534268
    assert bake_locus.carried(154723, scaffold, chromosome, "-") == 7379544


def test_a_fix_patch_gene_is_carried_to_its_chromosome(snapshot):
    aliases = tables()["chromAlias"]
    place, placed_from = bake_locus.placed(POLR2A, snapshot.placements("NW_016107299.1"), aliases)
    assert place["sequence"] == "NC_000017.11"
    # A fix patch is within a few bases of its region's length (this one
    # 154,723 against 154,725), so a position moves by the region's start,
    # give or take those two bases.
    assert abs(place["start"] - (7379544 + 104823 - 1)) <= 2
    assert abs(place["end"] - (7379544 + 135071 - 1)) <= 2
    assert place["strand"] == 1
    assert placed_from == {
        "sequence": "NW_016107299.1", "ucsc": "chr17_KV575245v1_fix", "scaffold": "HG2046_PATCH",
        "span": {"start": 104823, "end": 135071, "strand": 1},
        "placement": {"assembly": "GRCh38.p14", "orientation": "+",
                      "scaffold": {"start": 1, "end": 154723},
                      "chromosome": {"start": 7379544, "end": 7534268}},
        "method": "proportional",
    }


def test_a_scaffold_placed_reversed_turns_the_gene(snapshot):
    aliases = tables()["chromAlias"]
    place, placed_from = bake_locus.placed(KIR2DS4, snapshot.placements("NT_187669.1"), aliases)
    # KIR2DS4 reads on the scaffold's minus strand; the scaffold lies reversed
    # on chromosome 19, so the gene reads on its plus strand there, its two
    # ends mirrored from the region's far end.
    assert place["strand"] == 1
    assert place["sequence"] == "NC_000019.10"
    assert place["start"] == bake_locus.carried(36357, (1, 170665), (54698377, 54869056), "-")
    assert place["end"] == bake_locus.carried(20467, (1, 170665), (54698377, 54869056), "-")
    assert 54698377 <= place["start"] < place["end"] <= 54869056
    assert placed_from["ucsc"] == "chr19_KI270915v1_alt"
    assert placed_from["placement"]["orientation"] == "-"


@pytest.mark.parametrize("damage, finding", [
    (lambda rows: [], "0 times"),
    (lambda rows: rows * 2, "2 times"),
    (lambda rows: [dict(rows[0], parent_type="SCAFFOLD")], "not a chromosome"),
    (lambda rows: [dict(rows[0], ori="?")], "orientation"),
    (lambda rows: [dict(rows[0], alt_scaf_stop="100000")], "outside the part"),
])
def test_a_placement_that_cannot_be_read_stops_the_bake(damage, finding, snapshot):
    rows = damage(snapshot.placements("NW_016107299.1"))
    with pytest.raises(LookupError, match=finding):
        bake_locus.placed(POLR2A, rows, tables()["chromAlias"])


def test_a_patch_gene_is_baked_on_its_band(snapshot):
    track = _bake(POLR2A, snapshot)
    assert (track["chromosome"], track["sequence"]) == ("17", "NC_000017.11")
    assert track["locus"] == "17p13.1"
    assert track["placed_from"]["sequence"] == "NW_016107299.1"
    assert track["schema_version"] == 3
    # The record the resolver builds is cut from the scaffold, around the span
    # MANE gives there: the check holds it to the scaffold the track says.
    assert check_locus.problems_of(_target(POLR2A), track) == []


def test_a_reversed_alternate_locus_gene_is_baked_on_its_band(snapshot):
    track = _bake(KIR2DS4, snapshot)
    assert track["locus"] == "19q13.42"
    assert track["span"]["strand"] == 1
    assert check_locus.problems_of(_target(KIR2DS4), track) == []


@pytest.mark.parametrize("damage, finding", [
    (lambda t: t["span"].update(start=t["span"]["start"] + 50), "the placement carries"),
    (lambda t: t["span"].update(strand=-1), "placed + from 1"),
    (lambda t: t["placed_from"].update(method="by eye"), "placed by 'by eye'"),
    (lambda t: t["placed_from"]["placement"].update(orientation="-"), "the placement carries"),
    (lambda t: t["placed_from"].update(ucsc="chr3_KV575245v1_fix"), "not a scaffold of"),
    (lambda t: t.update(schema_version=2), "placed_from in schema 2"),
])
def test_the_check_catches_a_misplaced_gene(damage, finding, snapshot):
    track = copy.deepcopy(_bake(POLR2A, snapshot))
    damage(track)
    found = check_locus.problems_of(_target(POLR2A), track)
    assert any(finding in problem for problem in found), found


def test_the_check_holds_the_record_to_the_scaffold(snapshot):
    track = _bake(POLR2A, snapshot)
    other = _target(POLR2A, accession="NW_009646203.1")
    assert any("not placed from it" in p for p in check_locus.problems_of(other, track))
    narrow = _target(POLR2A, start=110000)
    assert any("outside the record's slice" in p
               for p in check_locus.problems_of(narrow, track))
    # A record cut from the chromosome is held to the chromosome, as before.
    on_chromosome = _target(POLR2A, accession="NC_000017.11", start=1, stop=8000000)
    assert check_locus.problems_of(on_chromosome, track) == []


# -- which Atlas entry ----------------------------------------------------------------


class _Atlas:
    """A snapshot as `atlas_match` reads one."""

    def __init__(self, mane=None, entries=(), hgnc=None, uniprot=None):
        self.mane = mane or {}
        self._entries = set(entries)
        self._hgnc = hgnc or {}
        self._uniprot = uniprot or {}

    def atlas_gene(self, ensembl):
        return {"Ensembl": ensembl} if ensembl in self._entries else None

    def hgnc_ensembl(self, hgnc_id):
        return self._hgnc.get(hgnc_id or "")

    def by_uniprot(self, accession):
        return tuple(self._uniprot.get(accession, ()))


ROW = {"ensembl_nuc": "ENST00000000001.3", "hgnc_id": "HGNC:1", "uniprot": "P00001"}


@pytest.mark.parametrize("atlas, found", [
    # MANE's Ensembl gene, where the Atlas has it.
    (_Atlas(mane={"ENST00000000001": {"Ensembl_Gene": "ENSG00000000001.7"}},
            entries=["ENSG00000000001"], hgnc={"HGNC:1": "ENSG00000000002"},
            uniprot={"P00001": ["ENSG00000000003"]}), ("ENSG00000000001", "mane")),
    # A patch gene's own Ensembl id, which the Atlas never saw: HGNC's.
    (_Atlas(mane={"ENST00000000001": {"Ensembl_Gene": "ENSG00000099999.1"}},
            entries=["ENSG00000000002"], hgnc={"HGNC:1": "ENSG00000000002"}),
     ("ENSG00000000002", "hgnc")),
    # Neither: the one entry that names the protein (C4B_2 is C4B's).
    (_Atlas(entries=["ENSG00000000003"], uniprot={"P00001": ["ENSG00000000003"]}),
     ("ENSG00000000003", "uniprot")),
    # Two entries naming it: none is chosen.
    (_Atlas(entries=["ENSG00000000003", "ENSG00000000004"],
            uniprot={"P00001": ["ENSG00000000003", "ENSG00000000004"]}), None),
    # Nothing at all: GSTT1, on an alternate haplotype the Atlas has no entry for.
    (_Atlas(), None),
])
def test_the_atlas_entry_is_found_by_mane_then_hgnc_then_uniprot(atlas, found):
    assert bake_locus.atlas_match(ROW, atlas) == found


def test_the_snapshot_finds_a_patch_gene_by_hgnc(snapshot):
    assert bake_locus.atlas_match(POLR2A, snapshot) == ("ENSG00000181222", "hgnc")
    track = _bake(POLR2A, snapshot)
    assert (track["expression"]["atlas_match"], track["expression"]["atlas_gene"]) == \
        ("hgnc", "POLR2A")
    assert track["expression"]["ensembl_gene"] == "ENSG00000181222"
    assert track["sources"]["hpa"]["url"].endswith("/ENSG00000181222.json")


def test_a_gene_the_atlas_has_no_entry_for_is_refused(snapshot):
    gstt1 = dict(POLR2A, gene="GSTT1", uniprot="P30711", hgnc_id="HGNC:4641",
                 ensembl_nuc="ENST00000612885.5")
    with pytest.raises(bake_locus.LocusRefused) as refused:
        _bake(gstt1, snapshot)
    assert str(refused.value) == "The Human Protein Atlas has no reading for this gene."


def test_an_entry_the_chain_did_not_choose_is_refused(snapshot):
    other = snapshot.atlas_gene("ENSG00000254647")
    cytoband, aliases, atlas = bake_locus._read(snapshot)
    with pytest.raises(ValueError, match="asked the Atlas for ENSG00000181222"):
        bake_locus.track_of(
            slug="polr2a", gene="POLR2A", uniprot="P24928", row=POLR2A, atlas_gene=other,
            match=("ENSG00000181222", "hgnc"), cytoband=cytoband, aliases=aliases, atlas=atlas,
            mane_release="v1.5", uniprot_source=RELEASE, retrieved="2026-11-13",
            placements=snapshot.placements("NW_016107299.1"))


def test_the_bake_is_the_same_bytes_every_time(snapshot):
    assert bake_locus.encode(_bake(POLR2A, snapshot)) == bake_locus.encode(_bake(POLR2A, snapshot))


# -- NCBI's place, read ----------------------------------------------------------------


def _x_track(locus: str = "Xp22.33") -> dict:
    """A gene in the pseudoautosomal tip of X's short arm."""
    return {"chromosome": "X", "locus": locus, "span": {"start": 1000, "end": 2000},
            "bands": [{"name": "p22.33", "start": 1, "end": 2781479},
                      {"name": "p22.32", "start": 2781480, "end": 4400000},
                      {"name": "q28", "start": 4400001, "end": 156040895}]}


@pytest.mark.parametrize("ncbi, verdict", [
    ("Xp22.33", "equal"),
    ("Xp22.3", "coarser"),
    ("Xp22.33-p22.32", "range"),
    ("Xp22.33|11p", "multi"),
    ("X;Y", "par"),
    ("Xp22.33 and Yp11.32", "par"),
    ("-", "none"),
    (None, "none"),
    ("Xp22.32", "wrong"),
    ("11p15.5", "wrong"),
])
def test_ncbi_places_are_read_as_ncbi_writes_them(ncbi, verdict):
    assert verify_locus.classify(_x_track(), ncbi) == verdict


def test_a_placed_gene_is_read_against_ncbi_as_placed(snapshot):
    track = _bake(POLR2A, snapshot)
    assert verify_locus.classify(track, "17p13.1") == "placed"
    assert verify_locus.classify(track, "17p13 alternate reference locus") == "placed"
    assert verify_locus.classify(track, "17q21") == "wrong"


def test_a_coarser_band_holds_its_finer_ones():
    track = _x_track()
    for coarser in ("X", "Xp", "Xp2", "Xp22", "Xp22.3", "Xp22.33"):
        assert verify_locus.within(track, coarser) is True, coarser
    assert verify_locus.within(track, "Xq") is False
    assert verify_locus.within(track, "Xp22.32") is False
