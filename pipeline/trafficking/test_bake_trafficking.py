"""The trafficking bake and its check, offline.

The entries under `fixtures/` are UniProt's own (release 2026_03), cut down to
what the bake reads. Nothing here fetches: the shared fetch is exercised
against a response this file builds.
"""

from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import fetch_tracks, uniprot, upload_tracks  # noqa: E402
from pipeline.targets import BY_SLUG  # noqa: E402
from pipeline.trafficking import bake_trafficking, check_trafficking  # noqa: E402
from pipeline.uniprot import Entry  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _entry(accession: str, retrieved: str = "2026-09-27") -> Entry:
    body = json.loads((FIXTURES / f"{accession}.json").read_text(encoding="utf-8"))
    return Entry(accession, body, "2026_03", "2026-09-02", retrieved)


def _track(slug: str) -> dict:
    target = BY_SLUG[slug]
    return bake_trafficking.payload(target, _entry(target.uniprot))


@pytest.fixture
def records(monkeypatch):
    """Every record as long as the table says, where the check looks for it."""
    monkeypatch.setattr(check_trafficking, "record_length", lambda target: target.aa)


# -- reading an entry ---------------------------------------------------------


def test_cftr_crosses_twelve_times():
    spans = _track("cftr")["transmembrane"]
    assert [(s["start"], s["end"]) for s in spans] == [
        (78, 98), (123, 146), (196, 216), (223, 243), (299, 319), (340, 358),
        (859, 879), (919, 939), (991, 1011), (1014, 1034), (1096, 1116), (1131, 1151),
    ]
    assert spans[0]["description"] == "Helical; Name=1"
    assert spans[0]["evidence"] == ["ECO:0000269"]


def test_tnf_is_anchored_by_its_one_span_and_says_which_way():
    track = _track("tnf")
    (span,) = track["transmembrane"]
    assert (span["start"], span["end"]) == (36, 56)
    assert "Signal-anchor for type II" in span["description"]
    # Predicted, not observed, and the payload says so.
    assert span["evidence"] == ["ECO:0000255"]
    # Myristoylation is a lipidation too, and not a GPI anchor.
    assert track["gpi_anchor"] is None
    entry_wide = [p for p in track["location"] if p["molecule"] is None]
    assert [(p["location"], p["topology"]) for p in entry_wide] == [
        ("Cell membrane", "Single-pass type II membrane protein")]
    soluble = [p for p in track["location"]
               if p["molecule"] == "Tumor necrosis factor, soluble form"]
    assert [p["location"] for p in soluble] == ["Secreted"]


def test_the_prion_protein_is_anchored_at_230_and_loses_the_rest():
    track = _track("prion")
    assert track["transmembrane"] == []
    assert track["gpi_anchor"] == {
        "site": 230,
        "residue": "GPI-anchor amidated serine",
        "signal": {"start": 231, "end": 253},
        "evidence": ["ECO:0000250"],
    }
    assert ("Cell membrane", "Lipid-anchor, GPI-anchor") in [
        (p["location"], p["topology"]) for p in track["location"]]


def test_insulin_crosses_nothing_and_is_secreted():
    track = _track("insulin")
    assert track["transmembrane"] == []
    assert track["gpi_anchor"] is None
    assert [p["location"] for p in track["location"]] == ["Secreted"]


def test_an_entry_with_no_location_says_none_rather_than_guessing():
    track = _track("hemoglobin")
    assert track["location"] == []
    assert track["transmembrane"] == []


def test_the_payload_names_the_release_and_the_day():
    track = _track("insulin")
    assert (track["release"], track["release_date"], track["retrieved"]) == (
        "2026_03", "2026-09-02", "2026-09-27")
    assert (track["entry_version"], track["sequence_version"]) == (284, 1)
    assert (track["slug"], track["gene"], track["uniprot"], track["residues"]) == (
        "insulin", "INS", "P01308", 110)


def test_two_bakes_of_one_entry_differ_only_on_the_day_fetched():
    target = BY_SLUG["cftr"]
    first = bake_trafficking.encode(bake_trafficking.payload(target, _entry("P13569")))
    again = bake_trafficking.encode(bake_trafficking.payload(target, _entry("P13569")))
    later = bake_trafficking.payload(target, _entry("P13569", retrieved="2026-10-01"))
    assert first == again
    assert json.loads(first) == {**later, "retrieved": "2026-09-27"}


def test_a_span_with_an_unknown_end_is_refused():
    entry = _entry("P01375")
    entry.body["features"][0]["location"]["end"] = {"modifier": "UNKNOWN"}
    with pytest.raises(ValueError):
        bake_trafficking.payload(BY_SLUG["tnf"], entry)


def test_an_entry_for_another_protein_is_refused():
    with pytest.raises(ValueError):
        bake_trafficking.payload(BY_SLUG["insulin"], _entry("P04156"))


def test_an_answer_that_names_no_release_is_refused():
    entry = _entry("P01308")
    unnamed = Entry(entry.accession, entry.body, None, None, entry.retrieved)
    with pytest.raises(ValueError):
        bake_trafficking.payload(BY_SLUG["insulin"], unnamed)


# -- the shared fetch ---------------------------------------------------------


class _Response(io.BytesIO):
    def __init__(self, body: dict, headers: dict):
        super().__init__(json.dumps(body).encode())
        self.headers = headers

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_the_record_builder_and_the_bake_share_one_fetch(monkeypatch):
    asked = []

    def urlopen(url, timeout):
        asked.append(url)
        return _Response(
            {"primaryAccession": "P01308", "sequence": {"value": "MALW", "length": 4}},
            {"X-UniProt-Release": "2026_03", "X-UniProt-Release-Date": "02-September-2026"},
        )

    monkeypatch.setattr(uniprot.urllib.request, "urlopen", urlopen)
    assert uniprot.canonical_sequence("P01308") == "MALW"
    entry = uniprot.fetch_entry("P01308")
    assert (entry.release, entry.release_date) == ("2026_03", "2026-09-02")
    assert asked == ["https://rest.uniprot.org/uniprotkb/P01308.json"] * 2

    from pipeline.mock import build_gene_record
    assert build_gene_record.canonical_sequence is uniprot.canonical_sequence


def test_the_release_date_is_read_in_english_whatever_the_locale():
    assert uniprot.release_day("02-September-2026") == "2026-09-02"
    assert uniprot.release_day("15-January-2027") == "2027-01-15"
    assert uniprot.release_day(None) is None


# -- the check ----------------------------------------------------------------


def test_every_fixture_passes_the_check(records):
    for slug in ("insulin", "hemoglobin", "tnf", "cftr", "prion"):
        assert check_trafficking.problems_of(BY_SLUG[slug], _track(slug)) == []


def test_a_span_past_the_last_residue_is_caught(records):
    track = _track("cftr")
    track["transmembrane"][-1]["end"] = 1481
    assert any("not within 1..1480" in p
               for p in check_trafficking.problems_of(BY_SLUG["cftr"], track))


def test_spans_out_of_order_are_caught(records):
    track = _track("cftr")
    track["transmembrane"][1]["start"] = 90
    assert any("overlaps" in p for p in check_trafficking.problems_of(BY_SLUG["cftr"], track))


def test_a_length_the_record_does_not_have_is_caught(monkeypatch):
    monkeypatch.setattr(check_trafficking, "record_length", lambda target: target.aa - 1)
    found = check_trafficking.problems_of(BY_SLUG["insulin"], _track("insulin"))
    assert found == ["insulin: 110 residues, the record translates 109"]


def test_no_record_to_measure_against_is_a_problem(monkeypatch):
    monkeypatch.setattr(check_trafficking, "record_length", lambda target: None)
    found = check_trafficking.problems_of(BY_SLUG["insulin"], _track("insulin"))
    assert len(found) == 1 and "fetch_tracks.py --kind record" in found[0]


def test_the_anchor_has_to_be_the_signal_the_table_names(records):
    track = _track("prion")
    track["gpi_anchor"] = None
    assert any("GPI-anchor signal" in p
               for p in check_trafficking.problems_of(BY_SLUG["prion"], track))
    # And one the table does not name is caught the other way round.
    stray = _track("insulin")
    stray["gpi_anchor"] = {"site": 100, "residue": "GPI-anchor amidated serine",
                           "signal": {"start": 101, "end": 110}, "evidence": []}
    assert any("GPI-anchor signal" in p
               for p in check_trafficking.problems_of(BY_SLUG["insulin"], stray))


def test_a_track_that_names_another_protein_is_caught(records):
    found = check_trafficking.problems_of(BY_SLUG["prion"], _track("insulin"))
    assert any("slug" in p for p in found) and any("gene" in p for p in found)


# -- the uploader and the fetcher know the kind -------------------------------


def test_the_uploader_refuses_what_the_check_refuses(records):
    target = BY_SLUG["cftr"]
    track = _track("cftr")
    provenance = upload_tracks.validate("trafficking", target, bake_trafficking.encode(track))
    assert provenance["release"] == "2026_03"
    assert provenance["retrieved"] == "2026-09-27"
    assert provenance["source"] == "UniProtKB"
    track["transmembrane"][0]["start"] = 0
    with pytest.raises(ValueError):
        upload_tracks.validate("trafficking", target, bake_trafficking.encode(track))


def test_the_uploader_and_the_fetcher_agree_on_where_the_file_is():
    target = BY_SLUG["tnf"]
    assert fetch_tracks.asset_path("trafficking", target) == \
        bake_trafficking.trafficking_asset(target)
    assert upload_tracks.asset_of("trafficking", target) == \
        fetch_tracks.DATA / bake_trafficking.trafficking_asset(target)
    assert "trafficking" in fetch_tracks.KINDS
