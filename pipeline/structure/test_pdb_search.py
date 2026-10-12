"""Finding a protein's PDB entries: what is believed, and what is tried again.

The network is never reached. Search pages are made here; RCSB's data API,
UniProt and InterPro answer from `fixtures/pdb/`, recorded on 2026-10-12.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.structure import pdb_search  # noqa: E402
from pipeline.structure.pdb_search import Incomplete, Missing, Unreachable  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "pdb"
WEDNESDAY = datetime(2026, 10, 7, 12, tzinfo=timezone.utc)


def recorded(name: str) -> bytes:
    return (FIXTURES / name).read_bytes()


def _now():
    return WEDNESDAY


def searching(ids, *, total=None, status=200, spoil=None):
    """A `post` that answers RCSB's search over `ids`, a page at a time, and
    notes each page asked for. `spoil(start, said)` may change an answer."""
    asked = []

    def post(url, body):
        assert url == pdb_search.SEARCH
        options = body["request_options"]
        assert options["results_content_type"] == ["experimental"]
        assert options["sort"] == [{"sort_by": "rcsb_id", "direction": "asc"}]
        assert body["return_type"] == "polymer_entity"
        start, rows = options["paginate"]["start"], options["paginate"]["rows"]
        asked.append(start)
        if status != 200:
            return status, b""
        if not ids:
            return 204, b""
        said = {"total_count": len(ids) if total is None else total,
                "result_set": [{"identifier": i} for i in ids[start:start + rows]]}
        if spoil is not None:
            said = spoil(start, said)
            if said is None:
                return 204, b""
        return 200, json.dumps(said).encode()

    post.asked = asked
    return post


def ids_of(count: int) -> list:
    return sorted(f"{1000 + i}_1" for i in range(count))


# ------------------------------------------------------------ discovery


@pytest.mark.parametrize("day, release", [
    (datetime(2026, 10, 7, 0, 0, tzinfo=timezone.utc), "2026-10-07"),    # a Wednesday, at the update
    (datetime(2026, 10, 6, 23, 59, tzinfo=timezone.utc), "2026-09-30"),  # the Tuesday before
    (datetime(2026, 10, 12, 9, tzinfo=timezone.utc), "2026-10-07"),      # the Monday after
    (datetime(2026, 10, 13, 23, tzinfo=timezone.utc), "2026-10-07"),
    (datetime(2026, 10, 14, 1, tzinfo=timezone.utc), "2026-10-14"),
])
def test_a_moment_falls_in_the_release_of_its_latest_wednesday(day, release):
    assert pdb_search.release_of(day) == release


def test_every_page_of_a_search_is_read_in_order():
    ids = ids_of(2 * pdb_search.PAGE + 7)
    post = searching(ids)
    found = pdb_search.search("P02766", post, _now)
    assert found.ids == tuple(ids)
    assert post.asked == [0, pdb_search.PAGE, 2 * pdb_search.PAGE]
    assert found.release == "2026-10-07"
    assert found.snapshot == {"query": "UniProt P02766", "release": "2026-10-07",
                              "entities": len(ids), "sha256": found.digest}
    # The digest names the list: another list, another digest.
    assert found.digest != pdb_search.search("P02766", searching(ids[:-1]), _now).digest


def test_the_search_names_the_accession_as_a_uniprot_one():
    seen = {}

    def post(url, body):
        seen.update(body)
        return 204, b""

    pdb_search.search("P61769", post, _now)
    nodes = seen["query"]["nodes"]
    assert [n["parameters"]["value"] for n in nodes] == ["P61769", "UniProt"]


def test_no_content_is_an_empty_answer_and_not_an_outage():
    found = pdb_search.search("Q6WN34", searching([]), _now)
    assert found.ids == () and found.snapshot["entities"] == 0


@pytest.mark.parametrize("status", [500, 502, 503, 429, 400])
def test_a_search_that_is_not_answered_is_unreachable(status):
    with pytest.raises(Unreachable):
        pdb_search.search("P02766", searching(ids_of(3), status=status), _now)


def test_a_total_that_moves_between_pages_is_incomplete():
    ids = ids_of(pdb_search.PAGE + 5)

    def spoil(start, said):
        return said if start == 0 else {**said, "total_count": said["total_count"] + 1}

    with pytest.raises(Incomplete, match="moved"):
        pdb_search.search("P02766", searching(ids, spoil=spoil), _now)


def test_a_search_that_runs_dry_before_its_total_is_incomplete():
    ids = ids_of(pdb_search.PAGE + 5)
    with pytest.raises(Incomplete, match="ran dry"):
        pdb_search.search("P02766", searching(ids, spoil=lambda s, said: said if s == 0 else None),
                          _now)
    with pytest.raises(Incomplete, match="empty page"):
        pdb_search.search("P02766", searching(
            ids, spoil=lambda s, said: said if s == 0 else {**said, "result_set": []}), _now)


def test_an_entity_twice_or_pages_out_of_order_are_incomplete():
    ids = ids_of(4)
    twice = lambda s, said: {**said, "result_set": said["result_set"][:3] + said["result_set"][:1]}  # noqa: E731
    with pytest.raises(Incomplete, match="twice"):
        pdb_search.search("P02766", searching(ids, spoil=twice), _now)
    turned = lambda s, said: {**said, "result_set": said["result_set"][::-1]}  # noqa: E731
    with pytest.raises(Incomplete, match="order"):
        pdb_search.search("P02766", searching(ids, spoil=turned), _now)


def test_an_answer_that_is_not_the_searchs_is_incomplete():
    def post(url, body):
        return 200, b"<html>maintenance</html>"

    with pytest.raises(Incomplete, match="something else"):
        pdb_search.search("P02766", post, _now)


def test_an_archive_updated_part_way_through_a_search_is_incomplete():
    ids = ids_of(pdb_search.PAGE + 5)
    moments = iter([datetime(2026, 10, 13, 23, 59, tzinfo=timezone.utc),
                    datetime(2026, 10, 14, 0, 1, tzinfo=timezone.utc)])
    with pytest.raises(Incomplete, match="updated part-way"):
        pdb_search.search("P02766", searching(ids), lambda: next(moments))


def test_a_sequence_too_short_for_the_service_is_not_asked_about():
    def post(url, body):
        raise AssertionError("the network was asked")

    assert pdb_search.sequence_search("CYIQNCPLG", 0.95, post, _now).ids == ()
    found = pdb_search.sequence_search("M" * pdb_search.SEQUENCE_MIN, 0.95,
                                       searching(ids_of(2)), _now)
    assert len(found.ids) == 2


# ------------------------------------------------------------ metadata


def answering(name: str, calls: list = None):
    """A `post` that answers the data API from a recorded answer, for the ids asked."""
    whole = json.loads(recorded(name))
    by_id = {e["rcsb_id"]: e for e in whole["data"]["polymer_entities"]}

    def post(url, body):
        assert url == pdb_search.GRAPHQL
        if calls is not None:
            calls.append(body["variables"]["ids"])
        return 200, json.dumps({"data": {"polymer_entities": [
            by_id[i] for i in body["variables"]["ids"] if i in by_id]}}).encode()

    return post


def test_an_entity_is_read_with_its_chains_and_its_entry():
    old, sharp, whole = pdb_search.metadata(["1BM7_1", "8C86_1", "1TTA_1"], answering("graphql_ttr.json"))
    assert (old.id, old.entity, old.entry.entry) == ("1BM7_1", "1", "1BM7")
    assert old.polymer == "Protein" and old.sources == (9606,) and old.hosts == (562,)
    assert len(old.sequence) == 127 and old.sequence.startswith("GPTGTGESK")
    (alignment,) = old.alignments
    assert (alignment.database, alignment.accession, alignment.isoform) == ("UniProt", "P02766", None)
    assert alignment.regions == ((1, 21, 127),) and alignment.provenance == "SIFTS"
    assert [c.asym for c in old.chains] == ["A", "B"]
    a = old.chains[0]
    assert a.auth == "A" and a.unobserved == frozenset(range(1, 10)) | frozenset(range(124, 128))
    assert a.strands and a.helices == ((75, 82),)
    # Density fit, a value a residue the entry places.
    assert set(a.rscc) == set(range(10, 124)) == set(a.rsrz)
    assert 0 < min(a.rscc.values()) and max(a.rscc.values()) <= 1
    assert a.clashes and not a.q_score
    assert old.entry.methods == ("X-RAY DIFFRACTION",) and old.entry.resolution == 2.0
    assert old.entry.status == "REL" and old.entry.released == "1998-08-05"
    assert old.entry.models == 1 and old.entry.nmr_conformer is None
    # Two chains in the file, two operators in the assembly: a tetramer of this
    # alone, with a tetramer's symmetry.
    assert (old.copies, old.alone, old.symmetry) == (4, True, "D2")
    assert sharp.entry.resolution == pytest.approx(1.1) and sharp.entry.rfree is not None
    assert whole.chains[0].unobserved == frozenset()


def test_copies_of_a_protein_are_told_from_its_partners():
    (b2m,) = pdb_search.metadata(["1K5N_2"], answering("graphql_b2m.json"))
    assert (b2m.copies, b2m.alone) == (1, False)          # one chain of an HLA complex
    (cftr,) = pdb_search.metadata(["9DW5_1"], answering("graphql_cftr.json"))
    assert (cftr.copies, cftr.alone) == (1, False)        # beside the fragment it was solved with
    assert cftr.entry.methods == ("ELECTRON MICROSCOPY",) and cftr.chains[0].q_score
    (sod1,) = pdb_search.metadata(["2C9V_1"], answering("graphql_sod1.json"))
    assert (sod1.copies, sod1.alone, sod1.symmetry) == (2, True, "C2")


def test_an_nmr_entry_names_its_representative_conformer():
    a_chain, b_chain = pdb_search.metadata(["2K6T_1", "2K6T_2"], answering("graphql_insl3.json"))
    assert a_chain.entry.methods == ("SOLUTION NMR",) and a_chain.entry.resolution is None
    assert a_chain.entry.nmr_conformer == "1" and a_chain.entry.models > 1
    assert {a.accession for a in a_chain.alignments + b_chain.alignments} == {"P51460"}


def test_entities_are_asked_about_a_batch_at_a_time_in_the_order_asked():
    calls = []
    found = pdb_search.metadata(["8C86_1", "1BM7_1", "1TTA_1"], answering("graphql_ttr.json", calls),
                                batch=2)
    assert calls == [["8C86_1", "1BM7_1"], ["1TTA_1"]]
    assert [e.id for e in found] == ["8C86_1", "1BM7_1", "1TTA_1"]


def test_an_entity_the_data_api_does_not_know_is_incomplete():
    with pytest.raises(Incomplete, match="does not know"):
        pdb_search.metadata(["1BM7_1", "9ZZZ_1"], answering("graphql_ttr.json"))


def test_the_data_apis_own_errors_are_incomplete_and_its_silence_unreachable():
    with pytest.raises(Incomplete, match="FieldUndefined"):
        pdb_search.metadata(["1BM7_1"], lambda url, body: (
            200, json.dumps({"errors": [{"message": "FieldUndefined"}]}).encode()))
    with pytest.raises(Incomplete, match="something else"):
        pdb_search.metadata(["1BM7_1"], lambda url, body: (200, b"not json"))
    with pytest.raises(Unreachable):
        pdb_search.metadata(["1BM7_1"], lambda url, body: (503, b""))


@pytest.mark.parametrize("expression, copies", [
    (None, 1), ("1", 1), ("1,2", 2), ("(1-60)", 60), ("1-3", 3), ("(1,2,5)", 3),
    ("(1-5)(6-10)", 25), ("(X0)(1-60)", 60), ("P", 1),
])
def test_an_assemblys_operators_are_counted(expression, copies):
    assert pdb_search.operators(expression) == copies


# ------------------------------------------------------------ transport and files


def test_no_content_reads_as_empty_and_not_found_as_missing(monkeypatch):
    answers = {"https://x/empty": (204, b""), "https://x/none": (404, b"no"),
               "https://x/there": (200, b"body"), "https://x/odd": (403, b"")}
    monkeypatch.setattr(pdb_search, "_request", lambda request: answers[request.full_url])
    assert pdb_search.get("https://x/empty") == b""
    assert pdb_search.get("https://x/there") == b"body"
    with pytest.raises(Missing):
        pdb_search.get("https://x/none")
    with pytest.raises(Unreachable):
        pdb_search.get("https://x/odd")


def test_an_entrys_file_is_asked_of_pdbe_by_its_lower_case_id():
    asked = []

    def get(url):
        asked.append(url)
        return b"data_1BM7\n"

    assert pdb_search.updated_cif("1BM7", get) == b"data_1BM7\n"
    assert asked == ["https://www.ebi.ac.uk/pdbe/entry-files/download/1bm7_updated.cif"]


# ------------------------------------------------------------ domains


def test_uniprots_domains_are_read_with_exact_ends_only():
    body = json.loads(recorded("uniprot_P00734.json"))
    found = pdb_search.domains_of(body)
    assert [(d.kind, d.name, d.start, d.end) for d in found] == [
        ("Domain", "Gla", 44, 89), ("Domain", "Kringle 1", 108, 186),
        ("Domain", "Kringle 2", 213, 291), ("Domain", "Peptidase S1", 364, 618)]
    assert {d.source for d in found} == {"UniProt"}
    body["features"].append({"type": "Domain", "description": "Loose",
                             "location": {"start": {"value": 5, "modifier": "OUTSIDE"},
                                          "end": {"value": 20, "modifier": "EXACT"}}})
    body["features"].append({"type": "Domain", "description": "Past the end",
                             "location": {"start": {"value": 600, "modifier": "EXACT"},
                                          "end": {"value": 700, "modifier": "EXACT"}}})
    assert len(pdb_search.domains_of(body)) == 4


def test_pfams_domains_are_its_folded_units_only():
    p53 = pdb_search.pfam_domains("P04637", lambda url: recorded("pfam_P04637.json"))
    # Its motifs and conserved sites are not domains.
    assert [(d.kind, d.name, d.start, d.end, d.source) for d in p53] == [
        ("Domain", "P53 DNA-binding domain", 100, 288, "Pfam PF00870")]
    prion = pdb_search.pfam_domains("P04156", lambda url: recorded("pfam_P04156.json"))
    # Nor is a family, or a region Pfam itself calls disordered.
    assert [(d.name, d.start, d.end) for d in prion] == [
        ("Prion/Doppel alpha-helical domain", 134, 251)]


def test_interpro_with_nothing_on_a_protein_answers_no_domains():
    assert pdb_search.pfam_domains("O00631", lambda url: b"") == []
    with pytest.raises(Incomplete, match="more than one page"):
        pdb_search.pfam_domains("P04637", lambda url: json.dumps(
            {"next": "https://www.ebi.ac.uk/…", "results": []}).encode())
    with pytest.raises(Incomplete, match="something else"):
        pdb_search.pfam_domains("P04637", lambda url: b"<html>")
