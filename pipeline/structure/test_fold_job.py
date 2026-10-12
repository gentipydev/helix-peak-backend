"""How a search for a protein's experimental structure ends, and what never ends it.

RCSB and PDBe answer from `fixtures/pdb/`; nothing reaches the network.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.structure import fold_job, mmcif, pdb_policy, pdb_search  # noqa: E402
from pipeline.structure.normalize import Unsuitable  # noqa: E402
from pipeline.structure.pdb_policy import Inputs, Product  # noqa: E402
from pipeline.structure.pdb_search import Domain, Incomplete, Missing, TooLarge, Unreachable  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "pdb"
NOW = datetime(2026, 10, 12, 9, tzinfo=timezone.utc)


def _body(accession: str) -> dict:
    return json.loads((FIXTURES / f"uniprot_{accession}.json").read_bytes())


RELAXIN = Inputs(
    slug="relaxin", gene="RLN2", accession="P04090", display="Relaxin",
    translation=_body("P04090")["sequence"]["value"],
    products=(Product(0, "B chain", 25, 53), Product(1, "A chain", 162, 185)),
    removed=((1, 24, "Signal"), (54, 161, "Connecting peptide")),
    disulfides=((35, 172), (47, 185), (171, 176)))


def services(metadata: str = "graphql_relaxin.json", ids=None, files=None, by_sequence=(),
             extra=()):
    """`post` and `get` that answer as RCSB and PDBe would, from the fixtures.

    `ids` are what the accession search finds (default: every recorded entity);
    `by_sequence` what the sequence search finds beside them; `extra` more
    recorded entities the data API knows; `files` maps an entry to its file's
    bytes, or to an exception to raise."""
    whole = json.loads((FIXTURES / metadata).read_bytes())["data"]["polymer_entities"]
    by_id = {e["rcsb_id"]: e for e in list(whole) + list(extra)}
    found = sorted(e["rcsb_id"] for e in whole) if ids is None else sorted(ids)
    files = files or {}
    asked = SimpleNamespace(search=0, sequence=0, graphql=[], files=[])

    def post(url, body):
        if url == pdb_search.SEARCH:
            by = body["query"].get("service") == "sequence"
            if by:
                asked.sequence += 1
            else:
                asked.search += 1
            listed = sorted(set(found) | set(by_sequence)) if by else found
            if not listed:
                return 204, b""
            start = body["request_options"]["paginate"]["start"]
            return 200, json.dumps({"total_count": len(listed), "result_set": [
                {"identifier": i} for i in listed[start:start + pdb_search.PAGE]]}).encode()
        assert url == pdb_search.GRAPHQL
        asked.graphql.append(body["variables"]["ids"])
        return 200, json.dumps({"data": {"polymer_entities": [
            by_id[i] for i in body["variables"]["ids"] if i in by_id]}}).encode()

    def get(url):
        entry = url.rsplit("/", 1)[-1].split("_")[0].upper()
        asked.files.append(entry)
        answer = files.get(entry, Missing(url))
        if isinstance(answer, Exception):
            raise answer
        return answer

    return post, get, asked


def choose(protein=RELAXIN, **said):
    kwargs = {k: said.pop(k) for k in ("limit", "skip", "recovery", "draw", "draws")
              if k in said}
    post, get, asked = services(**said)
    return fold_job.choose(protein, get=get, post=post, now=lambda: NOW, **kwargs), asked


RLX = (FIXTURES / "6rlx_cut.cif").read_bytes()
MV1 = (FIXTURES / "2mv1_cut.cif").read_bytes()
# 6RLX's file with its SIFTS map made another protein's: every chain of it is
# refused, by a documented rule.
RLX_ELSEWHERE = RLX.replace(b"P04090", b"Q99999")


def test_the_best_candidate_that_passes_its_own_file_is_chosen():
    choice, asked = choose(files={"6RLX": RLX, "2MV1": MV1})
    assert choice.outcome == "chosen"
    (only,) = choice.examinations
    assert only.candidate.entry == "6RLX" and only.examined is not None
    assert choice.winner is only
    assert [c.node for c in only.examined.chains] == ["chainA", "chainB"]
    # One search, and the one file: the runner-up's is never asked for.
    assert asked.files == ["6RLX"] and asked.search == 1
    assert choice.discovery.snapshot["entities"] == 4
    assert choice.entries["6RLX"].methods == ("X-RAY DIFFRACTION",)
    assert set(choice.entities) == {"2MV1_1", "2MV1_2", "6RLX_1", "6RLX_2"}
    assert fold_job.candidate_id(only.candidate) == "6RLX:B@26,A@162"


def test_a_candidate_its_file_refuses_gives_way_to_the_next():
    choice, asked = choose(files={"6RLX": RLX_ELSEWHERE, "2MV1": MV1})
    assert choice.outcome == "chosen" and choice.winner.candidate.entry == "2MV1"
    refused = [e for e in choice.examinations if e.refused is not None]
    assert refused[0].candidate.entry == "6RLX" and refused[0].refused.reason == "mapping"
    assert asked.files == ["6RLX", "2MV1"]
    assert choice.winner.examined.model == "1"       # 2MV1's representative conformer


def test_an_entry_pdbe_has_no_file_of_is_refused_for_that_and_not_asked_for_again():
    choice, asked = choose(files={"2MV1": MV1})
    first = choice.examinations[0]
    assert first.candidate.entry == "6RLX" and first.refused.reason == "file"
    assert choice.outcome == "chosen" and choice.winner.candidate.entry == "2MV1"
    nothing, asked = choose(files={})
    assert nothing.outcome == "none"
    # Several candidates in two entries (the two chains together, and each
    # alone): two downloads, and every candidate answered.
    assert asked.files == ["6RLX", "2MV1"]
    assert len(nothing.examinations) == len(nothing.ranking.candidates) > 4
    assert {e.refused.reason for e in nothing.examinations} == {"file"}


def test_a_search_that_finds_nothing_is_complete_and_nothing_passed():
    choice, asked = choose(ids=[])
    assert choice.outcome == "none" and choice.winner is None
    assert choice.ranking.candidates == [] and choice.examinations == []
    assert asked.files == [] and asked.graphql == []


@pytest.mark.parametrize("trouble", [Unreachable("PDBe did not answer"), Incomplete("half a page")])
def test_an_outage_is_raised_and_never_read_as_nothing_suitable(trouble):
    with pytest.raises(type(trouble)):
        choose(files={"6RLX": trouble, "2MV1": MV1})


def test_a_search_that_cannot_be_made_is_raised():
    def post(url, body):
        return 503, b""

    with pytest.raises(Unreachable):
        fold_job.choose(RELAXIN, post=post, get=lambda url: b"", now=lambda: NOW)


def test_a_file_that_cannot_be_read_is_no_verdict_on_its_entry():
    with pytest.raises(mmcif.CifError):
        choose(files={"6RLX": b"<html>try again later</html>", "2MV1": MV1})


def test_an_entry_too_large_to_read_is_refused_for_its_size():
    choice, _ = choose(files={"6RLX": TooLarge("big"), "2MV1": MV1})
    assert choice.examinations[0].refused.reason == "size"
    assert choice.outcome == "chosen" and choice.winner.candidate.entry == "2MV1"


def test_a_bound_on_downloads_leaves_the_search_incomplete_and_never_decided():
    choice, asked = choose(files={"6RLX": RLX_ELSEWHERE, "2MV1": MV1}, limit=1)
    # One file was read and its candidates refused; 2MV1's was not asked for.
    assert choice.outcome == "incomplete" and choice.winner is None
    assert asked.files == ["6RLX"]
    assert {e.candidate.entry for e in choice.examinations} == {"6RLX"}
    # The next attempt goes on from there, and finds it.
    done = [fold_job.candidate_id(e.candidate) for e in choice.examinations]
    again, asked = choose(files={"2MV1": MV1}, limit=1, skip=done)
    assert again.outcome == "chosen" and asked.files == ["2MV1"]


def test_every_candidate_refused_is_a_complete_search_with_nothing_suitable():
    choice, _ = choose(files={"6RLX": RLX_ELSEWHERE,
                              "2MV1": MV1.replace(b"P04090", b"Q99999")})
    assert choice.outcome == "none"
    assert {e.refused.reason for e in choice.examinations} == {"mapping"}


# ------------------------------------------------------------ the bake


def drawing(refuse=(), broken=None):
    """A `draw` that bakes nothing: it refuses the entries named, by a rule,
    breaks on `broken`, and says what it was handed."""
    seen = []

    def draw(candidate, examined, entry, entity):
        seen.append((fold_job.candidate_id(candidate), entry.entry,
                     None if entity is None else entity.id))
        if candidate.entry == broken:
            raise RuntimeError("PyMOL exited with status 1")
        if candidate.entry in refuse:
            raise Unsuitable("binding", "A ring of chainA lies on no stretch of its chain")
        return f"baked {candidate.entry}"

    return draw, seen


def test_a_candidate_is_baked_before_it_is_chosen_and_keeps_what_was_baked():
    draw, seen = drawing()
    choice, _ = choose(files={"6RLX": RLX, "2MV1": MV1}, draw=draw)
    assert choice.outcome == "chosen" and choice.winner.built == "baked 6RLX"
    # Two chains drawn, so no one entity is the one drawn.
    assert seen == [("6RLX:B@26,A@162", "6RLX", None)]


def test_a_candidate_the_bake_refuses_by_a_rule_gives_way_to_the_next():
    draw, seen = drawing(refuse={"6RLX"})
    choice, asked = choose(files={"6RLX": RLX, "2MV1": MV1}, draw=draw, draws=9)
    assert choice.outcome == "chosen" and choice.winner.candidate.entry == "2MV1"
    assert choice.winner.built == "baked 2MV1"
    turned_down = [e for e in choice.examinations if not e.passed]
    assert turned_down and {e.candidate.entry for e in turned_down} == {"6RLX"}
    # Its file passed; the bake is what refused it, and says so.
    assert all(e.examined is not None and e.refused.reason == "binding" for e in turned_down)
    assert choice.refused == [fold_job.candidate_id(e.candidate) for e in turned_down]


def test_a_bake_that_breaks_is_raised_and_never_a_refusal():
    draw, _ = drawing(broken="6RLX")
    with pytest.raises(RuntimeError, match="PyMOL exited"):
        choose(files={"6RLX": RLX, "2MV1": MV1}, draw=draw)


def test_a_bound_on_bakes_leaves_the_search_incomplete_and_the_next_one_untried():
    draw, seen = drawing(refuse={"6RLX"})
    choice, _ = choose(files={"6RLX": RLX, "2MV1": MV1}, draw=draw, draws=1)
    assert choice.outcome == "incomplete" and choice.winner is None
    assert len(seen) == 1
    # What was not tried is not recorded, so the next attempt tries it.
    assert [fold_job.candidate_id(e.candidate) for e in choice.examinations] == choice.refused
    again, _ = choose(files={"6RLX": RLX, "2MV1": MV1}, draw=draw, draws=9, skip=choice.refused)
    assert again.outcome == "chosen" and seen[0][0] not in [
        fold_job.candidate_id(e.candidate) for e in again.examinations]


def test_every_candidate_the_bake_refuses_is_a_complete_search_with_nothing_suitable():
    draw, seen = drawing(refuse={"6RLX", "2MV1"})
    choice, _ = choose(files={"6RLX": RLX, "2MV1": MV1}, draw=draw, draws=99)
    assert choice.outcome == "none" and choice.winner is None
    assert len(seen) == len(choice.ranking.candidates) == len(choice.refused)


# ------------------------------------------------------------ recovery


def unmapped_copy(of: str, name: str) -> dict:
    """A recorded entity under another name, as SIFTS has not mapped it yet."""
    whole = json.loads((FIXTURES / "graphql_relaxin.json").read_bytes())["data"]["polymer_entities"]
    made = json.loads(json.dumps(next(e for e in whole if e["rcsb_id"] == of)))
    entry = name.split("_")[0]
    made["rcsb_id"] = name
    made["rcsb_polymer_entity_align"] = None
    made["rcsb_polymer_entity_container_identifiers"]["entry_id"] = entry
    made["entry"]["rcsb_id"] = entry
    for assembly in made["entry"].get("assemblies") or []:
        assembly["rcsb_id"] = f"{entry}-1"
    return made


def test_an_entity_sifts_has_not_mapped_is_found_by_sequence_and_laid_on_the_protein():
    new_a = unmapped_copy("6RLX_1", "9NEW_1")
    choice, asked = choose(by_sequence=["9NEW_1"], extra=[new_a], files={})
    assert asked.sequence == 2                       # one search a product
    assert choice.recovered == ("9NEW_1",) and choice.recovery_refusals == {}
    assert asked.graphql == [["2MV1_1", "2MV1_2", "6RLX_1", "6RLX_2", "9NEW_1"]]
    found = [c for c in choice.ranking.candidates if c.entry == "9NEW"]
    assert found and all(p.mapped_by == "alignment" for c in found for p in c.picks)
    assert found[0].picks[0].first == 162            # the A chain, where it lies


def test_an_entity_found_by_sequence_that_sifts_maps_elsewhere_is_not_recovered():
    whole = json.loads((FIXTURES / "graphql_relaxin.json").read_bytes())["data"]["polymer_entities"]
    other = json.loads(json.dumps(next(e for e in whole if e["rcsb_id"] == "6RLX_1")))
    other["rcsb_id"] = "9OTH_1"
    choice, _ = choose(by_sequence=["9OTH_1"], extra=[other], files={})
    assert choice.recovery_refusals == {"9OTH_1": "mapped by SIFTS"}
    assert not [c for c in choice.ranking.candidates if c.entry == "9OTH"]


def test_recovery_can_be_left_out():
    choice, asked = choose(recovery=False, files={})
    assert asked.sequence == 0 and choice.recovered == ()


# ------------------------------------------------------------ inputs


def test_inputs_are_a_targets_products_and_its_annotated_domains():
    def region(label, start, end, kept=True):
        return SimpleNamespace(label=label, start=start, end=end, kept=kept)

    target = SimpleNamespace(
        slug="f2", gene="F2", uniprot="P00734", display="Prothrombin", cleaved=True,
        chain_label="Prothrombin", disulfides=((336, 482),),
        regions=(region("Signal peptide", 1, 24, False), region("Propeptide", 25, 43, False),
                 region("Prothrombin", 44, 622)))
    body = _body("P00734")
    pfam = [Domain("Domain", "Trypsin", 364, 613, "Pfam PF00089"),
            Domain("Domain", "Somewhere new", 300, 350, "Pfam PF99999")]
    made = fold_job.inputs_of(target, body["sequence"]["value"], body, pfam)
    assert (made.slug, made.gene, made.accession, made.display) == (
        "f2", "F2", "P00734", "Prothrombin")
    assert [(p.label, p.start, p.end) for p in made.products] == [("Prothrombin", 44, 622)]
    assert made.removed == ((1, 24, "Signal peptide"), (25, 43, "Propeptide"))
    assert made.disulfides == ((336, 482),) and made.kept == 579
    # UniProt's own names first; Pfam's only where UniProt names none.
    assert [(d.name, d.source) for d in made.domains] == [
        ("Gla", "UniProt"), ("Kringle 1", "UniProt"), ("Kringle 2", "UniProt"),
        ("Somewhere new", "Pfam PF99999"), ("Peptidase S1", "UniProt")]
