"""`scripts/draft_curated.py`'s own decisions, offline: what it writes of a
resolved protein and its baked fold. The network, the database and the bake are
not reached; their answers are the shapes the tool reads.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[1]
for path in (str(BACKEND), str(BACKEND / "scripts")):
    if path not in sys.path:
        sys.path.insert(0, path)

import draft_curated  # noqa: E402
from pipeline import check_assets  # noqa: E402
from pipeline.paths import CURATED  # noqa: E402
from pipeline.targets import BY_SLUG, Chain, Region, Source, Structure, Target  # noqa: E402


def test_a_gene_is_named_with_its_slug_and_perhaps_its_entry():
    assert draft_curated.parse_spec("PLN=phospholamban@1zll") == draft_curated.Spec(
        "PLN", "phospholamban", "1ZLL")
    assert draft_curated.parse_spec("s100b=s100b") == draft_curated.Spec("S100B", "s100b", None)
    for wrong in ("PLN", "PLN=Phospholamban", "PLN=phospholamban@1ZL", "=pln"):
        with pytest.raises(ValueError):
            draft_curated.parse_spec(wrong)


@pytest.mark.parametrize("slug", ["insulin", "p53", "amylase", "dystrophin", "tnf", "cftr"])
def test_the_row_written_is_the_row_read_back(slug):
    # Rows as the twenty are written, through every field the draft writes:
    # a slice, a transcript, declared variants, a span, a sampling, bridges.
    target = BY_SLUG[slug]
    namespace = {"Target": Target, "Source": Source, "Region": Region,
                 "Structure": Structure, "Chain": Chain}
    exec("ROWS = (\n" + draft_curated.row_source(target) + ")", namespace)
    assert namespace["ROWS"] == (target,)


def test_a_tube_is_written_as_one():
    target = dataclasses.replace(BY_SLUG["oxytocin"], proprotein=None)
    namespace = {"Target": Target, "Source": Source, "Region": Region,
                 "Structure": Structure, "Chain": Chain}
    exec("ROWS = (\n" + draft_curated.row_source(target) + ")", namespace)
    assert namespace["ROWS"] == (target,)
    assert "representation='tube'" in draft_curated.row_source(target)


def test_what_the_rules_leave_to_a_person_is_marked():
    text = draft_curated.row_source(BY_SLUG["insulin"])
    assert "# HAND: proprotein=" in text and "# HAND: True where" in text
    assert "structures/working/insulin.pdb, a working copy of 3I40" in text
    assert "# HAND: proprotein=" not in draft_curated.row_source(BY_SLUG["p53"])


def _provenance(scope: str, bridges: list) -> dict:
    return {"pdb": "6TQL", "scope": scope, "bridges": bridges, "sampling": 6,
            "representation": "cartoon",
            "chains": [{"node": "chainA", "first": 292, "last": 584}]}


def test_the_structure_is_the_working_copys_as_the_bake_drew_it():
    domain = draft_curated.structure_of(_provenance("domain", [[297, 306]]))
    assert domain == Structure(pdb="6TQL", chains=(Chain("chainA", "A"),), bonds=True,
                               residues=(292, 584), sampling=6)
    product = draft_curated.structure_of(_provenance("product", []))
    assert (product.residues, product.bonds) == (None, False)


def _outcome() -> dict:
    return {
        "chrome": {"pdb": "1ZLL", "modelled": None, "label": "the subunit", "count": 52,
                   "unit": "residues", "sentence": "NMR. One chain of five.",
                   "semantics": "The structure. Drag to turn it."},
        "chains": [{"node": "chainA", "tint": "mature1"}],
    }


def test_the_entry_drafts_every_word_and_derives_every_fact():
    target = dataclasses.replace(BY_SLUG["myoglobin"], slug="myo", display="Myoglobin")
    record = {"protein": {"translation": "M" * 154}, "exons": [{}, {}, {}], "peptides": []}
    entry = draft_curated.catalog_entry(target, record, 20, _outcome(), "Made by MB.")
    template = json.loads(CURATED.read_text())["proteins"][0]
    assert list(entry) == list(template)
    assert entry["facts"] == {"residues": 154, "exons": 3, "chains": 1, "bridges": 0}
    assert entry["chain"] == "DRAFT: myoglobin"
    assert entry["tracks"] == ["clinvar", "constraint", "impact", "structure"]
    row = {**entry, "chrome": entry["structure"]}
    assert check_assets.drafted(row) == ["display", "summary", "chain", "label", "sentence",
                                         "semantics"]
    # A record that names its peptides names no chain, and counts them.
    record["peptides"] = [{}, {}]
    entry = draft_curated.catalog_entry(target, record, 20, _outcome(), "Made by MB.")
    assert entry["chain"] is None and entry["facts"]["chains"] == 2


def test_the_sidecar_says_where_the_working_copy_came_from():
    record = json.dumps({"protein": {"translation": "MEKVQ"}}).encode()
    outcome = {
        "policy": 1, "candidate": "1ZLL:A@1", "measured": {"scene_bytes": 1},
        "sha256": {"glb": "g", "scene": "s", "folding": "f"},
        "provenance": {"entry": "1ZLL", "model": 1, "method": ["SOLUTION NMR"],
                       "resolution": None, "revision": "1.3", "file": {"sha256": "x"},
                       "scope": "product", "domains": [], "bridges": [],
                       "bridges_dropped": [], "secondary_structure": "DSSP",
                       "representation": "cartoon", "sampling": 8,
                       "residues": {"ordered": 52}, "sequence_differences": 0,
                       "chains": [{"node": "chainA", "auth": "A", "first": 1, "last": 52}]},
        "folding_provenance": {"chains": [{"node": "chainA", "pdb_chain": "A", "offset": 0}]},
    }
    said = draft_curated.sidecar(outcome, record, "REMARK   1 1ZLL\n")
    assert said["chains"] == [{"node": "chainA", "auth": "A", "first": 1, "last": 52,
                               "offset": 0}]
    assert said["sha256"]["translation"] == hashlib.sha256(b"MEKVQ").hexdigest()
    assert said["sha256"]["glb"] == "g" and said["entry"] == "1ZLL"
