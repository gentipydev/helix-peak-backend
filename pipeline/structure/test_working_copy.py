"""Which file a curated row's model and fold are cut from, offline.

The twenty are cut from their entries' own legacy files, as they always were;
a row added since is cut from a working copy of its entry
(`structures/working/<slug>.pdb`, what `normalize.working_copy` writes). The
working copy here is relaxin's, made from the cut-down 6RLX fixture the
normalising tests read, so nothing is fetched and no PyMOL is needed.
"""

from __future__ import annotations

import dataclasses
import sys
import urllib.request
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.folding import bake_folding, check_folding  # noqa: E402
from pipeline.structure import experimental, normalize, pdb  # noqa: E402
from pipeline.structure.test_normalize import FIXTURES, RELAXIN, candidate_of  # noqa: E402
from pipeline.targets import BY_SLUG, Chain, Structure  # noqa: E402

# The curated list as it stood on 2026-10-17, before any row was cut from a
# working copy: each keeps the file its stored model was baked from.
TWENTY = ("insulin", "hemoglobin", "myoglobin", "p53", "lysozyme", "relaxin", "oxytocin",
          "somatotropin", "ubiquitin", "dystrophin", "vasopressin", "glucagon", "app",
          "cftr", "erythropoietin", "leptin", "tnf", "sod1", "amylase", "prion")


def _copy_row(tmp_path, monkeypatch, text: str | None = None):
    """A row cut from relaxin's working copy, stored where a working copy is."""
    monkeypatch.setattr(pdb, "WORKING", tmp_path)
    if text is None:
        examined = normalize.examine(
            RELAXIN, candidate_of(RELAXIN, "graphql_relaxin.json", "6RLX", 2),
            (FIXTURES / "6rlx_cut.cif").read_bytes())
        text = normalize.working_copy(examined)
    target = dataclasses.replace(
        BY_SLUG["relaxin"], slug="relaxin-copy",
        structure=Structure(pdb="6RLX", chains=(Chain("chainA", "A"), Chain("chainB", "B")),
                            bonds=True))
    path = tmp_path / "relaxin-copy.pdb"
    path.write_text(text)
    return target, path


@pytest.mark.parametrize("slug", TWENTY)
def test_the_twenty_are_cut_from_their_entries_own_files(slug):
    target = BY_SLUG[slug]
    assert pdb.working_copy_of(target) is None
    assert pdb.source_of(target) == pdb.STRUCTURES / f"{target.structure.pdb}.pdb"
    # Exactly the call the structure bake made before working copies: no
    # source, PyMOL's dashes.
    assert pdb.export_source(target) == (None, True)
    assert not pdb.is_working_copy(pdb.source_of(target))


def test_a_rows_working_copy_is_what_it_is_cut_from(tmp_path, monkeypatch):
    target, path = _copy_row(tmp_path, monkeypatch)
    assert pdb.is_working_copy(path)
    assert pdb.working_copy_of(target) == path
    assert pdb.source_of(target) == path
    assert pdb.export_source(target) == (path, False)
    # And it is read as an entry is, in the record's numbering.
    assert [r.number for r in pdb.chain_residues(target.structure, "A", path)][:2] == [26, 27]


def test_a_working_copy_is_drawn_as_an_entry_chosen_by_rule_is():
    assert pdb.WORKING_COPY_DASHES == experimental.GAP_DASHES


def test_a_file_under_working_that_is_not_a_working_copy_is_refused(tmp_path, monkeypatch):
    legacy = (pdb.STRUCTURES / "6RLX.pdb").read_text()
    target, _ = _copy_row(tmp_path, monkeypatch, text=legacy)
    with pytest.raises(ValueError, match="not a working copy"):
        pdb.source_of(target)


def test_a_row_with_no_file_is_refused_and_nothing_is_fetched(tmp_path, monkeypatch):
    monkeypatch.setattr(pdb, "WORKING", tmp_path)

    def no_network(*_, **__):
        raise AssertionError("a curated row's file is never fetched")

    monkeypatch.setattr(urllib.request, "urlopen", no_network)
    target = dataclasses.replace(
        BY_SLUG["relaxin"], slug="nowhere",
        structure=dataclasses.replace(BY_SLUG["relaxin"].structure, pdb="9ZZZ"))
    with pytest.raises(FileNotFoundError, match="draft_curated"):
        pdb.source_of(target)
    with pytest.raises(FileNotFoundError):
        pdb.export_source(target)


def test_the_folding_check_reads_each_row_by_its_own_file(tmp_path, monkeypatch):
    insulin = BY_SLUG["insulin"]
    twenty = "HELIX and SHEET records of 3I40"
    assert check_folding.reading_of(insulin) == (pdb.STRUCTURES / "3I40.pdb", None, twenty)

    target, path = _copy_row(tmp_path, monkeypatch)
    assert check_folding.reading_of(target) == (
        path, normalize.BRIDGE_REACH,
        bake_folding.WORKING_COPY_SHAPE.format(pdb="6RLX"))
    assert check_folding.reading_of(target, reach=1.9)[1] == 1.9


def test_a_file_given_is_read_as_before_whatever_it_is(tmp_path, monkeypatch):
    # The on-demand path hands its working copy in, and checks the track before
    # it names the entry's own chains: the twenty's statement and its own reach.
    target, path = _copy_row(tmp_path, monkeypatch)
    assert check_folding.reading_of(target, path, 2.5) == (
        path, 2.5, "HELIX and SHEET records of 6RLX")
    assert check_folding.reading_of(target, path) == (
        path, None, "HELIX and SHEET records of 6RLX")
