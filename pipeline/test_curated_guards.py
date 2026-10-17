"""What keeps a curated row from going out unfinished or in the wrong place,
offline: its words a person's, its facts the record's, and its seed never
written over a protein built on demand. No database is reached: the seed's
query is answered by a stub.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import check_assets, seed_catalog  # noqa: E402
from pipeline.targets import BY_SLUG  # noqa: E402


class _Live:
    """A connection whose one query answers with these protein rows."""

    def __init__(self, rows):
        self.rows, self.asked = rows, []

    def execute(self, sql, parameters):
        self.asked.append(parameters)
        return self

    def fetchall(self):
        return self.rows


def _rows(*slugs):
    curated = seed_catalog.curated_rows()
    return [(seed_catalog.protein_row(BY_SLUG[s], curated[s]), [], []) for s in slugs]


def test_a_seed_over_its_own_curated_rows_is_no_conflict():
    live = _Live([("insulin", "INS", 0, 0), ("lysozyme", "LYZ", 3, 0)])
    assert seed_catalog.conflicts(_rows("insulin", "lysozyme"), live) == []
    assert live.asked == [(["insulin", "lysozyme"], ["INS", "LYZ"])]


def test_a_seed_never_writes_over_a_protein_built_on_demand():
    # A protein built on demand under the curated slug: the upsert would take it.
    found = seed_catalog.conflicts(_rows("tnf"), _Live([("tnf", "TNF", None, 1)]))
    assert found == ["tnf: the live row is a protein built on demand (resolver version 1); "
                     "a curated row would write over it"]
    # The gene served under another slug: the seed would be refused whole.
    found = seed_catalog.conflicts(_rows("insulin"), _Live([("ins", "INS", None, 1)]))
    assert found == ["insulin: INS is already served as 'ins'"]


def test_drafted_words_are_found_wherever_they_are():
    row = {"display": "Phospholamban", "summary": "DRAFT: Made by PLN.", "chain": None,
           "chrome": {"label": "the fold", "sentence": " DRAFT: NMR.", "semantics": "x"}}
    assert seed_catalog.drafted_fields(row) == ["summary", "sentence"]
    assert check_assets.drafted(row) == ["summary", "sentence"]
    row["summary"], row["chrome"]["sentence"] = "Holds back the calcium pump.", "One of five."
    assert seed_catalog.drafted_fields(row) == [] == check_assets.drafted(row)


def test_a_real_seed_refuses_drafted_words(monkeypatch, capsys):
    curated = seed_catalog.curated_rows()
    curated["insulin"] = {**curated["insulin"], "summary": "DRAFT: the hormone"}
    monkeypatch.setattr(seed_catalog, "curated_rows", lambda: curated)
    monkeypatch.setenv("DATABASE_URL", "postgresql://nowhere.invalid/none")
    monkeypatch.setattr(sys, "argv", ["seed_catalog.py", "--target", "insulin"])
    with pytest.raises(SystemExit, match="insulin.summary"):
        seed_catalog.main()
    # A dry run prints the rows and says what a real run would refuse.
    monkeypatch.setattr(sys, "argv", ["seed_catalog.py", "--dry-run", "--target", "insulin"])
    assert seed_catalog.main() == 0
    assert "words still drafted, which a real run refuses: insulin.summary" in capsys.readouterr().err


def test_an_unknown_target_is_refused(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["seed_catalog.py", "--dry-run", "--target", "nothing"])
    with pytest.raises(SystemExit, match="no such target"):
        seed_catalog.main()


def test_the_facts_are_the_records_and_a_chain_is_named_where_no_peptide_is(monkeypatch):
    found = []
    monkeypatch.setattr(check_assets, "problems", found)
    target = BY_SLUG["insulin"]
    record = {"protein": {"translation": "M" * 110}, "exons": [{}, {}, {}],
              "peptides": [{}, {}, {}]}
    row = {"facts": {"residues": 110, "exons": 3, "chains": 3, "bridges": 3}, "chain": None}
    check_assets.check_facts("insulin", target, record, row)
    assert found == []

    check_assets.check_facts("insulin", target, record,
                             {**row, "facts": {**row["facts"], "chains": 2}})
    assert len(found) == 1 and "the record gives" in found[0]

    found.clear()
    uncut = {**record, "peptides": []}
    check_assets.check_facts("insulin", target, uncut,
                             {"facts": {**row["facts"], "chains": 1}, "chain": None})
    assert len(found) == 1 and "names a chain None where the record names 0" in found[0]
