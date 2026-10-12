"""The mmCIF reader: the STAR syntax the PDB writes, and whole entries.

The syntax is held on small texts written here. Whole entries are the files
under `fixtures/pdb/` (PDBe's SIFTS-enriched mmCIF, cut down), and each is
held, item for item, to what Biopython's `MMCIF2Dict` reads from it, where
Biopython is installed (the service's `.venv`).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.structure import mmcif  # noqa: E402
from pipeline.structure.mmcif import CifError  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "pdb"
ENTRIES = ("6rlx_cut.cif", "2k6t_cut.cif", "1ppb_cut.cif")


def test_single_items_and_a_loop_are_categories():
    block = mmcif.parse(
        "data_TEST\n"
        "_entry.id TEST\n"
        "_cell.length_a 10.5\n"
        "_cell.length_b 20.5\n"
        "loop_\n"
        "_atom.id\n"
        "_atom.name\n"
        "1 CA\n"
        "2 CB 3 SG\n")
    assert block["entry"].value("id") == "TEST"
    assert block["cell"].items == ["length_a", "length_b"]
    assert block["cell"].value("length_b") == "20.5"
    assert len(block["atom"]) == 3
    assert list(block["atom"].rows("id", "name")) == [("1", "CA"), ("2", "CB"), ("3", "SG")]


def test_a_missing_item_reads_as_nothing_for_every_row():
    block = mmcif.parse("data_T\nloop_\n_a.x\n1\n2\n")
    assert block["a"].get("y") == [None, None]
    assert "y" not in block["a"] and "x" in block["a"]


def test_bare_dot_and_question_mark_are_nothing_and_quoted_ones_are_themselves():
    block = mmcif.parse("data_T\nloop_\n_a.x\n. ? '.' \"?\" 'x'\n")
    assert block["a"].get("x") == [None, None, ".", "?", "x"]


def test_a_quote_inside_a_value_ends_it_only_before_white_space():
    block = mmcif.parse("data_T\nloop_\n_a.x\n'the 5'-end' \"O5'\" 'a \"b\" c'\n")
    assert block["a"].get("x") == ["the 5'-end", "O5'", 'a "b" c']


def test_a_text_field_runs_between_semicolons_at_line_starts():
    block = mmcif.parse(
        "data_T\n"
        "_a.seq\n"
        ";MKV\n"
        "LLA; not the end\n"
        ";\n"
        "_a.after yes\n")
    assert block["a"].value("seq") == "MKV\nLLA; not the end"
    assert block["a"].value("after") == "yes"


def test_a_text_field_can_be_a_loops_value_and_more_can_follow_its_close():
    block = mmcif.parse("data_T\nloop_\n_a.id\n_a.text\n1\n;one\n; 2 two\n")
    assert list(block["a"].rows("id", "text")) == [("1", "one"), ("2", "two")]


def test_a_comment_runs_to_the_end_of_its_line_and_a_quoted_hash_is_a_value():
    block = mmcif.parse("data_T\n# a comment\n_a.x 1 # and another\n_a.y '#2'\n#\n")
    assert block["a"].value("x") == "1" and block["a"].value("y") == "#2"


def test_only_the_first_data_block_is_read():
    block = mmcif.parse("data_ONE\n_a.x 1\ndata_TWO\n_a.x 2\n_b.y 3\n")
    assert block["a"].value("x") == "1" and "b" not in block


@pytest.mark.parametrize("text, said", [
    ("_a.x 1\n", "before the first data block"),
    ("", "no data block"),
    ("data_T\nloop_\n_a.x\n_a.y\n1 2 3\n", "part-way through a row"),
    ("data_T\n_a.x\n;never closed\n", "never closed"),
    ("data_T\n_a.x 1\n_a.x 2\n", "written twice"),
    ("data_T\n_a.x\n_a.y 2\n", "has no value"),
    ("data_T\nloop_\n_a.x\n_b.y\n1 2\n", "several categories"),
    ("data_T\nstray\n", "unexpected token"),
])
def test_what_is_not_mmcif_is_refused(text, said):
    with pytest.raises(CifError, match=said):
        mmcif.parse(text)


@pytest.mark.parametrize("name", ENTRIES)
def test_a_whole_entry_is_read(name):
    block = mmcif.parse((FIXTURES / name).read_text())
    assert block["entry"].value("id") == name[:4].upper()
    atoms = block["atom_site"]
    assert len(atoms) > 500
    # Every row whole: each column as long as the table.
    assert {len(atoms.columns[item]) for item in atoms.items} == {len(atoms)}
    assert set(atoms.get("group_PDB")) <= {"ATOM", "HETATM"}
    assert all(float(x) == float(x) for x in atoms.get("Cartn_x"))


@pytest.mark.parametrize("name", ENTRIES)
def test_a_whole_entry_reads_as_biopython_reads_it(name):
    MMCIF2Dict = pytest.importorskip("Bio.PDB.MMCIF2Dict").MMCIF2Dict
    path = FIXTURES / name
    ours = {f"_{category}.{item}": table.columns[item]
            for category, table in mmcif.parse(path.read_text()).items()
            for item in table.items}
    theirs = {key: values for key, values in MMCIF2Dict(str(path)).items() if key != "data_"}
    assert sorted(ours) == sorted(theirs)

    def lines(value):
        # Biopython trims a text field's lines; this reader keeps them as written.
        return None if value is None else [line.rstrip() for line in value.strip().split("\n")]

    for key, values in theirs.items():
        # Biopython cannot tell a quoted `?` from a bare one, and keeps both.
        assert [lines(v) for v in ours[key]] == [
            None if v in ("?", ".") else lines(v) for v in values], key
