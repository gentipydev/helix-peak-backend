"""Assemblies: apart from the twenty, built from their entries, one frame for two states.

Offline: the entries are in `structure/structures/`. The baked morph is held
to its check where `bake_assembly.py` has written it, and skipped where not.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import seed_catalog  # noqa: E402
from pipeline.assemblies import check_assembly, seed_assemblies  # noqa: E402
from pipeline.assemblies.assemblies import ASSEMBLIES, BY_SLUG  # noqa: E402
from pipeline.assemblies.bake_assembly import morph_asset  # noqa: E402
from pipeline.assemblies.build import (  # noqa: E402
    build,
    chains_of,
    kabsch,
    operators,
    state_lines,
)
from pipeline.paths import DATA  # noqa: E402
from pipeline.structure.pdb import STRUCTURES  # noqa: E402
from pipeline.targets import TARGETS  # noqa: E402

HEMOGLOBIN = BY_SLUG["hemoglobin-a"]


def test_an_assembly_is_never_one_of_the_twenty():
    assert len(TARGETS) == 20
    assert sorted(seed_catalog.curated_rows()) == sorted(t.slug for t in TARGETS)
    for assembly in ASSEMBLIES:
        assert assembly.slug not in {t.slug for t in TARGETS}


def test_the_seeder_writes_only_an_assemblys_own_tables():
    for statement in (seed_assemblies._UPSERT, seed_assemblies._INSERT_TRACK):
        flat = " ".join(statement.split())
        assert flat.startswith(("insert into assembly ", "insert into assembly_track "))
        assert "protein" not in flat


def test_the_oxy_entry_holds_half_and_its_assembly_is_the_tetramer():
    entry = STRUCTURES / "2DN1.pdb"
    ops = operators(entry, 1)
    assert [op.chains for op in ops] == [("A", "B"), ("A", "B")]
    assert np.allclose(ops[0].matrix[:, :3], np.eye(3))
    relaxed = next(s for s in HEMOGLOBIN.states if s.name == "relaxed")
    assert chains_of(state_lines(relaxed, entry)) == ["A", "B", "C", "D"]


def test_the_deoxy_entry_already_holds_four_chains():
    tense = next(s for s in HEMOGLOBIN.states if s.name == "tense")
    assert chains_of(state_lines(tense, STRUCTURES / "2DN2.pdb")) == ["A", "B", "C", "D"]


def test_the_two_entries_do_not_share_a_frame_until_one_is_moved():
    report = build(HEMOGLOBIN).report
    # As deposited, the two crystals' coordinates are tens of angstroms apart.
    assert report["raw_rmsd_angstrom"] > 40
    # Held on alpha1-beta1, that dimer agrees to under an angstrom, and the
    # other turns by the textbook's ten to twenty degrees.
    assert report["held_rmsd_angstrom"] < 1.5
    assert 10 < report["turned_degrees"] < 20
    assert report["best_fit_rmsd_angstrom"] < report["tetramer_rmsd_angstrom"]


def test_kabsch_finds_a_known_turn():
    rng = np.random.default_rng(3)
    points = rng.normal(size=(40, 3))
    angle = np.radians(33)
    turn = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    moved = points @ turn.T + [4, -2, 7]
    rotation, shift = kabsch(points, moved)
    assert np.allclose(rotation, turn)
    assert np.allclose(shift, [4, -2, 7])


@pytest.mark.skipif(not (DATA / morph_asset(HEMOGLOBIN)).exists(),
                    reason="run bake_assembly.py first")
def test_the_baked_morph_checks_out():
    assert check_assembly.problems_of(HEMOGLOBIN) == []
