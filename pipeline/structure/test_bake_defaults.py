"""The structure bake's parameters for exports outside the twenty, and their defaults.

`write_pml` takes a `source` and an `origin`, and `verify_frame.py` a PDB file,
an output directory, an origin and chains, so that an assembly built from an
entry can be exported and audited. Unset, each is what the twenty are baked
and audited with; these pin that. The `write_pml` tests need the structure
bake's environment (trimesh) and skip without it.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.structure import verify_frame  # noqa: E402
from pipeline.structure.pdb import pdb_path  # noqa: E402
from pipeline.targets import BY_SLUG  # noqa: E402


def test_verify_frame_audits_insulin_unless_told_otherwise():
    pdb, out, origin, chains = verify_frame.arguments([])
    assert (pdb, out) == ("structures/3I40.pdb", "output/insulin")
    assert np.allclose(origin, [-18.0825, -1.3525, -9.326])
    assert chains == (("A", "chainA.obj"), ("B", "chainB.obj"))


def test_verify_frame_takes_another_export():
    pdb, out, origin, chains = verify_frame.arguments([
        "--pdb", "x.pdb", "--out", "o", "--origin", "1", "2", "3",
        "--chain", "A=alpha1.obj", "--chain", "B=beta1.obj",
    ])
    assert (pdb, out) == ("x.pdb", "o")
    assert np.allclose(origin, [1, 2, 3])
    assert chains == (("A", "alpha1.obj"), ("B", "beta1.obj"))


def test_the_default_script_centres_each_export_on_its_own_entry(tmp_path):
    pytest.importorskip("trimesh")
    from pipeline.structure import bake

    target = BY_SLUG["insulin"]
    script = bake.write_pml(target, tmp_path).read_text().splitlines()
    assert script[0] == f"load {pdb_path(target.structure)}, src"
    assert "ext = cmd.get_extent('chainA or chainB')" in script
    assert "cx = (ext[0][0] + ext[1][0]) / 2.0" in script


def test_a_source_and_an_origin_move_only_what_they_name(tmp_path):
    pytest.importorskip("trimesh")
    from pipeline.structure import bake

    target = BY_SLUG["insulin"]
    plain = bake.write_pml(target, tmp_path).read_text().splitlines()
    built = tmp_path / "built.pdb"
    given = bake.write_pml(
        target, tmp_path, source=built, origin=np.array([1.5, -2.0, 0.25]),
    ).read_text().splitlines()
    assert given[0] == f"load {built}, src"
    assert "cx, cy, cz = 1.5, -2.0, 0.25" in given
    assert not any(line.startswith("ext = ") for line in given)
    # Everything else, the view it pins and the files it saves, is the same.
    kept = [line for line in plain
            if not line.startswith(("load ", "ext = ", "cx = ", "cy = ", "cz = "))]
    assert [line for line in given
            if not line.startswith(("load ", "cx, cy, cz = "))] == kept
