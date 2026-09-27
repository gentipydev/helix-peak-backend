"""The AR bake: writing the USDZ, and where the track lives.

Offline. Reading the stored model and finding its size are the shared
`structure/glb.py` and `structure/frame.py`, tested beside them.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import fetch_tracks  # noqa: E402
from pipeline.structure_ar import bake_ar  # noqa: E402
from pipeline.structure.glb import Mesh  # noqa: E402
from pipeline.targets import BY_SLUG  # noqa: E402


def test_the_usdz_is_arkit_shaped_and_says_how_big_it_is(tmp_path):
    pytest.importorskip("pxr")
    from pxr import Usd, UsdGeom

    square = np.array([[-0.5, -0.4, 0], [0.5, -0.4, 0], [0.5, 0.4, 0.2], [-0.5, 0.4, 0.2]])
    mesh = Mesh(positions=square, normals=None, triangles=np.array([[0, 1, 2], [0, 2, 3]]))
    meta = {"pdb": "TEST", "nodes": ["chainA"], "bbox_angstrom": [20.0, 16.0, 4.0]}
    destination = tmp_path / "fold.usdz"
    bake_ar.write_usdz({"chainA": mesh}, 20.0, destination, meta)

    stage = Usd.Stage.Open(str(destination))
    assert UsdGeom.GetStageMetersPerUnit(stage) == pytest.approx(0.01)
    assert UsdGeom.GetStageUpAxis(stage) == UsdGeom.Tokens.y
    assert stage.GetDefaultPrim().GetPath() == "/Fold"
    assert json.loads(stage.GetRootLayer().customLayerData["helixpeek"]) == meta
    points = np.array(UsdGeom.Mesh(stage.GetPrimAtPath("/Fold/chainA")).GetPointsAttr().Get())
    assert np.allclose(points.max(axis=0) - points.min(axis=0), [20, 16, 4])
    assert points[:, 1].min() == pytest.approx(0)   # resting on the floor
    assert points[:, 0].min() == pytest.approx(-points[:, 0].max())

    again = tmp_path / "again.usdz"
    bake_ar.write_usdz({"chainA": mesh}, 20.0, again, meta)
    assert again.read_bytes() == destination.read_bytes()


def test_the_fetcher_and_the_uploader_agree_on_where_it_lives():
    insulin = BY_SLUG["insulin"]
    assert fetch_tracks.asset_path("structure_ar", insulin) == bake_ar.ar_asset(insulin)
    assert "structure_ar" in fetch_tracks.KINDS
