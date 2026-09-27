"""The AR bake: reading the stored model, finding its size, writing the USDZ.

Offline. The tests that read the stored `.glb`s skip until
`fetch_tracks.py --kind structure` has put them in `pipeline/data/`.
"""

from __future__ import annotations

import json
import struct
import sys
from pathlib import Path

import numpy as np
import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import fetch_tracks  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.structure_ar import bake_ar  # noqa: E402
from pipeline.structure_ar.glb import Mesh, read_glb  # noqa: E402
from pipeline.targets import BY_SLUG  # noqa: E402


def _glb(nodes: list[dict], meshes: list[tuple[np.ndarray, np.ndarray]]) -> bytes:
    """A minimal glTF binary: float32 positions and uint32 indices per mesh."""
    binary, views, accessors, gltf_meshes = b"", [], [], []
    for positions, triangles in meshes:
        for array, kind, component in ((positions.astype(np.float32), "VEC3", 5126),
                                       (triangles.astype(np.uint32).reshape(-1), "SCALAR", 5125)):
            views.append({"buffer": 0, "byteOffset": len(binary), "byteLength": array.nbytes})
            accessors.append({"bufferView": len(views) - 1, "componentType": component,
                              "count": len(array), "type": kind})
            binary += array.tobytes()
        gltf_meshes.append({"primitives": [{
            "attributes": {"POSITION": len(accessors) - 2}, "indices": len(accessors) - 1}]})
    document = json.dumps({
        "asset": {"version": "2.0"}, "scene": 0,
        "scenes": [{"nodes": list(range(len(nodes)))}], "nodes": nodes,
        "meshes": gltf_meshes, "accessors": accessors, "bufferViews": views,
        "buffers": [{"byteLength": len(binary)}],
    }).encode()
    document += b" " * (-len(document) % 4)
    binary += b"\0" * (-len(binary) % 4)
    body = (struct.pack("<II", len(document), 0x4E4F534A) + document
            + struct.pack("<II", len(binary), 0x004E4942) + binary)
    return struct.pack("<4sII", b"glTF", 2, 12 + len(body)) + body


def test_reads_each_named_node_with_its_transform():
    square = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float)
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    payload = _glb(
        [{"name": "chainA", "mesh": 0},
         {"name": "bonds", "mesh": 1, "translation": [0, 0, 2]}],
        [(square, faces), (square, faces)],
    )
    meshes = read_glb(payload)
    assert list(meshes) == ["chainA", "bonds"]
    assert np.allclose(meshes["chainA"].positions, square)
    assert np.allclose(meshes["bonds"].positions, square + [0, 0, 2])
    assert (meshes["bonds"].triangles == faces).all()
    assert meshes["chainA"].normals is None


def test_refuses_what_is_not_a_glb():
    with pytest.raises(ValueError):
        read_glb(b"PK\x03\x04" + b"\0" * 20)


def test_the_bridges_give_the_scale_exactly():
    rng = np.random.default_rng(7)
    length, centre = 31.7, np.array([4.0, -12.5, 30.25])
    bridges = [rng.uniform(-15, 15, size=(4, 3)) + centre for _ in range(3)]
    joint = rng.normal(size=(bake_ar._JOINT_VERTICES, 3))
    joint -= joint.mean(axis=0)
    blocks = []
    for atoms in bridges:
        model = (atoms - centre) / length
        rods = rng.normal(size=(5 * bake_ar._ROD_VERTICES, 3))
        blocks += [rods] + [m + joint * 0.01 for m in model]
    bonds = Mesh(positions=np.vstack(blocks), normals=None,
                 triangles=np.zeros((0, 3), int))
    fit = bake_ar.fit_bridges(bonds, bridges)
    assert fit.length == pytest.approx(length, rel=1e-9)
    assert np.allclose(fit.centre, centre)
    assert fit.rms < 1e-9
    # A mesh not laid out the bake's way is not read as one.
    assert bake_ar.fit_bridges(bonds, bridges[:2]) is None


def test_the_ca_fit_is_exact_where_the_ribbon_runs_through_the_atoms():
    # A helical CA trace and a model sampled along it, as a ribbon's spine is.
    # On a real cartoon the ribbon's width moves the answer a little: measured
    # on the thirteen stored models whose bridges give the size exactly, the
    # CA fit lands within 1.2% (see README.md, and the insulin test below).
    t = np.linspace(0, 6 * np.pi, 120)
    trace = np.stack([2.3 * np.cos(t), 1.5 * t, 2.3 * np.sin(t)], axis=1)
    spine = np.vstack([trace[:-1] + (trace[1:] - trace[:-1]) * f
                       for f in (0, 0.25, 0.5, 0.75)])
    length, centre = 30.0, np.array([3.0, 14.0, -2.0])
    fit = bake_ar.fit_size((spine - centre) / length, trace)
    assert fit.length == pytest.approx(length, rel=1e-3)
    assert np.allclose(fit.centre, centre, atol=0.05)
    assert fit.rms < 0.05


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


_STORED = DATA / BY_SLUG["insulin"].structure_asset


@pytest.mark.skipif(not _STORED.exists(), reason="run fetch_tracks.py --kind structure first")
def test_insulin_comes_back_at_the_size_its_bake_logged():
    # `structure/README.md`: "Bounding box before normalising: 24.31 x 18.66 x 20.33 A."
    structure = BY_SLUG["insulin"].structure
    meshes = read_glb(_STORED.read_bytes())
    vertices = np.vstack([m.positions for m in meshes.values()])
    fit = bake_ar.fit_bridges(meshes["bonds"], bake_ar.bridge_atoms(structure))
    extent = (vertices.max(axis=0) - vertices.min(axis=0)) * fit.length
    assert np.allclose(extent, [24.31, 18.66, 20.33], atol=0.005)
    assert fit.rms < 0.01
    # The ribbon fit, on its own, lands within about one per cent of it.
    ribbon = bake_ar.fit_size(vertices, bake_ar.ca_atoms(structure))
    assert ribbon.length == pytest.approx(fit.length, rel=0.015)
