"""Reading a stored structure `.glb`, and writing one back.

Offline. Moved with `glb.py` from `structure_ar/`, whose bake wrote them first.
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

from pipeline.structure.glb import Mesh, read_glb, write_glb  # noqa: E402


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


def test_the_room_glb_reads_back_as_written():
    rng = np.random.default_rng(3)
    meshes = {
        "chainA": Mesh(positions=rng.normal(size=(6, 3)), normals=None,
                       triangles=np.array([[0, 1, 2], [3, 4, 5]])),
        "bonds": Mesh(positions=rng.normal(size=(3, 3)),
                      normals=np.tile([0.0, 0.0, 1.0], (3, 1)),
                      triangles=np.array([[0, 1, 2]])),
    }
    payload = write_glb(meshes, scale=0.01)
    back = read_glb(payload)
    assert list(back) == ["chainA", "bonds"]
    for name, mesh in meshes.items():
        assert np.allclose(back[name].positions, mesh.positions * 0.01, atol=1e-7)
        assert (back[name].triangles == mesh.triangles).all()
    assert np.allclose(back["bonds"].normals, [0, 0, 1])
    assert write_glb(meshes, scale=0.01) == payload
