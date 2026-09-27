"""Read the meshes out of a stored structure `.glb`, node by node.

The structure bake writes one mesh per node (`chainA`, `chainB`, `bonds`)
through trimesh. This reads exactly that shape back -- positions, normals and
triangle indices, with each node's transform applied -- and nothing more of
glTF: no materials, no textures, no animation. It needs only numpy, so the AR
bake runs without the structure bake's PyMOL and trimesh.
"""

from __future__ import annotations

import json
import struct
from dataclasses import dataclass

import numpy as np

_JSON = 0x4E4F534A
_BIN = 0x004E4942
_COMPONENTS = {5121: np.uint8, 5123: np.uint16, 5125: np.uint32, 5126: np.float32}
_WIDTH = {"SCALAR": 1, "VEC2": 2, "VEC3": 3, "VEC4": 4, "MAT4": 16}


@dataclass(frozen=True)
class Mesh:
    """One node's triangles, in the model's own frame."""

    positions: np.ndarray        # (n, 3) float64
    normals: np.ndarray | None   # (n, 3) float64, where the file has them
    triangles: np.ndarray        # (m, 3) int64


def read_glb(payload: bytes) -> dict[str, Mesh]:
    """Every named node that carries a mesh, by node name, in file order."""
    magic, version, length = struct.unpack_from("<4sII", payload, 0)
    if magic != b"glTF" or version != 2:
        raise ValueError("not a glTF 2 binary")
    if length != len(payload):
        raise ValueError(f"header says {length} bytes, file is {len(payload)}")
    document, binary, offset = None, b"", 12
    while offset < length:
        size, kind = struct.unpack_from("<II", payload, offset)
        chunk = payload[offset + 8:offset + 8 + size]
        if kind == _JSON:
            document = json.loads(chunk)
        elif kind == _BIN:
            binary = chunk
        offset += 8 + size
    if document is None:
        raise ValueError("no JSON chunk")

    def accessor(index: int) -> np.ndarray:
        spec = document["accessors"][index]
        view = document["bufferViews"][spec["bufferView"]]
        dtype = np.dtype(_COMPONENTS[spec["componentType"]])
        width = _WIDTH[spec["type"]]
        start = view.get("byteOffset", 0) + spec.get("byteOffset", 0)
        stride = view.get("byteStride") or dtype.itemsize * width
        count = spec["count"]
        rows = np.ndarray(
            shape=(count, width), dtype=dtype, buffer=binary,
            offset=start, strides=(stride, dtype.itemsize),
        )
        return np.array(rows)

    def local(node: dict) -> np.ndarray:
        if "matrix" in node:
            return np.array(node["matrix"], dtype=float).reshape(4, 4).T
        matrix = np.eye(4)
        x, y, z, w = node.get("rotation", (0.0, 0.0, 0.0, 1.0))
        rotation = np.array([
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ])
        matrix[:3, :3] = rotation * np.array(node.get("scale", (1.0, 1.0, 1.0)))
        matrix[:3, 3] = node.get("translation", (0.0, 0.0, 0.0))
        return matrix

    nodes = document.get("nodes", [])
    scene = document.get("scenes", [{}])[document.get("scene", 0)]
    found: dict[str, Mesh] = {}

    def visit(index: int, parent: np.ndarray) -> None:
        node = nodes[index]
        world = parent @ local(node)
        if "mesh" in node:
            name = node.get("name") or f"node{index}"
            positions, normals, triangles = [], [], []
            base = 0
            for primitive in document["meshes"][node["mesh"]]["primitives"]:
                if primitive.get("mode", 4) != 4:
                    raise ValueError(f"{name}: only triangles are read")
                points = accessor(primitive["attributes"]["POSITION"]).astype(float)
                positions.append(points @ world[:3, :3].T + world[:3, 3])
                if "NORMAL" in primitive["attributes"]:
                    normal = accessor(primitive["attributes"]["NORMAL"]).astype(float)
                    normal = normal @ np.linalg.inv(world[:3, :3])
                    normal /= np.linalg.norm(normal, axis=1, keepdims=True).clip(1e-12)
                    normals.append(normal)
                indices = (accessor(primitive["indices"]).reshape(-1, 3).astype(np.int64)
                           if "indices" in primitive
                           else np.arange(len(points)).reshape(-1, 3))
                triangles.append(indices + base)
                base += len(points)
            if name in found:
                raise ValueError(f"two nodes are named {name!r}")
            found[name] = Mesh(
                positions=np.vstack(positions),
                normals=np.vstack(normals) if len(normals) == len(positions) else None,
                triangles=np.vstack(triangles),
            )
        for child in node.get("children", []):
            visit(child, world)

    for root in scene.get("nodes", range(len(nodes))):
        visit(root, np.eye(4))
    return found


def write_glb(meshes: dict[str, Mesh], scale: float = 1.0) -> bytes:
    """The meshes as a glTF binary, one named node each, every point times `scale`.

    Positions as float32, normals where there are any, indices as uint32; one
    buffer, no materials. The bytes depend on nothing but the meshes, so the
    same meshes always give the same file.
    """
    binary = bytearray()
    views, accessors, gltf_meshes, nodes = [], [], [], []

    def add(array: np.ndarray, kind: str, component: int, bounds: bool) -> int:
        data = np.ascontiguousarray(array)
        views.append({"buffer": 0, "byteOffset": len(binary), "byteLength": data.nbytes})
        accessor = {"bufferView": len(views) - 1, "componentType": component,
                    "count": len(data), "type": kind}
        if bounds:
            accessor["min"] = [float(v) for v in data.min(axis=0)]
            accessor["max"] = [float(v) for v in data.max(axis=0)]
        accessors.append(accessor)
        binary.extend(data.tobytes())
        binary.extend(b"\0" * (-len(binary) % 4))
        return len(accessors) - 1

    for name, mesh in meshes.items():
        attributes = {"POSITION": add((mesh.positions * scale).astype(np.float32),
                                      "VEC3", 5126, True)}
        if mesh.normals is not None:
            attributes["NORMAL"] = add(mesh.normals.astype(np.float32), "VEC3", 5126, False)
        indices = add(mesh.triangles.reshape(-1).astype(np.uint32), "SCALAR", 5125, False)
        gltf_meshes.append({"name": name, "primitives": [
            {"attributes": attributes, "indices": indices, "mode": 4}]})
        nodes.append({"name": name, "mesh": len(gltf_meshes) - 1})

    document = json.dumps({
        "asset": {"version": "2.0", "generator": "pipeline/structure_ar"},
        "scene": 0, "scenes": [{"nodes": list(range(len(nodes)))}],
        "nodes": nodes, "meshes": gltf_meshes, "accessors": accessors,
        "bufferViews": views, "buffers": [{"byteLength": len(binary)}],
    }, separators=(",", ":"), sort_keys=True).encode()
    document += b" " * (-len(document) % 4)
    body = (struct.pack("<II", len(document), _JSON) + document
            + struct.pack("<II", len(binary), _BIN) + bytes(binary))
    return struct.pack("<4sII", b"glTF", 2, 12 + len(body)) + body
