"""Bake `structure_ar`: each fold as a USDZ, at its real size, for AR.

    .venv/Scripts/python pipeline/structure_ar/bake_ar.py --all
    .venv/Scripts/python pipeline/structure_ar/bake_ar.py --target insulin

iOS AR Quick Look reads USDZ and nothing else. The walk's own `structure`
track is a `.glb` normalised so its longest axis is 1.0 -- the fold page needs
no scale constant that way -- and that normalisation throws the molecule's
real size away. This track is a second, separate object: the same meshes, back
at their size in angstroms, packaged as USDZ. The `structure` row, its bytes
and its provenance are not touched.

**Where the meshes come from.** The stored `.glb`, fetched and digest-checked
by `fetch_tracks.py --kind structure`, and never a new PyMOL run: the walk's
fold and the reader's AR fold are then the very same triangles, and this bake
runs anywhere numpy and usd-core do. The source `.glb`'s sha256 is recorded.

**How big it really is.** `structure/bake.py` exports PyMOL's cartoon as a pure
translation of the PDB frame -- no rotation -- then centres it and divides by
its longest axis, L. So each exported CA atom `p` sits on the ribbon at
`(p - c) / L` in the stored model, for one centre `c` and one length `L`. Both
are recovered by fitting the entry's own CA atoms (the chains and span the
bake exported) to the model's vertices: nearest vertex, then the least-squares
scale and shift, repeated until it settles. The real bounding box is then the
model's own extent times L. The fit's residual goes into the provenance; for
insulin it reproduces the bake's own logged 24.31 x 18.66 x 20.33 A.

**Scale.** One unit is one centimetre (`metersPerUnit` 0.01) and points are in
angstroms, so AR Quick Look shows 1 A as 1 cm: insulin stands 24 cm across.
The fold rests on its lowest point, centred over the anchor, because AR Quick
Look puts the model's origin on the floor it finds.

**Android.** Scene Viewer takes a `.glb`, in metres. The same placed meshes
are written beside the USDZ as `<slug>.glb`, every point times 0.01, so 1 A is
1 cm there too. It is uploaded beside the USDZ and named in the row's
provenance, as the walk's `.glb` is beside its `.fsceneb`.

**Colour.** None is baked. Which colour a chain takes is the app's decision
(`targets.Chain`: "a colour written down twice is a colour that will
disagree"), so every chain gets the walk's matte material -- roughness 0.65,
no metal -- at USD's default neutral base colour, and each chain stays its own
named mesh, as it is in the `.glb`.

Which proteins get the track: every target with a structure (`AR_TARGETS`).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.paths import DATA  # noqa: E402
from pipeline.structure_ar.glb import Mesh, read_glb, write_glb  # noqa: E402
from pipeline.targets import BY_SLUG, TARGETS, Structure, Target  # noqa: E402

STRUCTURES = BACKEND / "pipeline" / "structure" / "structures"

# Every protein with a fold gets an AR fold: recorded here, never in targets.py.
AR_TARGETS = tuple(t for t in TARGETS if t.structure is not None)

# One unit is a centimetre, and a point is written in angstroms.
METERS_PER_UNIT = 0.01

# The walk's material on the fold page: matte, so the form reads from shading.
ROUGHNESS = 0.65

# The fit stops when a round moves the length by less than this fraction.
_SETTLED = 1e-7
_ROUNDS = 60


def ar_asset(target: Target) -> str:
    return f"assets/models_ar/{target.slug}.usdz"


def ar_glb(target: Target) -> str:
    """The same fold as a `.glb` in metres, for Android's Scene Viewer."""
    return f"assets/models_ar/{target.slug}.glb"


def ca_atoms(structure: Structure, pdb: Path | None = None) -> np.ndarray:
    """The CA atoms of what the structure bake exported, in the PDB frame.

    The same selection `structure/bake.py` hands PyMOL: each exported chain,
    polymer only (ATOM records), clipped to the span where there is one. Of
    alternate locations, the first.
    """
    path = pdb or STRUCTURES / f"{structure.pdb}.pdb"
    wanted = {c.pdb_chain for c in structure.chains}
    span = structure.residues
    seen, out = set(), []
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("ENDMDL"):
            break   # the first model of an NMR or multi-model entry
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        chain, resi = line[21], int(line[22:26])
        if chain not in wanted or line[16] not in (" ", "A"):
            continue
        if span is not None and not span[0] <= resi <= span[1]:
            continue
        key = (chain, resi, line[26])
        if key in seen:
            continue
        seen.add(key)
        out.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if len(out) < 4:
        raise ValueError(f"{structure.pdb}: {len(out)} CA atoms in the exported chains")
    return np.array(out)


def _nearest(points: np.ndarray, vertices: np.ndarray) -> np.ndarray:
    """For each point, the nearest vertex. Brute force, in slices."""
    out = np.empty_like(points)
    squared = (vertices * vertices).sum(axis=1)
    step = max(1, 4_000_000 // max(1, len(vertices)))
    for start in range(0, len(points), step):
        chunk = points[start:start + step]
        distance = squared[None, :] - 2 * chunk @ vertices.T
        out[start:start + step] = vertices[distance.argmin(axis=1)]
    return out


@dataclass(frozen=True)
class Fit:
    """`(p - centre) / length` puts a PDB point into the stored model."""

    length: float           # angstroms per model unit: the bake's L
    centre: np.ndarray      # the bake's centre, in the PDB frame
    rms: float              # CA to its nearest vertex, in angstroms
    atoms: int


def fit_size(vertices: np.ndarray, atoms: np.ndarray) -> Fit:
    """Recover the bake's centre and length from the CA atoms."""
    low, high = atoms.min(axis=0), atoms.max(axis=0)
    model = vertices.max(axis=0) - vertices.min(axis=0)
    # A first guess from the two bounding boxes, then nearest-vertex rounds.
    scale = model.max() / (high - low).max()        # model units per angstrom
    shift = -scale * (low + high) / 2
    for _ in range(_ROUNDS):
        mapped = atoms * scale + shift
        target = _nearest(mapped, vertices)
        a_mean, t_mean = atoms.mean(axis=0), target.mean(axis=0)
        a, t = atoms - a_mean, target - t_mean
        new_scale = float((a * t).sum() / (a * a).sum())
        new_shift = t_mean - new_scale * a_mean
        settled = abs(new_scale - scale) / scale < _SETTLED
        scale, shift = new_scale, new_shift
        if settled:
            break
    mapped = atoms * scale + shift
    residual = np.linalg.norm(mapped - _nearest(mapped, vertices), axis=1)
    return Fit(
        length=1.0 / scale,
        centre=-shift / scale,
        rms=float(np.sqrt((residual ** 2).mean()) / scale),
        atoms=len(atoms),
    )


# How `structure/bake.py` builds one disulfide: five rods (a 12-section
# cylinder is 26 vertices) and then four joints (a once-subdivided icosphere
# is 42), centred on CB, SG, SG and CB in that order.
_ROD_VERTICES = 26
_JOINT_VERTICES = 42
_BRIDGE_VERTICES = 5 * _ROD_VERTICES + 4 * _JOINT_VERTICES


def bridge_atoms(structure: Structure, pdb: Path | None = None) -> list[np.ndarray]:
    """Each exported bridge's CB, SG, SG and CB, in the order the bake built them.

    The same reading of the file as `structure/bake.py`'s `ssbonds` and
    `atoms`: SSBOND records with both ends in the exported chains and span,
    in file order, and each atom's last ATOM line. Written again here rather
    than imported, because sharing code with that baker would need its twenty
    models re-baked byte for byte to prove the move changed nothing, and that
    needs PyMOL.
    """
    path = pdb or STRUCTURES / f"{structure.pdb}.pdb"
    text = path.read_text(errors="replace").splitlines()
    wanted = {c.pdb_chain for c in structure.chains}
    span = structure.residues
    at = {}
    for line in text:
        if line.startswith("ATOM") and line[12:16].strip() in ("CB", "SG"):
            at[(line[21], int(line[22:26]), line[12:16].strip())] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    out = []
    for line in text:
        if not line.startswith("SSBOND"):
            continue
        c1, r1, c2, r2 = line[15], int(line[17:21]), line[29], int(line[31:35])
        if c1 not in wanted or c2 not in wanted:
            continue
        if span is not None and not all(span[0] <= r <= span[1] for r in (r1, r2)):
            continue
        out.append(np.array([at[(c1, r1, "CB")], at[(c1, r1, "SG")],
                             at[(c2, r2, "SG")], at[(c2, r2, "CB")]]))
    return out


def fit_bridges(bonds: Mesh, bridges: list[np.ndarray]) -> Fit | None:
    """The bake's centre and length, exactly, from the bridges' joints.

    Each joint is a sphere centred on its atom, so the centroid of its
    vertices in the stored model is that atom's position there. Twelve points
    and more, known in both frames, fix the scale and shift exactly. None where
    the `bonds` mesh is not laid out as the bake lays it out.
    """
    if not bridges or len(bonds.positions) != _BRIDGE_VERTICES * len(bridges):
        return None
    model, pdb = [], []
    for i, atoms in enumerate(bridges):
        block = bonds.positions[i * _BRIDGE_VERTICES:(i + 1) * _BRIDGE_VERTICES]
        for j in range(4):
            start = 5 * _ROD_VERTICES + j * _JOINT_VERTICES
            model.append(block[start:start + _JOINT_VERTICES].mean(axis=0))
            pdb.append(atoms[j])
    model, pdb = np.array(model), np.array(pdb)
    p_mean, m_mean = pdb.mean(axis=0), model.mean(axis=0)
    p, m = pdb - p_mean, model - m_mean
    scale = float((p * m).sum() / (p * p).sum())
    shift = m_mean - scale * p_mean
    residual = np.linalg.norm(pdb * scale + shift - model, axis=1)
    return Fit(
        length=1.0 / scale,
        centre=-shift / scale,
        rms=float(np.sqrt((residual ** 2).mean()) / scale),
        atoms=len(pdb),
    )


def place(meshes: dict[str, Mesh], length: float) -> dict[str, Mesh]:
    """The meshes in angstroms, centred over the anchor, resting on their lowest point.

    Where AR puts a model is where its origin is: on the floor it finds.
    """
    stacked = np.vstack([m.positions for m in meshes.values()]) * length
    low, high = stacked.min(axis=0), stacked.max(axis=0)
    offset = np.array([-(low[0] + high[0]) / 2, -low[1], -(low[2] + high[2]) / 2])
    return {name: Mesh(positions=m.positions * length + offset, normals=m.normals,
                       triangles=m.triangles)
            for name, m in meshes.items()}


def write_usdz(
    meshes: dict[str, Mesh],
    length: float,
    destination: Path,
    metadata: dict,
) -> None:
    """The meshes at `length` angstroms per model unit, as an ARKit USDZ."""
    from pxr import Gf, Sdf, Usd, UsdGeom, UsdShade, UsdUtils, Vt

    placed = place(meshes, length)

    with tempfile.TemporaryDirectory() as work:
        layer = Path(work) / "fold.usdc"
        stage = Usd.Stage.CreateNew(str(layer))
        UsdGeom.SetStageMetersPerUnit(stage, METERS_PER_UNIT)
        UsdGeom.SetStageUpAxis(stage, UsdGeom.Tokens.y)
        root = UsdGeom.Xform.Define(stage, "/Fold")
        stage.SetDefaultPrim(root.GetPrim())
        stage.GetRootLayer().customLayerData = {"helixpeek": json.dumps(metadata, sort_keys=True)}

        material = UsdShade.Material.Define(stage, "/Fold/Looks/Matte")
        shader = UsdShade.Shader.Define(stage, "/Fold/Looks/Matte/Surface")
        shader.CreateIdAttr("UsdPreviewSurface")
        shader.CreateInput("roughness", Sdf.ValueTypeNames.Float).Set(ROUGHNESS)
        shader.CreateInput("metallic", Sdf.ValueTypeNames.Float).Set(0.0)
        material.CreateSurfaceOutput().ConnectToSource(shader.ConnectableAPI(), "surface")

        for name, mesh in meshes.items():
            points = placed[name].positions
            prim = UsdGeom.Mesh.Define(stage, f"/Fold/{name}")
            prim.CreatePointsAttr(Vt.Vec3fArray.FromNumpy(points.astype(np.float32)))
            prim.CreateFaceVertexCountsAttr(
                Vt.IntArray.FromNumpy(np.full(len(mesh.triangles), 3, dtype=np.int32)))
            prim.CreateFaceVertexIndicesAttr(
                Vt.IntArray.FromNumpy(mesh.triangles.reshape(-1).astype(np.int32)))
            if mesh.normals is not None:
                prim.CreateNormalsAttr(Vt.Vec3fArray.FromNumpy(mesh.normals.astype(np.float32)))
                prim.SetNormalsInterpolation(UsdGeom.Tokens.vertex)
            prim.CreateSubdivisionSchemeAttr(UsdGeom.Tokens.none)
            prim.CreateExtentAttr([
                Gf.Vec3f(*points.min(axis=0).astype(float)),
                Gf.Vec3f(*points.max(axis=0).astype(float)),
            ])
            UsdShade.MaterialBindingAPI.Apply(prim.GetPrim()).Bind(material)
        stage.GetRootLayer().Save()

        staged = Path(work) / "fold.usdz"
        if not UsdUtils.CreateNewARKitUsdzPackage(Sdf.AssetPath(str(layer)), str(staged)):
            raise RuntimeError("usd-core did not write the package")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(staged, destination)


def bake(target: Target) -> dict:
    """Write the target's USDZ under `pipeline/data/` and return its metadata."""
    structure = target.structure
    source = DATA / target.structure_asset
    if not source.exists():
        raise FileNotFoundError(
            f"{source}: run `fetch_tracks.py --kind structure` first")
    payload = source.read_bytes()
    meshes = read_glb(payload)
    vertices = np.vstack([m.positions for m in meshes.values()])
    # Exact where the model has bridges to read the scale from; fitted to the
    # ribbon's CA atoms where it has none.
    exact = fit_bridges(meshes["bonds"], bridge_atoms(structure)) if "bonds" in meshes else None
    fit = exact or fit_size(vertices, ca_atoms(structure))
    method = ("disulfide CB and SG atoms at the bridges' joints" if exact
              else "exported CA atoms fitted to the ribbon")
    extent = (vertices.max(axis=0) - vertices.min(axis=0)) * fit.length

    from pxr import Usd

    metadata = {
        "pdb": structure.pdb,
        "nodes": list(meshes),
        "bbox_angstrom": [round(float(v), 2) for v in extent],
        "scale": "1 A = 1 cm",
        "meters_per_unit": METERS_PER_UNIT,
        "size_fit": {
            "method": method,
            "atoms": fit.atoms,
            "rms_angstrom": round(fit.rms, 3),
        },
        "source_glb": {"sha256": hashlib.sha256(payload).hexdigest(), "bytes": len(payload)},
        "usd": f"usd-core {'.'.join(str(v) for v in Usd.GetVersion())}",
        "built_by": "pipeline/structure_ar/bake_ar.py",
    }
    write_usdz(meshes, fit.length, DATA / ar_asset(target), metadata)
    # Scene Viewer reads metres: angstroms times 0.01 shows 1 A as 1 cm there too.
    (DATA / ar_glb(target)).write_bytes(
        write_glb(place(meshes, fit.length), scale=METERS_PER_UNIT))
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if not args.all and not args.target:
        parser.error("name --target or --all")
    chosen = AR_TARGETS if args.all else tuple(BY_SLUG[s] for s in args.target)
    for target in chosen:
        if target.structure is None:
            print(f"{target.slug}: no structure, no AR fold", file=sys.stderr)
            continue
        metadata = bake(target)
        path = DATA / ar_asset(target)
        box = metadata["bbox_angstrom"]
        print(f"{target.slug:<16} {box[0]:7.2f} x {box[1]:7.2f} x {box[2]:7.2f} A  "
              f"fit {metadata['size_fit']['rms_angstrom']:.2f} A over "
              f"{metadata['size_fit']['atoms']} atoms  {path.stat().st_size:>9,} B")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
