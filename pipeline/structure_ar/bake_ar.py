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
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.paths import DATA  # noqa: E402
from pipeline.structure.frame import model_frame  # noqa: E402
from pipeline.structure.glb import Mesh, read_glb, write_glb  # noqa: E402
from pipeline.targets import BY_SLUG, TARGETS, Target  # noqa: E402

# Every protein with a fold gets an AR fold: recorded here, never in targets.py.
AR_TARGETS = tuple(t for t in TARGETS if t.structure is not None)

# One unit is a centimetre, and a point is written in angstroms.
METERS_PER_UNIT = 0.01

# The walk's material on the fold page: matte, so the form reads from shading.
ROUGHNESS = 0.65


def ar_asset(target: Target) -> str:
    return f"assets/models_ar/{target.slug}.usdz"


def ar_glb(target: Target) -> str:
    """The same fold as a `.glb` in metres, for Android's Scene Viewer."""
    return f"assets/models_ar/{target.slug}.glb"


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
    fit, method = model_frame(structure, meshes)
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
