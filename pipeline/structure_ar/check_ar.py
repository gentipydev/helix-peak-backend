"""Check every baked `structure_ar` USDZ against the model it was made from.

    .venv/Scripts/python pipeline/structure_ar/check_ar.py

For each AR target, the USDZ under `pipeline/data/assets/models_ar/`:

- passes every validator usd-core registers (the USDZ package and root-layer
  rules AR Quick Look needs, stage metadata, material bindings);
- is in centimetres (`metersPerUnit` 0.01), Y up, with `/Fold` its default
  prim, and rests on its lowest point, centred over the anchor;
- has a `.glb` beside it holding the same points in metres, for Scene Viewer;
- holds the same meshes as the stored `.glb` -- the same node names in the same
  order, the same vertex and triangle counts -- scaled by one factor, the same
  on every axis, and names that `.glb` by the sha256 of the file beside it;
- says the bounding box its own points have, in angstroms, and says how that
  size was found, with a residual that method can be trusted at.

Exits 1 with the problems listed. Reads `pipeline/data/` only.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.paths import DATA  # noqa: E402
from pipeline.structure_ar.bake_ar import AR_TARGETS, METERS_PER_UNIT, ar_asset, ar_glb  # noqa: E402
from pipeline.structure_ar.glb import read_glb  # noqa: E402

# How far a residual may run before the size it gave is not believed: the
# bridges fix the scale exactly, and a CA fit to a cartoon lands within about
# 1.2% of it (measured on the thirteen models that have both).
_EXACT_RMS = 0.01
_FITTED_RMS = 1.0


def metadata(stage) -> dict:
    return json.loads(stage.GetRootLayer().customLayerData["helixpeek"])


def problems_of(target) -> list[str]:
    from pxr import Usd, UsdGeom, UsdValidation

    out: list[str] = []
    path = DATA / ar_asset(target)
    source = DATA / target.structure_asset
    if not path.exists():
        return [f"{path}: not baked"]
    if not source.exists():
        return [f"{source}: fetch the structure track first"]

    stage = Usd.Stage.Open(str(path))
    validators = UsdValidation.ValidationRegistry().GetOrLoadAllValidators()
    for error in UsdValidation.ValidationContext(validators).Validate(stage):
        if error.GetType() == UsdValidation.ValidationErrorType.Error:
            out.append(f"validator: {error.GetMessage()}")

    meta = metadata(stage)
    if abs(UsdGeom.GetStageMetersPerUnit(stage) - METERS_PER_UNIT) > 1e-12:
        out.append(f"metersPerUnit is {UsdGeom.GetStageMetersPerUnit(stage)}")
    if UsdGeom.GetStageUpAxis(stage) != UsdGeom.Tokens.y:
        out.append(f"up axis is {UsdGeom.GetStageUpAxis(stage)}")
    root = stage.GetDefaultPrim()
    if not root or root.GetPath() != "/Fold":
        out.append("the default prim is not /Fold")
        return out

    glb_bytes = source.read_bytes()
    if meta["source_glb"]["sha256"] != hashlib.sha256(glb_bytes).hexdigest():
        out.append("made from another .glb than the stored one")
    model = read_glb(glb_bytes)
    meshes = [p for p in root.GetChildren() if p.IsA(UsdGeom.Mesh)]
    names = [p.GetName() for p in meshes]
    if names != list(model) or names != meta["nodes"]:
        out.append(f"nodes {names}, the .glb has {list(model)}")
        return out

    points, scaled = [], []
    for prim in meshes:
        mesh = UsdGeom.Mesh(prim)
        glb = model[prim.GetName()]
        these = np.array(mesh.GetPointsAttr().Get(), dtype=float)
        indices = np.array(mesh.GetFaceVertexIndicesAttr().Get()).reshape(-1, 3)
        if these.shape != glb.positions.shape or not np.array_equal(indices, glb.triangles):
            out.append(f"{prim.GetName()}: not the .glb's mesh")
        points.append(these)
        scaled.append(glb.positions)
    points, scaled = np.vstack(points), np.vstack(scaled)

    low, high = points.min(axis=0), points.max(axis=0)
    extent = high - low
    if not np.allclose(extent, meta["bbox_angstrom"], atol=0.006):
        out.append(f"points span {np.round(extent, 2)}, it says {meta['bbox_angstrom']}")
    if abs(low[1]) > 1e-3 or abs(low[0] + high[0]) > 1e-3 or abs(low[2] + high[2]) > 1e-3:
        out.append("does not rest centred on its lowest point")
    ratio = extent / (scaled.max(axis=0) - scaled.min(axis=0))
    if np.ptp(ratio) > 1e-4 * ratio.mean():
        out.append(f"scaled unevenly: {ratio}")

    room = DATA / ar_glb(target)
    if not room.exists():
        out.append(f"{room}: no .glb beside it for Scene Viewer")
    else:
        companion = read_glb(room.read_bytes())
        if list(companion) != names:
            out.append(f"the room .glb has nodes {list(companion)}")
        else:
            metres = np.vstack([companion[n].positions for n in names])
            if not np.allclose(metres, points * METERS_PER_UNIT, atol=1e-5):
                out.append("the room .glb is not the USDZ's points in metres")

    fit = meta["size_fit"]
    limit = _EXACT_RMS if fit["method"].startswith("disulfide") else _FITTED_RMS
    if fit["rms_angstrom"] > limit:
        out.append(f"size fit residual {fit['rms_angstrom']} A over {limit}")
    return out


def main() -> int:
    failed = 0
    for target in AR_TARGETS:
        found = problems_of(target)
        path = DATA / ar_asset(target)
        if found:
            failed += 1
            print(f"{target.slug}: FAILED", file=sys.stderr)
            for problem in found:
                print(f"  {problem}", file=sys.stderr)
        else:
            from pxr import Usd
            meta = metadata(Usd.Stage.Open(str(path)))
            box = meta["bbox_angstrom"]
            print(f"{target.slug:<16} ok  {box[0]:7.2f} x {box[1]:7.2f} x {box[2]:7.2f} A"
                  f"  {meta['size_fit']['method'].split()[0]}")
    print(f"\n{len(AR_TARGETS) - failed} of {len(AR_TARGETS)} AR folds check",
          file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
