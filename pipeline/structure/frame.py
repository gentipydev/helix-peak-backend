"""Where a PDB point lands in a stored structure model.

`bake.py` exports PyMOL's cartoon as a pure translation of the PDB frame -- no
rotation -- then centres it and divides by its longest axis, L. So a point `p`
of the entry sits at `(p - c) / L` in the stored `.glb`, for one centre `c` and
one length `L` per model. Neither is written down anywhere; this recovers both
from the stored model and the entry it was cut from, so a baker can place
something of its own in the frame the fold page draws, without PyMOL.

Moved unchanged from `structure_ar/bake_ar.py`, which found each model's real
size this way first, so that every track placed in the model's frame uses the
same one.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pipeline.structure.glb import Mesh
from pipeline.structure.pdb import bridge_atoms, ca_atoms
from pipeline.targets import Structure

# The fit stops when a round moves the length by less than this fraction.
_SETTLED = 1e-7
_ROUNDS = 60


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


def model_frame(structure: Structure, meshes: dict[str, Mesh]) -> tuple[Fit, str]:
    """The stored model's frame, and how it was found.

    Exact where the model has bridges to read the scale from; fitted to the
    ribbon's CA atoms where it has none.
    """
    vertices = np.vstack([m.positions for m in meshes.values()])
    exact = fit_bridges(meshes["bonds"], bridge_atoms(structure)) if "bonds" in meshes else None
    fit = exact or fit_size(vertices, ca_atoms(structure))
    method = ("disulfide CB and SG atoms at the bridges' joints" if exact
              else "exported CA atoms fitted to the ribbon")
    return fit, method
