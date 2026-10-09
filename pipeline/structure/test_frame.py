"""Finding where the entry's atoms sit in a stored structure model.

Offline. Moved with `frame.py` from `structure_ar/`, whose bake found each
model's size this way first. The tests that read the stored `.glb`s skip until
`fetch_tracks.py --kind structure` has put them in `pipeline/data/`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.paths import DATA  # noqa: E402
from pipeline.structure import frame  # noqa: E402
from pipeline.structure.glb import Mesh, read_glb  # noqa: E402
from pipeline.structure.pdb import bridge_atoms, ca_atoms  # noqa: E402
from pipeline.targets import BY_SLUG, TARGETS  # noqa: E402


def test_the_bridges_give_the_scale_exactly():
    rng = np.random.default_rng(7)
    length, centre = 31.7, np.array([4.0, -12.5, 30.25])
    bridges = [rng.uniform(-15, 15, size=(4, 3)) + centre for _ in range(3)]
    joint = rng.normal(size=(frame._JOINT_VERTICES, 3))
    joint -= joint.mean(axis=0)
    blocks = []
    for atoms in bridges:
        model = (atoms - centre) / length
        rods = rng.normal(size=(5 * frame._ROD_VERTICES, 3))
        blocks += [rods] + [m + joint * 0.01 for m in model]
    bonds = Mesh(positions=np.vstack(blocks), normals=None,
                 triangles=np.zeros((0, 3), int))
    fit = frame.fit_bridges(bonds, bridges)
    assert fit.length == pytest.approx(length, rel=1e-9)
    assert np.allclose(fit.centre, centre)
    assert fit.rms < 1e-9
    # A mesh not laid out the bake's way is not read as one.
    assert frame.fit_bridges(bonds, bridges[:2]) is None


def test_the_ca_fit_is_exact_where_the_ribbon_runs_through_the_atoms():
    # A helical CA trace and a model sampled along it, as a ribbon's spine is.
    # On a real cartoon the ribbon's width moves the answer a little: measured
    # on the thirteen stored models whose bridges give the size exactly, the
    # CA fit lands within 1.2% (see the insulin test below).
    t = np.linspace(0, 6 * np.pi, 120)
    trace = np.stack([2.3 * np.cos(t), 1.5 * t, 2.3 * np.sin(t)], axis=1)
    spine = np.vstack([trace[:-1] + (trace[1:] - trace[:-1]) * f
                       for f in (0, 0.25, 0.5, 0.75)])
    length, centre = 30.0, np.array([3.0, 14.0, -2.0])
    fit = frame.fit_size((spine - centre) / length, trace)
    assert fit.length == pytest.approx(length, rel=1e-3)
    assert np.allclose(fit.centre, centre, atol=0.05)
    assert fit.rms < 0.05


_STORED = DATA / BY_SLUG["insulin"].structure_asset


@pytest.mark.skipif(not _STORED.exists(), reason="run fetch_tracks.py --kind structure first")
def test_insulin_comes_back_at_the_size_its_bake_logged():
    # `structure/README.md`: "Bounding box before normalising: 24.31 x 18.66 x 20.33 A."
    structure = BY_SLUG["insulin"].structure
    meshes = read_glb(_STORED.read_bytes())
    vertices = np.vstack([m.positions for m in meshes.values()])
    fit = frame.fit_bridges(meshes["bonds"], bridge_atoms(structure))
    extent = (vertices.max(axis=0) - vertices.min(axis=0)) * fit.length
    assert np.allclose(extent, [24.31, 18.66, 20.33], atol=0.005)
    assert fit.rms < 0.01
    # The ribbon fit, on its own, lands within about one per cent of it.
    ribbon = frame.fit_size(vertices, ca_atoms(structure))
    assert ribbon.length == pytest.approx(fit.length, rel=0.015)


_BRIDGED = [t for t in TARGETS if t.structure.bonds]


@pytest.mark.skipif(not all((DATA / t.structure_asset).exists() for t in _BRIDGED),
                    reason="run fetch_tracks.py --kind structure first")
@pytest.mark.parametrize("target", _BRIDGED, ids=lambda t: t.slug)
def test_every_stored_bridge_is_read_where_the_bake_built_it(target):
    # `bridge_atoms` is the bake's own reading of the SSBOND records and atom
    # lines, so its atoms land on the joints of the stored `bonds` node to
    # float precision, and `model_frame` takes that exact answer.
    meshes = read_glb((DATA / target.structure_asset).read_bytes())
    fit, method = frame.model_frame(target.structure, meshes)
    assert method == "disulfide CB and SG atoms at the bridges' joints"
    assert fit.atoms == 4 * len(bridge_atoms(target.structure))
    assert fit.rms < 0.001


@pytest.mark.skipif(not (DATA / BY_SLUG["myoglobin"].structure_asset).exists(),
                    reason="run fetch_tracks.py --kind structure first")
def test_a_model_without_bridges_is_fitted_to_its_ribbon():
    target = BY_SLUG["myoglobin"]
    meshes = read_glb((DATA / target.structure_asset).read_bytes())
    fit, method = frame.model_frame(target.structure, meshes)
    assert method == "exported CA atoms fitted to the ribbon"
    assert fit.atoms == len(ca_atoms(target.structure))
    assert fit.rms < 0.5
