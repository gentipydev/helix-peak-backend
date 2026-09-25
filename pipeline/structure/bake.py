"""Bake `assets/models/<slug>.glb` — the fold on the last page of the walk.

Offline, run by hand, run once per protein. The app never touches this
directory at runtime; `hook/build.dart` compiles the committed `.glb` into a
`.fsceneb` and that is what ships.

This is `insulin_convert.pml` and `assemble_glb.py` made to take a target, and
every workaround they carried is still here because every one of them was found
by checking the output. Read `README.md` before changing any of it.

    pipeline/structure/venv/bin/python pipeline/structure/bake.py --target lysozyme
    pipeline/structure/venv/bin/python pipeline/structure/bake.py --all
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import shutil
import subprocess
import sys
import urllib.request

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

import numpy as np  # noqa: E402
import trimesh  # noqa: E402

from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import BY_SLUG, TARGETS, Structure, Target  # noqa: E402

STRUCTURES = HERE / "structures"
OUTPUT = HERE / "output"
PYMOL = shutil.which("pymol") or "/opt/homebrew/bin/pymol"

ROD_RADIUS = 0.5   # Angstrom. Comparable to a cartoon loop tube, so the
                   # bridges read as bonds rather than as more ribbon.
SECTIONS = 12
TUBE_RADIUS = 0.6  # For a peptide with no secondary structure to draw.

# What `AnatomyLayout`-sized pages can afford. Insulin's compiles to 173 kB.
FSCENEB_BUDGET_BYTES = 900_000


def pdb_path(structure: Structure) -> Path:
    path = STRUCTURES / f"{structure.pdb}.pdb"
    if not path.exists():
        url = f"https://files.rcsb.org/download/{structure.pdb}.pdb"
        print(f"  fetching {url}", flush=True)
        path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=120) as response:
            path.write_bytes(response.read())
    return path


def selection(chain: str, structure: Structure) -> str:
    """One chain, optionally clipped to the span the app claims to draw."""
    parts = [f"chain {chain}", "polymer"]
    if structure.residues is not None:
        parts.append(f"resi {structure.residues[0]}-{structure.residues[1]}")
    return " and ".join(parts)


def write_pml(target: Target, workspace: Path) -> Path:
    structure = target.structure
    lines = [
        f"load {pdb_path(structure)}, src",
        "remove solvent",
        "remove not polymer",
        "remove hydrogens",
        f"set cartoon_sampling, {structure.sampling}",
        f"set ribbon_sampling, {structure.sampling}",
    ]
    if structure.representation == "tube":
        # A nine-residue hormone has no helix and no sheet. Asked for a
        # cartoon, PyMOL draws a hairline; a tube down the backbone is the
        # shape the molecule actually has. Not `cartoon_trace_atoms`, which
        # threads the tube through the side chains as well.
        lines.append(f"set cartoon_tube_radius, {TUBE_RADIUS}")
    for chain in structure.chains:
        lines.append(f"create {chain.node}, src and {selection(chain.pdb_chain, structure)}")
        if structure.representation == "tube":
            lines.append(f"cartoon tube, {chain.node}")

    everything = " or ".join(c.node for c in structure.chains)
    lines += [
        "python",
        "from pymol import cmd",
        "# Export in model space, not camera space. PyMOL's OBJ exporter writes",
        "# whatever the current view transform produces, so pin the rotation to",
        "# identity and put the rotation origin at the molecule centre; the export",
        "# is then a pure translation of the PDB coordinates, which is the one",
        "# frame we can reason about.",
        f"ext = cmd.get_extent({everything!r})",
        "cx = (ext[0][0] + ext[1][0]) / 2.0",
        "cy = (ext[0][1] + ext[1][1]) / 2.0",
        "cz = (ext[0][2] + ext[1][2]) / 2.0",
        "cmd.set_view([",
        "    1.0, 0.0, 0.0,",
        "    0.0, 1.0, 0.0,",
        "    0.0, 0.0, 1.0,",
        "    0.0, 0.0, -100.0,",
        "    cx,  cy,  cz,",
        "    50.0, 150.0, -20.0,",
        "])",
        "print('EXPORT_ORIGIN %.6f %.6f %.6f' % (cx, cy, cz))",
        "python end",
    ]
    for chain in structure.chains:
        # The OBJ exporter ignores its selection argument and writes the whole
        # visible scene, so each object is exported by hiding everything else.
        lines += [
            "hide everything",
            f"show cartoon, {chain.node}",
            "refresh",
            f"save {workspace / (chain.node + '.obj')}",
        ]
    path = workspace / "convert.pml"
    path.write_text("\n".join(lines) + "\n")
    return path


def run_pymol(script: Path) -> np.ndarray:
    result = subprocess.run(
        [PYMOL, "-cq", str(script)], capture_output=True, text=True, cwd=str(HERE)
    )
    if result.returncode != 0:
        raise RuntimeError(f"pymol failed:\n{result.stdout}\n{result.stderr}")
    # PyMOL echoes the script's own source before running it, so the format
    # string itself comes past on stdout first. Match the line with numbers on
    # it, at the start of a line, which the echo never is.
    found = re.findall(
        r"^EXPORT_ORIGIN (-?\d+\.\d+) (-?\d+\.\d+) (-?\d+\.\d+)$",
        result.stdout,
        re.M,
    )
    if not found:
        raise RuntimeError(f"pymol printed no export origin:\n{result.stdout}\n{result.stderr}")
    return np.array([float(g) for g in found[-1]])


def atoms(path: Path) -> dict:
    """{(chain, resi, name): xyz} for the atoms the bridges need."""
    out = {}
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        if name not in ("SG", "CB", "CA"):
            continue
        out[(line[21], int(line[22:26]), name)] = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        )
    return out


def ssbonds(path: Path, structure: Structure) -> list[tuple]:
    """The SSBOND records with both ends inside what was exported.

    PyMOL cannot export sticks at all, so the bridges — the entire point of the
    page, where there are any — would be silently missing from a plain cartoon
    export. They are built from the crystallographic coordinates instead, with
    the pairs taken from the file's own records.
    """
    wanted = {c.pdb_chain for c in structure.chains}
    span = structure.residues
    out = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("SSBOND"):
            continue
        bond = (line[15], int(line[17:21]), line[29], int(line[31:35]))
        if bond[0] not in wanted or bond[2] not in wanted:
            continue
        if span is not None and not all(span[0] <= r <= span[1] for r in (bond[1], bond[3])):
            continue
        out.append(bond)
    return out


def rod(a, b):
    return trimesh.creation.cylinder(radius=ROD_RADIUS, segment=[a, b], sections=SECTIONS)


def joint(p):
    sphere = trimesh.creation.icosphere(subdivisions=1, radius=ROD_RADIUS)
    sphere.apply_translation(p)
    return sphere


def build_bonds(path: Path, structure: Structure, origin: np.ndarray):
    found = ssbonds(path, structure)
    if not structure.bonds:
        if found:
            print(f"  (ignoring {len(found)} SSBOND records; this page does not draw them)")
        return None
    if not found:
        raise ValueError(f"{structure.pdb}: bonds were asked for and the file has none")
    at, parts = atoms(path), []
    for c1, r1, c2, r2 in found:
        ca1, cb1, sg1 = (at[(c1, r1, n)] for n in ("CA", "CB", "SG"))
        ca2, cb2, sg2 = (at[(c2, r2, n)] for n in ("CA", "CB", "SG"))
        # CA->CB->SG->SG->CB->CA. The ribbon and tube run through CA, so a
        # bridge that starts at CB stops about 1.5 A short of the chain and
        # reads as floating beside it. Starting at CA buries both ends inside
        # the backbone. Spheres round the four bends where rods meet.
        parts += [
            rod(ca1, cb1), rod(cb1, sg1), rod(sg1, sg2), rod(sg2, cb2), rod(cb2, ca2),
            joint(cb1), joint(sg1), joint(sg2), joint(cb2),
        ]
        print(f"  bridge {c1}{r1}-{c2}{r2}: SG-SG {np.linalg.norm(sg1 - sg2):.2f} A")
    mesh = trimesh.util.concatenate(parts)
    mesh.apply_translation(-origin)   # PDB frame -> PyMOL export frame
    return mesh


def bake(target: Target) -> None:
    structure = target.structure
    workspace = OUTPUT / target.slug
    if workspace.exists():
        shutil.rmtree(workspace)
    workspace.mkdir(parents=True)

    print(f"{target.slug} <- {structure.pdb}", flush=True)
    origin = run_pymol(write_pml(target, workspace))

    meshes = {}
    for chain in structure.chains:
        mesh = trimesh.load(workspace / f"{chain.node}.obj", process=False)
        if len(mesh.vertices) == 0:
            raise ValueError(f"{target.slug}: {chain.node} exported no geometry")
        meshes[chain.node] = mesh
    bonds = build_bonds(pdb_path(structure), structure, origin)
    if bonds is not None:
        meshes["bonds"] = bonds

    # Nothing normalises the model but us: centre it and scale the longest axis
    # to exactly 1.0, so the Dart side carries no scale constant and the camera
    # can be framed once. See `_framingMargin` in structure_view.dart.
    stacked = np.vstack([np.asarray(m.vertices) for m in meshes.values()])
    low, high = stacked.min(axis=0), stacked.max(axis=0)
    centre, extent = (low + high) / 2, (high - low)
    scale = 1.0 / extent.max()
    print(f"  pre  {extent[0]:.2f} x {extent[1]:.2f} x {extent[2]:.2f} A", flush=True)

    scene = trimesh.Scene()
    for name, mesh in meshes.items():
        normals = np.asarray(mesh.vertex_normals).copy()  # keep PyMOL's smooth shading
        mesh.vertices = (np.asarray(mesh.vertices) - centre) * scale
        mesh.vertex_normals = normals
        scene.add_geometry(mesh, geom_name=name, node_name=name)
        print(f"  {name:8s} v={len(mesh.vertices):<7d} f={len(mesh.faces):<7d}")

    after = np.vstack([np.asarray(m.vertices) for m in meshes.values()])
    low2, high2 = after.min(axis=0), after.max(axis=0)
    spread, middle = high2 - low2, (low2 + high2) / 2
    if abs(spread.max() - 1.0) > 1e-6 or np.abs(middle).max() > 1e-6:
        raise ValueError(f"{target.slug}: normalisation left {spread} centred at {middle}")

    staged = workspace / f"{target.slug}.glb"
    scene.export(staged)
    destination = DATA / target.structure_asset
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(staged, destination)
    size = destination.stat().st_size
    print(f"  wrote {destination.relative_to(DATA)}  {size:,} B\n", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", action="append")
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    chosen = [BY_SLUG[s] for s in args.target] if args.target else list(TARGETS)
    if not args.target and not args.all:
        raise SystemExit("Pass --target <slug> or --all")
    for target in chosen:
        bake(target)
