"""Reading a PDB entry the way the structure bake exports it.

One reading of `structures/<entry>.pdb`, shared by every baker that starts
from the twenty entries: `bake.py` builds the bridges from it, `structure_ar`
recovers each stored model's size from it, and `verify_frame.py` audits the
export frame with it. Each function here moved unchanged from the baker that
wrote it first, so every caller reads the file exactly as it did before.

It needs only numpy, so a baker that reads the stored model instead of running
PyMOL can use it from the backend's own `.venv`.
"""

from __future__ import annotations

from pathlib import Path
import urllib.request

import numpy as np

from pipeline.targets import Structure

STRUCTURES = Path(__file__).resolve().parent / "structures"


def pdb_path(structure: Structure) -> Path:
    path = STRUCTURES / f"{structure.pdb}.pdb"
    if not path.exists():
        url = f"https://files.rcsb.org/download/{structure.pdb}.pdb"
        print(f"  fetching {url}", flush=True)
        path.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=120) as response:
            path.write_bytes(response.read())
    return path


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


def bridge_atoms(structure: Structure, pdb: Path | None = None) -> list[np.ndarray]:
    """Each exported bridge's CB, SG, SG and CB, in the order the bake built them.

    `ssbonds` and `atoms`, the bake's own reading of the file: SSBOND records
    with both ends in the exported chains and span, in file order, and each
    atom's last ATOM line. `bake.py` puts a joint sphere on each of these four
    atoms, in this order, which is what lets the stored model's size be read
    back from its `bonds` node.
    """
    path = pdb or STRUCTURES / f"{structure.pdb}.pdb"
    at = atoms(path)
    return [np.array([at[(c1, r1, "CB")], at[(c1, r1, "SG")],
                      at[(c2, r2, "SG")], at[(c2, r2, "CB")]])
            for c1, r1, c2, r2 in ssbonds(path, structure)]
