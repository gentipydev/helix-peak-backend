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

from dataclasses import dataclass
from pathlib import Path
import urllib.request

import numpy as np

from pipeline.targets import Structure

STRUCTURES = Path(__file__).resolve().parent / "structures"

# The residue names an entry uses, and the letter a record's translation uses.
ONE_LETTER = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O",
}


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


@dataclass(frozen=True)
class Residue:
    """One residue of an exported chain, as the entry numbers and names it."""

    number: int                 # the entry's own residue number
    letter: str                 # one letter; a modified residue reads as the one it modifies
    ca: np.ndarray | None       # its CA in the PDB frame, or None: listed, never located


def modified(path: Path) -> dict[tuple[str, int], str]:
    """{(chain, number): standard residue} from the entry's MODRES records.

    Relaxin's A chain starts on a pyroglutamate, written as HETATM `PCA`; its
    MODRES record says it is a modified glutamine. Only residues named here are
    read from HETATM lines, so a ligand -- a haem, a calcium ion, whose atom is
    also called CA -- is never taken for a residue.
    """
    out = {}
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("MODRES"):
            out[(line[16], int(line[18:22]))] = line[24:27].strip()
    return out


def inserted(path: Path) -> set[tuple[str, int]]:
    """(chain, number) of every residue the entry adds to the protein.

    SEQADV records with no database residue: expression tags, linkers and caps,
    which are the construct's and not the protein's. A substitution (an
    engineered mutation, a sequence conflict) names the database residue it
    replaces and stays: it is the protein's residue, carrying another side chain.
    """
    out = set()
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("SEQADV") and not line[39:42].strip():
            out.add((line[16], int(line[18:22])))
    return out


def substituted(path: Path) -> set[tuple[str, int]]:
    """(chain, number) of every residue the entry says differs from the protein's.

    SEQADV records that do name a database residue: engineered mutations,
    sequence conflicts, and the like. The entry's coordinates at such a residue
    carry its own side chain, not the protein's.
    """
    out = set()
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("SEQADV") and line[39:42].strip():
            out.add((line[16], int(line[18:22])))
    return out


def chain_residues(structure: Structure, chain: str, pdb: Path | None = None) -> list[Residue]:
    """Every residue the entry holds of one exported chain, in its numbering.

    A residue with any atom in the ATOM records, or in the HETATM records
    `modified` names, is listed, with its CA where it has one: of alternate
    locations the first, of models the first. A residue can be modelled without
    its CA (1HGU's Pro37 is one nitrogen), and then has none. Residues the
    experiment did not locate at all come from REMARK 465, with no CA.
    Residues `inserted` names are left out. Clipped to the span the structure
    exports, where it has one. Refuses an insertion code, which would make the
    entry's numbering something other than a count.
    """
    path = pdb or STRUCTURES / f"{structure.pdb}.pdb"
    span = structure.residues
    mods, extra = modified(path), inserted(path)
    names: dict[int, str] = {}
    cas: dict[int, np.ndarray] = {}
    listing = False
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("REMARK 465"):
            # The list starts after its own column header, not after the prose
            # above it that explains the columns.
            if line[11:].split() == ["M", "RES", "C", "SSSEQI"]:
                listing = True
                continue
            if listing and len(line) > 26 and line[19] == chain and line[15:18].strip():
                number = int(line[21:26])
                if line[26].strip():
                    raise ValueError(f"{structure.pdb} {chain}{number}: insertion code")
                names.setdefault(number, line[15:18].strip())
            continue
        if not line.startswith(("ATOM", "HETATM")) or line[21] != chain:
            continue
        number, name = int(line[22:26]), line[17:20].strip()
        if line[26].strip():
            raise ValueError(f"{structure.pdb} {chain}{number}: insertion code")
        if line.startswith("HETATM"):
            if (chain, number) not in mods:
                continue
            name = mods[(chain, number)]
        names[number] = name
        if line[12:16] == " CA " and line[16] in (" ", "A") and number not in cas:
            cas[number] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return [Residue(n, ONE_LETTER.get(names[n], "X"), cas.get(n)) for n in sorted(names)
            if (chain, n) not in extra
            and (span is None or span[0] <= n <= span[1])]


def secondary_structure(path: Path) -> dict[tuple[str, int], str]:
    """{(chain, number): "helix" or "strand"} from the entry's HELIX and SHEET records.

    Every helix class counts as a helix (alpha and 3-10 alike), and every strand
    of every sheet as a strand; a residue named by neither is coil. The records
    are the entry's own assignment, deposited with its coordinates.
    """
    out = {}
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("HELIX "):
            chain, first, last, label = line[19], int(line[21:25]), int(line[33:37]), "helix"
        elif line.startswith("SHEET "):
            chain, first, last, label = line[21], int(line[22:26]), int(line[33:37]), "strand"
        else:
            continue
        for number in range(first, last + 1):
            out[(chain, number)] = label
    return out
