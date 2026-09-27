"""Build each state of an assembly as one PDB file, both in one frame.

A state whose entry holds the whole molecule is that entry's chains. A state
whose entry holds less is built from its biological assembly: the REMARK 350
BIOMT operators applied to the chains they name, each copy lettered in turn.
Then every state after the first is superposed onto the first (Kabsch, on the
CA atoms of the subunits the assembly holds still), which is what makes the
two one frame: their entries come from different crystals, whose coordinates
agree about nothing.

Kept: the polymer chains' ATOM records, the hemes' and bound oxygen's HETATM
records (the bake draws neither, but a viewer of the pair can show where
oxygen binds), and the HELIX and SHEET records, so PyMOL draws the entry's own
secondary structure. Water and everything else is left out.
"""

from __future__ import annotations

import string
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pipeline.assemblies.assemblies import Assembly, State
from pipeline.structure.pdb import STRUCTURES

# The HETATM groups a state keeps: the haem, and oxygen bound to it.
KEPT_GROUPS = ("HEM", "OXY")


@dataclass(frozen=True)
class Operator:
    """One BIOMT transform, and the entry chains it applies to."""

    chains: tuple[str, ...]
    matrix: np.ndarray          # 3 x 4: rotation, then translation


def operators(path: Path, biomolecule: int) -> list[Operator]:
    """The entry's BIOMT operators for one biological assembly, in order."""
    found: dict[int, np.ndarray] = {}
    chains: dict[int, tuple[str, ...]] = {}
    current, applying = None, ()
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("REMARK 350"):
            continue
        words = line.split()
        if "BIOMOLECULE:" in words:
            current = int(words[words.index("BIOMOLECULE:") + 1])
        elif current == biomolecule and "CHAINS:" in words:
            applying = tuple(c.strip(",") for c in words[words.index("CHAINS:") + 1:] if c.strip(","))
        elif current == biomolecule and len(words) >= 8 and words[2] in ("BIOMT1", "BIOMT2", "BIOMT3"):
            number, row = int(words[3]), int(words[2][-1]) - 1
            found.setdefault(number, np.zeros((3, 4)))[row] = [float(v) for v in words[4:8]]
            chains[number] = applying
    if not found:
        raise ValueError(f"{path.name}: no biological assembly {biomolecule}")
    return [Operator(chains[n], found[n]) for n in sorted(found)]


def _moved(line: str, xyz: np.ndarray, chain: str | None = None) -> str:
    """[line] with its coordinates replaced and, given one, its chain letter."""
    out = f"{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"
    if chain is not None:
        out = out[:21] + chain + out[22:]
    return out


def _xyz(line: str) -> np.ndarray:
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def _kept(line: str) -> bool:
    if line.startswith("ATOM"):
        return True
    return line.startswith("HETATM") and line[17:20].strip() in KEPT_GROUPS


def _relettered(line: str, letters: dict[str, str]) -> str:
    """A HELIX or SHEET record for a copied chain, under its new letter."""
    if line.startswith("HELIX "):
        return line[:19] + letters[line[19]] + line[20:31] + letters[line[31]] + line[32:]
    return line[:21] + letters[line[21]] + line[22:32] + letters[line[32]] + line[33:]


def state_lines(state: State, entry: Path) -> list[str]:
    """The state's molecule, as PDB records: the entry's, or its assembly's."""
    text = entry.read_text(errors="replace").splitlines()
    structure = [l for l in text if l.startswith(("HELIX ", "SHEET "))]
    atoms = []
    for line in text:
        if line.startswith("ENDMDL"):
            break
        if _kept(line):
            atoms.append(line)
    if state.biomt is None:
        return structure + atoms
    out_structure, out_atoms = [], []
    used: list[str] = []
    for k, op in enumerate(operators(entry, state.biomt)):
        letters = {}
        for chain in op.chains:
            if k == 0:
                letters[chain] = chain
            else:
                letters[chain] = next(c for c in string.ascii_uppercase
                                      if c not in used and c not in letters.values())
            used.append(letters[chain])
        for line in structure:
            ends = (line[19], line[31]) if line.startswith("HELIX ") else (line[21], line[32])
            if all(c in letters for c in ends):
                out_structure.append(_relettered(line, letters))
        for line in atoms:
            if line[21] in letters:
                moved = op.matrix[:, :3] @ _xyz(line) + op.matrix[:, 3]
                out_atoms.append(_moved(line, moved, letters[line[21]]))
    return out_structure + out_atoms


def chains_of(lines: list[str]) -> list[str]:
    """The polymer chains the records hold, in order."""
    seen: list[str] = []
    for line in lines:
        if line.startswith("ATOM") and line[21] not in seen:
            seen.append(line[21])
    return seen


def ca_of(lines: list[str]) -> dict[tuple[str, int], np.ndarray]:
    """Each residue's CA, first alternate location."""
    out: dict[tuple[str, int], np.ndarray] = {}
    for line in lines:
        if line.startswith("ATOM") and line[12:16] == " CA " and line[16] in (" ", "A"):
            out.setdefault((line[21], int(line[22:26])), _xyz(line))
    return out


def kabsch(moving: np.ndarray, fixed: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The rotation R and shift t that best put [moving] on [fixed]: R p + t."""
    pm, fm = moving.mean(axis=0), fixed.mean(axis=0)
    u, _, vt = np.linalg.svd((moving - pm).T @ (fixed - fm))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rotation, fm - rotation @ pm


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def turn(rotation: np.ndarray) -> float:
    """How far a rotation turns, in degrees."""
    return float(np.degrees(np.arccos(np.clip((np.trace(rotation) - 1) / 2, -1, 1))))


@dataclass(frozen=True)
class Pair:
    """Both states' records in the first state's frame, and how they got there."""

    lines: dict[str, list[str]]         # by state name
    rotation: np.ndarray                # applied to the second state
    shift: np.ndarray
    report: dict


def build(assembly: Assembly) -> Pair:
    first, second = assembly.states
    lines = {s.name: state_lines(s, STRUCTURES / f"{s.pdb}.pdb") for s in assembly.states}
    wanted = [u.pdb_chain for u in assembly.subunits]
    for state in assembly.states:
        found = chains_of(lines[state.name])
        if sorted(found) != sorted(wanted):
            raise ValueError(f"{assembly.slug} {state.name} ({state.pdb}): chains "
                             f"{found}, the assembly is {wanted}")

    fixed, moving = ca_of(lines[first.name]), ca_of(lines[second.name])
    shared = sorted(set(fixed) & set(moving))
    held = {assembly.subunit(n).pdb_chain for n in assembly.superpose_on}
    on = [k for k in shared if k[0] in held]
    rest = [k for k in shared if k[0] not in held]
    p = np.array([moving[k] for k in shared])
    q = np.array([fixed[k] for k in shared])
    raw = rmsd(p, q)
    rotation, shift = kabsch(np.array([moving[k] for k in on]), np.array([fixed[k] for k in on]))
    moved = p @ rotation.T + shift
    index = {k: i for i, k in enumerate(shared)}
    best_r, best_t = kabsch(p, q)
    rest_r, _ = kabsch(moved[[index[k] for k in rest]], q[[index[k] for k in rest]])
    report = {
        "matched_ca": len(shared),
        "raw_rmsd_angstrom": round(raw, 2),
        "superposed_on": list(assembly.superpose_on),
        "held_rmsd_angstrom": round(rmsd(moved[[index[k] for k in on]], q[[index[k] for k in on]]), 2),
        "turned_rmsd_angstrom": round(rmsd(moved[[index[k] for k in rest]], q[[index[k] for k in rest]]), 2),
        "tetramer_rmsd_angstrom": round(rmsd(moved, q), 2),
        "best_fit_rmsd_angstrom": round(rmsd(p @ best_r.T + best_t, q), 2),
        "turned_degrees": round(turn(rest_r), 1),
        "entries_apart_degrees": round(turn(best_r), 1),
    }
    lines[second.name] = [
        _moved(l, rotation @ _xyz(l) + shift) if l.startswith(("ATOM", "HETATM")) else l
        for l in lines[second.name]
    ]
    return Pair(lines=lines, rotation=rotation, shift=shift, report=report)


def write(lines: list[str], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines + ["END"]) + "\n")
    return path
