"""Check every baked `folding` track against its record, its entry and its model.

    .venv/Scripts/python pipeline/fetch_tracks.py --kind structure --kind record
    .venv/Scripts/python pipeline/folding/check_folding.py

For each target, the payload under `pipeline/data/assets/folding/`:

- names its protein: the slug, gene, UniProt accession and PDB entry are the
  table's;
- has the structure's chains: one per exported chain, in the table's order,
  with its node name and PDB chain;
- covers the mature chain: each chain's residues run one by one from the first
  residue of the mature chain it is cut from to the last -- a kept region of a
  cleaved precursor, or the whole of an uncut one less what the table removes
  from its ends -- clipped to the span the structure exports. So the residue
  count is the mature chain's length. Where the record names that stretch as
  a mature peptide, the letters are the peptide's;
- reads the record's residues: each letter is the record's at that position,
  and where the entry carries another, the entry declares it (SEQADV);
- places the entry's atoms: each ordered residue's CA is the entry's CA for
  that residue, moved by the track's frame, to 1e-4 model units;
- is a chain: every CA 2.8 to 4.4 A from the next (3.8 for a trans peptide
  bond, 2.9 for a cis one, and a 2.5 A entry stretches it), and across
  residues with no place, no more than 3.8 A a residue;
- sits in the stored model: the frame names the stored `.glb` by its sha256
  and carries its bounding box, and every helix and coil CA lies on its own
  chain's ribbon, within 0.25 A of
  the nearest vertex (0.6 A for a tube), as `structure/verify_frame.py` holds
  insulin's. A strand's CA, which PyMOL's flat arrows smooth past, within 3 A;
- says only what it can: an ordered residue carries a CA and a helix, strand
  or coil label and the others neither, and an absent residue is never
  between two that are there.

Exits 1 with the problems listed. Reads `pipeline/data/` and the entries in
`pipeline/structure/structures/`.
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

from pipeline.folding.bake_folding import (  # noqa: E402
    FOLDING_TARGETS,
    SCHEMA_VERSION,
    folding_asset,
)
from pipeline.paths import DATA  # noqa: E402
from pipeline.structure.frame import _nearest  # noqa: E402
from pipeline.structure.glb import read_glb  # noqa: E402
from pipeline.structure.pdb import chain_residues, pdb_path, substituted  # noqa: E402
from pipeline.targets import Target  # noqa: E402

STATES = ("ordered", "disordered", "absent")
LABELS = ("helix", "strand", "coil")

# CA to the next CA, in angstroms: a cis peptide bond (amylase has two) closes
# it to 2.9, a trans one holds it at 3.8, and 1HGU, at 2.5 A, stretches one to 4.3.
STEP = (2.8, 4.4)
PER_RESIDUE = 3.8

# How far a CA may sit from its own chain's ribbon, in angstroms. PyMOL draws a
# helix as an oval 0.25 A thick and a loop as a tube of radius 0.2 through the
# CA atoms; a peptide drawn as a tube has radius 0.6 (`structure/bake.py`).
ON_RIBBON = 0.25
ON_TUBE = 0.6
ON_STRAND = 3.0
_SLACK = 0.005      # the coordinates are rounded to 1e-5 model units


def record_of(target: Target) -> dict | None:
    path = DATA / target.mock_asset
    return json.loads(path.read_bytes()) if path.exists() else None


def mature_chain(target: Target, ordered: list[int]) -> tuple[int, int]:
    """The mature chain the ordered positions belong to, from the region table."""
    if target.cleaved:
        holding = [r for r in target.regions
                   if r.kept and r.start <= min(ordered) and max(ordered) <= r.end]
        if len(holding) != 1:
            raise ValueError(f"{len(holding)} kept regions hold {min(ordered)}-{max(ordered)}")
        return holding[0].start, holding[0].end
    removed = sorted((r.start, r.end) for r in target.regions if not r.kept)
    start, end = 1, target.aa
    for low, high in removed:
        if low == start:
            start = high + 1
    for low, high in reversed(removed):
        if high == end:
            end = low - 1
    return start, end


def problems_of(target: Target, track: dict | None = None) -> list[str]:
    """Everything wrong with one protein's track: the file's, or [track]."""
    where = target.slug
    if track is None:
        path = DATA / folding_asset(target)
        if not path.exists():
            return [f"{where}: not baked ({path})"]
        track = json.loads(path.read_bytes())

    out: list[str] = []
    structure = target.structure
    for key, want in (("slug", target.slug), ("gene", target.gene),
                      ("uniprot", target.uniprot), ("pdb", structure.pdb),
                      ("schema_version", SCHEMA_VERSION)):
        if track.get(key) != want:
            out.append(f"{where}: {key} is {track.get(key)!r}, expected {want!r}")

    record = record_of(target)
    model_path = DATA / target.structure_asset
    if record is None or not model_path.exists():
        return out + [f"{where}: no stored record or model to hold it to; "
                      f"run fetch_tracks.py --kind record --kind structure first"]
    translation = record["protein"]["translation"]
    glb = model_path.read_bytes()
    meshes = read_glb(glb)

    frame = track.get("frame") or {}
    if frame.get("glb_sha256") != hashlib.sha256(glb).hexdigest():
        out.append(f"{where}: the frame names another .glb than the stored one")
    stacked = np.vstack([m.positions for m in meshes.values()])
    bounds = frame.get("bounds") or {}
    for side, want in (("min", stacked.min(axis=0)), ("max", stacked.max(axis=0))):
        got = bounds.get(side)
        if not isinstance(got, list) or len(got) != 3 or np.abs(np.array(got) - want).max() > 1e-6:
            out.append(f"{where}: the frame's bounds {side} is not the stored model's")
    try:
        centre = np.array(frame["centre_angstrom"], dtype=float)
        length = float(frame["length_angstrom"])
    except (KeyError, TypeError, ValueError):
        return out + [f"{where}: the frame has no centre or length"]

    chains = track.get("chains") or []
    wanted = [(c.node, c.pdb_chain) for c in structure.chains]
    if [(c.get("node"), c.get("pdb_chain")) for c in chains] != wanted:
        return out + [f"{where}: chains {[(c.get('node'), c.get('pdb_chain')) for c in chains]}, "
                      f"the structure exports {wanted}"]

    entry = pdb_path(structure)
    declared = substituted(entry)
    on_ribbon = ON_TUBE if structure.representation == "tube" else ON_RIBBON
    peptides = {(translation.find(p["translation"]) + 1,
                 translation.find(p["translation"]) + len(p["translation"])): p["translation"]
                for p in record.get("peptides") or [] if p.get("translation")}

    for chain in chains:
        at = f"{where} {chain['node']}"
        residues = chain.get("residues") or []
        states = [r.get("state") for r in residues]
        if any(s not in STATES for s in states):
            out.append(f"{at}: a state outside {STATES}")
            continue
        ordered = [r["n"] for r in residues if r["state"] == "ordered"]
        if not ordered:
            out.append(f"{at}: nothing ordered")
            continue

        # The mature chain, and every residue of it, one by one.
        try:
            start, end = mature_chain(target, ordered)
        except ValueError as problem:
            out.append(f"{at}: {problem}")
            continue
        offset = chain.get("offset")
        if structure.residues is not None and isinstance(offset, int):
            start = max(start, structure.residues[0] + offset)
            end = min(end, structure.residues[1] + offset)
        numbers = [r.get("n") for r in residues]
        if numbers != list(range(start, end + 1)):
            out.append(f"{at}: residues {numbers[:1]}..{numbers[-1:]} ({len(numbers)}), "
                       f"the mature chain is {start}-{end} ({end - start + 1})")
            continue
        if (chain.get("first"), chain.get("last")) != (start, end):
            out.append(f"{at}: says {chain.get('first')}-{chain.get('last')}, "
                       f"the residues run {start}-{end}")
        letters = "".join(r.get("aa", "") for r in residues)
        if letters != translation[start - 1:end]:
            out.append(f"{at}: the letters are not the record's {start}-{end}")
        if (start, end) in peptides and peptides[(start, end)] != letters:
            out.append(f"{at}: the letters are not the record's mature peptide")

        # Where the entry and the record part, the entry says so.
        for change in chain.get("entry_differs", []):
            number = change["n"] - (offset or 0)
            if (chain["pdb_chain"], number) not in declared:
                out.append(f"{at}: the entry differs at {change['n']} and declares no substitution")

        # Each state says only what it can, and nothing absent sits inside.
        for r in residues:
            has = ("ca" in r, "ss" in r)
            if r["state"] == "ordered":
                if has != (True, True) or r["ss"] not in LABELS or len(r["ca"]) != 3:
                    out.append(f"{at} {r['n']}: an ordered residue needs a CA and a label")
            elif has != (False, False):
                out.append(f"{at} {r['n']}: a {r['state']} residue has no place to give")
        inside = states[states.index("ordered"):len(states) - states[::-1].index("ordered")]
        if "absent" in inside:
            out.append(f"{at}: an absent residue between two that are there")

        # The CA atoms are the entry's, placed by the frame.
        if not isinstance(offset, int):
            out.append(f"{at}: no offset to read the entry by")
            continue
        entry_cas = {r.number + offset: r.ca for r in chain_residues(structure, chain["pdb_chain"])
                     if r.ca is not None}
        placed = {r["n"]: np.array(r["ca"], dtype=float) for r in residues if r["state"] == "ordered"}
        if sorted(entry_cas) != sorted(placed):
            out.append(f"{at}: ordered at {sorted(set(placed) ^ set(entry_cas))[:6]}, "
                       f"which the entry does not locate that way")
            continue
        for n, point in placed.items():
            if np.abs((entry_cas[n] - centre) / length - point).max() > 1e-4:
                out.append(f"{at} {n}: the CA is not the entry's in this frame")
                break

        # A chain: each CA a peptide bond from the next.
        previous = None
        for n in ordered:
            if previous is not None:
                gap = n - previous
                step = float(np.linalg.norm(placed[n] - placed[previous])) * length
                if gap == 1 and not STEP[0] <= step <= STEP[1]:
                    out.append(f"{at} {previous}-{n}: CA to CA {step:.2f} A")
                if gap > 1 and step > PER_RESIDUE * gap:
                    out.append(f"{at} {previous}-{n}: {step:.1f} A across {gap} residues")
            previous = n

        # On its own chain's ribbon in the stored model.
        mesh = meshes.get(chain["node"])
        if mesh is None:
            out.append(f"{at}: the stored model has no {chain['node']} node")
            continue
        points = np.array([placed[n] for n in ordered])
        distance = np.linalg.norm(points - _nearest(points, mesh.positions), axis=1) * length
        for n, far in zip(ordered, distance):
            label = next(r["ss"] for r in residues if r["n"] == n)
            allowed = ON_STRAND if label == "strand" else on_ribbon
            if far > allowed + _SLACK:
                out.append(f"{at} {n} ({label}): {far:.2f} A off its ribbon")

    if track.get("secondary_structure") != f"HELIX and SHEET records of {structure.pdb}":
        out.append(f"{where}: secondary_structure is {track.get('secondary_structure')!r}")
    if not track.get("built_by"):
        out.append(f"{where}: says nothing of what built it")
    return out


def main() -> int:
    found = [problem for target in FOLDING_TARGETS for problem in problems_of(target)]
    for problem in found:
        print(f"  {problem}", file=sys.stderr)
    total = sum((DATA / folding_asset(t)).stat().st_size
                for t in FOLDING_TARGETS if (DATA / folding_asset(t)).exists())
    print(f"{len(FOLDING_TARGETS)} folding tracks checked, {len(found)} problem(s), "
          f"{total:,} B", file=sys.stderr)
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
