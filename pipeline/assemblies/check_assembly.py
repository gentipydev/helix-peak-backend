"""Check every baked assembly morph against its entries, rebuilt.

    .venv/Scripts/python pipeline/assemblies/check_assembly.py

For each assembly, the morph under `pipeline/data/assets/assemblies/<slug>/`:

- names the assembly and its states' entries, as `assemblies.py` has them;
- holds four subunits in each state: `build.py` rebuilds both molecules from
  the entries in `structure/structures/` (the second from its biological
  assembly) and counts four chains in each;
- is one frame: each residue's CA, in each state, is the rebuilt molecule's
  CA placed by the morph's frame, to 1e-4 model units; the frame's
  superposition is the one the rebuild finds; and both states' exports were
  audited against one origin, every CA within 0.25 A of its own ribbon;
- is a chain in each state: every CA 2.8 to 4.4 A from the next;
- has one haem iron per subunit per state, and oxygen bound only where the
  state's entry is an oxy structure;
- names the models beside it by their sha256.

Exits 1 with the problems listed. Reads `pipeline/data/` and the entries only.
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

from pipeline.assemblies.assemblies import ASSEMBLIES, Assembly  # noqa: E402
from pipeline.assemblies.bake_assembly import SCHEMA_VERSION, morph_asset  # noqa: E402
from pipeline.assemblies.build import build, ca_of, chains_of  # noqa: E402
from pipeline.paths import DATA  # noqa: E402

STEP = (2.8, 4.4)
ON_RIBBON = 0.25


def problems_of(assembly: Assembly, track: dict | None = None) -> list[str]:
    where = assembly.slug
    if track is None:
        path = DATA / morph_asset(assembly)
        if not path.exists():
            return [f"{where}: not baked ({path})"]
        track = json.loads(path.read_bytes())
    out: list[str] = []
    for key, want in (("slug", assembly.slug), ("display", assembly.display),
                      ("schema_version", SCHEMA_VERSION)):
        if track.get(key) != want:
            out.append(f"{where}: {key} is {track.get(key)!r}, expected {want!r}")
    states = [(s.get("name"), s.get("pdb")) for s in track.get("states", [])]
    if states != [(s.name, s.pdb) for s in assembly.states]:
        return out + [f"{where}: states {states}"]

    pair = build(assembly)
    for state in assembly.states:
        if len(chains_of(pair.lines[state.name])) != 4:
            out.append(f"{where} {state.name}: {len(chains_of(pair.lines[state.name]))} chains")
    frame = track.get("frame") or {}
    report = frame.get("superposition") or {}
    for key, value in pair.report.items():
        if report.get(key) != value:
            out.append(f"{where}: superposition {key} is {report.get(key)!r}, the rebuild finds {value!r}")
    for name, audit in (frame.get("audit") or {}).items():
        if audit.get("max_angstrom", 99) > ON_RIBBON + 0.005:
            out.append(f"{where} {name}: a CA {audit.get('max_angstrom')} A off its ribbon")
    if sorted(frame.get("audit") or {}) != sorted(s.name for s in assembly.states):
        out.append(f"{where}: not every state's export was audited")

    centre = np.array(frame["centre_angstrom"])
    length = float(frame["length_angstrom"])
    low, high = np.array(frame["bounds"]["min"]), np.array(frame["bounds"]["max"])
    cas = {s.name: ca_of(pair.lines[s.name]) for s in assembly.states}
    if [c.get("node") for c in track.get("chains", [])] != [u.node for u in assembly.subunits]:
        return out + [f"{where}: subunits {[c.get('node') for c in track.get('chains', [])]}"]
    for chain, subunit in zip(track["chains"], assembly.subunits):
        at = f"{where} {subunit.node}"
        previous = None
        for residue in chain["residues"]:
            if residue["n"] != residue["entry"] + subunit.offset:
                out.append(f"{at} {residue['entry']}: numbered {residue['n']}")
            for state in assembly.states:
                placed = np.array(residue[state.name]["ca"])
                entry = cas[state.name].get((subunit.pdb_chain, residue["entry"]))
                if entry is None or np.abs((entry - centre) / length - placed).max() > 1e-4:
                    out.append(f"{at} {residue['entry']} {state.name}: not the entry's CA in this frame")
                    break
                if (placed < low - 1e-3).any() or (placed > high + 1e-3).any():
                    out.append(f"{at} {residue['entry']} {state.name}: outside the frame's box")
            if previous is not None and residue["entry"] == previous["entry"] + 1:
                for state in assembly.states:
                    step = np.linalg.norm(np.array(residue[state.name]["ca"])
                                          - np.array(previous[state.name]["ca"])) * length
                    if not STEP[0] <= step <= STEP[1]:
                        out.append(f"{at} {previous['entry']}-{residue['entry']} {state.name}: "
                                   f"CA to CA {step:.2f} A")
            previous = residue
        for state in assembly.states:
            if len(chain["haem"].get(state.name) or []) != 3:
                out.append(f"{at} {state.name}: no haem iron")
            bound = chain["oxygen"].get(state.name) or []
            if (state.ligand == "oxy") != bool(bound):
                out.append(f"{at} {state.name}: {len(bound)} oxygen atoms in a {state.ligand} state")

    for state in track["states"]:
        model = DATA / state["model"]["path"]
        if not model.exists() or hashlib.sha256(model.read_bytes()).hexdigest() != state["model"]["sha256"]:
            out.append(f"{where} {state['name']}: the model beside it is not the one it names")
    return out


def main() -> int:
    found = [p for a in ASSEMBLIES for p in problems_of(a)]
    for problem in found:
        print(f"  {problem}", file=sys.stderr)
    print(f"{len(ASSEMBLIES)} assembly morph(s) checked, {len(found)} problem(s)", file=sys.stderr)
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
