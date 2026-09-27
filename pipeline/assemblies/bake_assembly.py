"""Bake an assembly's morph pair: both states, one frame, as models and a CA morph.

    pipeline/structure/venv/bin/python pipeline/assemblies/bake_assembly.py --all   # PyMOL on PATH

For each assembly in `assemblies.py`:

1. `build.py` makes each state's molecule (the second state's tetramer from its
   biological assembly) and superposes the second onto the first.
2. Both are exported by the structure bake's own `export`, with one origin
   between them, so the two exports are one frame; `verify_frame.py` audits
   each against that origin, every CA on its own ribbon.
3. Both are cut to one box (`normalisation` over the two states' meshes
   together, longest axis 1.0) and written as `<state>.glb`.
4. The morph: per subunit, per residue both states place, its CA in each state
   in that same frame, and each state's secondary structure; each subunit's
   haem iron in each state, and the oxygen bound to it where there is any.

Writes under `pipeline/data/assets/assemblies/<slug>/`. Needs PyMOL and trimesh:
run it with the structure bake's environment.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.assemblies.assemblies import ASSEMBLIES, BY_SLUG, Assembly  # noqa: E402
from pipeline.assemblies.build import Pair, build, write  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.structure.glb import Mesh, write_glb  # noqa: E402
from pipeline.structure.pdb import ONE_LETTER, Residue  # noqa: E402
from pipeline.targets import BY_SLUG as TARGETS_BY_SLUG, Chain, Structure  # noqa: E402

OUTPUT = HERE / "output"
SCHEMA_VERSION = 1

# Model units, to five places, as the folding track has them.
_PLACES = 5

# Every CA within this of its own ribbon, in angstroms: `verify_frame.py`'s
# criterion. Hemoglobin is all helix and loop, which PyMOL draws through them.
_ON_RIBBON = 0.26


def morph_asset(assembly: Assembly) -> str:
    return f"assets/assemblies/{assembly.slug}/{assembly.slug}_morph.json"


def model_asset(assembly: Assembly, state: str) -> str:
    return f"assets/assemblies/{assembly.slug}/{state}.glb"


class _Export:
    """What the structure bake's `export` reads off a target: a slug and a structure."""

    def __init__(self, slug: str, structure: Structure):
        self.slug, self.structure = slug, structure


def _structure(assembly: Assembly, pdb: str) -> Structure:
    return Structure(
        pdb=pdb,
        chains=tuple(Chain(u.node, u.pdb_chain) for u in assembly.subunits),
        bonds=False,
        sampling=assembly.sampling,
    )


def _residues(lines: list[str], chain: str) -> dict[int, Residue]:
    """Each residue of a chain with a CA (first alternate location)."""
    out: dict[int, Residue] = {}
    for line in lines:
        if line.startswith("ATOM") and line[21] == chain and line[12:16] == " CA " and line[16] in (" ", "A"):
            number = int(line[22:26])
            out.setdefault(number, Residue(
                number, ONE_LETTER.get(line[17:20].strip(), "X"),
                np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return out


def _helices(lines: list[str]) -> dict[tuple[str, int], str]:
    out = {}
    for line in lines:
        if line.startswith("HELIX "):
            for n in range(int(line[21:25]), int(line[33:37]) + 1):
                out[(line[19], n)] = "helix"
        elif line.startswith("SHEET "):
            for n in range(int(line[22:26]), int(line[33:37]) + 1):
                out[(line[21], n)] = "strand"
    return out


def _group(lines: list[str], chain: str, group: str, atom: str | None = None) -> list[np.ndarray]:
    return [np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            for l in lines
            if l.startswith("HETATM") and l[21] == chain and l[17:20].strip() == group
            and (atom is None or l[12:16].strip() == atom) and l[16] in (" ", "A")]


def _checked_offsets(assembly: Assembly, pair: Pair) -> None:
    """The beta chains' numbering, against the catalog's own HBB record."""
    record_path = DATA / TARGETS_BY_SLUG["hemoglobin"].mock_asset
    if not record_path.exists():
        raise FileNotFoundError(f"{record_path}: run fetch_tracks.py --kind record first")
    translation = json.loads(record_path.read_bytes())["protein"]["translation"]
    for subunit in assembly.subunits:
        if subunit.gene != TARGETS_BY_SLUG["hemoglobin"].gene:
            continue
        for lines in pair.lines.values():
            for number, residue in _residues(lines, subunit.pdb_chain).items():
                if translation[number + subunit.offset - 1] != residue.letter:
                    raise ValueError(f"{subunit.node} {number}: {residue.letter}, the record "
                                     f"has {translation[number + subunit.offset - 1]}")


def bake(assembly: Assembly) -> dict:
    from pipeline.structure import bake as structure_bake
    from pipeline.structure import verify_frame

    first, second = assembly.states
    pair = build(assembly)
    _checked_offsets(assembly, pair)
    workspace = OUTPUT / assembly.slug
    sources = {s.name: write(pair.lines[s.name], workspace / f"{s.name}.pdb")
               for s in assembly.states}

    # One origin between the two exports: the middle of both states' atoms.
    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                    for s in assembly.states for l in pair.lines[s.name]
                    if l.startswith("ATOM") and l[76:78].strip() != "H"])
    origin = (xyz.min(axis=0) + xyz.max(axis=0)) / 2

    exported, audits = {}, {}
    for state in assembly.states:
        target = _Export(assembly.slug, _structure(assembly, state.pdb))
        state_dir = workspace / state.name
        state_dir.mkdir(parents=True, exist_ok=True)
        used, meshes = structure_bake.export(target, state_dir, sources[state.name], origin)
        if not np.allclose(used, origin, atol=1e-5):
            raise ValueError(f"{state.name}: PyMOL exported from {used}, not {origin}")
        exported[state.name] = meshes
        found = verify_frame.audit(
            pdb=str(sources[state.name]), out=str(state_dir), origin=origin,
            chains=tuple((u.pdb_chain, f"{u.node}.obj") for u in assembly.subunits))
        worst = max(row[3] for row in found)
        if worst > _ON_RIBBON:
            raise ValueError(f"{state.name}: a CA {worst:.2f} A off its ribbon: not this frame")
        audits[state.name] = {
            "chains": {chain: {"mean_angstrom": round(mean, 3), "max_angstrom": round(top, 3)}
                       for chain, _, mean, top in found},
            "max_angstrom": round(worst, 3),
        }

    # One box for both: normalise the two states' meshes together.
    together = {f"{s}/{n}": m for s, meshes in exported.items() for n, m in meshes.items()}
    centre, extent = structure_bake.normalisation(together)
    scale = 1.0 / extent.max()
    for state in assembly.states:
        meshes = {
            name: Mesh(positions=(np.asarray(m.vertices) - centre) * scale,
                       normals=np.asarray(m.vertex_normals).copy(),
                       triangles=np.asarray(m.faces))
            for name, m in exported[state.name].items()
        }
        destination = DATA / model_asset(assembly, state.name)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(write_glb(meshes))
    stacked = np.vstack([(np.asarray(m.vertices) - centre) * scale for m in together.values()])

    def place(p: np.ndarray) -> list[float]:
        return [round(float(v), _PLACES) for v in (p - origin - centre) * scale]

    helices = {s.name: _helices(pair.lines[s.name]) for s in assembly.states}
    chains, dropped = [], {}
    for subunit in assembly.subunits:
        by_state = {s.name: _residues(pair.lines[s.name], subunit.pdb_chain) for s in assembly.states}
        both = sorted(set(by_state[first.name]) & set(by_state[second.name]))
        dropped[subunit.node] = sorted(
            (set(by_state[first.name]) | set(by_state[second.name])) - set(both))
        residues = []
        for number in both:
            a, b = by_state[first.name][number], by_state[second.name][number]
            if a.letter != b.letter:
                raise ValueError(f"{subunit.node} {number}: {a.letter} in {first.pdb}, "
                                 f"{b.letter} in {second.pdb}")
            residues.append({
                "n": number + subunit.offset,
                "entry": number,
                "aa": a.letter,
                first.name: {"ca": place(a.ca),
                             "ss": helices[first.name].get((subunit.pdb_chain, number), "coil")},
                second.name: {"ca": place(b.ca),
                              "ss": helices[second.name].get((subunit.pdb_chain, number), "coil")},
            })
        haem, oxygen = {}, {}
        for state in assembly.states:
            irons = _group(pair.lines[state.name], subunit.pdb_chain, "HEM", "FE")
            if len(irons) != 1:
                raise ValueError(f"{subunit.node} {state.name}: {len(irons)} haem irons")
            haem[state.name] = place(irons[0])
            bound = _group(pair.lines[state.name], subunit.pdb_chain, "OXY")
            oxygen[state.name] = [place(p) for p in bound]
        chains.append({
            "node": subunit.node,
            "gene": subunit.gene,
            "uniprot": subunit.uniprot,
            "pdb_chain": subunit.pdb_chain,
            "residues": residues,
            "haem": haem,
            "oxygen": oxygen,
        })

    track = {
        "slug": assembly.slug,
        "display": assembly.display,
        "schema_version": SCHEMA_VERSION,
        "states": [
            {"name": s.name, "pdb": s.pdb, "ligand": s.ligand,
             "built": ("biological assembly %d of the entry" % s.biomt) if s.biomt
             else "the entry's own chains",
             "model": {"path": model_asset(assembly, s.name),
                       "sha256": hashlib.sha256((DATA / model_asset(assembly, s.name)).read_bytes()).hexdigest()}}
            for s in assembly.states
        ],
        "frame": {
            "space": "one frame for both states: centred, longest axis 1.0",
            "shared_with": first.name,
            "superposition": {**pair.report,
                              "rotation": [[round(float(v), 6) for v in row] for row in pair.rotation],
                              "shift_angstrom": [round(float(v), 4) for v in pair.shift]},
            "origin_angstrom": [round(float(v), 6) for v in origin],
            "centre_angstrom": [round(float(v), 6) for v in origin + centre],
            "length_angstrom": round(float(extent.max()), 6),
            "bounds": {"min": [round(float(v), 6) for v in stacked.min(axis=0)],
                       "max": [round(float(v), 6) for v in stacked.max(axis=0)]},
            "audit": audits,
        },
        "dropped": {node: numbers for node, numbers in dropped.items() if numbers},
        "chains": chains,
        "built_by": "pipeline/assemblies/bake_assembly.py",
    }
    destination = DATA / morph_asset(assembly)
    destination.write_bytes((json.dumps(track, indent=1) + "\n").encode())
    return track


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--assembly", action="append", default=None)
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if not args.all and not args.assembly:
        parser.error("name --assembly or --all")
    chosen = ASSEMBLIES if args.all else tuple(BY_SLUG[s] for s in args.assembly)
    for assembly in chosen:
        track = bake(assembly)
        report = track["frame"]["superposition"]
        sizes = {p: (DATA / p).stat().st_size for p in
                 [morph_asset(assembly)] + [model_asset(assembly, s.name) for s in assembly.states]}
        print(f"{assembly.slug}: entries {report['entries_apart_degrees']} deg apart "
              f"(raw CA RMSD {report['raw_rmsd_angstrom']} A); on "
              f"{'+'.join(report['superposed_on'])} {report['held_rmsd_angstrom']} A, the other "
              f"half {report['turned_rmsd_angstrom']} A, turned {report['turned_degrees']} deg; "
              f"tetramer {report['tetramer_rmsd_angstrom']} A (best fit {report['best_fit_rmsd_angstrom']} A)")
        for state, audit in track["frame"]["audit"].items():
            print(f"  {state}: every CA within {audit['max_angstrom']:.2f} A of its ribbon")
        for path, size in sizes.items():
            print(f"  {path}  {size:,} B")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
