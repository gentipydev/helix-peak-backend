"""Bake `folding`: each fold's CA trace, residue by residue, in the model's frame.

    .venv/Scripts/python pipeline/fetch_tracks.py --kind structure --kind record
    pipeline/structure/venv/bin/python pipeline/folding/bake_folding.py --all   # PyMOL on PATH
    pipeline/structure/venv/bin/python pipeline/folding/bake_folding.py --target prion

The `structure` track is a baked mesh: it cannot morph, and it carries no
residues. This track is what a fold animation needs instead, per residue of
each exported chain: where its CA sits in the finished fold, whether it sits
in a helix, a strand or a coil, and whether it has a place at all.

**Which residues.** Each chain of the structure, over the mature chain it is
cut from -- the region of the precursor it belongs to -- clipped to the span
the structure exports. Every residue in that stretch is listed once, in the
precursor's numbering, the numbering the walk uses:

- `ordered`: the entry locates its CA. It carries `ca` and `ss`.
- `disordered`: it was in the experiment and has no place in it -- listed in
  REMARK 465, or modelled without its CA -- between residues that do. The
  prion protein's flexible N-terminal half is ninety-three of these.
- `absent`: not in the entry at all, beyond the ends of what was crystallised
  (all of APP past its E1 domain). The track has nothing to say about it.

**Where.** In the stored structure model's own frame -- centred, longest axis
1.0 -- so the last frame of a fold lands on exactly the fold the walk's page
draws. That frame is the structure bake's: this re-runs its own `export` (the
PyMOL export origin and meshes) and `normalisation`, then holds the re-export
to the stored `.glb` vertex for vertex and refuses a frame that does not match
it. So the bake runs in the structure bake's environment, PyMOL and all. The
`.glb`'s sha256 and how closely it matched are recorded with the frame, and CA
coordinates are rounded to 1e-5 model units, a thousandth of an angstrom or
less. (A fit of the CA atoms to the stored ribbon, which is how the AR track
sizes a model with no bridges, is exact where there are bridges and up to
1.5 A out on glucagon's single helix, whose ribbon barely fixes the scale
along its axis.)

**Secondary structure.** The entry's own HELIX and SHEET records, which were
deposited with its coordinates: every helix class is `helix`, every strand of
every sheet `strand`, anything else `coil`. No new source.

**The residue letters** are the record's, what the gene makes. Where the entry
carries another residue at a position -- an engineered mutation, a sequence
conflict, all of them declared in its SEQADV records -- the coordinates are
the entry's and the letter the record's, and `entry_differs` says so.

**The bridges** (schema 2). Where the model draws its disulfides (`bonds`),
each one's atoms as the structure bake builds its rods -- CA, CB, SG, SG, CB,
CA, from the entry's SSBOND records and its own reading of the atoms -- placed
in the same frame. A fold drawn from the track can grow each bridge along the
path the model's rod takes, and end on exactly that rod.

**The ribbon** (schema 2). For each helix and strand residue, where the
model's ribbon passes it and which way the ribbon lies across: the centre and
the widest direction of the ribbon's cross-section nearest the CA, measured on
the stored model itself. A helix's ribbon runs within 0.25 A of its CA atoms,
but PyMOL flattens a sheet, so a strand's can run up to 3 A from them, and its
twist is PyMOL's own; an animation ending on these ends on the model's ribbon
rather than beside it.

**The cartoon** (schema 2). The half-width and half-thickness PyMOL gives a
helix's oval and a strand's slab, the radius of a loop, and the tube a peptide
with no secondary structure is drawn as. They are PyMOL 3.1.0's settings, held
to the running PyMOL by the bake, so a fold drawn from the track can end as
wide and as thick as the model it hands over to.

Which proteins get the track: every target with a structure (`FOLDING_TARGETS`),
recorded here, never in targets.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.paths import DATA  # noqa: E402
from pipeline.structure.glb import Mesh, read_glb  # noqa: E402
from pipeline.structure.pdb import (  # noqa: E402
    Residue,
    atoms,
    chain_residues,
    pdb_path,
    secondary_structure,
    ssbonds,
)
from pipeline.targets import BY_SLUG, TARGETS, Target  # noqa: E402

# Every protein with a fold gets the track.
FOLDING_TARGETS = tuple(t for t in TARGETS if t.structure is not None)

SCHEMA_VERSION = 2

# The cartoon PyMOL draws, in angstroms: each shape's half-width and
# half-thickness, as PyMOL's own settings name them. Measured on the stored
# models, a helix is 2.7 A wide and 0.5 thick and a loop 0.4 across, so these
# are half-extents. `bake` holds them to the PyMOL it runs.
PYMOL_CARTOON = {
    "cartoon_loop_radius": 0.2,
    "cartoon_oval_length": 1.35,
    "cartoon_oval_width": 0.25,
    "cartoon_rect_length": 1.4,
    "cartoon_rect_width": 0.4,
}

# `structure/bake.py`'s, which imports trimesh and cannot be read offline:
# the tube a peptide with no secondary structure is drawn as, and the radius
# of a bridge's rods. `bake` holds these to it too.
TUBE_RADIUS = 0.6
ROD_RADIUS = 0.5

# Measuring the ribbon, in angstroms: how far from a CA its ribbon may be
# looked for, how thick a slice of it is taken across, and how far out from
# its centre (an arrowhead's half-width is 2.4). A slice at least this much
# wider than it is thick is flat enough to have a direction.
_RIBBON_REACH = 3.2
_RIBBON_SLICE = 0.3
_RIBBON_RADIUS = 2.8
_RIBBON_FLAT = 2.5

# Model units, to five places: 1e-5 of CFTR's 103 A is a thousandth of an angstrom.
_PLACES = 5

# How far a re-exported vertex may sit from the stored one, in model units. A
# re-export on another machine's PyMOL build lands within 3e-6 (floating point);
# a different frame, or a different model, is off by orders of magnitude more.
_MATCHED = 1e-5

FRAME_METHOD = ("the structure bake's own export origin and normalisation, re-run "
                "and matched to the stored model vertex for vertex")


@dataclass(frozen=True)
class Frame:
    """`(p - centre) / length` puts a point of the entry into the stored model."""

    centre: np.ndarray      # angstroms, in the PDB frame
    length: float           # angstroms per model unit
    method: str
    matched: float | None   # the largest re-export-to-stored vertex distance, model units


def folding_asset(target: Target) -> str:
    return f"assets/folding/{target.slug}_folding.json"


def structure_frame(target: Target, stored: dict[str, Mesh]) -> Frame:
    """The frame the stored model was cut in, from the structure bake itself.

    Needs PyMOL and trimesh: run it with the structure bake's environment.
    """
    from pipeline.structure import bake as structure_bake

    with tempfile.TemporaryDirectory() as work:
        origin, meshes = structure_bake.export(target, Path(work))
    centre, extent = structure_bake.normalisation(meshes)
    scale = 1.0 / extent.max()
    if list(meshes) != list(stored):
        raise ValueError(f"{target.slug}: re-exported {list(meshes)}, stored {list(stored)}")
    matched = 0.0
    for name, mesh in meshes.items():
        placed = (np.asarray(mesh.vertices) - centre) * scale
        if placed.shape != stored[name].positions.shape:
            raise ValueError(f"{target.slug} {name}: re-exported {len(placed)} vertices, "
                             f"stored {len(stored[name].positions)}")
        matched = max(matched, float(np.abs(placed - stored[name].positions).max()))
    if matched > _MATCHED:
        raise ValueError(f"{target.slug}: the re-export is {matched:.1e} model units "
                         f"from the stored model, not the frame it was cut in")
    return Frame(centre=origin + centre, length=float(extent.max()),
                 method=FRAME_METHOD, matched=matched)


def numbering_offset(residues: list[Residue], translation: str) -> int:
    """The k that makes the entry's residue `number` the precursor's `number + k`.

    The one that matches the most residue letters; of equals, the first, so
    ubiquitin's one repeat is read as the precursor's first. The entry's DBREF
    records are not read for it: two of the twenty (2DN1, 2C9V) number their
    UniProt ends from before the initiator methionine was counted.
    """
    first, last = residues[0].number, residues[-1].number
    best, found = -1, None
    for k in range(1 - first, len(translation) - last + 1):
        matches = sum(1 for r in residues if translation[r.number + k - 1] == r.letter)
        if matches > best:
            best, found = matches, k
    if found is None:
        raise ValueError("the chain is longer than the precursor")
    return found


def mature_span(target: Target, located: list[int]) -> tuple[int, int]:
    """The mature chain the located residues belong to, in precursor numbering.

    A cleaved precursor is cut into its kept regions, so the chain is the kept
    region that holds every located residue. An uncut one is one chain, less
    whatever the table says is removed from its ends (the initiator methionine).
    """
    low, high = min(located), max(located)
    if target.cleaved:
        for region in target.regions:
            if region.kept and region.start <= low and high <= region.end:
                return region.start, region.end
        raise ValueError(f"{target.slug}: no kept region holds residues {low}-{high}")
    start, end = 1, target.aa
    for region in sorted(target.regions, key=lambda r: r.start):
        if not region.kept and region.start == start:
            start = region.end + 1
    for region in sorted(target.regions, key=lambda r: -r.end):
        if not region.kept and region.end == end:
            end = region.start - 1
    return start, end


def ribbon_at(ribbon: np.ndarray, ca: np.ndarray, before: np.ndarray, after: np.ndarray,
              length: float) -> list[float] | None:
    """The ribbon's centre and widest direction where it passes [ca].

    [ribbon] is the chain's vertices in model units, and [before] and [after]
    the CA atoms either side, which give the way the chain runs. The slice of
    ribbon across that way, through the vertex nearest the CA, is measured:
    its mean is where the ribbon is, and its first principal axis the way it
    lies across. None where no ribbon is near, or the slice is not flat.
    """
    along = after - before
    if np.linalg.norm(along) == 0:
        return None
    along = along / np.linalg.norm(along)
    distance = np.linalg.norm(ribbon - ca, axis=1)
    if distance.min() > _RIBBON_REACH / length:
        return None
    nearest = ribbon[distance.argmin()]
    offset = ribbon - nearest
    taken = offset[(np.abs(offset @ along) < _RIBBON_SLICE / length)
                   & (np.linalg.norm(offset, axis=1) < _RIBBON_RADIUS / length)]
    if len(taken) < 8:
        return None
    across = taken - np.outer(taken @ along, along)
    centre = across.mean(axis=0)
    _, _, axes = np.linalg.svd(across - centre)
    spread = (across - centre) @ axes.T
    extent = spread.max(axis=0) - spread.min(axis=0)
    if extent[0] < _RIBBON_FLAT * max(extent[1], 1e-12):
        return None
    widest = axes[0]
    # A direction has no sign; this one is given the sign of its largest part,
    # so that the same model always bakes the same bytes.
    if widest[np.abs(widest).argmax()] < 0:
        widest = -widest
    return ([round(float(v), _PLACES) for v in nearest + centre]
            + [round(float(v), 4) for v in widest])


def chain_track(target: Target, node: str, pdb_chain: str, translation: str,
                helices: dict, frame: Frame, ribbon: np.ndarray | None = None) -> dict:
    structure = target.structure
    residues = chain_residues(structure, pdb_chain)
    if not residues:
        raise ValueError(f"{target.slug}: the entry holds nothing of chain {pdb_chain}")
    k = numbering_offset(residues, translation)
    located = [r.number + k for r in residues if r.ca is not None]
    start, end = mature_span(target, located)
    if structure.residues is not None:
        start, end = max(start, structure.residues[0] + k), min(end, structure.residues[1] + k)
    # The entry's own ends: inside them a residue without a CA was in the
    # experiment and has no place; beyond them it was never there.
    first, last = residues[0].number + k, residues[-1].number + k
    by_position = {r.number + k: r for r in residues}

    out = []
    for n in range(start, end + 1):
        residue = by_position.get(n)
        entry = {"n": n, "aa": translation[n - 1]}
        if residue is not None and residue.ca is not None:
            placed = (residue.ca - frame.centre) / frame.length
            entry["state"] = "ordered"
            entry["ss"] = helices.get((pdb_chain, residue.number), "coil")
            entry["ca"] = [round(float(v), _PLACES) for v in placed]
        elif first <= n <= last:
            entry["state"] = "disordered"
        else:
            entry["state"] = "absent"
        out.append(entry)

    # Where the model's ribbon runs, by each helix and strand residue.
    if ribbon is not None:
        placed = {e["n"]: np.array(e["ca"]) for e in out if e["state"] == "ordered"}
        for e in out:
            if e.get("ss") not in ("helix", "strand"):
                continue
            n = e["n"]
            found = ribbon_at(ribbon, placed[n], placed.get(n - 1, placed[n]),
                              placed.get(n + 1, placed[n]), frame.length)
            if found is not None:
                e["ribbon"] = found

    differs = [
        {"n": r.number + k, "entry": r.letter, "record": translation[r.number + k - 1]}
        for r in residues
        if start <= r.number + k <= end and r.letter != translation[r.number + k - 1]
    ]
    return {
        "node": node,
        "pdb_chain": pdb_chain,
        "offset": k,
        "first": start,
        "last": end,
        "entry_differs": differs,
        "residues": out,
    }


def cartoon_of(target: Target) -> dict:
    """The shapes the model's chains are drawn with, in angstroms."""
    return {
        "representation": target.structure.representation,
        "loop_radius": PYMOL_CARTOON["cartoon_loop_radius"],
        "helix": {"half_width": PYMOL_CARTOON["cartoon_oval_length"],
                  "half_thickness": PYMOL_CARTOON["cartoon_oval_width"]},
        "strand": {"half_width": PYMOL_CARTOON["cartoon_rect_length"],
                   "half_thickness": PYMOL_CARTOON["cartoon_rect_width"]},
        "tube_radius": TUBE_RADIUS,
        "rod_radius": ROD_RADIUS,
    }


def bridge_track(target: Target, chains: list[dict], frame: Frame) -> list[dict]:
    """Each bridge the model draws, as the atoms its rods run through.

    The structure bake's own reading: `ssbonds` for the pairs and `atoms` for
    the coordinates, CA to CB to SG, across, and back. Numbered as the chains
    are, in the precursor, the lower number first. None where the model has
    no `bonds` node.
    """
    structure = target.structure
    if not structure.bonds:
        return []
    path = pdb_path(structure)
    at = atoms(path)
    by_chain = {c["pdb_chain"]: c for c in chains}
    out = []
    for c1, r1, c2, r2 in ssbonds(path, structure):
        ends = [(by_chain[c1], r1, c1), (by_chain[c2], r2, c2)]
        if ends[0][1] + ends[0][0]["offset"] > ends[1][1] + ends[1][0]["offset"]:
            ends.reverse()
        (one, n1, pdb1), (other, n2, pdb2) = ends
        atoms_along = [at[(pdb1, n1, "CA")], at[(pdb1, n1, "CB")], at[(pdb1, n1, "SG")],
                       at[(pdb2, n2, "SG")], at[(pdb2, n2, "CB")], at[(pdb2, n2, "CA")]]
        out.append({
            "a": n1 + one["offset"],
            "a_node": one["node"],
            "b": n2 + other["offset"],
            "b_node": other["node"],
            "path": [[round(float(v), _PLACES) for v in (p - frame.centre) / frame.length]
                     for p in atoms_along],
        })
    return out


def pymol_cartoon() -> dict[str, float]:
    """The cartoon settings of the PyMOL on PATH, read from a session of it."""
    from pipeline.structure import bake as structure_bake

    with tempfile.TemporaryDirectory() as work:
        script = Path(work) / "settings.pml"
        script.write_text("python\nfrom pymol import cmd\n" + "".join(
            f"print('CARTOON_SETTING {name} %.6f' % float(cmd.get({name!r})))\n"
            for name in PYMOL_CARTOON) + "python end\n")
        result = subprocess.run([structure_bake.PYMOL, "-cq", str(script)],
                                capture_output=True, text=True)
    found = dict(re.findall(r"^CARTOON_SETTING (\w+) (-?\d+\.\d+)$", result.stdout, re.M))
    return {name: float(value) for name, value in found.items()}


def hold_to_pymol() -> None:
    """Refuses to bake a cartoon other than the one PyMOL and the bake draw."""
    from pipeline.structure import bake as structure_bake

    running = pymol_cartoon()
    for name, want in PYMOL_CARTOON.items():
        if name not in running or abs(running[name] - want) > 1e-6:
            raise ValueError(f"PyMOL's {name} is {running.get(name)}, the track says {want}")
    for name, ours, theirs in (("TUBE_RADIUS", TUBE_RADIUS, structure_bake.TUBE_RADIUS),
                               ("ROD_RADIUS", ROD_RADIUS, structure_bake.ROD_RADIUS)):
        if ours != theirs:
            raise ValueError(f"{name} is {ours} here and {theirs} in structure/bake.py")


def payload(target: Target, record: dict, glb: bytes, frame: Frame) -> dict:
    """The track for one protein: its stored record and model, in `frame`."""
    structure = target.structure
    translation = record["protein"]["translation"]
    if record.get("gene") != target.gene:
        raise ValueError(f"{target.slug}: the record is {record.get('gene')}'s")
    if len(translation) != target.aa:
        raise ValueError(f"{target.slug}: the record translates {len(translation)} "
                         f"residues, the table says {target.aa}")
    helices = secondary_structure(pdb_path(structure))
    meshes = read_glb(glb)
    stacked = np.vstack([m.positions for m in meshes.values()])
    chains = [chain_track(target, c.node, c.pdb_chain, translation, helices, frame,
                          meshes[c.node].positions if structure.representation == "cartoon"
                          else None)
              for c in structure.chains]
    return {
        "slug": target.slug,
        "gene": target.gene,
        "uniprot": target.uniprot,
        "pdb": structure.pdb,
        "schema_version": SCHEMA_VERSION,
        "frame": {
            "space": "the stored structure model: centred, longest axis 1.0",
            "glb_sha256": hashlib.sha256(glb).hexdigest(),
            "centre_angstrom": [round(float(v), 6) for v in frame.centre],
            "length_angstrom": round(frame.length, 6),
            "method": frame.method,
            "matched_within": None if frame.matched is None else float(f"{frame.matched:.1e}"),
            # The stored model's bounding box, in its own units: what a viewer
            # frames the model by. A fold drawn from this track alone can be
            # framed exactly as the model is, without loading the model.
            "bounds": {"min": [round(float(v), 6) for v in stacked.min(axis=0)],
                       "max": [round(float(v), 6) for v in stacked.max(axis=0)]},
        },
        "secondary_structure": f"HELIX and SHEET records of {structure.pdb}",
        "cartoon": cartoon_of(target),
        "chains": chains,
        "bridges": bridge_track(target, chains, frame),
        "built_by": "pipeline/folding/bake_folding.py",
    }


def encode(track: dict) -> bytes:
    """Indented JSON with each residue on one line of its own."""
    marker = "\u0000residues\u0000"
    shell = {**track, "chains": [{**c, "residues": marker} for c in track["chains"]]}
    text = json.dumps(shell, indent=2, ensure_ascii=False)
    for chain in track["chains"]:
        rows = ",\n".join("        " + json.dumps(r, ensure_ascii=False)
                          for r in chain["residues"])
        text = text.replace(json.dumps(marker), "[\n" + rows + "\n      ]", 1)
    return (text + "\n").encode()


def bake(target: Target) -> dict:
    record_path, model_path = DATA / target.mock_asset, DATA / target.structure_asset
    for path in (record_path, model_path):
        if not path.exists():
            raise FileNotFoundError(
                f"{path}: run `fetch_tracks.py --kind structure --kind record` first")
    glb = model_path.read_bytes()
    frame = structure_frame(target, read_glb(glb))
    track = payload(target, json.loads(record_path.read_bytes()), glb, frame)
    destination = DATA / folding_asset(target)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(encode(track))
    return track


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if not args.all and not args.target:
        parser.error("name --target or --all")
    unknown = sorted(set(args.target or []) - set(BY_SLUG))
    if unknown:
        parser.error(f"no such target: {unknown}")
    chosen = FOLDING_TARGETS if args.all else tuple(BY_SLUG[s] for s in args.target)
    hold_to_pymol()
    total = 0
    for target in chosen:
        track = bake(target)
        size = (DATA / folding_asset(target)).stat().st_size
        total += size
        states = [r["state"] for c in track["chains"] for r in c["residues"]]
        ss = [r["ss"] for c in track["chains"] for r in c["residues"] if "ss" in r]
        print(f"{target.slug:<16} {len(states):>5} residues: {states.count('ordered'):>4} ordered "
              f"({ss.count('helix'):>3} helix {ss.count('strand'):>3} strand), "
              f"{states.count('disordered'):>3} disordered, {states.count('absent'):>3} absent, "
              f"{len(track['bridges'])} bridges  "
              f"frame matched to {track['frame']['matched_within']:.0e}  {size:>9,} B")
    print(f"{len(chosen)} tracks, {total:,} B", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
