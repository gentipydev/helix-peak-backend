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

Which proteins get the track: every target with a structure (`FOLDING_TARGETS`),
recorded here, never in targets.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
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
    chain_residues,
    pdb_path,
    secondary_structure,
)
from pipeline.targets import BY_SLUG, TARGETS, Target  # noqa: E402

# Every protein with a fold gets the track.
FOLDING_TARGETS = tuple(t for t in TARGETS if t.structure is not None)

SCHEMA_VERSION = 1

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


def chain_track(target: Target, node: str, pdb_chain: str, translation: str,
                helices: dict, frame: Frame) -> dict:
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
        },
        "secondary_structure": f"HELIX and SHEET records of {structure.pdb}",
        "chains": [chain_track(target, c.node, c.pdb_chain, translation, helices, frame)
                   for c in structure.chains],
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
    total = 0
    for target in chosen:
        track = bake(target)
        size = (DATA / folding_asset(target)).stat().st_size
        total += size
        states = [r["state"] for c in track["chains"] for r in c["residues"]]
        ss = [r["ss"] for c in track["chains"] for r in c["residues"] if "ss" in r]
        print(f"{target.slug:<16} {len(states):>5} residues: {states.count('ordered'):>4} ordered "
              f"({ss.count('helix'):>3} helix {ss.count('strand'):>3} strand), "
              f"{states.count('disordered'):>3} disordered, {states.count('absent'):>3} absent  "
              f"frame matched to {track['frame']['matched_within']:.0e}  {size:>9,} B")
    print(f"{len(chosen)} tracks, {total:,} B", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
