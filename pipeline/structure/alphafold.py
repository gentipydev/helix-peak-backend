"""An AlphaFold DB model as the fold of a protein built on demand.

The twenty are cut from experimental entries chosen by hand (`README.md`). A
protein resolved on demand has nobody to choose one, so its fold page draws
AlphaFold DB's model of its UniProt sequence, under rule R5.4: over the mature
span only, and not at all where the model is not confident.

    pipeline/structure/venv/bin/python pipeline/structure/alphafold.py \\
        --service https://helix-peak-backend.onrender.com --slug oca2 --out /tmp/oca2

bakes one by hand from the served row and record, writing nothing anywhere but
`--out`. The resolver's worker calls `build` with the same `Job`.

What it decides, and the rule for each:

- **Entry.** The API answers with the canonical entry and one for each isoform
  it has modelled. The canonical one, `AF-<accession>-F1`, is taken. A protein
  with none (those over 2,700 residues) is refused.
- **Sequence.** The model's residues are held to the record's protein: the
  same length, and no more residues different than the resolver allows the
  record to differ from UniProt by. Otherwise the two are different molecules.
- **Span.** The mature chain or chains: the first kept region to the last, so
  a signal peptide and the propeptides at either end are not drawn. Where a
  precursor is cut into several chains, what lies between them is drawn too,
  because the model is of the one precursor chain, and the caption says so.
- **Gate.** No model where the mean pLDDT over the span is under 50.
- **Colour.** The field's four bands of pLDDT, one node a band: the chain is
  exported whole, exactly as the twenty's are, and its triangles are then
  divided by the band of the residue each lies on. Shown a part at a time,
  PyMOL draws another cartoon, cut short at every break.
- **Bridges.** UniProt's pairs, each drawn only where the model puts the two
  sulfurs within 2.5 A of each other. The model file names none, so the bake
  is handed a copy of it that does.
- **Size.** Cartoon sampling by length, lowered until the compiled scene is
  inside the budget every model is held to, and refused if none is.

A refusal is a `Refused`, whose message is the sentence a reader is shown.
Anything else that goes wrong (the API not answering, PyMOL or the importer
breaking) is raised as it is, for the worker to try again.

The pure half of this (everything above `export_bands`) needs numpy alone, so
it is tested where the service's suite runs. The mesh half needs the structure
bake's environment, PyMOL, and the app's checkout for the scene importer.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Sequence

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

import numpy as np  # noqa: E402

from pipeline.paths import CLIENT  # noqa: E402
from pipeline.structure.pdb import ONE_LETTER, atoms  # noqa: E402

API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
SOURCE = "AlphaFold DB"
LICENCE = "CC BY 4.0"
CITATION = (
    "Jumper et al. 2021, Nature 596:583 (doi:10.1038/s41586-021-03819-2)",
    "Varadi et al. 2024, Nucleic Acids Res. 52:D368 (doi:10.1093/nar/gkad1011)",
)

# AlphaFold DB models no protein longer than this as one chain.
MODEL_LIMIT = 2_700

# R5.4's gate: the mean pLDDT over the span under which nothing is drawn.
GATE = 50.0

# How far apart two sulfurs may be and still be drawn as bonded. A disulfide is
# 2.05 A; the twenty's run 1.98 to 2.18.
BRIDGE_REACH = 2.5

# The four bands, best first, as the node and tint each is drawn as and the
# pLDDT it starts at. The names are the contract with the app's `ChainTint`.
BANDS = (
    ("plddtVeryHigh", 90.0),
    ("plddtConfident", 70.0),
    ("plddtLow", 50.0),
    ("plddtVeryLow", float("-inf")),
)
_SHARES = ("very_high", "confident", "low", "very_low")
BONDS = "bonds"

# The scene importer the app's own flutter_scene ships. A `.fsceneb` is a
# versioned container, so a scene is compiled by the version every stored
# scene was compiled by, and by no other.
FLUTTER_SCENE = "0.23.0"

# What a compiled scene weighs, about, for each residue at each unit of cartoon
# sampling, and what the twenty's scenes are held near.
BYTES_PER_RESIDUE_SAMPLE = 440
TARGET_BYTES = 500_000
MAX_SAMPLING = 8

# The export frame's audit, by `folding/check_folding.py`'s measures: a CA lies
# on its ribbon, within the 0.25 A the twenty's are held to (a tube's radius,
# where the chain is drawn as one), except in a strand, which PyMOL's flat
# arrows smooth past by up to 3 A. A turned or shifted frame puts every atom
# tens of angstroms off.
FRAME_ON = 0.25
FRAME_NEAR = 3.0
_FRAME_SLACK = 0.01

_TIMEOUT = 60
# Seconds waited before each try of a request: a blip of the API's is waited
# out here, rather than spending one of a bake's three tries on it.
_TRIES = (0, 2, 6)
_PYMOL_TIMEOUT = 600
_IMPORT_TIMEOUT = 300

Get = Callable[[str], bytes]


class Refused(Exception):
    """No model is drawn for this protein. The message is the reader's sentence."""


@dataclass(frozen=True)
class Job:
    """What a model is made for: the protein as its row and record have it."""

    slug: str
    accession: str
    display: str
    # The record's own translation, which is what the walk draws.
    protein: str
    # The kept regions of the precursor, 1-based inclusive.
    kept: tuple
    # UniProt's disulfide pairs, in precursor numbering.
    disulfides: tuple = ()
    # How many residues the model may differ from the record's protein at.
    allowed: int = 0


@dataclass(frozen=True)
class Entry:
    """The AlphaFold DB entry a model is read from."""

    accession: str
    entry_id: str
    version: int
    model_url: str
    created: Optional[str] = None


@dataclass(frozen=True)
class Model:
    """The model's one chain: a letter, a pLDDT and a CA for each residue."""

    letters: str
    plddt: np.ndarray
    ca: np.ndarray


@dataclass(frozen=True)
class Confidence:
    """How sure the model is over a span: the mean, and each band's share."""

    mean: float
    shares: tuple


@dataclass(frozen=True)
class Built:
    """One baked model: the two stored objects, and what the rows say of them."""

    glb: bytes
    scene: bytes
    chrome: dict
    chains: list
    provenance: dict


# ------------------------------------------------------------ the entry


def _get(url: str) -> bytes:
    """One address's bytes. An answer that is the server's own trouble (5xx,
    429) or the network's is asked for again, twice; any other answer, a 404
    among them, is the answer."""
    request = urllib.request.Request(url, headers={"User-Agent": "helix-peek-pipeline"})
    for wait in _TRIES:
        time.sleep(wait)
        try:
            with urllib.request.urlopen(request, timeout=_TIMEOUT) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            if exc.code < 500 and exc.code != 429:
                raise
            last = exc
        except OSError as exc:      # URLError, a timeout, a reset connection
            last = exc
    raise last


def no_model(accession: str, length: int) -> Refused:
    if length > MODEL_LIMIT:
        return Refused(f"AlphaFold DB has no model of proteins over {MODEL_LIMIT:,} residues; "
                       f"this one has {length:,}.")
    return Refused(f"AlphaFold DB holds no model of UniProt {accession}.")


def entry_of(accession: str, answer: object, length: int) -> Entry:
    """The canonical entry among those the API answered with.

    The answer lists the canonical sequence's model and one for each isoform
    AlphaFold DB has modelled (`AF-Q9H6X2-6-F1`), each under its own
    `uniprotAccession`.
    """
    wanted = f"AF-{accession}-F1"
    for entry in answer if isinstance(answer, list) else []:
        if not isinstance(entry, dict):
            continue
        if entry.get("entryId") == wanted and entry.get("uniprotAccession") == accession:
            url, version = entry.get("pdbUrl"), entry.get("latestVersion")
            if not isinstance(url, str) or not isinstance(version, int):
                raise ValueError(f"{wanted}: the API names no model file or version")
            return Entry(accession, wanted, version, url, entry.get("modelCreatedDate"))
    raise no_model(accession, length)


def fetch_entry(accession: str, length: int, get: Get = _get) -> Entry:
    """Ask AlphaFold DB for a protein's model. `Refused` where it has none."""
    try:
        answer = json.loads(get(API.format(accession=accession)))
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise no_model(accession, length) from exc
        raise
    return entry_of(accession, answer, length)


def read_model(text: str) -> Model:
    """The model file's chain, residue by residue.

    AlphaFold writes one chain, `A`, numbered from 1 with no gap, and each
    atom's B-factor is its residue's pLDDT. A file that is anything else is
    not read as one.
    """
    letters, plddt, ca = [], [], []
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM") or line[12:16] != " CA ":
            continue
        if line[21] != "A" or line[26].strip() or line[16] not in (" ", "A"):
            raise ValueError("the model is not one plain chain A")
        if int(line[22:26]) != len(letters) + 1:
            raise ValueError(f"the model's residue {len(letters) + 1} is numbered {int(line[22:26])}")
        letters.append(ONE_LETTER.get(line[17:20].strip(), "X"))
        plddt.append(float(line[60:66]))
        ca.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not letters:
        raise ValueError("the model file holds no residues")
    return Model("".join(letters), np.array(plddt), np.array(ca))


def check_sequence(accession: str, model: Model, protein: str, allowed: int) -> int:
    """Hold the model to the record's protein. Returns how many residues differ."""
    if len(model.letters) != len(protein):
        raise Refused(f"AlphaFold's model of {accession} is {len(model.letters):,} residues long; "
                      f"this record's protein is {len(protein):,}.")
    differ = sum(1 for a, b in zip(model.letters, protein) if a != b)
    if differ > allowed:
        raise Refused(f"AlphaFold's model of {accession} differs from this record's protein at "
                      f"{differ:,} residues, more than the {allowed:,} allowed.")
    return differ


# ------------------------------------------------------------ what is drawn


def span_of(kept: Sequence, length: int) -> tuple:
    """The first kept region to the last, or the whole chain where none is named."""
    if not kept:
        return 1, length
    first, last = min(start for start, _ in kept), max(end for _, end in kept)
    if not 1 <= first <= last <= length:
        raise ValueError(f"kept regions {first}-{last} lie outside a protein of {length}")
    return first, last


def band_of(plddt: np.ndarray) -> np.ndarray:
    """Each residue's band, as an index into `BANDS`."""
    band = np.full(len(plddt), len(BANDS) - 1)
    for index in range(len(BANDS) - 2, -1, -1):
        band[plddt >= BANDS[index][1]] = index
    return band


def confidence(model: Model, span: tuple) -> Confidence:
    scores = model.plddt[span[0] - 1:span[1]]
    counts = np.bincount(band_of(scores), minlength=len(BANDS))
    return Confidence(float(scores.mean()), tuple(float(c) / len(scores) for c in counts))


def gate(said: Confidence, span: tuple) -> None:
    if said.mean < GATE:
        raise Refused(f"AlphaFold's model is not confident here: mean pLDDT {said.mean:.1f} over "
                      f"residues {span[0]:,}–{span[1]:,}, under the {GATE:.0f} needed to draw it.")


def bridges_of(pairs: Sequence, at: dict, span: tuple) -> tuple:
    """The pairs the model bears out, and each it does not, with why.

    `at` is the bake's own reading of the model (`pdb.atoms`).
    """
    drawn, dropped = [], []
    for first, second in sorted(tuple(sorted(pair)) for pair in pairs):
        if not all(span[0] <= number <= span[1] for number in (first, second)):
            dropped.append({"pair": [first, second], "reason": "outside the span drawn"})
            continue
        ends = [at.get(("A", number, "SG")) for number in (first, second)]
        if ends[0] is None or ends[1] is None:
            dropped.append({"pair": [first, second], "reason": "no sulfur in the model"})
            continue
        apart = float(np.linalg.norm(ends[0] - ends[1]))
        if apart > BRIDGE_REACH:
            dropped.append({"pair": [first, second], "reason": "sulfurs apart in the model",
                            "distance": round(apart, 2)})
            continue
        drawn.append((first, second))
    return drawn, dropped


def working_copy(text: str, span: tuple, bridges: Sequence, path: Path) -> Path:
    """The model as the bake is handed it: the span's atoms, and its bridges
    as the `SSBOND` records `pdb.ssbonds` reads an entry's from."""
    lines = [f"SSBOND {n:3d} CYS A {first:4d}    CYS A {second:4d}"
             for n, (first, second) in enumerate(bridges, 1)]
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM") and span[0] <= int(line[22:26]) <= span[1]:
            lines.append(line)
    path.write_text("\n".join(lines + ["TER", "END"]) + "\n")
    return path


def sampling_for(residues: int) -> int:
    """PyMOL's cartoon sampling for a span this long: the twenty's 8 where it
    fits, less for a longer chain, to hold the scene near half a megabyte."""
    return max(1, min(MAX_SAMPLING, TARGET_BYTES // (BYTES_PER_RESIDUE_SAMPLE * residues)))


def percentages(shares: Sequence) -> list:
    """Shares as whole per cents that add to a hundred (largest remainders)."""
    exact = [share * 100 for share in shares]
    whole = [int(value) for value in exact]
    by_remainder = sorted(range(len(exact)), key=lambda i: (whole[i] - exact[i], i))
    for index in by_remainder[:100 - sum(whole)]:
        whole[index] += 1
    return whole


def describe(job: Job, entry: Entry, span: tuple, said: Confidence, bridges: Sequence,
             nodes: Sequence) -> tuple:
    """The fold page's words and its chains, from the model's own numbers.

    The seven keys the twenty's rows carry (`StructureChrome`), templated
    rather than written (R4.5): nothing here knows which protein it is.
    """
    first, last = span
    per_cents = percentages(said.shares)
    very_high, confident = per_cents[0], per_cents[1]
    # Only the bands the model has, as the legend shows only those.
    bands = ", ".join(
        f"{per_cent or 'under 1'}% {name}"
        for per_cent, name, share in zip(per_cents, ("very high", "confident", "low", "very low"),
                                         said.shares)
        if share > 0)
    if len(job.kept) > 1:
        sentence = (f"AlphaFold's precursor, mean pLDDT {said.mean:.1f}: "
                    f"every chain and what lies between them.")
    else:
        sentence = (f"AlphaFold prediction, mean pLDDT {said.mean:.1f}. "
                    f"{very_high + confident}% of residues at 70 or over.")
    drawn = ""
    if bridges:
        drawn = (" Its disulfide bridge is drawn." if len(bridges) == 1
                 else f" Its {len(bridges)} disulfide bridges are drawn.")
    semantics = (f"AlphaFold's predicted fold of {job.display}, residues {first:,} to {last:,}, "
                 f"coloured by confidence: {bands}.{drawn} Drag to turn it.")
    chrome = {
        "pdb": entry.entry_id,
        "modelled": None if (first, last) == (1, len(job.protein)) else [first, last],
        "label": "the predicted fold",
        "count": last - first + 1,
        "unit": "residues",
        "sentence": sentence,
        "semantics": semantics,
    }
    chains = [{"node": node, "tint": "cysteine" if node == BONDS else node} for node in nodes]
    return chrome, chains


# ------------------------------------------------------------ the mesh


def residue_of(points: np.ndarray, ca: np.ndarray) -> np.ndarray:
    """Which residue each point of a ribbon lies on, as an index into `ca`.

    The nearer end of the nearest stretch between two consecutive CA atoms: a
    residue's ribbon runs from half-way to the residue before it to half-way
    to the one after, which is where PyMOL changes a cartoon's colour. Nearest
    CA alone gives a strand's residues to their neighbours, because a flat
    arrow is smoothed past the atoms it is drawn for.
    """
    from scipy.spatial import cKDTree

    if len(ca) == 1:
        return np.zeros(len(points), dtype=int)
    _, near = cKDTree(ca).query(points, k=min(4, len(ca)))
    near = near.reshape(len(points), -1)
    best = np.full(len(points), np.inf)
    owner = np.zeros(len(points), dtype=int)
    for column in range(near.shape[1]):
        for shift in (-1, 0):
            start = np.clip(near[:, column] + shift, 0, len(ca) - 2)
            a, along = ca[start], ca[start + 1] - ca[start]
            t = np.clip(((points - a) * along).sum(axis=1) / (along * along).sum(axis=1), 0, 1)
            away = np.linalg.norm(points - (a + t[:, None] * along), axis=1)
            better = away < best
            best[better] = away[better]
            owner[better] = (start + (t >= 0.5))[better]
    return owner


def split_by_band(mesh, origin: np.ndarray, ca: np.ndarray, band: np.ndarray) -> dict:
    """A chain's cartoon as one mesh a band, in `BANDS`' order.

    Every triangle goes to the band of the residue its centre lies on, with
    its vertices and PyMOL's normals as they were: together the meshes are the
    cartoon, triangle for triangle. A band no residue is in has no mesh.
    """
    import trimesh

    vertices, faces = np.asarray(mesh.vertices), np.asarray(mesh.faces)
    normals = np.asarray(mesh.vertex_normals)
    owner = residue_of(vertices[faces].mean(axis=1) + origin, ca)   # export frame -> PDB frame
    of_face = band[owner]
    meshes = {}
    for index, (name, _) in enumerate(BANDS):
        chosen = faces[of_face == index]
        if not len(chosen):
            continue
        used, renumbered = np.unique(chosen, return_inverse=True)
        meshes[name] = trimesh.Trimesh(vertices=vertices[used], faces=renumbered.reshape(-1, 3),
                                       vertex_normals=normals[used], process=False)
    if sum(len(m.faces) for m in meshes.values()) != len(faces):
        raise ValueError("the bands do not add up to the cartoon")
    return meshes


def frame_audit(vertices: np.ndarray, origin: np.ndarray, ca: np.ndarray,
                on: float = FRAME_ON) -> dict:
    """`verify_frame.py`'s measure, each CA to the nearest point of its
    exported ribbon, for a model of any length: none further than
    `FRAME_NEAR`, and a quarter of them at least within `on`, which a chain
    drawn as a tube widens to the tube's radius."""
    from scipy.spatial import cKDTree

    away, _ = cKDTree(vertices + origin).query(ca)          # export frame -> PDB frame
    print(f"  chain A: {len(ca)} CA vs {len(vertices)} verts | min {away.min():.2f} A  "
          f"mean {away.mean():.2f} A  max {away.max():.2f} A", flush=True)
    on = float((away <= on + _FRAME_SLACK).mean())
    if away.max() > FRAME_NEAR or on < 0.25:
        raise ValueError(f"the export frame is off: CA atoms up to {away.max():.2f} A from the "
                         f"ribbon, {on:.0%} of them on it")
    return {"ca_mean": round(float(away.mean()), 3), "ca_max": round(float(away.max()), 3),
            "ca_on_ribbon": round(on, 3)}


def _target(job: Job, entry: Entry, sampling: int, bonds: bool, representation: str):
    from pipeline.targets import Chain, Source, Structure, Target

    return Target(
        slug=job.slug, gene=job.slug.upper(), uniprot=job.accession, display=job.display,
        source=Source(entry.entry_id),
        structure=Structure(pdb=entry.entry_id, chains=(Chain("chainA", "A"),), bonds=bonds,
                            sampling=sampling, representation=representation),
        cleaved=False, chain_label=job.display,
    )


def pymol_reading(source: Path) -> tuple:
    """PyMOL's version, and how many residues of `source` it draws as helix
    or strand. It assigns them itself, on loading: the model file names none.

    Asked in a script, as the bake asks for its export origin: on the command
    line PyMOL takes the per cent signs out of a format string.
    """
    from pipeline.structure import bake

    script = source.with_name("reading.pml")
    script.write_text("\n".join([
        f"load {source}, src",
        "python",
        "from pymol import cmd",
        "print('PYMOL_VERSION %s' % cmd.get_version()[0])",
        "print('PYMOL_ORDERED %d' % cmd.count_atoms('src and name CA and (ss h or ss s)'))",
        "python end",
    ]) + "\n")
    said = subprocess.run([bake.PYMOL, "-cq", str(script)], capture_output=True, text=True,
                          timeout=_PYMOL_TIMEOUT)
    version = re.search(r"^PYMOL_VERSION (\S+)$", said.stdout, re.M)
    ordered = re.search(r"^PYMOL_ORDERED (\d+)$", said.stdout, re.M)
    if said.returncode != 0 or not version or not ordered:
        raise RuntimeError(f"pymol could not read the model:\n{said.stdout}\n{said.stderr}")
    return version.group(1), int(ordered.group(1))


def export_bands(job: Job, entry: Entry, source: Path, workspace: Path, sampling: int,
                 ca: np.ndarray, band: np.ndarray, bonds: bool, representation: str) -> tuple:
    """The span's cartoon by band, and its bridges, in the bake's export frame.

    `bake.export` exactly as it runs for the twenty, with the model as its
    `source`. Returns the meshes by node, in the order the model holds them,
    and the frame's audit.
    """
    from pipeline.structure import bake

    origin, exported = bake.export(_target(job, entry, sampling, bonds, representation),
                                   workspace, source=source)
    chain = exported["chainA"]
    audit = frame_audit(np.asarray(chain.vertices), origin, ca,
                        bake.TUBE_RADIUS if representation == "tube" else FRAME_ON)
    meshes = split_by_band(chain, origin, ca, band)
    if BONDS in exported:
        meshes[BONDS] = exported[BONDS]
    return meshes, audit


# ------------------------------------------------------------ the scene


def importer_version(app: Path) -> Optional[str]:
    """The flutter_scene the app's lockfile pins, which is the importer it runs."""
    lock = app / "pubspec.lock"
    if not lock.exists():
        return None
    found = re.search(r'^  flutter_scene:\n(?:    .*\n)*?    version: "([^"]+)"',
                      lock.read_text(), re.M)
    return found.group(1) if found else None


def importer_unready(app: Path) -> Optional[str]:
    """Why no scene could be compiled from this checkout. None when one could."""
    version = importer_version(app)
    if version is None:
        return f"{app} is not the app's checkout, or its pubspec.lock names no flutter_scene."
    if version != FLUTTER_SCENE:
        return (f"The app's checkout pins flutter_scene {version}; every stored scene is "
                f"compiled by {FLUTTER_SCENE}.")
    return None


def compile_scene(glb: Path, out: Path, app: Path, dart: str) -> bytes:
    """The `.fsceneb` the app loads, compiled by the app's own importer."""
    unready = importer_unready(app)
    if unready:
        raise RuntimeError(unready)
    done = subprocess.run([dart, "run", "flutter_scene:import", "-i", str(glb), "-o", str(out)],
                          cwd=str(app), capture_output=True, text=True, timeout=_IMPORT_TIMEOUT)
    if done.returncode != 0 or not out.exists() or not out.stat().st_size:
        raise RuntimeError(f"the scene importer failed ({done.returncode}):\n"
                           f"{done.stdout}\n{done.stderr}")
    return out.read_bytes()


def build(job: Job, workspace: Path, app: Path, dart: str = "dart", get: Get = _get,
          budget: Optional[int] = None) -> Built:
    """Make one protein's model: the `.glb`, its compiled scene, and its words.

    Raises `Refused` with the reader's sentence where no model is drawn.
    """
    entry = fetch_entry(job.accession, len(job.protein), get)
    raw = get(entry.model_url)
    text = raw.decode("ascii", errors="replace")
    model = read_model(text)
    differ = check_sequence(job.accession, model, job.protein, job.allowed)
    span = span_of(job.kept, len(job.protein))
    said = confidence(model, span)
    gate(said, span)

    from pipeline.structure import bake

    budget = bake.FSCENEB_BUDGET_BYTES if budget is None else budget
    workspace.mkdir(parents=True, exist_ok=True)
    whole = workspace / f"{entry.entry_id}.pdb"
    whole.write_bytes(raw)
    bridges, dropped = bridges_of(job.disulfides, atoms(whole), span)
    source = working_copy(text, span, bridges, workspace / "span.pdb")
    pymol, ordered = pymol_reading(source)
    # A chain PyMOL finds no helix or strand in is drawn as the tube it is: a
    # cartoon of it is a hairline (README.md, on oxytocin's nine residues).
    representation = "cartoon" if ordered else "tube"
    ca = model.ca[span[0] - 1:span[1]]
    band = band_of(model.plddt[span[0] - 1:span[1]])
    residues = span[1] - span[0] + 1

    print(f"{job.slug} <- {entry.entry_id} v{entry.version}, residues {span[0]}-{span[1]}, "
          f"mean pLDDT {said.mean:.1f}", flush=True)
    for sampling in range(sampling_for(residues), 0, -1):
        meshes, audit = export_bands(job, entry, source, workspace, sampling, ca, band,
                                     bool(bridges), representation)
        nodes = list(meshes)
        staged = workspace / f"{job.slug}.glb"
        bake.assemble(job.slug, meshes).export(staged)
        scene = compile_scene(staged, workspace / f"{job.slug}.fsceneb", app, dart)
        print(f"  sampling {sampling}: {staged.stat().st_size:,} B, scene {len(scene):,} B",
              flush=True)
        if len(scene) <= budget:
            break
    else:
        raise Refused(f"AlphaFold's model of {residues:,} residues is too large to draw.")

    glb = staged.read_bytes()
    chrome, chains = describe(job, entry, span, said, bridges, nodes)
    provenance = {
        "pdb": entry.entry_id,
        "nodes": nodes,
        "source": SOURCE,
        "entry": entry.entry_id,
        "model_version": entry.version,
        "model_url": entry.model_url,
        "model_sha256": hashlib.sha256(raw).hexdigest(),
        "model_created": entry.created,
        "licence": LICENCE,
        "citation": list(CITATION),
        "span": [span[0], span[1]],
        "mean_plddt": round(said.mean, 2),
        "plddt_shares": {name: round(share, 4) for name, share in zip(_SHARES, said.shares)},
        "gate": GATE,
        "sequence_differences": differ,
        "bridges": [list(pair) for pair in bridges],
        "bridges_dropped": dropped,
        "representation": representation,
        "sampling": sampling,
        "frame": audit,
        "pymol": pymol,
        "importer": {"package": "flutter_scene", "version": FLUTTER_SCENE},
        "built_by": "pipeline/structure/alphafold.py",
    }
    return Built(glb, scene, chrome, chains, provenance)


# ------------------------------------------------------------ by hand


def job_from_service(service: str, slug: str, allowed: int, get: Get = _get) -> Job:
    """A `Job` read from what the service serves of a protein: its row, and
    the record its record track names."""
    base = service.rstrip("/")
    row = json.loads(get(f"{base}/protein/{slug}"))
    tracks = json.loads(get(f"{base}/protein/{slug}/tracks"))
    record = json.loads(get(tracks.get("tracks", tracks)["record"]["url"]))
    return Job(
        slug=slug, accession=row["uniprot"], display=row["display"],
        protein=record["protein"]["translation"],
        kept=tuple((r["start"], r["end"]) for r in row["regions"] if r.get("kept")),
        disulfides=tuple(tuple(pair) for pair in row["disulfides"]),
        allowed=allowed,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--service", required=True, help="the service to read the row from")
    parser.add_argument("--slug", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--app", type=Path, default=CLIENT,
                        help="the app's checkout, for its scene importer")
    parser.add_argument("--dart", default="dart")
    parser.add_argument("--allowed", type=int, default=0,
                        help="residues the model may differ from the record's protein at")
    args = parser.parse_args()
    try:
        built = build(job_from_service(args.service, args.slug, args.allowed), args.out,
                      args.app, args.dart)
    except Refused as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 3
    (args.out / "described.json").write_text(json.dumps(
        {"chrome": built.chrome, "chains": built.chains, "provenance": built.provenance},
        indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(built.chrome, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
