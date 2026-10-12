"""An experimental entry as the fold of a protein built on demand, with its fold.

The twenty are cut from entries chosen by hand, and each has a `folding` track
baked from the same entry, which the walk's fold page animates and ends on
the model (`folding/README.md`). This makes both for an entry chosen by rule
(`pdb_policy`, `normalize`), so that a protein built on demand folds on its
page exactly as the twenty do.

Both are made from **one working copy** of the entry
(`normalize.working_copy`): one model, one set of atoms, in the record's
numbering. And both are made by the twenty's own bakers, unchanged:

- **The model.** `bake.export` draws the working copy as it draws an entry,
  and `bake.assemble` normalises it; `alphafold.compile_scene` compiles it
  with the app's own importer; the cartoon's sampling is lowered until the
  scene fits what a model made on demand may weigh (`alphafold.BUDGET_BYTES`).
  One node a chain, `chainA` to `chainC`, and `bonds` where bridges are drawn.
  A gap the entry does not place is left open: PyMOL's dashed loop across it
  is not drawn (`GAP_DASHES`).
- **The fold.** `folding/bake_folding.payload`, in the frame of that same
  export: no second export, and nothing to match. It is then held to
  `folding/check_folding`'s whole rule, as every one of the twenty's is.

Then the pair is held to **the app's own rules**, before anything is stored:
`test/fold_bundle_check.dart` in the app's checkout reads the scene, the track
and the row as the walk's fold page will, binds every ring of the model to
its chain and every rod to its bridge, and plays the fold through to the
stored mesh, float for float (`app_check`).

What is refused here is `normalize.Unsuitable`, a documented rule, and the
next candidate is tried (`fold_job.choose`):

- `size`: no sampling of the cartoon fits what a model may weigh;
- `representation`: the cartoon PyMOL draws of it is not one the track can be
  laid on (a CA too far from its ribbon, by `folding/check_folding`'s own
  tolerances, or any other problem that check finds);
- `binding`: the app's check finds a rule of the animation broken, and says
  which in its own words.

PyMOL, the importer or the app's check breaking is raised as it is: no verdict
on the entry, and tried again.

By hand, for a protein the service serves, writing nothing anywhere but `--out`:

    pipeline/structure/venv/bin/python pipeline/structure/experimental.py \\
        --service https://helix-peak-backend.onrender.com --slug ttr --out /tmp/ttr
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

import numpy as np  # noqa: E402

from pipeline.structure import (  # noqa: E402
    alphafold, fold_job, normalize, pdb_policy, pdb_search, pdb_words)
from pipeline.structure.normalize import Examined, Unsuitable  # noqa: E402
from pipeline.structure.pdb_policy import Candidate, Inputs  # noqa: E402

SOURCE = "PDB"
LICENCE = "CC0 1.0"

# The PyMOL the twenty were baked by. A cartoon is PyMOL's own drawing, and
# another version's is another drawing.
PYMOL = "3.1.0"

FRAME_METHOD = ("the structure bake's own export origin and normalisation, of the export the "
                "model was made from")
BUILT_BY = "pipeline/structure/experimental.py"

# Whether PyMOL's dashed loop across a gap of up to ten unplaced residues is
# drawn (`bake.write_pml`). It is not: the loop is PyMOL's curve through space
# the entry has no atoms in, and the handoff's rule is that no connecting loop
# is made up. The ribbon stops at the gap and starts again after it, as
# CFTR's does among the twenty (its gaps are longer than ten). The dashes are
# also not a shape the app's binder lays on the chain where they curve off the
# straight line between the gap's ends: over the 52 proteins of the dry run, 10
# candidates were refused for them and somatotropin lost every one, where
# without them every candidate the rule chose binds (2026-10-12).
GAP_DASHES = False

# The app's check of a baked fold, in the app's checkout, and how long it may
# take: the first run after a change there compiles the app's test.
BUNDLE_CHECK = "test/fold_bundle_check.dart"
CHECK_TIMEOUT = 600


@dataclass(frozen=True)
class Built:
    """One baked fold: the model's two objects, its folding track, and what
    the rows say of them."""

    glb: bytes
    scene: bytes
    folding: bytes
    chrome: dict
    chains: list
    provenance: dict            # the structure track's
    folding_provenance: dict    # the folding track's
    working_copy: str
    audit: dict
    # What making it took, and what the app's check measured of it.
    measured: dict = field(default_factory=dict)


def assign_shape(source: Path, nodes: Sequence[str]) -> Dict[str, Tuple[Tuple[str, int, int], ...]]:
    """PyMOL's own helices and strands for a working copy that names none
    (`dss`), each node's as (kind, first, last)."""
    from pipeline.structure import bake

    script = source.with_name("shape.pml")
    script.write_text("\n".join([
        f"load {source}, src",
        "dss src",
        "python",
        "from pymol import cmd",
        "cmd.iterate('src and name CA', \"print('SHAPE %s %s %s' % (chain, resi, ss or 'L'))\")",
        "python end",
    ]) + "\n")
    said = subprocess.run([bake.PYMOL, "-cq", str(script)], capture_output=True, text=True,
                          timeout=alphafold._PYMOL_TIMEOUT)
    if said.returncode != 0:
        raise RuntimeError(f"pymol could not assign the shape:\n{said.stdout}\n{said.stderr}")
    letters = {letter: node for node, letter in normalize.LETTERS.items() if node in nodes}
    runs: Dict[str, List[Tuple[str, int, int]]] = {node: [] for node in nodes}
    for chain, number, shape in re.findall(r"^SHAPE (\S) (-?\d+) (\S)$", said.stdout, re.M):
        kind = {"H": "helix", "S": "strand"}.get(shape)
        node = letters.get(chain)
        if kind is None or node is None:
            continue
        n = int(number)
        held = runs[node]
        if held and held[-1][0] == kind and held[-1][2] == n - 1:
            held[-1] = (kind, held[-1][1], n)
        else:
            held.append((kind, n, n))
    return {node: tuple(spans) for node, spans in runs.items()}


def _unmeasured(track: dict, glb: bytes) -> int:
    """Take out of a track each ribbon measure its own check would turn down,
    and say how many.

    `bake_folding.ribbon_at` measures the ribbon by the vertex nearest a
    residue's CA, and `check_folding` holds what it found to where the
    twenty's ribbons run: 0.7 A from a helix's CA, 3.4 from a strand's, and
    inside the ribbon. In an entry nobody chose for its geometry a strand's
    first residue or a bulge can lie further than that from where PyMOL
    flattens its sheet, and there the nearest ribbon may be the next strand's.
    A measure that cannot be told from that is not kept: the residue carries
    none, as one `ribbon_at` could not measure carries none, and the app
    smooths the strand there as PyMOL does (`FoldBackbone`)."""
    from pipeline.folding import check_folding
    from pipeline.structure.frame import _nearest
    from pipeline.structure.glb import read_glb

    meshes = read_glb(glb)
    length = track["frame"]["length_angstrom"]
    dropped = 0
    for chain in track["chains"]:
        mesh = meshes[chain["node"]]
        for residue in chain["residues"]:
            ribbon = residue.get("ribbon")
            if ribbon is None:
                continue
            middle = np.array(ribbon[:3], dtype=float)
            from_ca = np.linalg.norm(middle - np.array(residue["ca"], dtype=float)) * length
            inside = np.linalg.norm(middle - _nearest(middle[None], mesh.positions)[0]) * length
            if (from_ca > check_folding.RIBBON_FROM_CA[residue["ss"]]
                    or inside > check_folding.RIBBON_INSIDE):
                del residue["ribbon"]
                dropped += 1
    return dropped


def author_offset(chain: normalize.Drawn) -> Optional[int]:
    """The k that makes the entry's own residue number `n` the precursor's
    `n + k`, where one k does for every residue the construct held; None
    where the author's numbering is not a count beside the precursor's (an
    insertion code, a jump). The whole map is in the audit either way."""
    found = set()
    for residue in chain.residues:
        if not re.fullmatch(r"-?\d+", residue.auth or ""):
            return None
        found.add(residue.number - int(residue.auth))
    return found.pop() if len(found) == 1 else None


def _bake_target(target, inputs: Inputs, candidate: Candidate, examined: Examined, sampling: int,
                 representation: str):
    """The `Target` the twenty's bakers are handed: the protein's own row, with
    the working copy's chains as its structure."""
    from pipeline.targets import Chain, Region, Structure

    regions = tuple(target.regions)
    if target.cleaved and not any(r.kept for r in regions):
        # Nothing named as kept: the products are what the removed regions
        # leave (`pdb_policy.products_of`), and the bakers read them as kept.
        regions += tuple(Region(p.label, "", p.start, p.end) for p in inputs.products)
    chain = examined.chains[0]
    span = (chain.first, chain.last) if candidate.scope == "domain" else None
    structure = Structure(
        pdb=examined.entry,
        chains=tuple(Chain(c.node, normalize.LETTERS[c.node]) for c in examined.chains),
        bonds=bool(examined.bridges), residues=span, sampling=sampling,
        representation=representation)
    return dataclasses.replace(target, structure=structure, regions=regions,
                               aa=len(inputs.translation))


def build(target, inputs: Inputs, candidate: Candidate, examined: Examined, record: bytes,
          entry: pdb_search.Entry, entity: Optional[pdb_search.Entity], workspace: Path,
          app: Path, dart: str = "dart", budget: Optional[int] = None) -> Built:
    """Bake one examined entry: the `.glb`, its compiled scene, the folding
    track, and the fold page's words.

    `target` is the protein's `Target` (the resolver's, rebuilt from its row)
    and `record` its stored record's bytes. Raises `Unsuitable` where the
    entry cannot be drawn by a documented rule."""
    from pipeline.folding import bake_folding, check_folding
    from pipeline.structure import bake

    budget = alphafold.BUDGET_BYTES if budget is None else budget
    workspace.mkdir(parents=True, exist_ok=True)
    source = workspace / f"{examined.entry}.pdb"
    nodes = [c.node for c in examined.chains]

    secondary, shape_by = examined.secondary, examined.secondary_method
    if shape_by is None or not any(secondary.values()):
        # The file assigns no helix or strand to what is drawn. PyMOL assigns
        # its own to a file that names none, as it loads it, so that
        # assignment is read back and written down: the model and the track
        # then share one, whatever it is.
        source.write_text(normalize.working_copy(examined, {}))
        secondary = assign_shape(source, nodes)
        shape_by = f"PyMOL {PYMOL} dss"
    text = normalize.working_copy(examined, secondary)
    source.write_text(text)

    pymol, ordered = alphafold.pymol_reading(source)
    if pymol != PYMOL:
        raise RuntimeError(f"PyMOL here is {pymol}; the twenty's models, and the cartoon the fold "
                           f"animation is measured against, are PyMOL {PYMOL}'s.")
    bake_folding.hold_to_pymol()
    # A chain PyMOL finds no helix or strand in is the tube it is (README.md).
    representation = "cartoon" if ordered else "tube"
    on_ribbon = bake.TUBE_RADIUS if representation == "tube" else alphafold.FRAME_ON
    placed = examined.placed
    record_said = json.loads(record)

    print(f"{target.slug} <- PDB {examined.entry}, {len(nodes)} chain(s), {placed} residues placed",
          flush=True)
    for sampling in range(alphafold.sampling_for(placed), 0, -1):
        drawn = _bake_target(target, inputs, candidate, examined, sampling, representation)
        (workspace / "export").mkdir(exist_ok=True)
        origin, meshes = bake.export(drawn, workspace / "export", source=source,
                                     dashes=GAP_DASHES)
        frames = {}
        for chain in examined.chains:
            ca = np.array([[a.x, a.y, a.z] for r in chain.residues for a in [r.atom("CA")] if a])
            try:
                frames[chain.node] = alphafold.frame_audit(
                    np.asarray(meshes[chain.node].vertices), origin, ca, on_ribbon)
            except ValueError as exc:
                # The export's frame is PyMOL's own origin, and the atoms are
                # the file's: what is off is the cartoon, which leaves some of
                # this chain's residues undrawn or strays from them.
                raise Unsuitable("representation", f"{chain.node}: {exc}") from exc
        # Before `assemble`, which normalises the meshes where they stand.
        centre, extent = bake.normalisation(meshes)
        staged = workspace / f"{target.slug}.glb"
        bake.assemble(target.slug, meshes).export(staged)
        scene = alphafold.compile_scene(staged, workspace / f"{target.slug}.fsceneb", app, dart)
        print(f"  sampling {sampling}: {staged.stat().st_size:,} B, scene {len(scene):,} B",
              flush=True)
        if len(scene) <= budget:
            break
    else:
        raise Unsuitable("size", f"its {placed:,} residues are too large to draw, at the "
                                 f"coarsest cartoon there is")
    glb = staged.read_bytes()

    # The fold, in the frame of the export the model was just made from.
    frame = bake_folding.Frame(centre=origin + centre, length=float(extent.max()),
                               method=FRAME_METHOD, matched=0.0)
    track = bake_folding.payload(drawn, record_said, glb, frame, entry=source, offset=0)
    unmeasured = _unmeasured(track, glb)
    problems = check_folding.problems_of(drawn, track, record=record_said, glb=glb, entry=source,
                                         reach=normalize.BRIDGE_REACH)
    if problems:
        # Kept beside the bake, for whoever asks why: a refusal leaves nothing else.
        (workspace / "folding.refused.json").write_bytes(bake_folding.encode(track))
        raise Unsuitable("representation", "; ".join(problems[:3])
                         + (f"; and {len(problems) - 3} more" if len(problems) > 3 else ""))
    # Checked as the twenty's are, against the working copy. Then each chain
    # is named as the entry names it, with the author's numbering offset
    # where there is one: a working copy's own (its letter, and 0) are no
    # fact about the entry.
    for said, chain in zip(track["chains"], examined.chains):
        if (said["node"], said["first"], said["last"]) != (chain.node, chain.first, chain.last):
            # The bakers read the span from the protein's regions, and the
            # choice from the same regions: two answers is a fault here, and
            # no verdict on the entry.
            raise RuntimeError(
                f"{target.slug}: the folding bake reads {said['node']} as residues "
                f"{said['first']}-{said['last']}, and the choice drew {chain.node} as "
                f"{chain.first}-{chain.last}")
        said["pdb_chain"] = chain.auth
        said["offset"] = author_offset(chain)
    track["secondary_structure"] = f"{shape_by}, {examined.entry}"
    track["built_by"] = BUILT_BY
    folding = bake_folding.encode(track)

    chrome, chains = pdb_words.describe(inputs, candidate, examined, entry, entity)
    states = [r["state"] for c in track["chains"] for r in c["residues"]]
    provenance = {
        "pdb": examined.entry,
        "nodes": [c["node"] for c in chains],
        "source": SOURCE,
        "entry": examined.entry,
        "method": list(examined.methods),
        "resolution": entry.resolution,
        "model": examined.model,
        "revision": entry.revision,
        "revised": entry.revised,
        "file": {"url": pdb_search.UPDATED_CIF.format(entry=examined.entry.lower()),
                 "sha256": examined.sha256},
        "licence": LICENCE,
        "scope": candidate.scope,
        "domains": list(candidate.domains),
        "chains": [{"node": c.node, "product": inputs.products[c.product].label,
                    "entity": c.entity, "asym": c.asym, "auth": c.auth,
                    "first": c.first, "last": c.last,
                    "placed": sum(1 for r in c.residues if r.placed)} for c in examined.chains],
        "residues": {state: states.count(state) for state in ("ordered", "disordered", "absent")},
        "sequence_differences": len(examined.audit.get("mismatches", [])),
        "bridges": [[b.a, b.b] for b in examined.bridges],
        "bridges_dropped": examined.audit.get("bridges_dropped", []),
        "secondary_structure": shape_by,
        "representation": representation,
        "sampling": sampling,
        "gap_dashes": GAP_DASHES,
        "ribbon_unmeasured": unmeasured,
        "frame": frames,
        "policy": pdb_policy.POLICY_VERSION,
        "pymol": pymol,
        "importer": {"package": "flutter_scene", "version": alphafold.FLUTTER_SCENE},
        "built_by": BUILT_BY,
    }
    folding_provenance = {
        "pdb": track["pdb"],
        "chains": [{key: chain[key] for key in ("node", "pdb_chain", "offset", "first", "last")}
                   for chain in track["chains"]],
        "frame": track["frame"],
        "secondary_structure": track["secondary_structure"],
        "schema_version": track["schema_version"],
        "built_by": track["built_by"],
    }
    audit = {
        "entry": examined.entry, "file_sha256": examined.sha256, "model": examined.model,
        "record_sha256": hashlib.sha256(record).hexdigest(),
        "policy": pdb_policy.POLICY_VERSION,
        "candidate": candidate.as_json(),
        "examined": normalize.summary(examined),
        "map": [{"node": c.node, "n": r.number, "label": r.label, "auth": r.auth,
                 "residue": r.name, "letter": r.letter, "placed": r.placed}
                for c in examined.chains for r in c.residues],
        "working_copy_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "tools": {"pymol": pymol, "flutter_scene": alphafold.FLUTTER_SCENE},
    }
    return Built(glb, scene, folding, chrome, chains, provenance, folding_provenance, text, audit)


# ------------------------------------------------------------ the app's check


def served_row(protein: dict) -> dict:
    """A `protein` row as `/protein/{slug}` serves it, as far as the app reads
    it to open a walk (`app/catalog.py`, `ProteinTarget.fromJson`): what the
    app's check is handed, with the fold's own `chains` and `structure` put
    in by `bundle`."""
    structure = protein.get("structure") or {}
    return {
        "slug": protein["slug"], "display": protein["display"], "gene": protein["gene"],
        "uniprot": protein["uniprot"], "accession": protein["accession"],
        "summary": protein["summary"],
        "facts": {key: protein[key] for key in ("residues", "exons", "chains", "bridges")},
        "catalog_order": protein.get("catalog_order"),
        "chain": protein.get("chain_name"),
        "chains": structure.get("chains", []),
        "structure": structure.get("chrome"),
        "tracks": {},
    }


def bundle(built: Built, row: dict, folder: Path) -> Path:
    """What the app's check reads of a baked fold, written into `folder`: the
    protein's row as the service would serve it with this fold published, the
    model, the scene compiled from it, and the track."""
    folder.mkdir(parents=True, exist_ok=True)
    served = {**row, "chains": built.chains, "structure": built.chrome}
    (folder / "row.json").write_text(json.dumps(served, ensure_ascii=False), encoding="utf-8")
    (folder / "model.glb").write_bytes(built.glb)
    (folder / "model.fsceneb").write_bytes(built.scene)
    (folder / "folding.json").write_bytes(built.folding)
    return folder


def app_check(folder: Path, app: Path, flutter: str, timeout: float = CHECK_TIMEOUT) -> dict:
    """Hold a bundle to the app's own rules, and return what the check measured.

    `test/fold_bundle_check.dart` is run from the app's checkout, and writes
    how it ended beside the bundle (`check.json`); that verdict is what is
    read, with the process's exit held to it. `unsuitable` is a rule of the
    animation this entry breaks, raised as `Unsuitable("binding")` in the
    check's own words. A `fault`, no verdict, or a verdict the exit
    contradicts is a `RuntimeError`: no verdict on the entry."""
    verdict = folder / "check.json"
    if verdict.exists():
        verdict.unlink()
    try:
        ran = subprocess.run(
            [flutter, "test", "--no-pub", BUNDLE_CHECK], cwd=str(app),
            env={**os.environ, "FOLD_BUNDLE": str(folder)}, stdin=subprocess.DEVNULL,
            capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"the app's check of the fold had not finished after "
                           f"{timeout:.0f} s") from exc
    except OSError as exc:
        raise RuntimeError(f"the app's check could not be started ({flutter}): {exc}") from exc
    last = " ".join((ran.stdout + " " + ran.stderr).split())[-300:]
    try:
        said = json.loads(verdict.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"the app's check gave no verdict (exit {ran.returncode}): "
                           f"{last}") from exc
    how, failed = said.get("verdict"), said.get("failed") or {}
    why = f"{failed.get('step')}: {' '.join(str(failed.get('detail', '')).split())[:240]}"
    if how == "ok" and ran.returncode == 0:
        return said
    if how == "unsuitable" and ran.returncode != 0:
        raise Unsuitable("binding", why)
    raise RuntimeError(f"the app's check of the fold broke (verdict {how!r}, exit "
                       f"{ran.returncode}): {why if failed else last}")


class Drawer:
    """`fold_job.Draw` for one protein: each candidate whose file passed is
    baked (`build`) and held to the app's check (`app_check`), in a folder of
    its own under `workspace`."""

    def __init__(self, target, inputs: Inputs, record: bytes, row: dict, workspace: Path,
                 app: Path, dart: str, flutter: str, budget: Optional[int] = None):
        self.target, self.inputs, self.record, self.row = target, inputs, record, row
        self.workspace, self.app, self.dart, self.flutter = workspace, app, dart, flutter
        self.budget = budget
        self.drawn = 0

    def __call__(self, candidate: Candidate, examined: Examined, entry: pdb_search.Entry,
                 entity: Optional[pdb_search.Entity]) -> Built:
        self.drawn += 1
        folder = self.workspace / f"{self.drawn:02d}-{examined.entry}"
        started = time.monotonic()
        built = build(self.target, self.inputs, candidate, examined, self.record, entry, entity,
                      folder / "bake", self.app, self.dart, self.budget)
        baked = time.monotonic()
        checked = app_check(bundle(built, self.row, folder / "bundle"), self.app, self.flutter)
        return dataclasses.replace(built, measured={
            "bake_seconds": round(baked - started, 2),
            "check_seconds": round(time.monotonic() - baked, 2),
            "glb_bytes": len(built.glb), "scene_bytes": len(built.scene),
            "folding_bytes": len(built.folding),
            **{key: checked[key] for key in ("vertices", "residues", "bridges", "rings",
                                             "bind_ms", "frames_ms", "slowest_frame_us",
                                             "largest_step", "largest_nudge")
               if key in checked}})


@dataclass
class Made:
    """How making one protein's experimental fold ended: the search, and the
    winner's bake where there is one."""

    inputs: Inputs
    choice: fold_job.Choice
    uniprot_release: Optional[str]

    @property
    def built(self) -> Optional[Built]:
        winner = self.choice.winner
        return None if winner is None else winner.built


def inputs_for(target, record: bytes, *, fetch_entry: Optional[Callable] = None,
               pfam: Optional[Callable] = None) -> Tuple[Inputs, Optional[str]]:
    """What a protein's structure is chosen for (`fold_job.inputs_of`), and the
    UniProt release its domains were read from. `fetch_entry` and `pfam` stand
    in for UniProt and InterPro in a test."""
    from pipeline import uniprot

    entry = (fetch_entry or uniprot.fetch_entry)(target.uniprot)
    translation = json.loads(record)["protein"]["translation"]
    inputs = fold_job.inputs_of(target, translation, entry.body,
                                (pfam or pdb_search.pfam_domains)(target.uniprot))
    return inputs, entry.release


def make(target, record: bytes, row: dict, workspace: Path, app: Path, dart: str, flutter: str,
         *, inputs: Optional[Inputs] = None, release: Optional[str] = None,
         skip: Sequence[str] = (), resume: Optional[Tuple[str, Sequence[str]]] = None,
         budget: Optional[int] = None, limit: int = fold_job.EXAMINE_LIMIT,
         draws: int = fold_job.DRAW_LIMIT, fetch_entry: Optional[Callable] = None,
         pfam: Optional[Callable] = None, **services) -> Made:
    """Search, choose, bake and check one protein's experimental fold.

    `target` is its `Target`, `record` its stored record's bytes, and `row`
    the protein as the service serves it. `inputs` (and the `release` they
    were read from) are `inputs_for`'s, read here where not given. `skip` and
    `resume` are what earlier attempts refused (`fold_job.choose`).
    `services` are `fold_job.choose`'s (`get`, `post`, `now`), for a test to
    stand in for the network; `fetch_entry` and `pfam` likewise, for UniProt
    and InterPro."""
    if inputs is None:
        inputs, release = inputs_for(target, record, fetch_entry=fetch_entry, pfam=pfam)
    drawer = Drawer(target, inputs, record, row, workspace, app, dart, flutter, budget)
    choice = fold_job.choose(inputs, limit=limit, skip=skip, resume=resume, draw=drawer,
                             draws=draws, **services)
    return Made(inputs, choice, release)


# ------------------------------------------------------------ a protein, as a job


def target_json(target) -> dict:
    """A `Target` as a job file carries it, without the structure it has yet
    to be given."""
    said = dataclasses.asdict(dataclasses.replace(target, structure=None))
    return json.loads(json.dumps(said))


def target_of_json(said: dict):
    """The `Target` a job file carries (`target_json`)."""
    from pipeline.targets import Region, Source, Target

    return Target(**{
        **said,
        "source": Source(**said["source"]),
        "structure": None,
        "regions": tuple(Region(**region) for region in said["regions"]),
        "disulfides": tuple(tuple(pair) for pair in said["disulfides"]),
        "uniprot_variants": tuple(tuple(v) for v in said.get("uniprot_variants") or ()),
    })


def target_of_served(protein: dict):
    """The `Target` of a protein built on demand, from what the service serves
    of it: what `resolve.target_of` rebuilds from its row, less the slice of
    the genome its record was cut from, which no bake of a fold reads."""
    from pipeline.targets import Region, Source, Target

    provenance = protein.get("provenance") or {}
    if "regions_from" not in provenance:
        raise ValueError(f"{protein.get('slug')} is not a protein built on demand: its row "
                         f"carries no table of regions")
    return Target(
        slug=protein["slug"], gene=protein["gene"], uniprot=protein["uniprot"],
        display=protein["display"],
        source=Source(protein["accession"], protein_id=protein.get("protein_id"),
                      transcript_id=protein.get("transcript_id")),
        structure=None,
        cleaved=provenance["cleaved"],
        chain_label=provenance["chain_label"],
        regions=tuple(Region(**region) for region in provenance["regions_from"]),
        disulfides=tuple(tuple(pair) for pair in protein.get("disulfides") or ()),
        aa=protein["facts"]["residues"],
        uniprot_variants=tuple(tuple(v) for v in provenance.get("uniprot_variants") or ()),
        mature_peptides=protein.get("mature_peptides", False),
    )


def _read(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def write(made: Made, out: Path) -> dict:
    """Leave what was made in `out`, and return what `outcome.json` says: how
    the search ended, and of a fold that was baked, its files, its words, its
    provenance and what it measured."""
    from pipeline.structure import fold_report

    out.mkdir(parents=True, exist_ok=True)
    choice, built = made.choice, made.built
    said = {
        "outcome": choice.outcome,
        "policy": pdb_policy.POLICY_VERSION,
        "selection": fold_report.selection(made.inputs, choice, made.uniprot_release),
    }
    if built is not None:
        (out / "model.glb").write_bytes(built.glb)
        (out / "model.fsceneb").write_bytes(built.scene)
        (out / "folding.json").write_bytes(built.folding)
        (out / "working_copy.pdb").write_text(built.working_copy, encoding="utf-8")
        audit = {**built.audit, "selection": said["selection"], "measured": built.measured}
        (out / "audit.json").write_text(json.dumps(audit, ensure_ascii=False), encoding="utf-8")
        said.update(
            candidate=fold_job.candidate_id(choice.winner.candidate),
            chrome=built.chrome, chains=built.chains, provenance=built.provenance,
            folding_provenance=built.folding_provenance, measured=built.measured,
            sha256={"glb": hashlib.sha256(built.glb).hexdigest(),
                    "scene": hashlib.sha256(built.scene).hexdigest(),
                    "folding": hashlib.sha256(built.folding).hexdigest()})
    (out / "outcome.json").write_text(json.dumps(said, ensure_ascii=False, indent=1),
                                      encoding="utf-8")
    return said


def main(arguments: List[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--service", help="a protein the service serves, with --slug")
    parser.add_argument("--slug")
    parser.add_argument("--job", type=Path,
                        help="or a job file: a target (`target_json`), the path of its record, "
                             "the row the app is to read, and candidates to skip")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--app", type=Path, default=BACKEND.parent / "helix-peek",
                        help="the app's checkout, for its scene importer and its check")
    parser.add_argument("--dart", default=None)
    parser.add_argument("--flutter", default=None)
    parser.add_argument("--limit", type=int, default=fold_job.EXAMINE_LIMIT)
    parser.add_argument("--draws", type=int, default=fold_job.DRAW_LIMIT)
    parser.add_argument("--describe", nargs=2, metavar=("NAME", "HOW"), default=None,
                        help="also write report.json, this protein's part of a dry run's "
                             "report (scripts/dry_fold.py)")
    args = parser.parse_args(arguments)

    sdk = Path.home() / "flutter" / "bin"
    dart = args.dart or str(sdk / "cache" / "dart-sdk" / "bin" / "dart"
                            if (sdk / "cache" / "dart-sdk" / "bin" / "dart").exists()
                            else shutil.which("dart") or "dart")
    flutter = args.flutter or str(sdk / "flutter" if (sdk / "flutter").exists()
                                  else shutil.which("flutter") or "flutter")
    skip: Sequence[str] = ()
    job: dict = {}
    if args.job:
        job = json.loads(args.job.read_text(encoding="utf-8"))
        target, row, skip = target_of_json(job["target"]), job["row"], job.get("skip") or ()
        record = (args.job.parent / job["record"]).read_bytes()
    elif args.service and args.slug:
        base = args.service.rstrip("/")
        row = json.loads(_read(f"{base}/protein/{args.slug}"))
        tracks = json.loads(_read(f"{base}/protein/{args.slug}/tracks"))
        stored = (tracks.get("tracks") or tracks)["record"]
        record = _read(stored["url"])
        if hashlib.sha256(record).hexdigest() != stored["sha256"]:
            raise SystemExit("The record served is not the one its row names.")
        target = target_of_served(row)
    else:
        parser.error("--service with --slug, or --job")

    started = time.monotonic()
    made = make(target, record, row, args.out / "work", args.app, dart, flutter, skip=skip,
                limit=args.limit, draws=args.draws)
    said = write(made, args.out)
    winner = made.choice.winner
    if args.describe:
        from pipeline.structure import fold_report

        name, how = args.describe
        described = fold_report.describe(name, how, target, made.inputs, made.choice,
                                         made.uniprot_release)
        if job.get("hand_pick"):
            described["historical"] = fold_report._historical(job["hand_pick"],
                                                              made.choice.ranking)
        described["seconds"] = round(time.monotonic() - started, 1)
        (args.out / "report.json").write_text(
            json.dumps(described, ensure_ascii=False, default=str), encoding="utf-8")
    print(f"{target.slug}: {said['outcome']}"
          + (f" {said['candidate']}: {said['chrome']['sentence']}" if winner else ""))
    for examination in made.choice.examinations:
        if not examination.passed:
            print(f"  passed over {fold_job.candidate_id(examination.candidate)}: "
                  f"{examination.refused.reason}, {examination.refused.detail}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
