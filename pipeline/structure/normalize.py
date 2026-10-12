"""An entry chosen by rule, read once from its own file and held to the protein.

`pdb_policy.rank` chooses from what RCSB's data API says of each entry. This
holds the best of those choices to the entry's own SIFTS-enriched mmCIF
(`pdb_search.updated_cif`) before anything is drawn from it, and reads, once,
everything both bakes need:

- **One model.** An NMR entry's representative conformer where it names one,
  else its first; never an average. A crystal or EM entry's first.
- **One set of atoms.** Of alternate locations, per residue, the one with the
  highest mean occupancy, then the first label: a residue is never half one
  conformation and half another. Hydrogens are dropped. A position the entry
  holds two residue types at (microheterogeneity) is refused.
- **One map to the protein.** Each residue by its label number, through
  SIFTS's own mapping in the file, to UniProt, which numbers the record's
  protein. The stretch drawn is one contiguous segment of that mapping, so an
  insertion, a deletion or a fusion partner is never stitched across.
- **The letters.** Each placed residue's, read as the gene would make it (a
  modified residue as its parent), against the record's.
- **One chain, unbroken.** Each pair of placed residues no further apart than
  the residues between them could reach: a chain break inside the stretch is
  refused, never bridged.
- **Chains that belong together.** Several products are drawn from one entry
  only where they touch there.
- **Bridges and shape.** The entry's disulfides between residues drawn, each
  only where its sulfurs are within 2.5 A of each other; and its helices and
  strands as the file assigns them, which in PDBe's file is PDBe's own DSSP
  run and not the depositors' records (a helix can be a residue shorter at
  either end than the legacy file's).

A refusal is `Unsuitable`, a documented rule with its reason, and the next
candidate is examined. A file that cannot be read is not a verdict on the
entry; that is raised as it is.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from pipeline.structure import mmcif
from pipeline.structure.pdb import ONE_LETTER
from pipeline.structure.pdb_policy import Candidate, Inputs, MISMATCH_SHARE, Pick

# Two cysteines are bridged where their sulfurs are this close: a disulfide is
# 2.05 A, and the twenty's run 1.98 to 2.18 (`alphafold.BRIDGE_REACH`).
BRIDGE_REACH = 2.5

# CA to CA, in angstroms, as `folding/check_folding.py` holds the twenty: 3.8
# for a trans peptide bond, 2.9 for a cis one, and room for a poorly resolved
# one; across residues with no place, no more than 3.8 a residue.
STEP = (2.8, 4.4)
PER_RESIDUE = 3.8
_SLACK = 1e-6       # a chain drawn fully extended is at the limit, not past it

# Two chains touch where any of their heavy atoms are this close.
CONTACT = 4.5

# What PDB format, which the bakes read the entry in, can number.
MAX_NUMBER = 9_999
MAX_ATOMS = 99_999

NODES = ("chainA", "chainB", "chainC")


class Unsuitable(ValueError):
    """The entry cannot be drawn by a documented rule; `reason` is its code."""

    def __init__(self, reason: str, detail: str):
        super().__init__(f"{reason}: {detail}")
        self.reason = reason
        self.detail = detail


@dataclass(frozen=True)
class Atom:
    name: str
    element: str
    x: float
    y: float
    z: float
    occupancy: float
    b: float
    hetero: bool


@dataclass(frozen=True)
class Residue:
    """One residue of a chain drawn, numbered as the precursor is."""

    number: int                 # the precursor's
    label: int                  # the entry's label_seq_id
    auth: str                   # the entry's author number and insertion code
    name: str                   # as the entry names it (MSE, SEP, ...)
    letter: str                 # what the entry holds, a modified residue as its parent
    parent: Optional[str]       # the standard residue a modified one is, where it is one
    atoms: Tuple[Atom, ...]
    alternate: Optional[str] = None
    record: str = ""            # the record's letter here: what the gene makes

    def atom(self, name: str) -> Optional[Atom]:
        return next((a for a in self.atoms if a.name == name), None)

    @property
    def placed(self) -> bool:
        return self.atom("CA") is not None


@dataclass(frozen=True)
class Drawn:
    """One chain of the model: a product, cut from one chain of the entry.

    `span` is what the fold page counts and the folding track lists, as the
    twenty's do (R5.5): a mature product whole, whatever of it the construct
    held, or the stretch drawn around a domain. `residues` are the
    construct's residues inside it, placed or not; a residue of the span the
    construct never held is in neither, and the track calls it absent."""

    node: str
    product: int
    entity: str
    asym: str
    auth: str
    span: Tuple[int, int]
    residues: Tuple[Residue, ...]

    @property
    def first(self) -> int:
        return self.span[0]

    @property
    def last(self) -> int:
        return self.span[1]


@dataclass(frozen=True)
class Bridge:
    a: int
    a_node: str
    b: int
    b_node: str
    distance: float


@dataclass
class Examined:
    """An entry as both bakes read it: one model, the chains drawn, their
    bridges and shape, and what was chosen and left out on the way."""

    entry: str
    sha256: str
    model: str
    methods: Tuple[str, ...]
    chains: Tuple[Drawn, ...]
    bridges: Tuple[Bridge, ...]
    secondary: Dict[str, Tuple[Tuple[str, int, int], ...]]
    secondary_method: Optional[str]
    audit: Dict[str, object] = field(default_factory=dict)

    @property
    def placed(self) -> int:
        return sum(1 for c in self.chains for r in c.residues if r.placed)

    @property
    def spanned(self) -> int:
        return sum(c.last - c.first + 1 for c in self.chains)


# ------------------------------------------------------------ reading the file


def _int(value: Optional[str]) -> Optional[int]:
    try:
        return int(value) if value is not None else None
    except ValueError:
        return None


def _float(value: Optional[str], default: float = 0.0) -> float:
    try:
        return float(value) if value is not None else default
    except ValueError:
        return default


def _model(block: mmcif.Block, nmr: bool) -> str:
    """The model drawn: an NMR entry's representative conformer, else the first."""
    models = sorted({m for m in block["atom_site"].get("pdbx_PDB_model_num") if m is not None},
                    key=lambda m: (_int(m) is None, _int(m) or 0, m))
    if not models:
        raise Unsuitable("model", "the file holds no atoms")
    if nmr and "pdbx_nmr_representative" in block:
        chosen = block["pdbx_nmr_representative"].value("conformer_id")
        if chosen in models:
            return chosen
    return models[0]


def _sifts(block: mmcif.Block, asym: str) -> Dict[int, Tuple[str, int, str]]:
    """label_seq_id -> (UniProt accession, number, residue), from SIFTS's rows
    in the file. A residue mapped twice must be mapped the same way twice."""
    if "pdbx_sifts_xref_db" not in block:
        return {}
    table = block["pdbx_sifts_xref_db"]
    out: Dict[int, Tuple[str, int, str]] = {}
    for row_asym, seq, acc, num, res in table.rows("asym_id", "seq_id", "unp_acc", "unp_num",
                                                      "unp_res"):
        if row_asym != asym or acc is None:
            continue
        label, number = _int(seq), _int(num)
        if label is None or number is None:
            continue
        mapped = (acc, number, res or "")
        held = out.get(label)
        if held is not None and held != mapped:
            raise Unsuitable("mapping", f"SIFTS maps label {label} of chain {asym} two ways")
        out[label] = mapped
    return out


def _scheme(block: mmcif.Block, asym: str) -> Dict[int, Tuple[str, str]]:
    """label_seq_id -> (residue name, author number + insertion code), for every
    residue of the chain's sequence. Two residue types at one position is
    microheterogeneity, which no single model can draw."""
    table = block.get("pdbx_poly_seq_scheme")
    if table is None:
        raise Unsuitable("file", "the file has no polymer sequence scheme")
    out: Dict[int, Tuple[str, str]] = {}
    for row_asym, seq, mon, auth, ins in table.rows("asym_id", "seq_id", "mon_id",
                                                    "pdb_seq_num", "pdb_ins_code"):
        if row_asym != asym:
            continue
        label = _int(seq)
        if label is None:
            continue
        if label in out and out[label][0] != mon:
            raise Unsuitable("microheterogeneity",
                             f"chain {asym} holds {out[label][0]} and {mon} at label {label}")
        out[label] = (mon or "", f"{auth or ''}{ins or ''}")
    return out


def _parents(block: mmcif.Block, asym: str) -> Dict[int, str]:
    """label_seq_id -> the standard residue a modified one is, where the entry says."""
    table = block.get("pdbx_struct_mod_residue")
    if table is None:
        return {}
    return {_int(seq): parent for row_asym, seq, parent in
            table.rows("label_asym_id", "label_seq_id", "parent_comp_id")
            if row_asym == asym and _int(seq) is not None and parent}


def _atoms(block: mmcif.Block, asym: str, model: str) -> Dict[int, List[Tuple[str, Atom]]]:
    """label_seq_id -> [(alt id, atom)] of one chain in one model, hydrogens dropped."""
    table = block["atom_site"]
    out: Dict[int, List[Tuple[str, Atom]]] = {}
    for (group, element, name, alt, row_asym, seq, x, y, z, occupancy, b,
         row_model) in table.rows(
            "group_PDB", "type_symbol", "label_atom_id", "label_alt_id", "label_asym_id",
            "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv",
            "pdbx_PDB_model_num"):
        if row_asym != asym or row_model != model or element in ("H", "D"):
            continue
        label = _int(seq)
        if label is None:
            continue
        out.setdefault(label, []).append((alt or "", Atom(
            name=name or "", element=element or "", x=_float(x), y=_float(y), z=_float(z),
            occupancy=_float(occupancy, 1.0), b=_float(b), hetero=group == "HETATM")))
    return out


def _one_conformation(found: List[Tuple[str, Atom]]) -> Tuple[Tuple[Atom, ...], Optional[str]]:
    """One residue's atoms in one conformation: those with no alternate, and of
    the alternates the one with the highest mean occupancy, then the first."""
    alternates: Dict[str, List[Atom]] = {}
    order: List[str] = []
    for alt, atom in found:
        if alt:
            if alt not in alternates:
                alternates[alt] = []
                order.append(alt)
            alternates[alt].append(atom)
    if not alternates:
        return tuple(atom for _, atom in found), None
    chosen = max(order, key=lambda a: (sum(x.occupancy for x in alternates[a])
                                       / len(alternates[a]), -order.index(a)))
    return tuple(atom for alt, atom in found if not alt or alt == chosen), chosen


def _distance(a: Atom, b: Atom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


# ------------------------------------------------------------ the chains


def _chain(block: mmcif.Block, inputs: Inputs, candidate: Candidate, pick: Pick, node: str,
           model: str, sequence: str, audit: Dict[str, list]) -> Drawn:
    asym = pick.asym
    product = inputs.products[pick.product]
    entity_beg, ref_beg, span = pick.segment
    sifts = _sifts(block, asym)
    if pick.mapped_by == "alignment":
        # SIFTS mapped nothing of the entity, and its one gapless place on the
        # protein stands in for the map; the file must not say otherwise.
        named = sorted({acc for acc, _, _ in sifts.values()})
        if named:
            raise Unsuitable("mapping", f"SIFTS maps chain {asym} to {', '.join(named)}")
    else:
        if not sifts:
            raise Unsuitable("mapping", f"PDBe's file carries no SIFTS map of chain {asym}")
        # The segment, residue by residue, as the file maps it.
        for k in range(span):
            label = entity_beg + k
            mapped = sifts.get(label)
            if mapped is None or mapped[0] != inputs.accession or mapped[1] != ref_beg + k:
                raise Unsuitable("mapping", f"the file's SIFTS map of chain {asym} at label "
                                            f"{label} is {mapped}, RCSB's segment says "
                                            f"{inputs.accession} {ref_beg + k}")
    scheme = _scheme(block, asym)
    parents = _parents(block, asym)
    atoms = _atoms(block, asym, model)

    # What is listed: of a product, the whole of the construct's part of it,
    # placed or not; of a domain, what the rule chose to draw around it, from
    # the first residue the entry places there to the last.
    low, high = max(product.start, ref_beg), min(product.end, ref_beg + span - 1)
    if candidate.scope == "domain":
        low, high = max(low, pick.first), min(high, pick.last)
    residues: List[Residue] = []
    for number in range(low, high + 1):
        label = entity_beg + (number - ref_beg)
        name, auth = scheme.get(label, ("", ""))
        held, alternate = _one_conformation(atoms.get(label, []))
        if alternate is not None:
            audit["alternates"].append({"node": node, "n": number, "kept": alternate})
        letter = sequence[label - 1] if 0 < label <= len(sequence) else "X"
        parent = parents.get(label)
        if parent is None and name and name not in ONE_LETTER:
            parent = next((standard for standard, one in ONE_LETTER.items() if one == letter),
                          None)
        if parent is not None:
            audit["modified"].append({"node": node, "n": number, "residue": name,
                                      "parent": parent})
        residues.append(Residue(number, label, auth, name, letter, parent, held, alternate,
                                inputs.translation[number - 1]))
    placed = [r for r in residues if r.placed]
    if not placed:
        raise Unsuitable("coverage", f"the file places no residue of chain {asym} in the product")
    reach = (product.start, product.end)
    if candidate.scope == "domain":
        reach = (placed[0].number, placed[-1].number)
        residues = [r for r in residues if reach[0] <= r.number <= reach[1]]

    # Unbroken: each placed residue a peptide bond, or the residues between,
    # from the last.
    for before, after in zip(placed, placed[1:]):
        gap = after.number - before.number
        step = _distance(before.atom("CA"), after.atom("CA"))
        if ((gap == 1 and not STEP[0] <= step <= STEP[1])
                or (gap > 1 and step > PER_RESIDUE * gap + _SLACK)):
            raise Unsuitable("chain break", f"chain {asym}: CA {before.number} to CA "
                                            f"{after.number} is {step:.1f} A across {gap} "
                                            f"residue(s)")

    for residue in placed:
        if residue.letter != inputs.translation[residue.number - 1]:
            audit["mismatches"].append({"node": node, "n": residue.number,
                                        "entry": residue.letter,
                                        "record": inputs.translation[residue.number - 1]})
    if residues[-1].number > MAX_NUMBER:
        raise Unsuitable("numbering", f"residue {residues[-1].number} cannot be written in PDB "
                                      f"format, which numbers to {MAX_NUMBER:,}")
    return Drawn(node, pick.product, pick.entity, asym, pick.auth, reach, tuple(residues))


def _touch(chains: Sequence[Drawn], bonded: set) -> bool:
    """Whether every chain is joined to the others, by a disulfide or by atoms
    within `CONTACT` of each other."""
    import numpy as np

    if len(chains) < 2:
        return True
    points = [np.array([[a.x, a.y, a.z] for r in c.residues for a in r.atoms]) for c in chains]
    joined = {0}
    changed = True
    while changed:
        changed = False
        for i in range(len(chains)):
            if i in joined:
                continue
            for j in list(joined):
                pair = frozenset((chains[i].node, chains[j].node))
                if pair in bonded or _near(points[i], points[j]):
                    joined.add(i)
                    changed = True
                    break
    return len(joined) == len(chains)


def _near(a, b) -> bool:
    import numpy as np

    if not len(a) or not len(b):
        return False
    for start in range(0, len(a), 512):
        chunk = a[start:start + 512]
        gaps = ((chunk[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)
        if float(gaps.min()) <= CONTACT * CONTACT:
            return True
    return False


def _bridges(block: mmcif.Block, chains: Sequence[Drawn], audit: Dict[str, list]
             ) -> Tuple[Bridge, ...]:
    """The entry's disulfides between residues drawn, where its sulfurs meet."""
    table = block.get("struct_conn")
    if table is None:
        return ()
    where: Dict[Tuple[str, int], Tuple[Drawn, Residue]] = {
        (c.asym, r.label): (c, r) for c in chains for r in c.residues}
    out: List[Bridge] = []
    seen = set()
    for kind, asym1, seq1, sym1, asym2, seq2, sym2 in table.rows(
            "conn_type_id", "ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_symmetry",
            "ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_symmetry"):
        if kind != "disulf":
            continue
        one, other = where.get((asym1, _int(seq1))), where.get((asym2, _int(seq2)))
        if one is None or other is None:
            continue
        ends = sorted([(one[0].node, one[1].number, one[1]), (other[0].node, other[1].number,
                                                              other[1])],
                      key=lambda e: (e[1], e[0]))
        pair = (ends[0][0], ends[0][1], ends[1][0], ends[1][1])
        if pair in seen:
            continue
        seen.add(pair)
        said = {"a": pair[1], "a_node": pair[0], "b": pair[3], "b_node": pair[2]}
        if any(end[2].record and end[2].record != "C" for end in ends):
            # A cysteine the entry has and the gene does not make: not a
            # bridge of this protein's (R5.6).
            audit["bridges_dropped"].append({**said, "reason": "not a cysteine in the record"})
            continue
        if (sym1 or "1_555") != "1_555" or (sym2 or "1_555") != "1_555":
            audit["bridges_dropped"].append({**said, "reason": "to a copy in the crystal"})
            continue
        sulfurs = [ends[0][2].atom("SG"), ends[1][2].atom("SG")]
        if sulfurs[0] is None or sulfurs[1] is None:
            audit["bridges_dropped"].append({**said, "reason": "no sulfur in the entry"})
            continue
        # A bridge is drawn as rods from CA through CB to the sulfur and back.
        if any(end[2].atom(name) is None for end in ends for name in ("CA", "CB")):
            audit["bridges_dropped"].append({**said, "reason": "a cysteine without its CA or CB"})
            continue
        apart = _distance(sulfurs[0], sulfurs[1])
        if apart > BRIDGE_REACH:
            audit["bridges_dropped"].append({**said, "reason": "sulfurs apart",
                                             "distance": round(apart, 2)})
            continue
        out.append(Bridge(pair[1], pair[0], pair[3], pair[2], round(apart, 3)))
    return tuple(sorted(out, key=lambda b: (b.a, b.b, b.a_node, b.b_node)))


def _secondary(block: mmcif.Block, chains: Sequence[Drawn]
               ) -> Tuple[Dict[str, Tuple[Tuple[str, int, int], ...]], Optional[str]]:
    """Each chain's helices and strands as the entry assigns them, in the
    precursor's numbering and clipped to what is drawn; and None where the
    entry assigns none, for the bake to ask PyMOL."""
    by_asym = {c.asym: c for c in chains}
    found: Dict[str, List[Tuple[str, int, int]]] = {c.node: [] for c in chains}
    any_records = False
    spans: List[Tuple[str, Optional[str], Optional[str], Optional[str]]] = []
    helices = block.get("struct_conf")
    if helices is not None:
        # Every helix class is a helix (alpha, 3-10, pi), as the twenty's are
        # read; a polyproline stretch, which DSSP also files here, is not one
        # a cartoon draws.
        spans += [("helix", asym, first, last) for kind, asym, first, last in helices.rows(
            "conf_type_id", "beg_label_asym_id", "beg_label_seq_id", "end_label_seq_id")
            if (kind or "").startswith("HELX") and "_PP_" not in (kind or "")]
    strands = block.get("struct_sheet_range")
    if strands is not None:
        spans += [("strand", asym, first, last) for asym, first, last in strands.rows(
            "beg_label_asym_id", "beg_label_seq_id", "end_label_seq_id")]
    for kind, asym, first, last in spans:
        chain = by_asym.get(asym)
        if chain is None or _int(first) is None or _int(last) is None:
            continue
        any_records = True
        numbers = {r.label: r.number for r in chain.residues}
        inside = [numbers[label] for label in range(_int(first), _int(last) + 1)
                  if label in numbers]
        if inside:
            found[chain.node].append((kind, min(inside), max(inside)))
    if not any_records:
        return {node: () for node in found}, None
    # PDBe's file carries PDBe's own assignment, and says by what (DSSP), where
    # the legacy file carries the depositors' records.
    kinds = block.get("struct_conf_type")
    by = sorted({c for c in (kinds.get("criteria") if kinds is not None else []) if c})
    return ({node: tuple(sorted(ranges, key=lambda r: (r[1], r[2], r[0])))
             for node, ranges in found.items()},
            f"{', '.join(by) or 'the assignment'} in the entry's file from PDBe "
            f"(struct_conf, struct_sheet_range)")


# ------------------------------------------------------------ the whole


def examine(inputs: Inputs, candidate: Candidate, text: bytes) -> Examined:
    """Hold `candidate` to its entry's own file, and read what both bakes need.

    Raises `Unsuitable` where the file shows a documented rule broken, and
    `mmcif.CifError` where it cannot be read at all, which is no verdict."""
    block = mmcif.parse(text.decode("utf-8", errors="replace"))
    named = (block.get("entry") or mmcif.Table("entry")).value("id")
    if (named or "").upper() != candidate.entry.upper():
        raise mmcif.CifError(f"the file is {named}'s, not {candidate.entry}'s")
    if "atom_site" not in block:
        raise mmcif.CifError("the file has no atoms")
    # Upper case, as RCSB's data API writes them (the file may not).
    methods = tuple(sorted({m.upper() for m in
                            (block.get("exptl") or mmcif.Table("exptl")).get("method") if m}))
    nmr = any("NMR" in m for m in methods)
    model = _model(block, nmr)
    sequences = {}
    poly = block.get("entity_poly")
    if poly is not None:
        sequences = {entity: "".join((seq or "").split())
                     for entity, seq in poly.rows("entity_id", "pdbx_seq_one_letter_code_can")}

    audit: Dict[str, list] = {"alternates": [], "modified": [], "mismatches": [],
                              "bridges_dropped": []}
    chains = []
    for node, pick in zip(NODES, sorted(candidate.picks, key=lambda p: p.product)):
        entity_id = pick.entity.rsplit("_", 1)[-1]
        sequence = sequences.get(entity_id)
        if not sequence:
            raise mmcif.CifError(f"the file has no sequence for entity {entity_id}")
        chains.append(_chain(block, inputs, candidate, pick, node, model, sequence, audit))

    allowed = int(MISMATCH_SHARE * sum(1 for c in chains for r in c.residues if r.placed))
    if len(audit["mismatches"]) > allowed:
        raise Unsuitable("fidelity", f"{len(audit['mismatches'])} residue(s) placed are not the "
                                     f"record's, {allowed} allowed")
    if sum(len(r.atoms) for c in chains for r in c.residues) > MAX_ATOMS:
        raise Unsuitable("numbering", f"more atoms than PDB format numbers ({MAX_ATOMS:,})")

    bridges = _bridges(block, chains, audit)
    bonded = {frozenset((b.a_node, b.b_node)) for b in bridges if b.a_node != b.b_node}
    if not _touch(chains, bonded):
        raise Unsuitable("apart", "the chains drawn do not touch in this entry")
    secondary, method = _secondary(block, chains)

    unbridged = [list(pair) for pair in _uniprot_pairs(inputs, chains)
                 if not any((b.a, b.b) == tuple(pair) for b in bridges)]
    audit["uniprot_bridges_not_drawn"] = unbridged
    return Examined(
        entry=candidate.entry.upper(), sha256=hashlib.sha256(text).hexdigest(), model=model,
        methods=methods, chains=tuple(chains), bridges=bridges, secondary=secondary,
        secondary_method=method, audit=audit)


def _uniprot_pairs(inputs: Inputs, chains: Sequence[Drawn]) -> List[Tuple[int, int]]:
    """UniProt's disulfides with both residues among those drawn (R5.6)."""
    drawn = {r.number for c in chains for r in c.residues}
    return [pair for pair in inputs.disulfides if set(pair) <= drawn]


def summary(examined: Examined) -> dict:
    """What the report and the audit say of an examined entry, without its atoms."""
    return {
        "entry": examined.entry, "sha256": examined.sha256, "model": examined.model,
        "methods": list(examined.methods),
        "chains": [{"node": c.node, "product": c.product, "entity": c.entity, "asym": c.asym,
                    "auth": c.auth, "first": c.first, "last": c.last,
                    "placed": sum(1 for r in c.residues if r.placed),
                    "listed": c.last - c.first + 1,
                    "construct": [c.residues[0].number, c.residues[-1].number]}
                   for c in examined.chains],
        "bridges": [asdict(b) for b in examined.bridges],
        "secondary_method": examined.secondary_method,
        "helix": sum(1 for spans in examined.secondary.values() for s in spans if s[0] == "helix"),
        "strand": sum(1 for spans in examined.secondary.values() for s in spans
                      if s[0] == "strand"),
        **examined.audit,
    }


# ------------------------------------------------------------ the working copy

# Chains of the working copy, by node: what `bake.export` and the readers in
# `pdb.py` select by.
LETTERS = {"chainA": "A", "chainB": "B", "chainC": "C"}


def _atom_name(atom: Atom) -> str:
    """An atom's name in PDB format's four columns: a one-letter element's name
    starts in the second (` CA ` is a carbon, `CA  ` a calcium)."""
    if len(atom.name) >= 4 or len(atom.element) == 2:
        return atom.name[:4].ljust(4)
    return (" " + atom.name).ljust(4)


def working_copy(examined: Examined, secondary: Optional[Dict[str, Sequence]] = None) -> str:
    """The entry as both bakes read it: one PDB-format text, in the record's
    numbering, with one chain a node (A, B, C), holding nothing but what is drawn.

    It is what `alphafold.working_copy` hands the bake of a predicted fold, for
    an experimental one: the twenty's own readers (`pdb.atoms`, `ssbonds`,
    `chain_residues`, `secondary_structure`) and PyMOL read it as they read an
    entry's legacy file, so nothing downstream knows the difference.

    - `ATOM` records for every residue placed, under the name of the residue
      the gene makes: a modified residue (a selenomethionine, a pyroglutamate)
      is drawn as its parent, with its own atoms, and the audit names it.
      Coordinates are the entry's, to the digit.
    - `REMARK 465` for residues the construct held and the entry never placed.
    - `SEQADV` where a placed residue is not the record's.
    - `HELIX` and `SHEET` from `secondary` (each node's (kind, first, last)),
      which defaults to the entry's own assignment; and `SSBOND` for the
      bridges drawn.
    """
    secondary = examined.secondary if secondary is None else secondary
    lines = [f"REMARK   1 {examined.entry} MODEL {examined.model}, AS DRAWN: THE RECORD'S NUMBERING"]
    three = {one: name for name, one in ONE_LETTER.items()}

    def name_of(residue: Residue) -> str:
        return residue.parent or residue.name

    for chain in examined.chains:
        for residue in chain.residues:
            if residue.placed and residue.letter not in three:
                raise Unsuitable("residue", f"{residue.name} at {residue.number} has no "
                                            f"standard residue to be drawn as")
    for chain in examined.chains:
        letter = LETTERS[chain.node]
        # Declared for every residue the construct held, placed or not: the
        # folding track lists each, and its check holds each difference to a
        # declaration, as it holds the twenty's to their SEQADV records.
        for residue in chain.residues:
            if residue.record and residue.letter != residue.record:
                lines.append(f"SEQADV {examined.entry[:4]:<4} {three.get(residue.letter, 'UNK'):>3} "
                             f"{letter} {residue.number:4d}  UNP              "
                             f"{three.get(residue.record, 'UNK'):>3}       CONFLICT")
    unplaced = [(LETTERS[c.node], r) for c in examined.chains for r in c.residues if not r.placed]
    if unplaced:
        lines.append("REMARK 465   M RES C SSSEQI")
        for letter, residue in unplaced:
            lines.append(f"REMARK 465     {three.get(residue.letter, 'UNK'):>3} {letter} "
                         f"{residue.number:5d} ")
    serial = 0
    for chain in examined.chains:
        by_number = {r.number: r for r in chain.residues}
        for index, (kind, first, last) in enumerate(
                s for s in secondary.get(chain.node, ()) if s[0] == "helix"):
            serial += 1
            a, b = three.get(by_number[first].letter, "UNK"), three.get(by_number[last].letter, "UNK")
            lines.append(f"HELIX  {serial:3d} {serial:3d} {a:>3} {LETTERS[chain.node]} {first:4d}  "
                         f"{b:>3} {LETTERS[chain.node]} {last:4d}  1"
                         f"{'':30}{last - first + 1:5d}")
    serial = 0
    for chain in examined.chains:
        by_number = {r.number: r for r in chain.residues}
        for kind, first, last in (s for s in secondary.get(chain.node, ()) if s[0] == "strand"):
            serial += 1
            a, b = three.get(by_number[first].letter, "UNK"), three.get(by_number[last].letter, "UNK")
            sheet = f"S{serial:02d}"[-3:]
            lines.append(f"SHEET    1 {sheet:>3} 1 {a:>3} {LETTERS[chain.node]}{first:4d}  "
                         f"{b:>3} {LETTERS[chain.node]}{last:4d}   0")
    for n, bridge in enumerate(examined.bridges, 1):
        lines.append(f"SSBOND {n:3d} CYS {LETTERS[bridge.a_node]} {bridge.a:4d}    "
                     f"CYS {LETTERS[bridge.b_node]} {bridge.b:4d}")
    serial = 0
    for chain in examined.chains:
        letter = LETTERS[chain.node]
        for residue in chain.residues:
            for atom in residue.atoms:
                serial += 1
                lines.append(
                    f"ATOM  {serial:5d} {_atom_name(atom)} {name_of(residue):>3} {letter}"
                    f"{residue.number:4d}    {atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
                    f"{atom.occupancy:6.2f}{min(atom.b, 999.99):6.2f}          "
                    f"{atom.element.upper():>2}")
        serial += 1
        last = [r for r in chain.residues if r.placed][-1]
        lines.append(f"TER   {serial:5d}      {name_of(last):>3} {letter}{last.number:4d}")
    lines.append("END")
    return "\n".join(lines) + "\n"
