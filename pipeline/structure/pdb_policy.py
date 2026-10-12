"""Choosing a protein's experimental structure, by rule: the suitability policy.

A protein built on demand has nobody to choose an entry for it. This is the
rule that does, over what `pdb_search.py` found, with every number in it named
here and versioned (`POLICY_VERSION`): a change to any of them is a new policy,
and what an older one decided (a cached "nothing suitable") is not reused.

Two kinds of thing are decided, kept apart:

- **Eligibility**, which rejects, each time with its reason: the method, the
  source organism (the protein's own, not the expression host's), the mapping
  (SIFTS's, to the canonical sequence), one contiguous stretch of it, the
  letters, the coverage, and the experiment's support for the region drawn,
  judged by what the method can measure: density fit for a crystal, map fit for
  EM, geometry for NMR.
- **Rank**, among what is eligible: the mature product before a domain; more
  of a precursor's products drawn together before fewer; how much of what is
  drawn the entry places (in tenths below the best, so a residue more never
  beats a much better experiment); a mature construct before a precursor's;
  fewer residues that are not the gene's; then the experiment (a method with
  a resolution before one without, a region measured residue by residue
  before one that is not, the finer resolution, the better fit); and then the
  ids, so the same inputs choose the same chains every time.

Coverage counts residues the entry places, each once, mapped to the protein:
not the construct's length, not its copies.

Nothing here touches the network or a file. Policy 1 was frozen on 2026-10-12,
once the dry-run report over 52 proteins (`scripts/dry_fold.py`) had been read
with the user: every threshold below, and the three orderings the report put
to them, each marked where it is. `README.md` has what was weighed.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from pipeline.structure.pdb_search import Alignment, Chain, Domain, Entity, Entry

POLICY_VERSION = 1

HUMAN = 9606

# The experimental methods a fold may come from. Computed and integrative
# models are neither.
CRYSTAL = frozenset({"X-RAY DIFFRACTION", "NEUTRON DIFFRACTION", "ELECTRON CRYSTALLOGRAPHY"})
EM = frozenset({"ELECTRON MICROSCOPY"})
NMR = frozenset({"SOLUTION NMR", "SOLID-STATE NMR"})
METHODS = CRYSTAL | EM | NMR

# --- What makes a stretch of protein worth drawing.
# A mature product: this share of its residues placed by the entry.
PRODUCT_COVERAGE = 0.70
# Otherwise an annotated domain: this share of its residues, and at least
# this many residues long.
DOMAIN_COVERAGE = 0.90
DOMAIN_MIN = 30
# What is drawn of a domain: the annotated stretch, and no more than this many
# residues either side of it. An annotation's ends are approximate (p53's core
# is solved as 94-293 where Pfam says 100-288), but a domain's label must not
# be put on a long tail beside it (the prion protein's 90-133, which NMR
# places and which has no fixed shape).
DOMAIN_FLANK = 10
# A model has one node a product, and three colours for them.
MAX_PRODUCTS = 3

# --- Fidelity. Residues placed that are not the record's, as a
# share of those drawn, rounded down: a peptide under fifty residues may have
# none, where a precursor-wide allowance would have let it differ throughout.
MISMATCH_SHARE = 0.02

# --- The experiment's support.
CRYSTAL_RESOLUTION = 3.0        # angstroms
CRYSTAL_RFREE = 0.30
EM_RESOLUTION = 4.0             # angstroms
NMR_RAMACHANDRAN = 5.0          # per cent of residues, the entry's
# Over the residues drawn: no more than this share fitting their density
# poorly (RSRZ over 2; RSCC under 0.8), or their map (Q-score more than the
# margin under what the resolution leads one to expect), or in a clash.
RSRZ_OUTLIER = 2.0
RSRZ_SHARE = 0.10
RSCC_POOR = 0.80
RSCC_SHARE = 0.20
Q_MARGIN = 0.10
# Half, and not a third (settled 2026-10-12): a hormone on its receptor or a
# loop in a 3-4 angstrom map sits under the map's average, and at 0.30
# oxytocin's only entries go, and mTOR's fullest.
Q_SHARE = 0.50
CLASH_SHARE = 0.30
# A residue measure counts only where it covers most of what is drawn.
MEASURED = 0.5

# --- Rank. Coverage is compared in tenths below the best of its scope: a few
# flexible residues at a terminus are not worth a worse experiment.
COVERAGE_BIN = 0.10
# An entry is a precursor's structure where it places this many residues the
# precursor loses (of a signal or transit peptide, a propeptide): a construct
# that runs a residue or two past a cut is not one.
PRECURSOR_MIN = 5


def expected_q(resolution: float) -> float:
    """The Q-score a residue's atoms reach, on average, at an EM map's
    resolution (Pintilie et al. 2020, Nature Methods 17:328)."""
    return 1.119 - 0.178 * resolution


# ------------------------------------------------------------ the protein


@dataclass(frozen=True)
class Product:
    """A mature product: a stretch of the precursor that is kept, and drawn as
    one chain of the model. Numbered in the precursor, 1-based inclusive."""

    index: int
    label: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class Inputs:
    """What a structure is chosen for: the protein as its record has it."""

    slug: str
    gene: str
    accession: str
    display: str
    translation: str
    products: Tuple[Product, ...]
    # What the precursor loses, as (start, end, label).
    removed: Tuple[Tuple[int, int, str], ...] = ()
    domains: Tuple[Domain, ...] = ()
    # UniProt's disulfide pairs, in precursor numbering.
    disulfides: Tuple[Tuple[int, int], ...] = ()

    @property
    def kept(self) -> int:
        return sum(p.length for p in self.products)


def products_of(regions: Sequence, cleaved: bool, length: int, chain_label: str
                ) -> Tuple[Tuple[Product, ...], Tuple[Tuple[int, int, str], ...]]:
    """The mature products, and what is removed, from a protein's regions.

    A cleaved precursor is cut into its kept regions, each a product. An uncut
    one is one chain, less what is removed from its ends (an initiator
    methionine): its regions are a table of domains, not of products, which is
    how the twenty's are written (`folding/bake_folding.mature_span`).
    `regions` are `targets.Region`s, or anything with their fields.
    """
    removed = tuple((r.start, r.end, r.label) for r in regions if not r.kept)
    if cleaved:
        kept = sorted((r for r in regions if r.kept), key=lambda r: r.start)
        if kept:
            return tuple(Product(i, r.label, r.start, r.end) for i, r in enumerate(kept)), removed
        # A precursor with nothing named as kept (apelin's overlapping
        # peptides leave the resolver only its signal and propeptide): the
        # products are what the removed regions leave.
        left, cursor = [], 1
        for low, high, _ in sorted(removed):
            if low > cursor:
                left.append((cursor, low - 1))
            cursor = max(cursor, high + 1)
        if cursor <= length:
            left.append((cursor, length))
        return tuple(Product(i, chain_label, a, b) for i, (a, b) in enumerate(left)), removed
    start, end = 1, length
    for low, high, _ in sorted(removed):
        if low == start:
            start = high + 1
    for low, high, _ in sorted(removed, reverse=True):
        if high == end:
            end = low - 1
    return (Product(0, chain_label, start, end),), removed


def domains_within(domains: Iterable[Domain], products: Sequence[Product]) -> Tuple[Domain, ...]:
    """The annotated domains of the mature products, each clipped to the
    product it lies in most, where at least half of it and `DOMAIN_MIN`
    residues are left (Pfam's prion domain runs past the mature chain into
    its GPI signal). A Pfam domain overlapping one UniProt names, by half of
    the shorter, is UniProt's."""
    clipped = []
    for d in domains:
        best = None
        for p in products:
            low, high = max(d.start, p.start), min(d.end, p.end)
            if high >= low and (best is None or high - low > best[1] - best[0]):
                best = (low, high)
        if best is None:
            continue
        length = best[1] - best[0] + 1
        if length < DOMAIN_MIN and length < d.end - d.start + 1:
            continue
        if 2 * length < d.end - d.start + 1:
            continue
        clipped.append(Domain(d.kind, d.name, best[0], best[1], d.source))
    named = [d for d in clipped if d.source == "UniProt"]
    out = list(named)
    for d in clipped:
        if d.source == "UniProt":
            continue
        if any(2 * (min(d.end, u.end) - max(d.start, u.start) + 1)
               >= min(d.end - d.start + 1, u.end - u.start + 1) for u in named):
            continue
        out.append(d)
    return tuple(sorted(out, key=lambda d: (d.start, d.end, d.name)))


# ------------------------------------------------------------ candidates


@dataclass(frozen=True)
class Pick:
    """One chain of an entry, drawn as one product: which chain, by its label
    and author ids, which stretch of SIFTS's map, and what of the product it
    places."""

    product: int
    entity: str             # "1BM7_1"
    asym: str
    auth: str
    segment: Tuple[int, int, int]     # entity_beg, ref_beg, length
    first: int                        # the product's residues placed, first to last
    last: int
    placed: int
    mismatches: Tuple[int, ...]       # precursor numbers where the entry's letter is not the record's
    mapped_by: str = "SIFTS"          # or "alignment", for an entity SIFTS has not mapped


@dataclass
class Candidate:
    """An entry, the chains it would draw, and how it measures up."""

    entry: str
    scope: str                        # "product" or "domain"
    picks: Tuple[Pick, ...]
    domains: Tuple[str, ...]
    method: str
    resolution: Optional[float]
    useful: int
    coverage: float                   # of the mature protein
    share: float                      # of the product(s) drawn
    precursor: bool                   # the entry places residues the precursor loses
    mismatches: int
    allowed: int
    tier: int                         # 0: the region's own measures; 1: only the entry's
    support: Optional[float]
    measures: Dict[str, Optional[float]] = field(default_factory=dict)
    # Tenths below the best share any candidate of its scope and chain count
    # reaches: 0 within a tenth of the best. Set by `rank`.
    shortfall: int = 0

    def key(self) -> tuple:
        """What it is ranked by, best first: see the module's docstring.

        Two orderings in it were the user's to settle (2026-10-12). Fewer
        residues that are not the gene's come before the experiment: R5.4's
        "prefer wild-type", so a coarser wild type is drawn before a sharper
        mutant. And a method with a resolution comes before one without, ahead
        of whether the region was measured residue by residue: a crystal whose
        data nobody deposited still ranks before an NMR model."""
        return (
            0 if self.scope == "product" else 1,
            -len(self.picks),
            self.shortfall,
            1 if self.precursor else 0,
            self.mismatches,
            0 if self.resolution is not None else 1,
            self.tier,
            self.resolution if self.resolution is not None else 99.0,
            -(round(self.support, 3) if self.support is not None else 0.0),
            self.entry,
            tuple((p.product, _entity_number(p.entity), p.asym, p.segment[0]) for p in self.picks),
        )

    def as_json(self) -> dict:
        said = asdict(self)
        said["picks"] = [asdict(p) for p in self.picks]
        return said


_KEY_NAMES = ("scope", "chains together", "coverage", "a mature construct", "fidelity",
              "a resolution", "the region's own measures", "resolution", "support", "entry id",
              "chains")


def _entity_number(entity: str) -> int:
    tail = entity.rsplit("_", 1)[-1]
    return int(tail) if tail.isdigit() else 0


@dataclass(frozen=True)
class Rejection:
    """Why one entity, chain or entry is not drawn."""

    entity: str
    chain: Optional[str]
    reason: str              # a code, counted in reports
    detail: str


@dataclass
class Ranking:
    """Every eligible candidate, best first, and every rejection."""

    candidates: List[Candidate]
    rejections: List[Rejection]

    def reasons(self) -> Dict[str, int]:
        counted: Dict[str, int] = {}
        for rejection in self.rejections:
            counted[rejection.reason] = counted.get(rejection.reason, 0) + 1
        return dict(sorted(counted.items(), key=lambda kv: (-kv[1], kv[0])))


def why(winner: Candidate, runner_up: Optional[Candidate]) -> str:
    """Why the winner ranks above the runner-up, by the first thing that differs."""
    if runner_up is None:
        return "the only eligible candidate"
    for name, ours, theirs in zip(_KEY_NAMES, winner.key(), runner_up.key()):
        if ours == theirs:
            continue
        if name == "chains together":
            return (f"it draws {len(winner.picks)} of the products together; the runner-up "
                    f"{len(runner_up.picks)}")
        if name == "coverage":
            return (f"coverage: {winner.share:.0%} of what it draws is placed, against "
                    f"{runner_up.share:.0%}")
        if name == "scope":
            return "it places the mature product; the runner-up only a domain"
        if name == "a mature construct":
            return "it places nothing the precursor loses; the runner-up does"
        if name == "fidelity":
            return (f"fidelity: {winner.mismatches} residue(s) not the gene's against "
                    f"{runner_up.mismatches}")
        if name == "a resolution":
            return (f"its method has a resolution ({_angstroms(winner.resolution)}); the "
                    f"runner-up's has none")
        if name == "the region's own measures":
            return "the region drawn is measured residue by residue; the runner-up's is not"
        if name == "resolution":
            return (f"resolution: {_angstroms(winner.resolution)} against "
                    f"{_angstroms(runner_up.resolution)}")
        if name == "support":
            return (f"support: {winner.support:.0%} of residues fit the data, against "
                    f"{(runner_up.support or 0):.0%}")
        return f"{name}: {winner.entry} before {runner_up.entry}"
    return "identical"


def _angstroms(value: Optional[float]) -> str:
    return "none" if value is None else f"{value:.2f} Å"


# ------------------------------------------------------------ eligibility


def _method(entity: Entity) -> Tuple[Optional[str], Optional[str]]:
    """The method class a fold is judged by, or why there is none."""
    methods = set(entity.entry.methods)
    if not methods:
        return None, "the entry names no method"
    outside = methods - METHODS
    if outside:
        return None, f"method {', '.join(sorted(outside))}"
    if methods & CRYSTAL:
        return "crystal", None
    if methods & EM:
        return "em", None
    return "nmr", None


def _share(values: Dict[int, float], residues: Sequence[int], bad) -> Optional[float]:
    """The share of `residues` whose value is bad, or None where too few are measured."""
    measured = [values[r] for r in residues if r in values]
    if not residues or len(measured) < MEASURED * len(residues):
        return None
    return sum(1 for v in measured if bad(v)) / len(measured)


def _quality(entry: Entry, chains: Dict[str, Chain], kind: str, drawn: Dict[str, List[int]]
             ) -> Tuple[Optional[str], int, Optional[float], Dict[str, Optional[float]]]:
    """Whether the experiment supports the residues drawn: a reason if not, and
    otherwise the tier, the support and what was measured.

    `drawn` is, for each chain (asym), the label numbers of the residues it draws.
    """
    measures: Dict[str, Optional[float]] = {
        "resolution": entry.resolution, "rfree": entry.rfree,
        "clashscore": entry.clashscore, "ramachandran_outliers": entry.ramachandran_outliers}

    def over(name: str, bad) -> Optional[float]:
        """A residue measure's bad share over every chain drawn together."""
        total, poor = 0, 0
        for asym, residues in drawn.items():
            values = getattr(chains[asym], name)
            share = _share(values, residues, bad)
            if share is None:
                return None
            measured = sum(1 for r in residues if r in values)
            total += measured
            poor += round(share * measured)
        return poor / total if total else None

    clashed = sum(sum(1 for r in residues if r in chains[asym].clashes)
                  for asym, residues in drawn.items())
    count = sum(len(residues) for residues in drawn.values())
    clash_share = clashed / count if count else 0.0
    measures["clash_share"] = round(clash_share, 4)
    if clash_share > CLASH_SHARE:
        return f"{clash_share:.0%} of the residues drawn clash", 0, None, measures

    if kind == "crystal":
        if entry.resolution is None:
            return "a crystal structure with no resolution", 0, None, measures
        if entry.resolution > CRYSTAL_RESOLUTION:
            return f"resolution {entry.resolution:.2f} Å", 0, None, measures
        if entry.rfree is not None and entry.rfree > CRYSTAL_RFREE:
            return f"R-free {entry.rfree:.3f}", 0, None, measures
        rsrz = over("rsrz", lambda v: v > RSRZ_OUTLIER)
        rscc = over("rscc", lambda v: v < RSCC_POOR)
        measures.update(rsrz_outliers=rsrz, rscc_poor=rscc)
        if rsrz is not None and rsrz > RSRZ_SHARE:
            return f"{rsrz:.0%} of the residues drawn are RSRZ outliers", 0, None, measures
        if rscc is not None and rscc > RSCC_SHARE:
            return f"{rscc:.0%} of the residues drawn fit their density poorly", 0, None, measures
        tier = 0 if (rsrz is not None and rscc is not None and entry.rfree is not None) else 1
        support = None if rscc is None else 1 - rscc
        return None, tier, support, measures

    if kind == "em":
        if entry.resolution is None:
            return "an EM structure with no resolution", 0, None, measures
        if entry.resolution > EM_RESOLUTION:
            return f"resolution {entry.resolution:.2f} Å", 0, None, measures
        floor = expected_q(entry.resolution) - Q_MARGIN
        q = over("q_score", lambda v: v < floor)
        measures["q_under"] = q
        if q is not None and q > Q_SHARE:
            return f"{q:.0%} of the residues drawn fit their map poorly", 0, None, measures
        return None, (0 if q is not None else 1), (None if q is None else 1 - q), measures

    # NMR: no density and no map; the model's own geometry is what there is.
    rama = entry.ramachandran_outliers
    if rama is not None and rama > NMR_RAMACHANDRAN:
        return f"{rama:.1f}% Ramachandran outliers", 0, None, measures
    return None, (0 if rama is not None else 1), 1 - clash_share, measures


# ------------------------------------------------------------ recovery

# An entity SIFTS has not mapped is drawn only where it lies on the protein in
# exactly one place, without gaps: at least this identical, over at least this
# many residues. A repeat, or a paralogue's sequence, lies in two and is refused.
RECOVERY_IDENTITY = 0.95
RECOVERY_MIN = 20
_MISMATCH_COST = 3


def recover(inputs: Inputs, entity: Entity) -> Tuple[Optional[Alignment], str]:
    """A gapless alignment of an entity SIFTS has not mapped to any UniProt
    entry, where exactly one place on the protein takes it; else why not.

    Each diagonal's best stretch is found by its score (a match +1, a mismatch
    -3), which trims a tag or a linker off either end."""
    if any(a.database == "UniProt" for a in entity.alignments):
        return None, "mapped by SIFTS"
    s, t = entity.sequence, inputs.translation
    found = []
    for offset in range(-(len(s) - RECOVERY_MIN), len(t) - RECOVERY_MIN + 1):
        low, high = max(0, -offset), min(len(s), len(t) - offset)
        best, start, score, begin = 0, low, 0, low
        best_end = low
        for i in range(low, high):
            score += 1 if s[i] == t[i + offset] else -_MISMATCH_COST
            if score <= 0:
                score, begin = 0, i + 1
            elif score > best:
                best, start, best_end = score, begin, i + 1
        length = best_end - start
        if length < RECOVERY_MIN:
            continue
        same = sum(1 for i in range(start, best_end) if s[i] == t[i + offset])
        if same >= RECOVERY_IDENTITY * length:
            found.append((offset, start, length))
    if not found:
        return None, "no stretch of it lies on the protein"
    if len(found) > 1:
        return None, f"it lies on the protein in {len(found)} places"
    offset, start, length = found[0]
    return Alignment("UniProt", inputs.accession, None,
                     ((start + 1, start + offset + 1, length),), "alignment"), ""


# ------------------------------------------------------------ the rule


@dataclass(frozen=True)
class _Stretch:
    """One chain's one SIFTS segment, measured against each product."""

    entity: Entity
    chain: str
    auth: str
    segment: Tuple[int, int, int]
    placed: Dict[int, int]            # precursor number -> label number, residues with coordinates
    letters: Dict[int, str]           # precursor number -> the entry's letter
    mapped_by: str = "SIFTS"


def _stretches(inputs: Inputs, entity: Entity, rejections: List[Rejection]) -> List[_Stretch]:
    mapped = [a for a in entity.alignments
              if a.database == "UniProt" and a.accession == inputs.accession]
    if not mapped:
        others = sorted({f"{a.database} {a.accession}" for a in entity.alignments})
        rejections.append(Rejection(entity.id, None, "mapping",
                                    f"mapped to {', '.join(others) or 'nothing'}, not "
                                    f"UniProt {inputs.accession}"))
        return []
    canonical = [a for a in mapped if a.isoform in (None, inputs.accession)]
    if not canonical:
        rejections.append(Rejection(entity.id, None, "isoform",
                                    f"mapped to isoform {mapped[0].isoform}"))
        return []
    out = []
    length = len(inputs.translation)
    for alignment in canonical:
        for segment in alignment.regions:
            entity_beg, ref_beg, span = segment
            if ref_beg < 1 or ref_beg + span - 1 > length:
                rejections.append(Rejection(entity.id, None, "mapping",
                                            f"maps to residues {ref_beg}-{ref_beg + span - 1}, "
                                            f"beyond the protein's {length}"))
                continue
            if entity_beg < 1 or entity_beg + span - 1 > len(entity.sequence):
                rejections.append(Rejection(entity.id, None, "mapping",
                                            "maps residues beyond the entity's sequence"))
                continue
            for chain in entity.chains:
                placed, letters = {}, {}
                for k in range(span):
                    label = entity_beg + k
                    if label in chain.unobserved:
                        continue
                    placed[ref_beg + k] = label
                    letters[ref_beg + k] = entity.sequence[label - 1]
                if placed:
                    out.append(_Stretch(entity, chain.asym, chain.auth, segment, placed, letters,
                                        "SIFTS" if alignment.provenance != "alignment"
                                        else "alignment"))
    return out


def _pick(inputs: Inputs, stretch: _Stretch, product: Product,
          within: Optional[Tuple[int, int]] = None) -> Optional[Pick]:
    """What a stretch places of a product, or of the part of it `within`."""
    low, high = product.start, product.end
    if within is not None:
        low, high = max(low, within[0]), min(high, within[1])
    drawn = sorted(n for n in stretch.placed if low <= n <= high)
    if not drawn:
        return None
    mismatches = tuple(n for n in drawn if stretch.letters[n] != inputs.translation[n - 1])
    return Pick(product.index, stretch.entity.id, stretch.chain, stretch.auth, stretch.segment,
                drawn[0], drawn[-1], len(drawn), mismatches, stretch.mapped_by)


def _precursor(inputs: Inputs, stretches: Sequence[_Stretch]) -> bool:
    """Whether the chains drawn place residues the precursor loses (a signal or
    transit peptide, a propeptide), which makes the entry a precursor's
    structure, where it places `PRECURSOR_MIN` of them. An initiator
    methionine is not counted: constructs keep it."""
    lost = [(low, high) for low, high, label in inputs.removed
            if label != "Initiator methionine"]
    return sum(1 for s in stretches for n in s.placed
               if any(low <= n <= high for low, high in lost)) >= PRECURSOR_MIN


def _candidate(inputs: Inputs, scope: str, entry: Entry, label: str,
               chosen: Sequence[Tuple[Pick, _Stretch]], domains: Sequence[str], kind: str,
               rejections: List[Rejection]) -> Optional[Candidate]:
    """One way of drawing the protein from `entry`, if the experiment supports it.
    `label` names it in a rejection."""
    picks = [pick for pick, _ in chosen]
    useful = sum(p.placed for p in picks)
    mismatches = sum(len(p.mismatches) for p in picks)
    allowed = int(MISMATCH_SHARE * useful)
    chains = ", ".join(p.asym for p in picks)
    if mismatches > allowed:
        rejections.append(Rejection(label, chains, "fidelity",
                                    f"{mismatches} residue(s) drawn are not the record's, "
                                    f"{allowed} allowed over {useful}"))
        return None
    drawn: Dict[str, List[int]] = {}
    held: Dict[str, Chain] = {}
    for pick, stretch in chosen:
        drawn[pick.asym] = sorted(stretch.placed[n] for n in stretch.placed
                                  if pick.first <= n <= pick.last)
        held[pick.asym] = next(c for c in stretch.entity.chains if c.asym == pick.asym)
    reason, tier, support, measures = _quality(entry, held, kind, drawn)
    if reason is not None:
        rejections.append(Rejection(label, chains, "quality", reason))
        return None
    spanned = sum(inputs.products[p.product].length for p in picks)
    return Candidate(
        entry=entry.entry, scope=scope, picks=tuple(picks), domains=tuple(domains),
        method=", ".join(entry.methods), resolution=entry.resolution,
        useful=useful, coverage=useful / inputs.kept if inputs.kept else 0.0,
        share=useful / spanned if spanned else 0.0,
        precursor=_precursor(inputs, [stretch for _, stretch in chosen]),
        mismatches=mismatches, allowed=allowed, tier=tier, support=support, measures=measures)


def rank(inputs: Inputs, entities: Sequence[Entity]) -> Ranking:
    """Every eligible way of drawing the protein from these entities, best first."""
    rejections: List[Rejection] = []
    candidates: List[Candidate] = []
    by_entry: Dict[str, List[Tuple[Entity, _Stretch]]] = {}

    for entity in entities:
        if entity.entry.status not in (None, "REL"):
            rejections.append(Rejection(entity.id, None, "status",
                                        f"the entry is {entity.entry.status}"))
            continue
        kind, reason = _method(entity)
        if kind is None:
            rejections.append(Rejection(entity.id, None, "method", reason))
            continue
        if entity.sources != (HUMAN,):
            rejections.append(Rejection(entity.id, None, "organism",
                                        f"source organism(s) {list(entity.sources) or 'unnamed'}"))
            continue
        for stretch in _stretches(inputs, entity, rejections):
            by_entry.setdefault(entity.entry.entry, []).append((entity, stretch))
            qualified = False
            for product in inputs.products:
                pick = _pick(inputs, stretch, product)
                if pick is None:
                    continue
                if pick.placed >= PRODUCT_COVERAGE * product.length:
                    qualified = True
                    found = _candidate(inputs, "product", entity.entry, entity.id,
                                       [(pick, stretch)], [], kind, rejections)
                    if found is not None:
                        candidates.append(found)
                    continue
                covered = [d for d in inputs.domains
                           if product.start <= d.start and d.end <= product.end
                           and d.end - d.start + 1 >= DOMAIN_MIN
                           and sum(1 for n in stretch.placed if d.start <= n <= d.end)
                           >= DOMAIN_COVERAGE * (d.end - d.start + 1)]
                if covered:
                    qualified = True
                    around = (min(d.start for d in covered) - DOMAIN_FLANK,
                              max(d.end for d in covered) + DOMAIN_FLANK)
                    found = _candidate(inputs, "domain", entity.entry, entity.id,
                                       [(_pick(inputs, stretch, product, around), stretch)],
                                       [d.name for d in covered], kind, rejections)
                    if found is not None:
                        candidates.append(found)
            if not qualified:
                best = max((sum(1 for n in stretch.placed if p.start <= n <= p.end) / p.length
                            for p in inputs.products), default=0.0)
                rejections.append(Rejection(entity.id, stretch.chain, "coverage",
                                            f"places {best:.0%} of a mature product and no "
                                            f"whole annotated domain"))

    if len(inputs.products) > 1:
        candidates += _together(inputs, by_entry, rejections)

    # Coverage is the share of the product(s) drawn that the entry places,
    # compared against the best of its scope and chain count in tenths of
    # that best. So a few flexible residues never beat a better experiment; on
    # a long protein known only in pieces the pieces are told apart; and of a
    # precursor's several hormones the longest does not win for its length.
    groups: Dict[Tuple[str, int], float] = {}
    for c in candidates:
        group = (c.scope, len(c.picks))
        groups[group] = max(groups.get(group, 0.0), c.share)
    for c in candidates:
        best = groups[(c.scope, len(c.picks))]
        if best:
            c.shortfall = int((1 - c.share / best) / COVERAGE_BIN + 1e-9)
    candidates.sort(key=Candidate.key)
    return Ranking(candidates, rejections)


def _together(inputs: Inputs, by_entry: Dict[str, List[Tuple[Entity, _Stretch]]],
              rejections: List[Rejection]) -> List[Candidate]:
    """Entries that place several mature products at once, each from an entity
    of its own: a precursor's chains as they are found together after the cut.
    Whether those chains touch is for the file to say (`normalize`)."""
    out = []
    for entry, found in sorted(by_entry.items()):
        # Every qualifying chain for each product, best first.
        options: Dict[int, List[Tuple[tuple, Pick, _Stretch]]] = {}
        for entity, stretch in found:
            # Each chain's colour is its product's place among the kept regions,
            # and there are three of them.
            for product in inputs.products[:MAX_PRODUCTS]:
                pick = _pick(inputs, stretch, product)
                if pick is None or pick.placed < PRODUCT_COVERAGE * product.length:
                    continue
                order = (-pick.placed, len(pick.mismatches), _entity_number(entity.id),
                         pick.asym, pick.segment[0])
                options.setdefault(product.index, []).append((order, pick, stretch))
        if len(options) < 2:
            continue
        # A mature chain is a polymer of its own after the cut, so each product
        # comes from an entity of its own. The product placed best is served
        # first; the most products drawn is three.
        chosen: List[Tuple[Pick, _Stretch]] = []
        used = set()
        for index in sorted(options, key=lambda i: (min(o[0] for o in options[i]), i)):
            for _, pick, stretch in sorted(options[index], key=lambda o: o[0]):
                if pick.entity not in used:
                    used.add(pick.entity)
                    chosen.append((pick, stretch))
                    break
        chosen = sorted(sorted(chosen, key=lambda c: -c[0].placed)[:MAX_PRODUCTS],
                        key=lambda c: c[0].product)
        if len(chosen) < 2:
            continue
        # Copies of one product (a precursor's identical repeats) are not
        # chains that belong together: two ubiquitins are two molecules.
        sequences = [inputs.translation[inputs.products[c[0].product].start - 1:
                                        inputs.products[c[0].product].end] for c in chosen]
        if len(set(sequences)) < len(sequences):
            rejections.append(Rejection(entry, ", ".join(c[0].asym for c in chosen), "copies",
                                        "the products drawn together are copies of one sequence"))
            continue
        first = chosen[0][1].entity
        kind, _ = _method(first)
        candidate = _candidate(inputs, "product", first.entry, entry, chosen, [], kind,
                               rejections)
        if candidate is not None:
            out.append(candidate)
    return out
