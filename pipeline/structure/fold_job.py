"""The fold of a protein built on demand: an experimental entry first.

`choose` finds every experimental entry of the protein (`pdb_search`), ranks
every way of drawing it (`pdb_policy`), and holds the best to its own file in
turn (`normalize.examine`), and then to the bake (`draw`, where one is given:
`experimental.Drawer`), until one passes. What it returns says how the search
ended, which is what decides what happens next:

- **chosen**: an entry passed, and is drawn with the twenty's fold animation;
- **none**: the search was complete and nothing passed, so AlphaFold DB's
  model may be drawn instead, as before, without an animation;
- **incomplete**: candidates were left unexamined when the bound on downloads
  or on bakes was reached; the search resumes, and nothing is decided yet.

A candidate is passed over, and the next one tried, only for a documented
rule it fails, raised as `normalize.Unsuitable` with its reason. An outage is
never that, nor one of the three endings: `pdb_search.Unreachable` and
`Incomplete`, a file that cannot be read (`mmcif.CifError`), and a bake or a
check that breaks, are raised as they are.
"""

from __future__ import annotations

import dataclasses
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from pipeline.structure import normalize, pdb_policy, pdb_search
from pipeline.structure.normalize import Examined, Unsuitable
from pipeline.structure.pdb_policy import Candidate, Inputs, Ranking

# The most entry files examined in one attempt. A candidate that fails its own
# file usually fails for its entry's sake, and the next is another entry.
EXAMINE_LIMIT = 8

# The most candidates baked in one attempt. A bake is a PyMOL export, a scene
# compiled and the app's own check: seconds for a peptide, a minute or two for
# two thousand residues.
DRAW_LIMIT = 3

# Whether entities SIFTS has not mapped are looked for by sequence.
RECOVERY = True

# What bakes an examined candidate, and holds what it baked to the app's own
# rules: `(candidate, examined, entry, entity) -> built`, the entity being the
# one drawn where one chain is. It raises `Unsuitable` where the candidate
# cannot be drawn by a documented rule, and anything else where it broke.
Draw = Callable[[Candidate, Examined, pdb_search.Entry, Optional[pdb_search.Entity]], object]


@dataclass
class Examination:
    """One candidate held to its file, and to the bake where there is one:
    what came of it. `examined` is set once its file passed, and `refused`
    where a rule turned it down, at either step; a winner has `examined` and
    no `refused`, and `built` where it was baked."""

    candidate: Candidate
    examined: Optional[Examined] = None
    refused: Optional[Unsuitable] = None
    built: object = None

    @property
    def passed(self) -> bool:
        return self.examined is not None and self.refused is None


@dataclass
class Choice:
    """How the search for a protein's experimental structure ended."""

    inputs: Inputs
    outcome: str                                  # chosen | none | incomplete
    discovery: pdb_search.Discovery
    recovered: Tuple[str, ...]
    recovery_refusals: Dict[str, str]
    ranking: Ranking
    examinations: List[Examination]
    entries: Dict[str, pdb_search.Entry]
    entities: Dict[str, pdb_search.Entity] = field(default_factory=dict)
    seconds: Dict[str, float] = field(default_factory=dict)

    @property
    def winner(self) -> Optional[Examination]:
        return next((e for e in self.examinations if e.passed), None)

    @property
    def refused(self) -> List[str]:
        """The candidates a rule turned down, by `candidate_id`: what a
        resumed search is told to skip."""
        return [candidate_id(e.candidate) for e in self.examinations if not e.passed]


def inputs_of(target, translation: str, uniprot_body: dict,
              pfam: Sequence[pdb_search.Domain] = ()) -> Inputs:
    """What a structure is chosen for, from a protein's `Target` (as `targets.py`
    writes the twenty, or as the resolver rebuilds one from its row), the
    record's translation, its UniProt entry, whose domains are read here (the
    resolver keeps only its processing), and Pfam's domains on it
    (`pdb_search.pfam_domains`), for those UniProt names none of."""
    products, removed = pdb_policy.products_of(target.regions, target.cleaved,
                                               len(translation), target.chain_label)
    domains = pdb_policy.domains_within(list(pdb_search.domains_of(uniprot_body)) + list(pfam),
                                        products)
    return Inputs(slug=target.slug, gene=target.gene, accession=target.uniprot,
                  display=target.display, translation=translation, products=products,
                  removed=removed, domains=domains,
                  disulfides=tuple(tuple(pair) for pair in target.disulfides))


def _recovered(inputs: Inputs, known: Sequence[str], post: pdb_search.Post,
               now: Callable[[], datetime]) -> List[str]:
    """Entities a sequence search finds that the accession search did not."""
    found = set()
    for product in inputs.products:
        sequence = inputs.translation[product.start - 1:product.end]
        found |= set(pdb_search.sequence_search(sequence, pdb_policy.RECOVERY_IDENTITY, post,
                                                now).ids)
    return sorted(found - set(known))


def choose(inputs: Inputs, *, get: pdb_search.Get = pdb_search.get_bounded,
           post: pdb_search.Post = pdb_search.post, limit: int = EXAMINE_LIMIT,
           recovery: bool = RECOVERY, skip: Sequence[str] = (),
           draw: Optional[Draw] = None, draws: int = DRAW_LIMIT,
           now: Callable[[], datetime] = lambda: datetime.now(timezone.utc)) -> Choice:
    """Search, rank, and examine in rank order until an entry passes.

    `skip` names candidates already examined and refused in an earlier
    attempt, by `candidate_id`, so a resumed search goes on from where it was.
    `draw` bakes each candidate whose file passes, at most `draws` of them an
    attempt; without it the first whose file passes is chosen, unbaked."""
    seconds: Dict[str, float] = {}
    started = time.monotonic()
    discovery = pdb_search.search(inputs.accession, post, now)
    recovered = _recovered(inputs, discovery.ids, post, now) if recovery else []
    seconds["search"] = round(time.monotonic() - started, 2)

    started = time.monotonic()
    entities = pdb_search.metadata(list(discovery.ids) + recovered, post)
    refusals: Dict[str, str] = {}
    usable = []
    for entity in entities:
        if entity.id in recovered:
            alignment, why = pdb_policy.recover(inputs, entity)
            if alignment is None:
                refusals[entity.id] = why
                continue
            entity = dataclasses.replace(entity, alignments=(alignment,))
        usable.append(entity)
    seconds["metadata"] = round(time.monotonic() - started, 2)

    started = time.monotonic()
    ranking = pdb_policy.rank(inputs, usable)
    entries = {e.entry.entry: e.entry for e in usable}
    seconds["rank"] = round(time.monotonic() - started, 2)

    started = time.monotonic()
    examinations: List[Examination] = []
    files: Dict[str, bytes] = {}
    refused_entries: Dict[str, Unsuitable] = {}
    skipped = set(skip)
    by_id = {e.id: e for e in usable}
    drawn = 0
    outcome = "none"
    for candidate in ranking.candidates:
        if candidate_id(candidate) in skipped:
            continue
        if candidate.entry in refused_entries:
            # The entry's file itself was refused: every chain of it is too.
            examinations.append(Examination(candidate, refused=refused_entries[candidate.entry]))
            continue
        if candidate.entry not in files and len(files) >= limit:
            outcome = "incomplete"
            break
        if candidate.entry not in files:
            try:
                files[candidate.entry] = pdb_search.updated_cif(candidate.entry, get)
            except pdb_search.Missing:
                refused_entries[candidate.entry] = Unsuitable(
                    "file", "PDBe has no SIFTS-enriched file of the entry")
            except pdb_search.TooLarge:
                refused_entries[candidate.entry] = Unsuitable(
                    "size", f"the entry's file is over {pdb_search.DOWNLOAD_LIMIT:,} bytes")
            if candidate.entry in refused_entries:
                examinations.append(Examination(candidate,
                                                refused=refused_entries[candidate.entry]))
                continue
        try:
            examined = normalize.examine(inputs, candidate, files[candidate.entry])
        except Unsuitable as refused:
            examinations.append(Examination(candidate, refused=refused))
            continue
        built = None
        if draw is not None:
            if drawn >= draws:
                # Left for the next attempt, unrecorded: it has not been tried.
                outcome = "incomplete"
                break
            drawn += 1
            lone = candidate.picks[0].entity if len(candidate.picks) == 1 else None
            try:
                built = draw(candidate, examined, entries[candidate.entry], by_id.get(lone))
            except Unsuitable as refused:
                examinations.append(Examination(candidate, examined=examined, refused=refused))
                continue
        examinations.append(Examination(candidate, examined=examined, built=built))
        outcome = "chosen"
        break
    seconds["examine"] = round(time.monotonic() - started, 2)
    return Choice(inputs, outcome, discovery, tuple(recovered), refusals, ranking, examinations,
                  entries, by_id, seconds)


def candidate_id(candidate: Candidate) -> str:
    """A candidate's name, stable across attempts: its entry and chains."""
    return candidate.entry + ":" + ",".join(
        f"{p.asym}@{p.segment[1]}" for p in sorted(candidate.picks, key=lambda p: p.product)) + (
        ":domain" if candidate.scope == "domain" else "")
