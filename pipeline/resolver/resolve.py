"""A protein the catalog does not list, resolved from its index entry.

A buildable `protein_index` entry already names where its gene is -- its
RefSeqGene, or the stretch of chromosome MANE puts its transcript on -- and
which transcript and protein are MANE Select. Its UniProt entry names the rest:
the processing that cuts the precursor, and its disulfides. This turns the two
into the `Target` row the twenty were written by hand as, and hands that row to
the bakers that already exist, unchanged: `mock/build_gene_record.build` makes
the record track here, and `constraint/score_protein.score_protein` makes the
ESM-2 track later, on a GPU. Nothing here re-implements a baker, so a protein
resolved on demand is drawn from a record made exactly as the twenty's were.

What resolver version 1 decides, and the rule for each:

- Regions are UniProt's processing features in precursor numbering. An
  initiator methionine, a signal or transit peptide and a propeptide are
  removed; a chain or a peptide is kept. Domains are not read: UniProt's
  overlap, and a protein's regions have to tile it. Where processing features
  overlap each other, the peptides inside a chain are left out first; if that
  is not enough, the protein is its leader and one chain.
- A precursor is `cleaved` when it loses more than its initiator methionine,
  the line `targets.py` draws for hemoglobin and SOD1.
- Disulfides are UniProt's bonds within this chain whose two residues are
  cysteines in the record's own protein. One that is not is left out, and the
  provenance names it.
- `mature_peptides` is False, so the walk stops at the precursor. Which of a
  record's peptides are a cleavage series worth a page was decided by hand for
  each of the twenty, and an automatic answer would sometimes draw one that is
  not there.
- The record's protein may differ from UniProt's by common alleles, as
  dystrophin's does at three of 3,685 residues. Up to 1% of residues, and
  never fewer than three, are declared, as `targets.py` declares them, and
  recorded in the provenance. More, or a different length, is refused: the two
  would be describing different molecules.
- The prose is templated from the index and the record, never written, and the
  provenance says so.

A refusal is a `Refused` whose message is the sentence the reader is shown.
Anything else that goes wrong (NCBI or UniProt not answering) is raised as it
is, for the worker to retry.
"""

from __future__ import annotations

import re
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Mapping

from app.genbank_parser import GeneNotFound
from pipeline.mock import build_gene_record as builder
from pipeline.seed_catalog import _aliases, protein_row
from pipeline.targets import GENE_PAGE_BUDGET_BP, Region, Source, Target
from pipeline.uniprot import Entry

# Written into every row this makes, and onto the request that made it. The
# twenty carry 0: they were written by hand.
RESOLVER_VERSION = 1

# How far the record's protein may differ from UniProt's and still be read as
# the same molecule carrying other common alleles: a share of its length, with
# a floor so a short peptide is not held to a fraction of one residue.
VARIANT_SHARE = 0.01
VARIANT_FLOOR = 3

_LEADERS = ("Initiator methionine", "Signal", "Transit peptide")
_REMOVED = _LEADERS + ("Propeptide",)
_KEPT = ("Chain", "Peptide")

# "Insulin B chain" is the B chain; "C peptide" is the C-peptide. The letter is
# what the field calls the piece, and what `partition` names residues by.
_LETTER = re.compile(r"(?:^|\s)([A-Z])[\s-](?:chain|peptide)$")


class Refused(Exception):
    """The resolver declined this protein. The message is the reader's sentence."""


@dataclass(frozen=True)
class IndexRow:
    """The `protein_index` columns the resolver reads."""

    uniprot: str
    gene: str
    name: str
    length: int
    refseq_nuc: str | None
    refseq_prot: str | None
    chrom_acc: str | None
    chrom_start: int | None
    chrom_end: int | None
    refseqgene: str | None
    buildable: bool
    unavailable_reason: str | None = None


@dataclass(frozen=True)
class Feature:
    """One UniProt feature, in precursor numbering."""

    kind: str
    start: int
    end: int
    description: str


@dataclass(frozen=True)
class Resolution:
    """Everything one resolved protein writes: its row, its aliases, its record."""

    target: Target
    record: bytes
    protein: dict
    aliases: list[tuple[str, str]]


def slug_for(gene: str) -> str:
    """The slug a protein resolved on demand takes: lower(gene), as 0001 set."""
    return gene.lower()


def source_for(row: IndexRow) -> Source:
    """Where the record comes from: the RefSeqGene, or the MANE transcript's span.

    The MANE Select protein and transcript are named either way, so a record
    holding several transcripts of the gene draws the one MANE chose.
    """
    if row.refseqgene:
        return Source(row.refseqgene, protein_id=row.refseq_prot, transcript_id=row.refseq_nuc)
    if row.chrom_acc and row.chrom_start and row.chrom_end:
        return Source(row.chrom_acc, row.chrom_start, row.chrom_end,
                      protein_id=row.refseq_prot, transcript_id=row.refseq_nuc)
    raise Refused(f"MANE places no transcript of {row.gene} on the genome.")


def features_of(body: dict, length: int, kinds: tuple[str, ...]) -> tuple[list[Feature], list[dict]]:
    """UniProt's features of these kinds, and the ones whose ends are not exact.

    A boundary UniProt calls unknown or outside the sequence cannot be drawn,
    so the feature is left out and named in what this returns second.
    """
    found: list[Feature] = []
    skipped: list[dict] = []
    for feature in body.get("features") or []:
        kind = feature.get("type")
        if kind not in kinds:
            continue
        location = feature.get("location") or {}
        start, end = location.get("start") or {}, location.get("end") or {}
        first, last = start.get("value"), end.get("value")
        description = (feature.get("description") or "").strip()
        exact = start.get("modifier") == "EXACT" and end.get("modifier") == "EXACT"
        if (not exact or not isinstance(first, int) or not isinstance(last, int)
                or not 1 <= first <= last <= length):
            skipped.append({"type": kind, "start": first, "end": last, "description": description})
            continue
        found.append(Feature(kind, first, last, description))
    return found, skipped


def _overlap(features: list[Feature]) -> bool:
    ordered = sorted(features, key=lambda f: (f.start, f.end))
    return any(after.start <= before.end for before, after in zip(ordered, ordered[1:]))


def _leader(features: list[Feature]) -> list[Feature]:
    """The removed run at the front: an initiator methionine, a signal peptide."""
    run: list[Feature] = []
    cursor = 1
    for feature in sorted(features, key=lambda f: (f.start, -f.end)):
        if feature.kind in _LEADERS and feature.start == cursor:
            run.append(feature)
            cursor = feature.end + 1
    return run


def _label(feature: Feature, display: str) -> str:
    text = feature.description
    if feature.kind == "Initiator methionine":
        return "Initiator methionine"
    if feature.kind == "Signal":
        return "Signal peptide"
    if feature.kind == "Transit peptide":
        return f"Transit peptide ({text[0].lower()}{text[1:]})" if text else "Transit peptide"
    if feature.kind == "Propeptide":
        # "Removed in mature form" says what happens to it, not what it is.
        return text if text and not text.startswith("Removed") else "Propeptide"
    return text or display


def _shorts(features: list[Feature]) -> list[str]:
    """What each region's residues are named by, the way `targets.py` names them.

    A leader is `S`, `T` or `Met`; one kept chain is named by position alone;
    pieces of a cut precursor take the letter their names give them ("Insulin B
    chain" is `B`), and are numbered in order where the names give no unique
    letters.
    """
    fixed = {"Initiator methionine": "Met", "Signal": "S", "Transit peptide": "T"}
    kept = [f for f in features if f.kind in _KEPT]
    pieces = (kept if len(kept) > 1 else []) + [f for f in features if f.kind == "Propeptide"]

    def shorts_from(names: dict[Feature, str]) -> list[str]:
        return [fixed.get(f.kind) or names.get(f, "") for f in features]

    letters = {}
    for piece in pieces:
        match = _LETTER.search(piece.description)
        if match:
            letters[piece] = match.group(1)
    shorts = shorts_from(letters)
    named = [s for s in shorts if s]
    if len(letters) == len(pieces) and len(set(named)) == len(named):
        return shorts
    # Some piece has no letter, or two share one: number them, which cannot
    # collide with each other or with a leader's S, T or Met.
    numbered = {f: str(n) for n, f in enumerate([p for p in pieces if p.kind in _KEPT], 1)}
    numbered.update({f: f"Pro{n}" for n, f in enumerate(
        [p for p in pieces if p.kind == "Propeptide"], 1)})
    return shorts_from(numbered)


def _named_gaps(regions: tuple[Region, ...]) -> tuple[Region, ...]:
    """Make a region of each gap in a cut precursor that a lone chain borders.

    `partition` names a gap between two pieces by their short names ("B / C
    cleavage site"), and a lone chain has none -- its residues are cited as
    "mature 157", which a short name would break -- so the gap would read
    "S /  cleavage site". It also numbers a lone chain from the end of the
    removed run at residue 1, which a gap after the leader breaks. A region of
    its own fixes both, and says what `partition` assumes of any gap in a cut
    precursor: it is the residues a protease cuts at, and it goes with the cut.
    """
    filled: list[Region] = []
    for region in regions:
        if filled and region.start > filled[-1].end + 1 and not (filled[-1].short and region.short):
            filled.append(Region("Cleavage site", "", filled[-1].end + 1, region.start - 1,
                                 kept=False))
        filled.append(region)
    return tuple(filled)


def regions_of(features: list[Feature], length: int, display: str) -> tuple[
        tuple[Region, ...], str, bool, str]:
    """The regions, the rule that chose them, whether the precursor is cut, and
    what its one chain is called."""
    chosen = sorted(features, key=lambda f: (f.start, f.end))
    rule = "processing"
    if _overlap(chosen):
        chosen = [f for f in chosen if f.kind != "Peptide"]
        rule = "processing, without the peptides inside its chains"
    if _overlap(chosen):
        leader = _leader(chosen)
        after = leader[-1].end + 1 if leader else 1
        chosen = leader + ([Feature("Chain", after, length, "")] if after <= length else [])
        rule = "leader and one chain: UniProt's processing features overlap"

    shorts = _shorts(chosen)
    regions = tuple(
        Region(_label(f, display), short, f.start, f.end, kept=f.kind in _KEPT)
        for f, short in zip(chosen, shorts)
    )
    cleaved = any(f.kind in _REMOVED and f.kind != "Initiator methionine" for f in chosen)
    if cleaved:
        regions = _named_gaps(regions)
    kept = [r for r in regions if r.kept]
    chain_label = kept[0].label if len(kept) == 1 else display
    return regions, rule, cleaved, chain_label


def disulfides_of(body: dict, protein: str) -> tuple[tuple[tuple[int, int], ...], list[list[int]]]:
    """UniProt's bonds within this chain, and those the record's residues refuse.

    A bond to another molecule is one position long and is not a pair here.
    """
    pairs: set[tuple[int, int]] = set()
    dropped: list[list[int]] = []
    found, _ = features_of(body, len(protein), ("Disulfide bond",))
    for bond in found:
        if bond.start == bond.end:
            continue
        if protein[bond.start - 1] == "C" and protein[bond.end - 1] == "C":
            pairs.add((bond.start, bond.end))
        else:
            dropped.append([bond.start, bond.end])
    return tuple(sorted(pairs)), dropped


def chromosome_of(accession: str | None) -> str | None:
    """`NC_000011.10` is chromosome 11; 23 and 24 are X and Y."""
    match = re.fullmatch(r"NC_0000(\d\d)\.\d+", accession or "")
    if not match or not 1 <= int(match.group(1)) <= 24:
        return None
    return {23: "X", 24: "Y"}.get(int(match.group(1)), str(int(match.group(1))))


def summary_for(row: IndexRow) -> str:
    """The one line the search card shows, from the index alone."""
    where = chromosome_of(row.chrom_acc)
    made = f"Made by {row.gene} on chromosome {where}." if where else f"Made by {row.gene}."
    return (f"{made} Built on demand from UniProt {row.uniprot} "
            f"and MANE Select {row.refseq_nuc}, not curated by hand.")


def _record_protein(target: Target) -> str:
    """The record's protein, read the way `build` will check it, before the
    table knows which differences from UniProt to declare."""
    try:
        record = builder.fetch(target)
        payload = builder.tidy(builder.extract_gene(record, target.gene), target)
        return builder.check_frame(payload, f"{target.slug} (as fetched)")
    except (ValueError, GeneNotFound) as exc:
        raise Refused(str(exc)) from exc


def _encoded(payload: dict) -> bytes:
    """The record's bytes exactly as the builder writes them to disk."""
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / "record.json"
        builder.write(payload, path)
        return path.read_bytes()


def resolve(row: IndexRow, entry: Entry, *, mane_release: str | None = None,
            budget: int = GENE_PAGE_BUDGET_BP) -> Resolution:
    """Resolve one buildable index entry into its row, aliases and record."""
    if not row.buildable:
        raise Refused(row.unavailable_reason or "The index does not call this protein buildable.")
    expected = (entry.body.get("sequence") or {}).get("value") or ""
    if not expected:
        raise Refused(f"UniProt {row.uniprot} carries no sequence.")

    display = row.name
    features, skipped = features_of(entry.body, len(expected), _REMOVED + _KEPT)
    regions, rule, cleaved, chain_label = regions_of(features, len(expected), display)
    target = Target(
        slug=slug_for(row.gene), gene=row.gene, uniprot=row.uniprot, display=display,
        source=source_for(row), structure=None, cleaved=cleaved, chain_label=chain_label,
        regions=regions, aa=len(expected), mature_peptides=False,
    )

    stated = _record_protein(target)
    if len(stated) != len(expected):
        raise Refused(f"The record's protein is {len(stated):,} residues long; "
                      f"UniProt {row.uniprot} is {len(expected):,}.")
    found = tuple((i, a, b) for i, (a, b) in enumerate(zip(stated, expected), 1) if a != b)
    allowed = max(VARIANT_FLOOR, int(len(expected) * VARIANT_SHARE))
    if len(found) > allowed:
        raise Refused(f"The record's protein differs from UniProt {row.uniprot} at "
                      f"{len(found):,} residues, more than the {allowed:,} common alleles allowed.")
    disulfides, dropped = disulfides_of(entry.body, stated)
    target = replace(target, uniprot_variants=found, disulfides=disulfides)

    try:
        payload = builder.build(target, budget)
    except (ValueError, GeneNotFound) as exc:
        raise Refused(str(exc)) from exc

    kept = sum(1 for region in regions if region.kept)
    curated = {
        "order": None,
        "display": display,
        "summary": summary_for(row),
        # Only an uncut protein's coding sequence is named as its one chain.
        "chain": None if cleaved else display,
        "facts": {
            "residues": len(payload["protein"]["translation"]),
            "exons": len(payload["exons"]),
            "chains": max(1, kept),
            "bridges": len(disulfides),
        },
        "chains": [],
        "chrome": None,
        "impact_explanations": False,
    }
    try:
        protein = protein_row(target, curated)
    except ValueError as exc:
        raise Refused(str(exc)) from exc
    protein.update(
        provenance={
            "source": "pipeline/resolver",
            "prose": "templated",
            "resolver_version": RESOLVER_VERSION,
            "uniprot": {"accession": entry.accession, "release": entry.release,
                        "release_date": entry.release_date, "retrieved": entry.retrieved},
            "mane": {"release": mane_release, "transcript": row.refseq_nuc,
                     "protein": row.refseq_prot},
            "regions_rule": rule,
            # The table this protein was resolved as, so a later bake rebuilds
            # the same `Target` from the row (`target_of`).
            "regions_from": [
                {"label": r.label, "short": r.short, "start": r.start, "end": r.end,
                 "kept": r.kept} for r in regions
            ],
            "chain_label": chain_label,
            "cleaved": cleaved,
            "uniprot_variants": [list(v) for v in found],
            "skipped_features": skipped,
            "disulfides_dropped": dropped,
        },
        resolver_version=RESOLVER_VERSION,
        catalog_order=None,
    )
    return Resolution(target, _encoded(payload), protein, _aliases(target, curated))


def target_of(protein: Mapping) -> Target:
    """The `Target` a resolved row was made as, rebuilt from the row itself."""
    provenance = protein["provenance"]
    return Target(
        slug=protein["slug"], gene=protein["gene"], uniprot=protein["uniprot"],
        display=protein["display"],
        source=Source(protein["accession"], protein["slice_start"], protein["slice_end"],
                      protein_id=protein["protein_id"], transcript_id=protein["transcript_id"]),
        structure=None,
        cleaved=provenance["cleaved"],
        chain_label=provenance["chain_label"],
        regions=tuple(Region(**region) for region in provenance["regions_from"]),
        disulfides=tuple(tuple(pair) for pair in protein["disulfides"]),
        aa=protein["residues"],
        uniprot_variants=tuple(tuple(v) for v in provenance["uniprot_variants"]),
        mature_peptides=protein["mature_peptides"],
    )
