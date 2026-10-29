from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

HERE = Path(__file__).resolve().parent
if str(HERE.parents[1]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[1]))

from pipeline.denovo.alphabet import AMINO_ACIDS, MAX_LENGTH, MIN_LENGTH
from pipeline.paths import DATA
from pipeline.targets import TARGETS, Target

CANONICAL = frozenset(AMINO_ACIDS)
MATURE = "peptides"
LEADER = "signal_peptide"
PRODUCTS = (MATURE, "proprotein", "protein")
FIELDS = (MATURE, LEADER, "proprotein", "protein")

Entry = Tuple[str, str, str]
Source = Tuple[str, str, str]


def modelable(sequence: str) -> bool:
    return (
        MIN_LENGTH <= len(sequence) <= MAX_LENGTH
        and CANONICAL.issuperset(sequence)
    )


def record_path(target: Target) -> Path:
    return DATA / target.mock_asset


def missing(targets: Sequence[Target] = TARGETS) -> List[Target]:
    return [t for t in targets if not record_path(t).is_file()]


def read_record(target: Target) -> Dict:
    path = record_path(target)
    if not path.is_file():
        raise FileNotFoundError(
            "%s: no stored record at %s. Run pipeline/fetch_tracks.py --kind record"
            % (target.slug, path)
        )
    with open(path) as handle:
        return json.load(handle)


def entries(payload: Dict) -> List[Entry]:
    found: List[Entry] = []
    for peptide in payload.get("peptides") or ():
        found.append((MATURE, peptide.get("product") or "?", peptide["translation"]))
    for field in FIELDS[1:]:
        piece = payload.get(field)
        if piece and piece.get("translation"):
            found.append((field, piece.get("product") or field, piece["translation"]))
    return found


def of(target: Target, with_leaders: bool = False) -> List[Entry]:
    found = entries(read_record(target))
    if not with_leaders:
        found = [e for e in found if e[0] in PRODUCTS]
    return [e for e in found if modelable(e[2])]


def leaders(targets: Sequence[Target] = TARGETS) -> List[Tuple[Target, Entry]]:
    found = []
    for target in targets:
        for entry in entries(read_record(target)):
            if entry[0] == LEADER and modelable(entry[2]):
                found.append((target, entry))
    return found


def collect(
    targets: Sequence[Target] = TARGETS, with_leaders: bool = False
) -> Dict[str, List[Source]]:
    found: Dict[str, List[Source]] = {}
    for target in targets:
        for field, label, sequence in of(target, with_leaders):
            found.setdefault(sequence, []).append((target.slug, field, label))
    return found


def by_protein(
    targets: Sequence[Target] = TARGETS, with_leaders: bool = False
) -> List[Tuple[Target, List[Entry]]]:
    rows = [(t, of(t, with_leaders)) for t in targets]
    return [row for row in rows if row[1]]


def withheld(targets: Sequence[Target] = TARGETS) -> List[Tuple[Target, int]]:
    found = []
    for target in targets:
        if read_record(target).get("peptides"):
            continue
        sized = [
            r for r in target.regions
            if MIN_LENGTH <= (r.end - r.start + 1) <= MAX_LENGTH and r.kept
        ]
        if sized:
            found.append((target, len(sized)))
    return found


def under_floor(targets: Sequence[Target] = TARGETS) -> List[Tuple[Target, List[Entry]]]:
    found = []
    for target in targets:
        small = [
            e for e in entries(read_record(target))
            if len(e[2]) < MIN_LENGTH and CANONICAL.issuperset(e[2])
        ]
        if small:
            found.append((target, small))
    return found


def report(
    targets: Sequence[Target] = TARGETS,
    with_leaders: bool = False,
    markdown: bool = False,
) -> str:
    absent = missing(targets)
    if absent:
        return ("%d of %d records are not fetched (%s ...). Run "
                "pipeline/fetch_tracks.py --kind record"
                % (len(absent), len(targets), absent[0].slug))

    unique = collect(targets, with_leaders)
    rows = by_protein(targets, with_leaders)
    lines: List[str] = []

    if markdown:
        lines.append("| protein | peptides of %d-%d residues |" % (MIN_LENGTH, MAX_LENGTH))
        lines.append("|---|---|")
        for target, found in rows:
            lines.append("| %s | %s |" % (
                target.slug,
                ", ".join("%s %d" % (label, len(s)) for _, label, s in found),
            ))
        return "\n".join(lines)

    counted: Dict[str, int] = {}
    for sources in unique.values():
        counted[sources[0][1]] = counted.get(sources[0][1], 0) + 1
    residues = sum(len(s) for s in unique)

    lines.append("%d curated proteins, model range %d-%d residues"
                 % (len(targets), MIN_LENGTH, MAX_LENGTH))
    lines.append("%d unique peptides, %d residues, from %d proteins"
                 % (len(unique), residues, len(rows)))
    lines.append("by kind: %s"
                 % ", ".join("%s %d" % kv for kv in sorted(counted.items())))
    lines.append("")
    for target, found in rows:
        lines.append("  %-14s %s" % (
            target.slug,
            ", ".join("%s %d" % (label, len(s)) for _, label, s in found),
        ))

    held = withheld(targets)
    if held:
        lines.append("")
        lines.append("no peptide on the record, though targets.py sizes some in range")
        for target, count in held:
            lines.append("  %-14s %d region(s)" % (target.slug, count))

    if not with_leaders:
        held_back = leaders(targets)
        if held_back:
            lines.append("")
            lines.append("%d signal peptides set aside, as fetch_peptides.py sets them "
                         "aside: a leader is cleaved and destroyed, not a product, and "
                         "the corpus takes Peptide and Propeptide features only"
                         % len(held_back))

    small = under_floor(targets)
    if small:
        lines.append("")
        lines.append("on the record but under %d residues" % MIN_LENGTH)
        for target, found in small:
            lines.append("  %-14s %s"
                         % (target.slug,
                            ", ".join("%s %d" % (l, len(s)) for _, l, s in found)))
    return "\n".join(lines)


def write(path: Path, targets: Sequence[Target] = TARGETS,
          with_leaders: bool = False) -> int:
    unique = collect(targets, with_leaders)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(sorted(unique)) + "\n")
    return len(unique)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--with-leaders", action="store_true")
    parser.add_argument("--markdown", action="store_true")
    parser.add_argument("--write", type=Path, default=None)
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    print(report(with_leaders=args.with_leaders, markdown=args.markdown))
    if args.write is not None:
        count = write(args.write, with_leaders=args.with_leaders)
        print("\nwrote %d peptides to %s" % (count, args.write))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
