from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

HERE = Path(__file__).resolve().parent
if str(HERE.parents[1]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[1]))

from pipeline.denovo.alphabet import MAX_LENGTH, MIN_LENGTH
from pipeline.targets import TARGETS, Target

Piece = Tuple[str, int]


def length(region) -> int:
    return region.end - region.start + 1


def in_range(size: int) -> bool:
    return MIN_LENGTH <= size <= MAX_LENGTH


def pieces(target: Target, kept_only: bool = True) -> List[Piece]:
    return [
        (region.label, length(region))
        for region in target.regions
        if in_range(length(region)) and (region.kept or not kept_only)
    ]


def offered(target: Target) -> List[Piece]:
    return pieces(target) if target.mature_peptides else []


def sources(targets: Sequence[Target] = TARGETS) -> List[Tuple[Target, List[Piece]]]:
    found = [(t, offered(t)) for t in targets]
    return sorted(
        [(t, p) for t, p in found if p],
        key=lambda row: (-len(row[1]), row[0].slug),
    )


def too_short(targets: Sequence[Target] = TARGETS) -> List[Tuple[Target, List[Piece]]]:
    found = []
    for target in targets:
        small = [
            (region.label, length(region))
            for region in target.regions
            if region.kept and length(region) < MIN_LENGTH
        ]
        if small:
            found.append((target, small))
    return found


def withheld(targets: Sequence[Target] = TARGETS) -> List[Tuple[Target, List[Piece]]]:
    return [
        (t, pieces(t)) for t in targets
        if not t.mature_peptides and pieces(t)
    ]


def whole(targets: Sequence[Target] = TARGETS) -> List[Target]:
    return [t for t in targets if in_range(t.aa)]


def _named(found: Sequence[Piece]) -> str:
    return ", ".join("%s %d" % piece for piece in found)


def report(targets: Sequence[Target] = TARGETS, markdown: bool = False) -> str:
    lines: List[str] = []
    rows = sources(targets)
    trainable = [row for row in rows if len(row[1]) >= 2]

    if markdown:
        lines.append("| protein | pieces of %d-%d the record carries |"
                     % (MIN_LENGTH, MAX_LENGTH))
        lines.append("|---|---|")
        for target, found in rows:
            lines.append("| %s | %s |" % (target.slug, _named(found)))
        return "\n".join(lines)

    lines.append("%d curated proteins, model range %d-%d residues"
                 % (len(targets), MIN_LENGTH, MAX_LENGTH))
    smallest = min(targets, key=lambda t: t.aa)
    lines.append("whole proteins in range: %d (smallest is %s at %d)"
                 % (len(whole(targets)), smallest.slug, smallest.aa))
    lines.append("")
    lines.append("what GeneRecord.peptides can offer the model")
    for target, found in rows:
        lines.append("  %-12s %d piece(s): %s"
                     % (target.slug, len(found), _named(found)))
    lines.append("")
    lines.append("%d can offer a set, %d can offer two or more"
                 % (len(rows), len(trainable)))

    held = withheld(targets)
    if held:
        lines.append("")
        lines.append("in range but mature_peptides is False, so no record carries them")
        for target, found in held:
            lines.append("  %-12s %d piece(s): %s"
                         % (target.slug, len(found), _named(found)))

    small = too_short(targets)
    if small:
        lines.append("")
        lines.append("kept regions under %d residues, which the model cannot take"
                     % MIN_LENGTH)
        for target, found in small:
            lines.append("  %-12s %s" % (target.slug, _named(found)))
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--markdown", action="store_true")
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])
    print(report(markdown=args.markdown))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
