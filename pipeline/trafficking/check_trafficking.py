"""Check every baked `trafficking` track against its record and its table row.

    .venv/Scripts/python pipeline/fetch_tracks.py --kind record   # the records, digest-checked
    .venv/Scripts/python pipeline/trafficking/check_trafficking.py

For each target, the payload under `pipeline/data/assets/trafficking/`:

- names its protein: the slug, gene and UniProt accession are the table's;
- is as long as the record: `residues` is the table's `aa`, and the length of
  the protein the stored record translates, so a span numbered on UniProt's
  sequence lands on the residue the walk draws;
- keeps every span inside that chain: each transmembrane span, and the GPI
  anchor's site and signal, fall within 1..residues, the spans in order and
  apart, the signal running from the residue after the site to the last;
- agrees with the region table about the GPI-anchor signal: where the table
  names one the anchor's signal is that region, and where it names none there
  is no anchor;
- says where it came from: the release, the release's date and the day it
  was fetched.

Exits 1 with the problems listed. Reads `pipeline/data/` only.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import Target  # noqa: E402
from pipeline.trafficking.bake_trafficking import (  # noqa: E402
    SCHEMA_VERSION,
    TRAFFICKING_TARGETS,
    trafficking_asset,
)

# What `targets.py` calls the stretch a GPI anchor takes the place of.
GPI_SIGNAL_LABEL = "GPI-anchor signal"

_RELEASE = re.compile(r"^\d{4}_\d{2}$")
_DAY = re.compile(r"^\d{4}-\d{2}-\d{2}$")


def record_length(target: Target) -> int | None:
    """The length of the protein the stored record translates, or None."""
    path = DATA / target.mock_asset
    if not path.exists():
        return None
    protein = json.loads(path.read_bytes()).get("protein") or {}
    return len(protein.get("translation", ""))


def problems_of(target: Target, track: dict | None = None) -> list[str]:
    """Everything wrong with one protein's track: the file's, or [track]."""
    where = target.slug
    if track is None:
        path = DATA / trafficking_asset(target)
        if not path.exists():
            return [f"{where}: not baked ({path})"]
        track = json.loads(path.read_bytes())

    out: list[str] = []
    for key, want in (("slug", target.slug), ("gene", target.gene),
                      ("uniprot", target.uniprot), ("schema_version", SCHEMA_VERSION)):
        if track.get(key) != want:
            out.append(f"{where}: {key} is {track.get(key)!r}, expected {want!r}")

    residues = track.get("residues")
    if residues != target.aa:
        out.append(f"{where}: {residues} residues, the table says {target.aa}")
    recorded = record_length(target)
    if recorded is None:
        out.append(f"{where}: no record to hold the length to; "
                   f"run fetch_tracks.py --kind record first")
    elif recorded != residues:
        out.append(f"{where}: {residues} residues, the record translates {recorded}")
    if not isinstance(residues, int) or residues < 1:
        return out

    def inside(start: object, end: object, what: str) -> bool:
        if not (isinstance(start, int) and isinstance(end, int)
                and 1 <= start <= end <= residues):
            out.append(f"{where}: {what} {start}-{end} is not within 1..{residues}")
            return False
        return True

    previous = 0
    for span in track.get("transmembrane", []):
        start, end = span.get("start"), span.get("end")
        if inside(start, end, "transmembrane span"):
            if start <= previous:
                out.append(f"{where}: transmembrane span {start}-{end} is out of order "
                           f"or overlaps the one before")
            previous = end

    gpi = track.get("gpi_anchor")
    if gpi is not None:
        site = gpi.get("site")
        signal = gpi.get("signal") or {}
        site_inside = inside(site, site, "GPI anchor site")
        signal_inside = inside(signal.get("start"), signal.get("end"), "GPI-anchor signal")
        if site_inside and signal_inside and (
                signal["start"] != site + 1 or signal["end"] != residues):
            out.append(f"{where}: the GPI-anchor signal {signal['start']}-{signal['end']} "
                       f"is not everything after residue {site}")
    named = [(r.start, r.end) for r in target.regions if r.label == GPI_SIGNAL_LABEL]
    carried = [] if gpi is None else [((gpi.get("signal") or {}).get("start"),
                                       (gpi.get("signal") or {}).get("end"))]
    if named != carried:
        out.append(f"{where}: the table's GPI-anchor signal is {named or 'none'}, "
                   f"the track's is {carried or 'none'}")

    for place in track.get("location", []):
        if not place.get("location"):
            out.append(f"{where}: a location with no place named: {place}")

    if track.get("source") != "UniProtKB":
        out.append(f"{where}: source is {track.get('source')!r}")
    if not _RELEASE.match(str(track.get("release"))):
        out.append(f"{where}: release is {track.get('release')!r}")
    for key in ("release_date", "retrieved"):
        if not _DAY.match(str(track.get(key))):
            out.append(f"{where}: {key} is {track.get(key)!r}")
    return out


def main() -> int:
    found = [problem for target in TRAFFICKING_TARGETS for problem in problems_of(target)]
    for problem in found:
        print(f"  {problem}", file=sys.stderr)
    print(f"{len(TRAFFICKING_TARGETS)} trafficking tracks checked, "
          f"{len(found)} problem(s)", file=sys.stderr)
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
