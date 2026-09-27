"""Hold the baked loci to HGNC's, which are curated apart from them.

    .venv/Scripts/python pipeline/locus/verify_locus.py

For each baked track, the location HGNC gives its gene symbol
(rest.genenames.org), and whether the baked place lies inside it. HGNC may
write a coarser place than the bake (21q22 where the bake reads 21q22.11), but
never a different one: the gene's span has to fall within the bands HGNC's
place covers, read off the track's own copy of the chromosome.

Nothing is written down on either side. The places the plan names to check
against (insulin 11p15.5, hemoglobin beta 11p15.4, p53 17p13.1, CFTR 7q31.2,
dystrophin Xp21, SOD1 21q22) are what this reads from HGNC, and what the bake
must agree with. Exits 1 on any disagreement, or on a place it cannot read.
"""

from __future__ import annotations

import json
import re
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.locus.bake_locus import LOCUS_TARGETS, locus_asset  # noqa: E402
from pipeline.paths import DATA  # noqa: E402

HGNC = "https://rest.genenames.org/fetch/symbol/{}"

_PLACE = re.compile(r"^(\d{1,2}|X|Y)([pq][\d.]*)(?:-([pq][\d.]*))?$")


def hgnc_location(symbol: str) -> str | None:
    request = urllib.request.Request(
        HGNC.format(symbol),
        headers={"Accept": "application/json", "User-Agent": "helixpeek-locus-verify"})
    with urllib.request.urlopen(request, timeout=60) as response:
        docs = json.load(response)["response"]["docs"]
    return docs[0].get("location") if docs else None


def _covers(prefix: str, name: str) -> bool:
    """Whether band [name] is [prefix] or one of its sub-bands: p15 covers
    p15.5, and p15.5 covers only itself."""
    return name == prefix or name.startswith(prefix + ".")


def within(track: dict, place: str) -> bool | None:
    """Whether [track]'s gene lies inside [place], a location as HGNC writes
    it; None where the place is not one this can read."""
    found = _PLACE.match(place.strip())
    if not found:
        return None
    chromosome, first, last = found.group(1), found.group(2), found.group(3) or found.group(2)
    if chromosome != track["chromosome"]:
        return False
    bands = track["bands"]
    starts = [i for i, b in enumerate(bands) if _covers(first, b["name"])]
    ends = [i for i, b in enumerate(bands) if _covers(last, b["name"])]
    if not starts or not ends:
        return None
    lo, hi = min(starts + ends), max(starts + ends)
    span = track["span"]
    return bands[lo]["start"] <= span["start"] and span["end"] <= bands[hi]["end"]


def main() -> int:
    wrong = 0
    for target in LOCUS_TARGETS:
        path = DATA / locus_asset(target)
        if not path.exists():
            print(f"{target.slug:<16} not baked")
            wrong += 1
            continue
        track = json.loads(path.read_bytes())
        place = hgnc_location(target.gene)
        verdict = None if place is None else within(track, place)
        label = {True: "agrees", False: "DISAGREES", None: "CANNOT READ"}[verdict]
        print(f"{target.slug:<16} {target.gene:<6} baked {track['locus']:<14} "
              f"HGNC {place or '-':<14} {label}")
        if verdict is not True:
            wrong += 1
    print(f"\n{len(LOCUS_TARGETS)} loci held to HGNC's: {wrong} problem(s)")
    return 1 if wrong else 0


if __name__ == "__main__":
    raise SystemExit(main())
