"""Hold the baked loci to NCBI Gene's, which is curated apart from them.

    .venv/bin/python -I pipeline/locus/verify_locus.py                  # the curated, in pipeline/data
    .venv/bin/python -I pipeline/locus/verify_locus.py --dir OUT/locus  # a survey's tracks
        [--gene-info Homo_sapiens.gene_info.gz] [--hgnc hgnc_complete_set.txt]

Each track's place is held to the `map_location` NCBI Gene gives its gene,
read in bulk from `Homo_sapiens.gene_info.gz` (downloaded into a temporary
directory unless `--gene-info` names a copy). The gene is NCBI's by the
GeneID MANE gives the track's transcript, as the pinned snapshot holds MANE
v1.5 (`inputs.py`), never by its symbol, which moves. Every difference is
one of:

- `equal`: written the same;
- `coarser`: NCBI writes a band the track's place lies in (21q22 for
  21q22.11);
- `range`: NCBI writes a range of bands the track's place lies in. NCBI
  ranges over all of a gene's transcripts, the track over MANE Select's one
  (AIM2: 1q23.1-q23.2 for 1q23.1);
- `multi`: NCBI names more than one place (`19q13.42|11p`), one of them
  holding the track's;
- `par`: NCBI names the gene on X and on Y, a pseudoautosomal gene
  (`Xp22.33 and Yp11.32`), the X or Y place holding the track's;
- `placed`: the gene MANE puts on a patch or alternate locus scaffold,
  carried to its chromosome (`placed_from`), and NCBI's place holds it;
- `none`: NCBI gives no place (`-`), or MANE names no GeneID;
- `wrong`: NCBI's place does not hold the track's. There must be none.

HGNC is a second opinion, written by hand and sometimes coarser or older
than NCBI's (`--hgnc`, its complete set): it is counted, never held to.
Exits 1 on any `wrong`. Writes nothing.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
import re
import sys
import tempfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from app.protein_index import versionless  # noqa: E402
from pipeline.locus import inputs as locus_inputs  # noqa: E402
from pipeline.locus.bake_locus import LOCUS_TARGETS, locus_asset  # noqa: E402
from pipeline.paths import DATA  # noqa: E402

GENE_INFO = ("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/"
             "Homo_sapiens.gene_info.gz")

CLASSES = ("equal", "coarser", "range", "multi", "par", "placed", "none", "wrong")

_PLACE = re.compile(r"^(\d{1,2}|X|Y)(?:([pq][\d.]*)(?:-([pq][\d.]*))?)?$")

# How NCBI and HGNC end the place of a gene on an alternate locus.
_ALTERNATE = re.compile(r"\s+alternate reference locus$")


def _covers(prefix: str, name: str) -> bool:
    """Whether band [name] is [prefix] or lies in it. Each digit of a band's
    name narrows it (region, band, then sub-bands after the point), so a
    coarser name is a prefix of a finer one, the point aside: p covers the
    arm, p1 covers p13.1, p15 covers p15.5, p21.3 covers p21.33 (the 850-band
    sub-band of the 400-band p21.3), and p15.5 covers only itself."""
    return name.replace(".", "").startswith(prefix.replace(".", ""))


def within(track: dict, place: str) -> bool | None:
    """Whether [track]'s gene lies inside [place], a location as NCBI or HGNC
    writes it (a band, a range of bands, an arm, or a whole chromosome);
    None where the place is not one this can read."""
    found = _PLACE.match(_ALTERNATE.sub("", place.strip()))
    if not found:
        return None
    chromosome, first = found.group(1), found.group(2)
    if chromosome != track["chromosome"]:
        return False
    if first is None:
        return True
    last = found.group(3) or first
    bands = track["bands"]
    starts = [i for i, b in enumerate(bands) if _covers(first, b["name"])]
    ends = [i for i, b in enumerate(bands) if _covers(last, b["name"])]
    if not starts or not ends:
        return None
    lo, hi = min(starts + ends), max(starts + ends)
    span = track["span"]
    return bands[lo]["start"] <= span["start"] and span["end"] <= bands[hi]["end"]


def _any_within(track: dict, places: list[str]) -> bool:
    return any(within(track, place) for place in places)


def classify(track: dict, ncbi: str | None) -> str:
    """How [track]'s place reads against NCBI's [ncbi] (`CLASSES`)."""
    if not ncbi or ncbi.strip() in ("-", ""):
        return "none"
    if ncbi == track["locus"]:
        return "placed" if "placed_from" in track else "equal"
    if "|" in ncbi:
        return "multi" if _any_within(track, ncbi.split("|")) else "wrong"
    parts = re.split(r"\s+and\s+|;", ncbi)
    if len(parts) > 1:
        return "par" if _any_within(track, parts) else "wrong"
    if not within(track, ncbi):
        return "wrong"
    if "placed_from" in track:
        return "placed"
    return "range" if "-" in ncbi else "coarser"


def read_gene_info(path: Path) -> dict[int, str]:
    """NCBI Gene's map_location, by GeneID."""
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return {int(row["GeneID"]): row["map_location"]
                for row in csv.DictReader(stream, delimiter="\t")}


def read_hgnc(path: Path) -> dict[str, str]:
    """HGNC's location, by HGNC id."""
    with path.open(encoding="utf-8") as stream:
        return {row["hgnc_id"]: row["location"]
                for row in csv.DictReader(stream, delimiter="\t")}


def verify(tracks: list[dict], snapshot: locus_inputs.Inputs, places: dict[int, str],
           hgnc: dict[str, str] | None = None) -> dict:
    """Every track's class against NCBI, and HGNC's opinion where given."""
    out = {"classes": collections.Counter(), "rows": [], "hgnc": collections.Counter()}
    for track in tracks:
        mane = snapshot.mane.get(versionless(track["transcript"]["ensembl"])) or {}
        gene_id = int(mane["NCBI_GeneID"].split(":")[-1]) if mane.get("NCBI_GeneID") else None
        ncbi = places.get(gene_id) if gene_id is not None else None
        verdict = classify(track, ncbi)
        out["classes"][verdict] += 1
        if verdict != "equal":
            out["rows"].append((track["slug"], track["locus"], ncbi, verdict))
        if hgnc is not None:
            place = hgnc.get(mane.get("HGNC_ID") or "")
            opinion = None if not place else within(track, place)
            out["hgnc"][{True: "holds it", None: "cannot be read",
                         False: "stale (NCBI agrees)" if verdict != "wrong"
                         else "disagrees"}[opinion]] += 1
    return out


def _fetch(url: str, folder: Path) -> Path:
    path = folder / url.rsplit("/", 1)[-1]
    request = urllib.request.Request(url, headers={"User-Agent": "helixpeek-locus-verify"})
    with urllib.request.urlopen(request, timeout=300) as response:
        path.write_bytes(response.read())
    return path


def main(arguments: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--dir", type=Path, default=None,
                        help="a directory of `*_locus.json` tracks (default: the curated)")
    parser.add_argument("--gene-info", type=Path, default=None)
    parser.add_argument("--hgnc", type=Path, default=None)
    args = parser.parse_args(arguments)

    if args.dir is not None:
        paths = sorted(args.dir.glob("*_locus.json"))
    else:
        paths = [DATA / locus_asset(target) for target in LOCUS_TARGETS]
    missing = [path for path in paths if not path.exists()]
    tracks = [json.loads(path.read_bytes()) for path in paths if path.exists()]

    snapshot = locus_inputs.load()
    with tempfile.TemporaryDirectory() as folder:
        gene_info = args.gene_info or _fetch(GENE_INFO, Path(folder))
        digest = hashlib.sha256(gene_info.read_bytes()).hexdigest()
        found = verify(tracks, snapshot, read_gene_info(gene_info),
                       read_hgnc(args.hgnc) if args.hgnc else None)
    for slug, locus, ncbi, verdict in found["rows"]:
        print(f"{verdict:<7} {slug:<16} baked {locus:<16} NCBI {ncbi or '-'}")
    classes = found["classes"]
    print(f"\n{len(tracks)} loci held to NCBI Gene (gene_info sha256 {digest}): "
          + ", ".join(f"{name} {classes[name]:,}" for name in CLASSES))
    if found["hgnc"]:
        print("HGNC, a second opinion: "
              + ", ".join(f"{name} {count:,}" for name, count in sorted(found["hgnc"].items())))
    for path in missing:
        print(f"not baked: {path}")
    return 1 if classes["wrong"] or missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
