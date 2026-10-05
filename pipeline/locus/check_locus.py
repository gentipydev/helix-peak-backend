"""Check every baked `locus` track against itself and its table row.

    .venv/Scripts/python pipeline/locus/check_locus.py

For each target, the payload under `pipeline/data/assets/locus/`:

- names its protein: the slug, gene and UniProt accession are the table's;
- is on GRCh38, read from UCSC's hg38 tables, on a chromosome 1 to 22, X or
  Y named by its RefSeq accession;
- draws the whole chromosome: the bands run from base 1 to its last base
  without a gap or an overlap, the short arm's before the long arm's, each
  with a stain cytoBand uses;
- places the gene inside it, and names its band truly: the band or bands
  are exactly those the span overlaps, begin and end where they do, and the
  locus is written from them;
- agrees with the record where the record says where it is: a record cut
  from a chromosome is cut from this one, around this span, and a table row
  that names its transcript names the one the span is;
- reads the Atlas as it reads itself: its categories for how specific the
  gene's RNA is, and the tissues and cell types named, highest first, only
  where the gene is specific to something, each by a name the Atlas lists;
  its tissue cell type pairs, each a tissue it lists and a cell type; where
  in a cell and where it is secreted to, in the Atlas's own words;
- takes the zoom's path the rule takes (`path_of`), recomputed from the
  readings, with a class for its cell;
- says where it came from: each UCSC table's genome, update time and
  digest, the MANE release, the UniProt release, and the Protein Atlas
  version, release date, Ensembl version and licence.

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

from pipeline.locus.bake_locus import (  # noqa: E402
    ASSEMBLY,
    CELL_SPECIFICITY,
    GENOME,
    HPA_LICENCE,
    LOCUS_TARGETS,
    SCHEMA_VERSION,
    STAINS,
    TISSUE_SPECIFICITY,
    UNSPECIFIC,
    locus_asset,
    locus_of,
    path_of,
)
from pipeline.locus.cell_types import CLASSES, SECRETOME, SUBCELLULAR, single  # noqa: E402
from pipeline.locus.tissues import pair_tissue, tissue_name  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import Target  # noqa: E402

CHROMOSOMES = tuple(str(n) for n in range(1, 23)) + ("X", "Y")

_BAND = re.compile(r"^[pq]\d+(\.\d+)?$")
_SEQUENCE = re.compile(r"^NC_0000\d\d\.\d+$")
_UPDATED = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}$")
_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_MANE = re.compile(r"^v\d+(\.\d+)*$")
_RELEASE = re.compile(r"^\d{4}_\d{2}$")
_DAY = re.compile(r"^\d{4}-\d{2}-\d{2}$")
_ENSEMBL = re.compile(r"^ENSG\d{11}$")
_VERSION = re.compile(r"^\d+\.\d+$")


def problems_of(target: Target, track: dict | None = None) -> list[str]:
    """Everything wrong with one protein's track: the file's, or [track]."""
    where = target.slug
    if track is None:
        path = DATA / locus_asset(target)
        if not path.exists():
            return [f"{where}: not baked ({path})"]
        track = json.loads(path.read_bytes())

    out: list[str] = []
    for key, want in (("slug", target.slug), ("gene", target.gene),
                      ("uniprot", target.uniprot), ("schema_version", SCHEMA_VERSION),
                      ("assembly", ASSEMBLY), ("genome", GENOME)):
        if track.get(key) != want:
            out.append(f"{where}: {key} is {track.get(key)!r}, expected {want!r}")
    chromosome = track.get("chromosome")
    if chromosome not in CHROMOSOMES:
        out.append(f"{where}: chromosome {chromosome!r} is not 1 to 22, X or Y")
    if not _SEQUENCE.match(str(track.get("sequence"))):
        out.append(f"{where}: sequence {track.get('sequence')!r} is not a chromosome")

    # The chromosome, whole: its bands from base 1 to the last, in order.
    bands = track.get("bands") or []
    length = track.get("length")
    if not bands:
        return out + [f"{where}: no bands"]
    at = 1
    seen_q = False
    names: set[str] = set()
    for band in bands:
        name = str(band.get("name"))
        if not _BAND.match(name):
            out.append(f"{where}: band {name!r} is not a band name")
        if name in names:
            out.append(f"{where}: band {name} twice")
        names.add(name)
        if name.startswith("q"):
            seen_q = True
        elif seen_q:
            out.append(f"{where}: short-arm band {name} after the long arm began")
        if band.get("stain") not in STAINS:
            out.append(f"{where}: band {name} has stain {band.get('stain')!r}")
        if band.get("start") != at or not isinstance(band.get("end"), int) \
                or band["end"] < band["start"]:
            out.append(f"{where}: band {name} runs {band.get('start')}-{band.get('end')}, "
                       f"expected to start at {at}")
            return out
        at = band["end"] + 1
    if length != at - 1:
        out.append(f"{where}: length {length}, the bands end at {at - 1}")

    # The gene inside it, at the bands it overlaps.
    span = track.get("span") or {}
    start, end, strand = span.get("start"), span.get("end"), span.get("strand")
    if not (isinstance(start, int) and isinstance(end, int) and 1 <= start <= end <= at - 1):
        return out + [f"{where}: span {start}-{end} is not on the chromosome"]
    if strand not in (1, -1):
        out.append(f"{where}: strand {strand!r}")
    inside = [b for b in bands if b["start"] <= end and b["end"] >= start]
    band = track.get("band") or {}
    want_names = [b["name"] for b in inside]
    if band.get("names") != want_names:
        out.append(f"{where}: band names {band.get('names')}, the span overlaps {want_names}")
    if (band.get("start"), band.get("end")) != (inside[0]["start"], inside[-1]["end"]):
        out.append(f"{where}: band runs {band.get('start')}-{band.get('end')}, "
                   f"its bands {inside[0]['start']}-{inside[-1]['end']}")
    if track.get("locus") != locus_of(str(chromosome), want_names):
        out.append(f"{where}: locus {track.get('locus')!r}, "
                   f"expected {locus_of(str(chromosome), want_names)!r}")

    # The record, where it says where it is.
    source = target.source
    if source.accession.startswith("NC_"):
        if track.get("sequence") != source.accession:
            out.append(f"{where}: on {track.get('sequence')}, the record is cut from "
                       f"{source.accession}")
        if source.seq_start is not None and source.seq_stop is not None and not (
                source.seq_start <= start and end <= source.seq_stop):
            out.append(f"{where}: span {start}-{end} is outside the record's slice "
                       f"{source.seq_start}-{source.seq_stop}")
    transcript = (track.get("transcript") or {}).get("refseq")
    if source.transcript_id and transcript != source.transcript_id:
        out.append(f"{where}: transcript {transcript}, the table names {source.transcript_id}")

    # Where in the body, as the Atlas reads it.
    expression = track.get("expression") or {}
    if not _ENSEMBL.match(str(expression.get("ensembl_gene"))):
        out.append(f"{where}: Ensembl gene {expression.get('ensembl_gene')!r}")
    for key, categories, unit in (("tissue", TISSUE_SPECIFICITY, "ntpm"),
                                  ("cell_type", CELL_SPECIFICITY, "ncpm")):
        section = expression.get(key) or {}
        specificity = section.get("specificity")
        if specificity not in categories:
            out.append(f"{where}: {key} specificity {specificity!r}")
        specific = section.get("specific")
        if not isinstance(specific, list):
            out.append(f"{where}: {key} names no list")
            continue
        values = [found.get(unit) for found in specific]
        if any(not isinstance(found.get("name"), str) or not found["name"]
               for found in specific):
            out.append(f"{where}: {key} names one without a name")
        if any(not isinstance(value, (int, float)) or value <= 0 for value in values):
            out.append(f"{where}: {key} has a level that is not above zero")
        elif values != sorted(values, reverse=True):
            out.append(f"{where}: {key} is not highest first")
        if specificity in UNSPECIFIC and specific:
            out.append(f"{where}: {key} is {specificity}, yet names {len(specific)}")
        if specificity not in UNSPECIFIC and specificity in categories and not specific:
            out.append(f"{where}: {key} is {specificity}, yet names none")
        for found in specific:
            name = found.get("name")
            if not isinstance(name, str):
                continue
            known = tissue_name(name) if key == "tissue" else single(name)
            if known is None:
                out.append(f"{where}: {key} {name!r} is not a name the Atlas lists")
    pairs = expression.get("tissue_cell_type")
    if not isinstance(pairs, list):
        out.append(f"{where}: tissue_cell_type names no list")
        pairs = []
    for pair in pairs:
        tissue, cell = (pair or {}).get("tissue"), (pair or {}).get("cell_type")
        if not isinstance(tissue, str) or pair_tissue(tissue) is None:
            out.append(f"{where}: tissue cell type tissue {tissue!r} is not one the Atlas lists")
        if not isinstance(cell, str) or not cell:
            out.append(f"{where}: tissue cell type in {tissue!r} names no cell type")
    subcellular = expression.get("subcellular")
    if not isinstance(subcellular, dict):
        out.append(f"{where}: subcellular names no locations")
    else:
        for part in ("main", "additional"):
            locations = subcellular.get(part)
            if not isinstance(locations, list):
                out.append(f"{where}: subcellular {part} is not a list")
                continue
            for location in locations:
                if location not in SUBCELLULAR:
                    out.append(f"{where}: subcellular {part} {location!r} is not the Atlas's word")
    secretome = expression.get("secretome", "absent")
    if secretome is not None and secretome not in SECRETOME:
        out.append(f"{where}: secretome {secretome!r} is not the Atlas's word")

    # The zoom's path, as the rule takes it from those readings.
    path = track.get("path")
    if not isinstance(path, dict):
        out.append(f"{where}: no path")
    else:
        try:
            want = path_of(expression)
        except (KeyError, LookupError, TypeError, ValueError) as exc:
            out.append(f"{where}: the path's readings do not resolve: {exc}")
        else:
            if path != want:
                out.append(f"{where}: path is {path!r}, the rule takes {want!r}")
            if path.get("cell_type") is not None and path.get("cell_class") not in CLASSES:
                out.append(f"{where}: path cell class {path.get('cell_class')!r}")

    # Where it came from.
    sources = track.get("sources") or {}
    for key, table in (("cytoband", "cytoBand"), ("chrom_alias", "chromAlias")):
        found = sources.get(key) or {}
        if found.get("table") != table or found.get("genome") != GENOME:
            out.append(f"{where}: {key} is {found.get('table')!r} on {found.get('genome')!r}")
        if not _UPDATED.match(str(found.get("updated"))):
            out.append(f"{where}: {key} updated {found.get('updated')!r}")
        if not _DIGEST.match(str(found.get("sha256"))):
            out.append(f"{where}: {key} digest {found.get('sha256')!r}")
    if not _MANE.match(str((sources.get("mane") or {}).get("release"))):
        out.append(f"{where}: MANE release {(sources.get('mane') or {}).get('release')!r}")
    uniprot = sources.get("uniprot") or {}
    if not _RELEASE.match(str(uniprot.get("release"))):
        out.append(f"{where}: UniProt release {uniprot.get('release')!r}")
    hpa = sources.get("hpa") or {}
    if not _VERSION.match(str(hpa.get("version"))):
        out.append(f"{where}: Protein Atlas version {hpa.get('version')!r}")
    if not _DAY.match(str(hpa.get("release_date"))):
        out.append(f"{where}: Protein Atlas release date {hpa.get('release_date')!r}")
    if not str(hpa.get("ensembl", "")).isdigit():
        out.append(f"{where}: Protein Atlas Ensembl version {hpa.get('ensembl')!r}")
    if hpa.get("licence") != HPA_LICENCE:
        out.append(f"{where}: Protein Atlas licence {hpa.get('licence')!r}")
    if not _DAY.match(str(track.get("retrieved"))):
        out.append(f"{where}: retrieved {track.get('retrieved')!r}")
    return out


def main() -> int:
    found = [problem for target in LOCUS_TARGETS for problem in problems_of(target)]
    for problem in found:
        print(f"  {problem}", file=sys.stderr)
    print(f"{len(LOCUS_TARGETS)} locus tracks checked, {len(found)} problem(s)",
          file=sys.stderr)
    return 1 if found else 0


if __name__ == "__main__":
    raise SystemExit(main())
