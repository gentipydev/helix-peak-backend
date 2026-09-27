"""Bake `locus`: where on its chromosome each protein's gene lies, by band.

    .venv/Scripts/python pipeline/locus/bake_locus.py --all
    .venv/Scripts/python pipeline/locus/bake_locus.py --target insulin

The app's zoom ends on a chromosome with the gene marked at its band. This
track carries what that needs, per protein:

- `span`: where the gene lies on GRCh38, 1-based and inclusive, and its
  strand. It is the MANE Select transcript's span, found the way the protein
  index finds it: the protein's UniProt entry joined to MANE by
  `rows_for_entry` (`app/protein_index.py`), the row for the table's gene.
- `bands`: every cytogenetic band of that chromosome, in order from the end
  of the short arm, each with its Giemsa stain, from the UCSC Genome Browser's
  cytoBand table for hg38. UCSC's 0-based starts are made 1-based, so every
  position in the track is counted the way MANE and GenBank count.
- `band`: the band or bands the span lies in, and where they begin and end;
  and `locus`, the place as cytogenetics writes it: 11p15.5, or Xp21.2-p21.1
  for a gene that runs across two bands.

A band is a stain pattern seen down a microscope at low resolution, millions
of bases long. It says where a gene lies, not that the gene can be seen there.

Each payload names its sources and their versions: the cytoBand and
chromAlias tables by genome and by the time UCSC's API says each was last
updated, with a digest of the rows read; the MANE release; the UniProt release
and entry version. UCSC's download server does not answer from every network,
so the tables are read through its REST API, which serves the same tables.

Which proteins get the track: all twenty (`LOCUS_TARGETS`), recorded here,
never in targets.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from app.protein_index import rows_for_entry  # noqa: E402
from pipeline.mane import current_summary  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import BY_SLUG, TARGETS, Target  # noqa: E402
from pipeline.uniprot import Entry, fetch_entry  # noqa: E402

# Every protein gets the track: each has a chromosome to zoom to.
LOCUS_TARGETS = TARGETS

SCHEMA_VERSION = 1

ASSEMBLY = "GRCh38"
GENOME = "hg38"
UCSC_API = "https://api.genome.ucsc.edu/getData/track?genome={};track={}"

_AGENT = "helixpeek-locus-bake"

# The stains cytoBand gives a band: G-negative, four depths of G-positive,
# the centromere, the variable heterochromatin, and the stalks of the
# acrocentric chromosomes.
STAINS = ("gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar", "stalk")


@dataclass(frozen=True)
class Table:
    """One UCSC table as its API served it."""

    track: str
    genome: str
    # When UCSC last updated the table, as its API reports it: the table's
    # version, since UCSC tables carry no other.
    updated: str
    rows: tuple[dict, ...]
    url: str

    @property
    def sha256(self) -> str:
        """A digest of the rows, in a canonical order: the same table gives the
        same digest however the API happens to order its answer."""
        canonical = sorted(json.dumps(row, sort_keys=True) for row in self.rows)
        return hashlib.sha256("\n".join(canonical).encode()).hexdigest()

    def provenance(self) -> dict:
        return {
            "table": self.track,
            "genome": self.genome,
            "updated": self.updated,
            "rows": len(self.rows),
            "sha256": self.sha256,
            "url": self.url,
        }


def locus_asset(target: Target) -> str:
    return f"assets/locus/{target.slug}_locus.json"


def _get(url: str) -> tuple[bytes, dict[str, str]]:
    request = urllib.request.Request(url, headers={"User-Agent": _AGENT})
    with urllib.request.urlopen(request, timeout=300) as response:
        return response.read(), {k.lower(): v for k, v in response.headers.items()}


def table_of(answer: dict, url: str) -> Table:
    """A UCSC API answer as a table: its rows, flattened out of the per-
    chromosome lists the API groups a whole-genome answer into."""
    track = answer["track"]
    found = answer[track]
    rows = found if isinstance(found, list) else [
        row for chrom in found.values() for row in chrom]
    return Table(
        track=track,
        genome=answer["genome"],
        updated=answer["dataTime"],
        rows=tuple(rows),
        url=url,
    )


def fetch_table(track: str) -> Table:
    url = UCSC_API.format(GENOME, track)
    body, _ = _get(url)
    return table_of(json.loads(body), url)


def chromosome_of(accession: str, aliases: Table) -> str:
    """UCSC's name for a RefSeq sequence ("NC_000011.10" -> "chr11"), as its
    own chromAlias table gives it."""
    for row in aliases.rows:
        if row["alias"] == accession and row.get("source") == "refseq":
            return row["chrom"]
    raise LookupError(f"UCSC's chromAlias names no chromosome {accession}")


def mane_row(target: Target, entry: Entry, mane: dict, release: str) -> dict:
    """The protein index's row for [target]: its entry joined to MANE, the
    buildable row for the table's gene."""
    rows = [row for row in rows_for_entry(entry.body, mane, {}, release)
            if row["buildable"] and row["gene"] == target.gene]
    if len(rows) != 1:
        raise LookupError(
            f"{target.slug}: {len(rows)} buildable MANE Select rows for {target.gene} "
            f"in {entry.accession}")
    return rows[0]


def bands_of(chrom: str, cytoband: Table) -> list[dict]:
    """Every band of [chrom], from the end of the short arm, 1-based."""
    found = sorted(
        (row for row in cytoband.rows if row["chrom"] == chrom),
        key=lambda row: row["chromStart"])
    return [{
        "name": row["name"],
        "start": row["chromStart"] + 1,
        "end": row["chromEnd"],
        "stain": row["gieStain"],
    } for row in found]


def locus_of(chromosome: str, names: list[str]) -> str:
    """How cytogenetics writes the place: 11p15.5; Xp21.2-p21.1 across two."""
    return chromosome + names[0] + ("" if len(names) == 1 else "-" + names[-1])


def payload(target: Target, entry: Entry, mane: dict, mane_release: str,
            cytoband: Table, aliases: Table, retrieved: str) -> dict:
    """The track for one protein, from the sources as they were read."""
    body = entry.body
    if body.get("primaryAccession") != target.uniprot:
        raise ValueError(f"{target.slug}: asked for {target.uniprot}, "
                         f"UniProt answered {body.get('primaryAccession')}")
    row = mane_row(target, entry, mane, mane_release)
    chrom = chromosome_of(row["chrom_acc"], aliases)
    bands = bands_of(chrom, cytoband)
    if not bands:
        raise LookupError(f"{target.slug}: cytoBand has no bands for {chrom}")
    start, end = row["chrom_start"], row["chrom_end"]
    inside = [b for b in bands if b["start"] <= end and b["end"] >= start]
    if not inside:
        raise LookupError(f"{target.slug}: {chrom}:{start}-{end} falls in no band")
    chromosome = chrom[3:] if chrom.startswith("chr") else chrom
    names = [b["name"] for b in inside]
    audit = body.get("entryAudit", {})
    return {
        "slug": target.slug,
        "gene": target.gene,
        "uniprot": target.uniprot,
        "assembly": ASSEMBLY,
        "genome": cytoband.genome,
        "chromosome": chromosome,
        "sequence": row["chrom_acc"],
        "length": bands[-1]["end"],
        "span": {"start": start, "end": end, "strand": row["chrom_strand"]},
        "transcript": {"refseq": row["refseq_nuc"], "ensembl": row["ensembl_nuc"]},
        "band": {"names": names, "start": inside[0]["start"], "end": inside[-1]["end"]},
        "locus": locus_of(chromosome, names),
        "bands": bands,
        "sources": {
            "cytoband": cytoband.provenance(),
            "chrom_alias": aliases.provenance(),
            "mane": {"release": mane_release},
            "uniprot": {
                "release": entry.release,
                "release_date": entry.release_date,
                "entry_version": audit.get("entryVersion"),
            },
        },
        "retrieved": retrieved,
        "schema_version": SCHEMA_VERSION,
        "built_by": "pipeline/locus/bake_locus.py",
    }


def encode(track: dict) -> bytes:
    return (json.dumps(track, indent=2, ensure_ascii=False) + "\n").encode()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if not args.all and not args.target:
        parser.error("name --target or --all")
    unknown = sorted(set(args.target or []) - set(BY_SLUG))
    if unknown:
        parser.error(f"no such target: {unknown}")
    chosen = LOCUS_TARGETS if args.all else tuple(BY_SLUG[s] for s in args.target)

    cytoband = fetch_table("cytoBand")
    aliases = fetch_table("chromAlias")
    mane, mane_release = current_summary(_get)
    retrieved = datetime.now(timezone.utc).date().isoformat()
    print(f"cytoBand {cytoband.genome}, updated {cytoband.updated}, {len(cytoband.rows)} rows; "
          f"chromAlias updated {aliases.updated}; MANE {mane_release}")

    failed = []
    for target in chosen:
        try:
            track = payload(target, fetch_entry(target.uniprot), mane, mane_release,
                            cytoband, aliases, retrieved)
        except (LookupError, ValueError) as exc:
            failed.append(f"{target.slug}: {exc}")
            print(f"{target.slug:<16} FAILED: {exc}")
            continue
        destination = DATA / locus_asset(target)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(encode(track))
        span = track["span"]
        print(f"{target.slug:<16} {track['gene']:<6} {track['locus']:<14} "
              f"{track['sequence']}:{span['start']}-{span['end']}")
    if failed:
        print(f"\n{len(failed)} protein(s) did not resolve to a band:")
        for line in failed:
            print("  " + line)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
