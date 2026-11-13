"""The locus bake's inputs, pinned: one snapshot of every source it reads.

Every `locus` track is made from this one file and from nothing fetched at
bake time: a curated protein's re-bake, a protein built on demand, the census
of every protein the search offers. The Human Protein Atlas's site serves its
newest release (26.0 since November 2026), and a bake that read it gene by
gene would move a curated zoom to another release without anyone deciding to.
A move to another release is a new snapshot, pinned here on purpose.

The snapshot (`snapshot.py` builds it) holds, from pinned downloads:

- the Atlas 25.1's reading of every gene it has an entry for: the fields
  `bake_locus.expression_of` reads, the gene's consensus RNA level in each
  tissue the Atlas measured, and its level in each single cell type;
- UCSC's hg38 cytoBand, chromAlias, altLocations and fixLocations tables, as
  its API serves them;
- NCBI's placement of each GRCh38.p14 patch and alternate locus scaffold on
  its chromosome;
- MANE v1.5's Select rows;
- HGNC's Ensembl gene for each HGNC id, from a dated archive.

It is gzipped text: a header line of JSON, then one line per Atlas gene,
`["ENSG…",{…}]`, in Ensembl order, of which only the genes asked for are
ever parsed. Its name carries the first twelve hex digits of the sha256 of
its text, and `INPUTS_SHA256` pins the whole text: `load` refuses any other
file. It lives under `pipeline/data/`, which git ignores; the digest here is
what makes it reproducible, since `snapshot.py` builds the same text from the
same sources on the same date (`INPUTS_CREATED`).

Where it is: `HELIXPEEK_LOCUS_INPUTS` (the file, or the folder it is in), else
`pipeline/data/locus_inputs/` beside this code, never `paths.DATA`, which a
worker points at a directory of its own.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

# The snapshot every locus is made from: the sha256 of its uncompressed text,
# and the date it was made on, which every track made from it carries as
# `retrieved`.
INPUTS_SHA256 = "bb475d41dcf93379c8c1af63026bef9af35495abbb658e4416d7c7d05c9ed4b9"
INPUTS_CREATED = "2026-11-13"

ENV = "HELIXPEEK_LOCUS_INPUTS"
FOLDER = Path(__file__).resolve().parents[1] / "data" / "locus_inputs"

KIND = "helixpeek-locus-inputs"
FORMAT = 1

_GENE_LINE = re.compile(rb'^\["(ENSG\d{11})",')


class InputsUnavailable(RuntimeError):
    """No snapshot, or not the pinned one: the message says which, and where."""


def file_name(digest: str) -> str:
    return f"locus-inputs.{digest[:12]}.jsonl.gz"


def default_path(digest: str = INPUTS_SHA256) -> Path:
    """Where the snapshot pinned as [digest] is looked for."""
    named = os.environ.get(ENV)
    if named:
        place = Path(named)
        return place / file_name(digest) if place.is_dir() else place
    return FOLDER / file_name(digest)


@dataclass(frozen=True)
class Inputs:
    """One snapshot, read: its header, and its genes' lines, unparsed."""

    path: Path
    sha256: str
    header: dict
    genes: dict = field(repr=False)

    @property
    def created(self) -> str:
        return self.header["created"]

    @property
    def sources(self) -> dict:
        return self.header["sources"]

    def table(self, track: str) -> dict:
        """One UCSC table as its API served it: `genome`, `updated`, `url` and
        `rows`."""
        return self.header["tables"][track]

    @property
    def mane(self) -> dict:
        """MANE's Select rows, keyed by versionless ENST (`parse_mane`)."""
        return self.header["mane"]

    @property
    def mane_release(self) -> str:
        return self.sources["mane"]["release"]

    @property
    def atlas(self) -> dict:
        """The Atlas release: `version`, `release_date`, `ensembl`, `licence`
        and `site`."""
        return self.sources["hpa"]

    def atlas_entry(self, ensembl: str) -> dict | None:
        """The Atlas's whole line for a gene: `atlas` (its fields, as its own
        answer names them), `tissues` and `cells` (levels, in the order
        `consensus_tissues` and `cell_types` name them). None where the Atlas
        has no entry for it."""
        line = self.genes.get(ensembl)
        if line is None:
            return None
        found, entry = json.loads(line)
        if found != ensembl:
            raise InputsUnavailable(f"{self.path}: the line for {ensembl} names {found}")
        return entry

    def atlas_gene(self, ensembl: str) -> dict | None:
        """The Atlas's answer for a gene, as its own JSON names the fields
        the bake reads; None where it has no entry."""
        entry = self.atlas_entry(ensembl)
        return None if entry is None else entry["atlas"]

    def levels(self, ensembl: str) -> tuple[dict, dict] | None:
        """A gene's consensus RNA level in each tissue the Atlas measured
        (nTPM) and in each single cell type (nCPM), by the names its tables
        write; None where it has no entry."""
        entry = self.atlas_entry(ensembl)
        if entry is None:
            return None
        tissues = dict(zip(self.header["consensus_tissues"], entry["tissues"]))
        cells = dict(zip(self.header["cell_types"], entry["cells"]))
        return tissues, cells

    def by_uniprot(self, accession: str) -> tuple[str, ...]:
        """The Atlas genes whose entries name this UniProt accession."""
        return tuple(self.header["uniprot"].get(accession, ()))

    def hgnc_ensembl(self, hgnc_id: str | None) -> str | None:
        """The Ensembl gene HGNC gives an HGNC id (`HGNC:5`)."""
        return self.header["hgnc"].get(hgnc_id or "")

    def placements(self, accession: str) -> list[dict]:
        """NCBI's placements of a GRCh38.p14 scaffold on its chromosome, by the
        scaffold's RefSeq accession (`NW_016107299.1`)."""
        return list(self.header["placements"].get(accession, ()))


def read(path: Path, expected: str | None = INPUTS_SHA256) -> Inputs:
    """The snapshot at [path], held to [expected] (None takes any)."""
    if not path.is_file():
        raise InputsUnavailable(
            f"no locus inputs snapshot at {path}: build it with "
            f"`python pipeline/locus/snapshot.py` (it is not in git)")
    digest = hashlib.sha256()
    header = None
    genes: dict[str, bytes] = {}
    with gzip.open(path, "rb") as stream:
        for number, line in enumerate(stream):
            digest.update(line)
            body = line.rstrip(b"\n")
            if number == 0:
                header = json.loads(body)
                continue
            found = _GENE_LINE.match(body)
            if found is None:
                raise InputsUnavailable(f"{path}: line {number + 1} is not a gene's")
            genes[found.group(1).decode()] = body
    sha256 = digest.hexdigest()
    if expected is not None and sha256 != expected:
        raise InputsUnavailable(
            f"{path} is the locus inputs snapshot {sha256[:12]}, not the pinned "
            f"{expected[:12]}; rebuild it with `python pipeline/locus/snapshot.py`")
    if not isinstance(header, dict) or header.get("kind") != KIND \
            or header.get("format") != FORMAT:
        raise InputsUnavailable(f"{path} is not a locus inputs snapshot of format {FORMAT}")
    return Inputs(path=path, sha256=sha256, header=header, genes=genes)


_LOADED: dict[tuple, Inputs] = {}


def load(path: Path | None = None, expected: str | None = INPUTS_SHA256) -> Inputs:
    """The pinned snapshot, read once per process. Raises `InputsUnavailable`
    where it is missing or is not the pinned one."""
    where = Path(path) if path is not None else default_path(expected or INPUTS_SHA256)
    key = (str(where), expected)
    if key not in _LOADED:
        _LOADED[key] = read(where, expected)
    return _LOADED[key]
