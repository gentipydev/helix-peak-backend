"""Response models for the structured gene endpoint.

Coordinates are 1-based inclusive throughout, matching what NCBI shows for the
record rather than Biopython's 0-based half-open internals.

The payload carries only what the app draws. Lengths a client can compute from
what is already here -- a span's ``end - start + 1``, a translation's
``len()`` -- are not sent, so there is no second copy of a number to disagree
with the first.

Annotations stay in ``typing`` form (``List``/``Optional``, not ``list[...]`` or
``X | Y``): Pydantic evaluates them at runtime and this runs on Python 3.9.
"""

from typing import Dict, List, Optional

from pydantic import BaseModel


class Segment(BaseModel):
    """One contiguous stretch of a feature, i.e. one part of a ``join(...)``."""

    start: int
    end: int


class Location(BaseModel):
    start: int
    end: int
    strand: Optional[int] = None


class Exon(BaseModel):
    number: Optional[int] = None
    start: int
    end: int


class Transcript(BaseModel):
    segments: List[Segment] = []


class Protein(BaseModel):
    product: Optional[str] = None
    translation: str
    segments: List[Segment] = []


class Peptide(BaseModel):
    product: Optional[str] = None
    segments: List[Segment] = []
    translation: str


class GeneResponse(BaseModel):
    """One gene lifted out of a GenBank record."""

    gene: str

    location: Location
    sequence: str

    transcript: Optional[Transcript] = None
    protein: Optional[Protein] = None
    exons: List[Exon] = []
    signal_peptide: Optional[Peptide] = None
    proprotein: Optional[Peptide] = None
    peptides: List[Peptide] = []


# ---------------------------------------------------------------------------
# The catalog.
#
# These mirror the Dart entities in `protein_target.dart` field for field, so a
# `ProteinDetail` deserialises straight into a `ProteinTarget` and
# `tool/check_assets.py --against <base-url>` can diff the two against
# `tool/targets.py` without a translation layer in between.
# ---------------------------------------------------------------------------


class Track(BaseModel):
    """Where one track for one protein is, and what state it is in.

    ``url`` points at Supabase storage, which the client fetches directly --
    this service names the bytes and never carries them.

    ``provenance`` is free-form on purpose. A constraint track names its model,
    revision, method and vocabulary; an impact track names its scorer and
    units; a ClinVar snapshot names its retrieval date and what it excluded; a
    structure names its entry, resolution and deviations. The About sheet reads
    whichever applies rather than holding a constant per source, because two
    constraint tracks can now come from two different models.
    """

    state: str
    reason: Optional[str] = None

    url: Optional[str] = None
    format: Optional[str] = None
    bytes: Optional[int] = None
    sha256: Optional[str] = None
    content_encoding: Optional[str] = None
    provenance: dict = {}


class ProteinFacts(BaseModel):
    """The four numbers on a search card. Counts, not lengths of anything sent."""

    residues: int
    exons: int
    chains: int
    bridges: int


class Region(BaseModel):
    """A named stretch of the precursor, 1-based inclusive, UniProt numbering.

    ``kept`` is whether the stretch is still there once the precursor has been
    cut: tolerance in a piece that gets thrown away is not the same fact about
    the molecule as tolerance in the part that goes on to work.
    """

    label: str
    short: str
    start: int
    end: int
    origin: int
    kept: bool = True


class StructureChain(BaseModel):
    """One node of the baked model, and the colour it takes.

    ``node`` is the contract with ``structure_view.dart``, which looks the name
    up to give each node a material. The tints run in the order the record
    lists the mature peptides, which is why insulin's B chain is the first and
    its A chain the third.
    """

    node: str
    tint: str


class StructureChrome(BaseModel):
    """What the fold page says about the entry it is drawn from."""

    pdb: str
    modelled: Optional[List[int]] = None
    label: str
    count: int
    unit: str
    sentence: str
    semantics: str


class ProteinSummary(BaseModel):
    """What a search result needs, and nothing more.

    ``tracks`` carries state strings only. A card shows what is ready; it does
    not need a URL for a track nobody has opened.
    """

    slug: str
    display: str
    gene: str
    uniprot: str
    accession: str
    summary: str
    facts: ProteinFacts

    # Where this protein sits in the catalog's reading order. It used to be the
    # index of a const list; a list the client no longer holds cannot carry it,
    # so it travels per row. Null for a protein resolved on demand, which has
    # no place in a sequence someone chose and sorts after all of them.
    catalog_order: Optional[int] = None

    chain: Optional[str] = None
    chains: List[StructureChain] = []
    structure: Optional[StructureChrome] = None
    tracks: Dict[str, str] = {}


class ProteinDetail(ProteinSummary):
    """One resolved protein, whole: a `targets.py` row and a catalog row at once."""

    mature_peptides: bool = True
    transcript_id: Optional[str] = None
    protein_id: Optional[str] = None

    regions: List[Region] = []
    disulfides: List[List[int]] = []

    # How this row was arrived at: per-field source, the positions where the
    # record and UniProt disagree, and whether the prose was written or
    # templated. A resolved protein must be able to say how it was resolved.
    provenance: dict = {}
    resolver_version: int = 0


class CatalogPage(BaseModel):
    proteins: List[ProteinSummary] = []
    # The slug to pass back as ``cursor``; null on the last page.
    next: Optional[str] = None


class Candidate(BaseModel):
    """A protein the catalog does not hold yet, offered for resolution.

    Kept separate from ``proteins`` in the search response so the screen can
    keep saying only what it knows: these have not been resolved, and nothing
    about them has been checked.
    """

    gene: str
    name: str
    taxon: int
    source: str


class SearchResponse(BaseModel):
    proteins: List[ProteinSummary] = []
    candidates: List[Candidate] = []


class Suggestion(BaseModel):
    """One protein a reader might mean, and what the app can do with it.

    ``status`` is ``listed`` (one of the twenty), ``ready`` (built earlier),
    ``buildable`` or ``unavailable``; ``reason`` says why for the last. ``slug``
    is where a listed or ready protein opens, and the slug a build would give a
    buildable one.
    """

    uniprot: str
    gene: Optional[str] = None
    name: str
    # The app's own name for a listed protein ("Hemoglobin (beta chain)").
    display: Optional[str] = None
    length: int
    slug: Optional[str] = None
    status: str
    reason: Optional[str] = None


class SuggestResponse(BaseModel):
    q: str
    # "UniProt 2026_03 · MANE v1.5": what the index was built from.
    release: Optional[str] = None
    suggestions: List[Suggestion] = []


class ResolveRequest(BaseModel):
    gene: str
    taxon: int = 9606


class ResolveResponse(BaseModel):
    slug: Optional[str] = None
    state: str
    reason: Optional[str] = None
