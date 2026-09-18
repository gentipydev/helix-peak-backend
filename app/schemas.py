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

from typing import List, Optional

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
