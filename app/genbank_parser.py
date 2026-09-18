"""Extract one gene's features from a parsed GenBank record."""

from typing import Any, Dict, List, Optional

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


class GeneNotFound(Exception):
    """Raised when no feature in the record carries the requested gene name."""


def _qualifier(feature: SeqFeature, key: str) -> Optional[str]:
    """Return the first value of ``key``, or None.

    Every qualifier Biopython parses is a list, even single-valued ones:
    ``feature.qualifiers["gene"]`` is ``["INS"]``, never ``"INS"``.
    """
    values = feature.qualifiers.get(key)
    if not values:
        return None
    return values[0]


def _span(feature: SeqFeature) -> Dict[str, int]:
    """Return the feature's outer span as a 1-based inclusive range.

    Biopython locations are 0-based half-open, so 5224..5410 in the flat file
    arrives as start=5223, end=5410. GenBank's own numbering is what a reader
    sees on NCBI, so that is what this service reports.
    """
    return {
        "start": int(feature.location.start) + 1,
        "end": int(feature.location.end),
    }


def _segments(feature: SeqFeature) -> List[Dict[str, int]]:
    """Return the feature's parts, so a ``join(...)`` keeps its introns.

    ``int(location.start)``/``end`` collapse a CompoundLocation to its outer
    span, which would silently splice the introns out of the INS CDS.
    """
    return [
        {"start": int(part.start) + 1, "end": int(part.end)}
        for part in feature.location.parts
    ]


def _translation(feature: SeqFeature, record: SeqRecord) -> str:
    """Return the feature's protein sequence.

    Prefers the record's own ``/translation`` (authoritative, and already
    accounts for ``/codon_start``); mat_peptide and sig_peptide carry none, so
    those are translated from the extracted nucleotides.
    """
    stated = _qualifier(feature, "translation")
    if stated is not None:
        return stated
    return str(feature.extract(record.seq).translate())


def _peptide(feature: SeqFeature, record: SeqRecord) -> Dict[str, Any]:
    """Return one sig_peptide, proprotein or mat_peptide.

    Only the parts are reported, never the outer span: a peptide that spans an
    intron is drawn from its segments, and the span its segments already imply
    would be a second copy of the same fact.
    """
    return {
        "product": _qualifier(feature, "product"),
        "segments": _segments(feature),
        "translation": _translation(feature, record),
    }


def _exon(feature: SeqFeature) -> Dict[str, Any]:
    number = _qualifier(feature, "number")
    exon = _span(feature)
    exon["number"] = int(number) if number is not None and number.isdigit() else None
    return exon


def extract_gene(record: SeqRecord, gene: str) -> Dict[str, Any]:
    """Return the structured payload for ``gene`` within ``record``.

    Matching is an exact comparison against the ``/gene`` qualifier, never a
    prefix or substring test. A RefSeqGene region carries several genes and they
    overlap: in NG_007114 the INS-IGF2 readthrough starts at the very same base
    as INS and shares its signal peptide, so ``"INS-IGF2".startswith("INS")``
    would fold the wrong features into the answer.

    Pure: takes an already-parsed record and does no I/O, so it can be tested
    straight off a saved fixture.
    """
    features = [f for f in record.features if _qualifier(f, "gene") == gene]
    if not features:
        raise GeneNotFound(gene)

    by_type = {}  # type: Dict[str, List[SeqFeature]]
    for feature in features:
        by_type.setdefault(feature.type, []).append(feature)

    anchor = by_type.get("gene", features)[0]
    location = _span(anchor)
    location["strand"] = anchor.location.strand

    transcript = None
    for feature in by_type.get("mRNA", []):
        transcript = {"segments": _segments(feature)}
        break

    protein = None
    for feature in by_type.get("CDS", []):
        protein = {
            "product": _qualifier(feature, "product"),
            "translation": _translation(feature, record),
            "segments": _segments(feature),
        }
        break

    exons = [_exon(f) for f in by_type.get("exon", [])]
    # Exon 3 of INS trails the mat_peptides in the flat file, so file order is
    # not transcript order.
    exons.sort(key=lambda exon: (exon["number"] is None, exon["number"] or 0, exon["start"]))

    signal_peptide = None
    for feature in by_type.get("sig_peptide", []):
        signal_peptide = _peptide(feature, record)
        break

    proprotein = None
    for feature in by_type.get("proprotein", []):
        proprotein = _peptide(feature, record)
        break

    return {
        "gene": gene,
        "location": location,
        "sequence": str(anchor.extract(record.seq)),
        "transcript": transcript,
        "protein": protein,
        "exons": exons,
        "signal_peptide": signal_peptide,
        "proprotein": proprotein,
        "peptides": [_peptide(f, record) for f in by_type.get("mat_peptide", [])],
    }
