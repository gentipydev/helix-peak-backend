"""Bake `assets/mock/gene_*.json` -- the gene records the walk opens on.

The records are stored as the `record` track (Phase 2 of HANDOFF-ONDEMAND.md):
this writes them under `pipeline/data/`, and `upload_tracks.py --kind record`
sends them. The file names keep the `mock` they were born with, when the app
bundled them for a mode with no backend.

Each file is exactly what `GET /gene/{id}/{gene}` would answer, because the
parsing is the backend's own: this imports `app.genbank_parser.extract_gene`
from this repo rather than reimplementing it, so a record cannot drift from the
contract it stands in for.

Six things happen here that the live service does not do yet, all noted at
their call sites: the transcript is chosen where a record holds several, exons
are synthesised from the mRNA feature for a record that annotates no `exon`
features, a connecting peptide, a proprotein and a lone chain the record leaves
out are filled in from the table, and a gene over the gene page's budget has its
introns scaled. All are recorded in the payload, and all are things the backend
should learn to do before these records are replaced by live calls.

    .venv/bin/python pipeline/mock/build_gene_record.py --all
    .venv/bin/python pipeline/mock/build_gene_record.py --target dystrophin
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import sys
import tempfile
import time
import urllib.request

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import GENE_PAGE_BUDGET_BP, TARGETS, Target  # noqa: E402

from Bio import Entrez, SeqIO  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from app.genbank_parser import extract_gene  # noqa: E402

Entrez.tool = "helixpeek-mockbake"
Entrez.email = os.environ.get("NCBI_EMAIL", "")

MIN_INTRON_BP = 60


# ----------------------------------------------------------------- fetching


CACHE = Path(
    os.environ.get("HELIXPEEK_GB_CACHE", Path(tempfile.gettempdir()) / "helixpeek-genbank")
)


def fetch(target: Target):
    """The GenBank record, whole or sliced, through a flat-file cache.

    A slice is how the two genes with no RefSeqGene are reached. NCBI re-bases a
    sliced record to 1, so every coordinate below is already relative to the
    slice and nothing else has to know it was one.

    The cache is not an optimisation for its own sake: dystrophin's RefSeqGene
    is 2.2 Mb of flat file, and re-fetching it on every run would make this
    script painful to iterate on and rude to NCBI.
    """
    source = target.source
    stem = source.accession.replace("/", "_")
    if source.seq_start is not None:
        stem += f"_{source.seq_start}_{source.seq_stop}"
    cached = CACHE / f"{stem}.gb"
    if not cached.exists():
        kwargs = dict(db="nucleotide", id=source.accession, rettype="gb", retmode="text")
        if source.seq_start is not None:
            kwargs["seq_start"] = source.seq_start
            kwargs["seq_stop"] = source.seq_stop
        handle = Entrez.efetch(**kwargs)
        try:
            body = handle.read()
        finally:
            handle.close()
        cached.parent.mkdir(parents=True, exist_ok=True)
        cached.write_text(body)
    with cached.open() as handle:
        record = SeqIO.read(handle, "genbank")
    return select_isoform(record, target)


def _covers(outer: list, inner: list) -> bool:
    """Whether every part of `inner` sits inside some part of `outer`."""
    return all(
        any(int(o.start) <= int(i.start) and int(i.end) <= int(o.end) for o in outer)
        for i in inner
    )


def select_isoform(record, target: Target):
    """Drop every transcript of this gene but the one the table names.

    A RefSeqGene holds one transcript and this does nothing. A chromosome slice
    holds all of them — RLN2 has six CDS features, four of them predicted `XP_`
    models — and `extract_gene` takes whichever comes first in the file, which
    is a coin toss it was never asked to win. Choosing here rather than in the
    parser keeps the parser the backend's, and puts "which isoform" where the
    rest of the per-protein decisions live.

    The mRNA usually need not be named in the table: the one that contains every
    segment of the chosen CDS is the one that codes for it. Not when two
    transcripts differ only in what they never translate. AMY1A's slice holds
    two mRNAs around one set of CDS coordinates, and the shortest of them is not
    the MANE Select one, so the table names it (`Source.transcript_id`).
    """
    wanted = target.source.protein_id
    named = target.source.transcript_id
    if named is not None and wanted is None:
        raise ValueError(f"{target.slug}: a transcript_id needs the protein_id it codes for")
    mine = [f for f in record.features if f.qualifiers.get("gene", [None])[0] == target.gene]
    coding = [f for f in mine if f.type == "CDS"]
    if wanted is None or len(coding) <= 1:
        return record

    chosen = [f for f in coding if wanted in f.qualifiers.get("protein_id", [])]
    if len(chosen) != 1:
        available = sorted(p for f in coding for p in f.qualifiers.get("protein_id", []))
        raise ValueError(f"{target.slug}: {wanted!r} is not one of {available}")
    cds = chosen[0]

    transcripts = [
        f for f in mine
        if f.type == "mRNA" and _covers(f.location.parts, cds.location.parts)
    ]
    if named is not None:
        candidates = sorted(p for f in transcripts for p in f.qualifiers.get("transcript_id", []))
        transcripts = [f for f in transcripts if named in f.qualifiers.get("transcript_id", [])]
        if len(transcripts) != 1:
            raise ValueError(
                f"{target.slug}: {named!r} is not one of the mRNAs holding {wanted!r}: {candidates}"
            )
    transcript = min(transcripts, key=lambda f: len(f.location), default=None)

    def keep(feature) -> bool:
        if feature.qualifiers.get("gene", [None])[0] != target.gene:
            return True
        if feature.type == "CDS":
            return feature is cds
        if feature.type == "mRNA":
            return feature is transcript
        if feature.type in ("exon", "mat_peptide", "sig_peptide", "proprotein"):
            return transcript is None or _covers(transcript.location.parts, feature.location.parts)
        return True

    record.features = [f for f in record.features if keep(f)]
    return record


def canonical_sequence(uniprot: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot}.json"
    with urllib.request.urlopen(url, timeout=60) as response:
        return json.load(response)["sequence"]["value"]


# --------------------------------------------------------------- tidying

# NCBI writes a UniProt feature's whole note into `/product` on some records:
# oxytocin's arrives as "Oxytocin. /evidence=ECO:0000269|PubMed:13591312.
# /id=PRO_0000020495". The app draws `/product` as a label on the molecule, so
# the trailing qualifiers have to go. The live backend will need the same cut
# before it can serve these ten; it has never met a record that needed it.
_QUALIFIER_BLEED = re.compile(r"\s*\.?\s*/[A-Za-z_]+=.*$", re.S)


def clean_product(value):
    if value is None:
        return None
    cleaned = _QUALIFIER_BLEED.sub("", value).strip().rstrip(".").strip()
    return cleaned or None


def tidy(payload: dict, target: Target) -> dict:
    """Clean the labels and drop the peptides that are not a cleavage series."""
    for key in ("protein", "signal_peptide", "proprotein"):
        if payload.get(key):
            payload[key]["product"] = clean_product(payload[key]["product"])

    peptides = []
    seen = set()
    for peptide in payload["peptides"]:
        peptide["product"] = clean_product(peptide["product"])
        # RLN2's B chain is annotated twice, identically. One page cannot show
        # the same piece of the molecule cut off twice.
        fingerprint = (peptide["translation"], tuple(
            (s["start"], s["end"]) for s in peptide["segments"]
        ))
        if fingerprint in seen:
            continue
        seen.add(fingerprint)
        peptides.append(peptide)

    payload["peptides"] = peptides if target.mature_peptides else []
    if target.mature_peptides:
        fill_removed_peptides(payload, target)
        fill_chain(payload, target)
    fill_proprotein(payload, target)

    spans = sorted((s["start"], s["end"]) for p in payload["peptides"] for s in p["segments"])
    for before, after in zip(spans, spans[1:]):
        if after[0] <= before[1]:
            raise ValueError(
                f"{target.slug}: mature peptides overlap at {before} and {after}. "
                "The stage cuts the precursor into pieces and cannot draw one base twice."
            )
    return payload


def fill_removed_peptides(payload: dict, target: Target) -> None:
    """Add the excised middle piece a record does not annotate, from the table.

    NC_000009.12 annotates relaxin's B and A chains as `mat_peptide` but not the
    C-peptide between them, so the mature-peptide page read the whole 108
    residues from the end of B to the start of A as one cut site: "One dibasic
    cut releases two chains", for a precursor cut twice, exactly as insulin's
    is. The table already names that stretch, so this draws it from there.

    Only a removed region that lies between two annotated peptides and overlaps
    none of them is filled, which is the shape of a connecting peptide and not
    of a signal peptide. Insulin's record annotates its C-peptide and this does
    nothing there. The live backend has the same gap to close.
    """
    protein = payload.get("protein")
    if not protein or len(payload["peptides"]) < 2:
        return
    cds = cds_positions(payload)
    codon = {position: i // 3 for i, position in enumerate(cds)}

    def residues(peptide: dict) -> tuple[int, int]:
        """The peptide's 1-based inclusive residue span in the precursor."""
        found = [codon[p] for s in peptide["segments"] for p in (s["start"], s["end"]) if p in codon]
        return min(found) + 1, max(found) + 1

    spans = [residues(p) for p in payload["peptides"]]
    for region in target.regions:
        if region.kept:
            continue
        if not any(end < region.start for _, end in spans):
            continue
        if not any(region.end < start for start, _ in spans):
            continue
        if any(start <= region.end and region.start <= end for start, end in spans):
            continue

        segments = segments_over(payload, cds[(region.start - 1) * 3 : region.end * 3])
        translation = protein["translation"][region.start - 1 : region.end]
        translated = str(Seq(extract(payload, segments)).translate())
        if translated != translation:
            raise ValueError(
                f"{target.slug}: {region.label} {region.start}-{region.end} translates to "
                f"{translated!r}, the protein says {translation!r}"
            )
        payload["peptides"].append(
            {"product": region.label, "segments": segments, "translation": translation}
        )
        spans.append((region.start, region.end))

    payload["peptides"].sort(key=residues)


def fill_chain(payload: dict, target: Target) -> None:
    """Add the chains a record names no peptide for, where the table names them.

    NG_021471 annotates erythropoietin's 27-residue signal peptide and nothing
    after it. Walked as it came, the protein page called the other 166 residues
    'proprotein' and the walk stopped there, where lysozyme's and leptin's go on
    to the chain they become. UniProt names that chain (P01588, Chain 28-193)
    and the table pins it, so this draws it from there.

    A cut a record misses need not be a signal peptide. NG_007462 annotates
    neither of TNF's forms: it is a type II membrane protein, so there is no
    leader to cleave, and ADAM17 sheds the soluble part off the part that stays
    in the membrane (P01375, `Site 76-77 Cleavage; by ADAM17`). So the leader
    here is whatever comes off the front — the record's signal peptide where it
    has one, and nothing at all where the table's chains tile the precursor from
    its first residue.

    Only a record with no peptide at all is filled, and only where the table
    says the precursor is cut and its kept regions tile what the leader leaves,
    exactly. A record that names any chain keeps its own. An uncut protein is
    never cut here, which is what stops p53's and dystrophin's tables of
    *domains* being read as chains — every target with a signal peptide is
    already `cleaved`, so asking for it narrows this to the cut a record missed.
    A table whose regions do not tile stops the bake rather than cutting the
    protein into them. The live backend has the same gap to close.
    """
    protein = payload.get("protein")
    signal = payload.get("signal_peptide")
    if not protein or payload["peptides"] or not target.cleaved:
        return
    leader = len(signal["translation"]) if signal else 0
    length = len(protein["translation"])
    chains = sorted(
        (r for r in target.regions if r.kept and r.start > leader),
        key=lambda r: r.start,
    )
    if not chains:
        return

    cursor = leader + 1
    for region in chains:
        if region.start != cursor:
            raise ValueError(
                f"{target.slug}: {region.label} starts at {region.start}, but the table's chains "
                f"have to tile the precursor from {cursor}"
            )
        cursor = region.end + 1
    if cursor != length + 1:
        raise ValueError(
            f"{target.slug}: the table's chains end at {cursor - 1}, "
            f"the precursor at {length}"
        )

    cds = cds_positions(payload)
    for region in chains:
        segments = segments_over(payload, cds[(region.start - 1) * 3 : region.end * 3])
        translation = protein["translation"][region.start - 1 : region.end]
        translated = str(Seq(extract(payload, segments)).translate())
        if translated != translation:
            raise ValueError(
                f"{target.slug}: {region.label} {region.start}-{region.end} translates to "
                f"{translated!r}, the protein says {translation!r}"
            )
        payload["peptides"].append(
            {"product": region.label, "segments": segments, "translation": translation}
        )


def fill_proprotein(payload: dict, target: Target) -> None:
    """Name what the signal peptide leaves, where the record does not.

    The protein page draws a precursor as its leader and the proprotein behind
    it, and names that block from the record's `proprotein` feature. Insulin's
    RefSeqGene carries one; the chromosome slices oxytocin and relaxin come
    from do not, so the table supplies the name and this supplies the rest:
    every residue after the signal peptide, translated and checked against the
    protein. A record that does annotate a proprotein keeps its own.
    """
    protein = payload.get("protein")
    signal = payload.get("signal_peptide")
    if target.proprotein is None or payload.get("proprotein") or not protein or not signal:
        return
    leader = len(signal["translation"])
    translation = protein["translation"][leader:]
    cds = cds_positions(payload)
    segments = segments_over(payload, cds[leader * 3 : len(protein["translation"]) * 3])
    translated = str(Seq(extract(payload, segments)).translate())
    if translated != translation:
        raise ValueError(
            f"{target.slug}: the proprotein after a {leader}-residue signal peptide "
            f"translates to {translated!r}, the protein says {translation!r}"
        )
    payload["proprotein"] = {
        "product": target.proprotein,
        "segments": segments,
        "translation": translation,
    }


def cds_positions(payload: dict) -> list[int]:
    """Every coding position, in transcript order."""
    reverse = payload["location"].get("strand") == -1
    cds: list[int] = []
    for segment in sorted(payload["protein"]["segments"], key=lambda s: s["start"], reverse=reverse):
        positions = range(segment["start"], segment["end"] + 1)
        cds.extend(reversed(positions) if reverse else positions)
    return cds


def segments_over(payload: dict, picked: list[int]) -> list[dict]:
    """Positions in transcript order, as the feature segments that cover them."""
    step = -1 if payload["location"].get("strand") == -1 else 1
    segments: list[dict] = []
    run = [picked[0]]
    for position in picked[1:]:
        if position == run[-1] + step:
            run.append(position)
        else:
            segments.append({"start": min(run), "end": max(run)})
            run = [position]
    segments.append({"start": min(run), "end": max(run)})
    return segments


# ------------------------------------------------------------- coordinates


def segments_of(payload: dict) -> list[list[dict]]:
    """Every segment list in the payload, so a remap can visit them all once."""
    found = [payload["exons"]]
    for key in ("transcript", "protein", "signal_peptide", "proprotein"):
        if payload.get(key):
            found.append(payload[key]["segments"])
    for peptide in payload["peptides"]:
        found.append(peptide["segments"])
    return found


def transcript_spans(payload: dict) -> list[tuple[int, int]]:
    """The exonic stretches, merged, ascending. What compression must keep."""
    raw = []
    if payload.get("transcript"):
        raw = [(s["start"], s["end"]) for s in payload["transcript"]["segments"]]
    if not raw:
        raw = [(e["start"], e["end"]) for e in payload["exons"]]
    if not raw:
        return []
    raw.sort()
    merged = [list(raw[0])]
    for start, end in raw[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(a, b) for a, b in merged]


def extract(payload: dict, spans: list[dict]) -> str:
    """Pull `spans` out of the payload's own sequence, in transcript order.

    The same arithmetic `AnatomyModel.baseAt` does in Dart: `sequence` is
    already reverse-complemented for a minus-strand gene, so the index runs
    from the far end there.
    """
    start, end = payload["location"]["start"], payload["location"]["end"]
    reverse = payload["location"].get("strand") == -1
    sequence = payload["sequence"]
    ordered = sorted(spans, key=lambda s: s["start"], reverse=reverse)
    out = []
    for span in ordered:
        positions = range(span["start"], span["end"] + 1)
        if reverse:
            positions = reversed(positions)
        out.append("".join(sequence[end - p if reverse else p - start] for p in positions))
    return "".join(out)


def check_frame(payload: dict, where: str) -> str:
    """Translate the CDS out of the payload and demand its own `/translation` back.

    This is the check that proves a remap did not quietly shift a frame: the
    coordinates, the sequence and the `/translation` all have to agree, and
    compression rebuilds two of the three from scratch. Returns the protein,
    which is what everything downstream is scored and drawn from.
    """
    protein = payload.get("protein")
    if not protein:
        raise ValueError(f"{where}: no CDS in the record")
    stated = protein["translation"]
    translated = str(Seq(extract(payload, protein["segments"])).translate()).rstrip("*")
    if translated != stated:
        for i, (a, b) in enumerate(zip(translated, stated), 1):
            if a != b:
                raise ValueError(f"{where}: translated CDS diverges at residue {i}: {a} vs {b}")
        raise ValueError(
            f"{where}: translated CDS is {len(translated)} residues, /translation is {len(stated)}"
        )
    return stated


def check_uniprot(stated: str, target: Target) -> None:
    """Compare the record's protein with UniProt, and name every difference.

    A genomic reference and a UniProt entry can hold different common alleles —
    NG_012232 and P11532 disagree at three of dystrophin's 3,685 residues — so
    this cannot be a plain equality. It is not a tolerance either: the
    differences have to be exactly the ones the table declares, or the bake
    stops. Anything else is NCBI having revised the record under us.
    """
    expected = canonical_sequence(target.uniprot)
    if len(expected) != target.aa:
        raise ValueError(f"{target.slug}: UniProt is {len(expected)} aa, table says {target.aa}")
    if len(stated) != len(expected):
        raise ValueError(
            f"{target.slug}: record protein is {len(stated)} aa, {target.uniprot} is {len(expected)}"
        )
    found = tuple(
        (i, a, b) for i, (a, b) in enumerate(zip(stated, expected), 1) if a != b
    )
    if found != target.uniprot_variants:
        raise ValueError(
            f"{target.slug}: differs from {target.uniprot} at {found or 'nothing'}; "
            f"the table declares {target.uniprot_variants or 'no differences'}"
        )


# ------------------------------------------------------------- compression


def clip_to_transcript(payload: dict) -> dict:
    """Narrow the record to the transcript it carries, where the gene is wider.

    A `gene` feature spans every variant RefSeq annotates at that locus, and on
    a chromosome slice that can be most of the page: RLN2's chosen mRNA spans
    4,853 bases inside a 39,463-base gene, so 88% of its gene page came out as
    "outside the transcript" — true, and nothing anyone came to see.

    On a RefSeqGene this does nothing, because the gene feature and the
    transcript already agree there — insulin's are the same 1,431 bases. So the
    ten records end up consistent with each other rather than one of them being
    a picture of a locus while the rest are pictures of a gene.
    """
    spans = transcript_spans(payload)
    if not spans:
        return payload
    low, high = spans[0][0], spans[-1][1]
    if low <= payload["location"]["start"] and high >= payload["location"]["end"]:
        return payload
    payload["sequence"] = extract(payload, [{"start": low, "end": high}])
    payload["location"] = {**payload["location"], "start": low, "end": high}
    return payload


def intron_scale(payload: dict, budget: int) -> float:
    """The largest scale that fits the gene page's budget, or 1.0 if it already does."""
    spans = transcript_spans(payload)
    kept = sum(b - a + 1 for a, b in spans)
    gaps = gap_lengths(payload, spans)
    if kept + sum(gaps) <= budget:
        return 1.0
    if kept >= budget:
        raise ValueError(
            f"Exons alone are {kept:,} bp, over the {budget:,} bp budget. "
            "Compressing introns cannot help; this gene needs a different page."
        )
    low, high = 0.0, 1.0
    for _ in range(60):
        middle = (low + high) / 2
        total = kept + sum(max(MIN_INTRON_BP, round(g * middle)) for g in gaps)
        low, high = (middle, high) if total <= budget else (low, middle)
    return low


def gap_lengths(payload: dict, spans: list[tuple[int, int]]) -> list[int]:
    start, end = payload["location"]["start"], payload["location"]["end"]
    edges = [(start - 1, start - 1), *spans, (end + 1, end + 1)]
    return [
        b[0] - a[1] - 1
        for a, b in zip(edges, edges[1:])
        if b[0] - a[1] - 1 > 0
    ]


def compress(payload: dict, scale: float) -> dict:
    """Shorten every intron by `scale`, keeping exons and the frame intact.

    A shortened intron keeps its own first and last bases rather than a middle
    slice, so the donor `GT` and acceptor `AG` that make it an intron at all are
    still the bases on screen. Nothing exonic is touched: the mRNA page, the
    protein page and the constraint track are the real, whole molecule.
    """
    start, end = payload["location"]["start"], payload["location"]["end"]
    spans = transcript_spans(payload)

    # Ascending runs of genomic positions to keep, exons and clipped introns
    # alike, as (source_start, source_end) pairs.
    keep: list[tuple[int, int]] = []
    cursor = start
    for span_start, span_end in [*spans, (end + 1, end + 1)]:
        gap = span_start - cursor
        if gap > 0:
            shortened = min(gap, max(MIN_INTRON_BP, round(gap * scale)))
            head = shortened // 2
            tail = shortened - head
            if head:
                keep.append((cursor, cursor + head - 1))
            if tail:
                keep.append((span_start - tail, span_start - 1))
        if span_start <= end:
            keep.append((span_start, min(span_end, end)))
        cursor = span_end + 1

    mapping: dict[int, int] = {}
    kept_bases: list[str] = []
    reverse = payload["location"].get("strand") == -1
    sequence = payload["sequence"]
    position = start
    for source_start, source_end in keep:
        for source in range(source_start, source_end + 1):
            mapping[source] = position
            kept_bases.append(sequence[end - source if reverse else source - start])
            position += 1
    new_end = position - 1

    def remap(coordinate: int, what: str) -> int:
        if coordinate not in mapping:
            raise ValueError(f"compression dropped {what} at base {coordinate:,}")
        return mapping[coordinate]

    out = json.loads(json.dumps(payload))
    for group in segments_of(out):
        for segment in group:
            segment["start"] = remap(segment["start"], "a feature boundary")
            segment["end"] = remap(segment["end"], "a feature boundary")
    out["location"] = {"start": start, "end": new_end, "strand": payload["location"].get("strand")}
    # `kept_bases` is ascending by genomic position, and each base was already
    # read out of a sequence the backend reverse-complemented. So a minus-strand
    # gene needs its bases put back into transcript order and nothing more —
    # complementing them a second time would hand the protein page the wrong
    # strand, which is exactly what `check_frame` catches below.
    out["sequence"] = "".join(reversed(kept_bases)) if reverse else "".join(kept_bases)
    out["intron_scale"] = round(scale, 6)
    out["real_span_bp"] = end - start + 1
    # Each intron's length before it was shortened, in transcript order, which
    # is how the gene page numbers them. The scale alone cannot give them back:
    # an intron shortened to the MIN_INTRON_BP floor was shortened by less.
    introns = [b[0] - a[1] - 1 for a, b in zip(spans, spans[1:]) if b[0] - a[1] - 1 > 0]
    out["real_intron_bp"] = list(reversed(introns)) if reverse else introns
    return out


# ------------------------------------------------------------------- build


def build(target: Target, budget: int) -> dict:
    record = fetch(target)
    payload = extract_gene(record, target.gene)

    if not payload["exons"] and payload.get("transcript"):
        # A chromosome slice annotates mRNA and CDS but no `exon` features,
        # and the gene page's caption counts exons. A transcript's segments are
        # its exons by definition, so this names them rather than inventing
        # them. The backend should do the same before it serves these genes.
        payload["exons"] = [
            {"number": n, "start": s["start"], "end": s["end"]}
            for n, s in enumerate(
                sorted(payload["transcript"]["segments"],
                       key=lambda s: s["start"],
                       reverse=payload["location"].get("strand") == -1),
                1,
            )
        ]

    payload = tidy(payload, target)
    check_uniprot(check_frame(payload, f"{target.slug} (as fetched)"), target)

    payload = clip_to_transcript(payload)
    check_frame(payload, f"{target.slug} (clipped)")

    span = payload["location"]["end"] - payload["location"]["start"] + 1
    if len(payload["sequence"]) != span:
        raise ValueError(f"{target.slug}: sequence is {len(payload['sequence'])} bp, span is {span}")

    scale = intron_scale(payload, budget)
    if scale < 1.0:
        payload = compress(payload, scale)
        # The whole point of the gate: compression rebuilt both the sequence
        # and every coordinate, so the frame is re-proved from the result.
        check_frame(payload, f"{target.slug} (compressed)")
        new_span = payload["location"]["end"] - payload["location"]["start"] + 1
        if len(payload["sequence"]) != new_span:
            raise ValueError(f"{target.slug}: compressed sequence and span disagree")

    return payload


def write(payload: dict, destination: Path) -> int:
    encoded = (json.dumps(payload, indent=4) + "\n").encode()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as temp:
        temp.write(encoded)
        temporary = temp.name
    try:
        os.chmod(temporary, 0o644)
        os.replace(temporary, destination)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return len(encoded)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", action="append", help="slug; repeatable. Default: all.")
    parser.add_argument("--all", action="store_true", help="Every target. The default anyway.")
    parser.add_argument("--budget", type=int, default=GENE_PAGE_BUDGET_BP)
    args = parser.parse_args()

    if not Entrez.email:
        raise SystemExit("Set NCBI_EMAIL, the way the backend's .env does.")

    chosen = [t for t in TARGETS if not args.target or t.slug in args.target]
    if args.target and len(chosen) != len(set(args.target)):
        raise SystemExit(f"Unknown slug in {args.target}")

    for target in chosen:
        payload = build(target, args.budget)
        size = write(payload, DATA / target.mock_asset)
        span = payload["location"]["end"] - payload["location"]["start"] + 1
        note = ""
        if payload.get("intron_scale"):
            note = f"  introns 1:{round(1 / payload['intron_scale']):,} of {payload['real_span_bp']:,} bp"
        print(
            f"{target.slug:14} {target.gene:5} {span:>7,} bp  "
            f"{len(payload['exons']):>3} exons  {len(payload['protein']['translation']):>5} aa  "
            f"{size:>8,} B{note}",
            flush=True,
        )
        time.sleep(0.4)  # Entrez, unauthenticated.
