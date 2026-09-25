"""AlphaGenome Variant Impact, per base, for every gene the app draws.

One asset per gene, `assets/impact/<slug>_avi.json`: the three substitutions of
every base the gene page can reach, scored by the Atlas and filed under the
coordinate the app already uses. The app never calls anything at tap time.

The whole problem here is coordinates. AVI is indexed by GRCh38 `chr:pos`; the
app has no chromosome anywhere. Every position it holds — in `assets/mock/`, in
`AnatomyStage.positions`, in a tracer — is a 1-based offset inside the `NG_` or
`NC_` record the gene was lifted from, and three genes (DMD, APP, CFTR) have
had their introns compressed to fit the gene page's budget on top of that.

The map back is recoverable without a per-gene table, because of two things the
rest of the pipeline already guarantees. The record is clipped to its transcript
(R1.4), so every drawn base is exonic or intronic and nothing flanks. And
`compress()` keeps a shortened intron's own first `d // 2` and last `d - d // 2`
bases rather than a middle slice, so every drawn base is still a real base. So:

  * pair the record's exons with GENCODE's MANE Select exons in transcript
    order — by position and strand, never by the record's `number` field, which
    is null four times in TP53 and skips 2, 4 and 9;
  * anchor each exon on the end it shares with an intron, and let the first
    exon's 5' end and the last exon's 3' end float, because RefSeqGene and
    GENCODE disagree there (CFTR's first exon is 185 bases to GENCODE's 124);
  * run a compressed intron's head forward from the exon before it and its tail
    backward from the exon after it.

None of that is trusted. Every score the Atlas returns carries the reference
base it was called against, and `check_sequence` demands that base equal the one
the app draws at that cell, complemented where the gene is on the chromosome's
minus strand. A wrong offset, a wrong strand or a missed complement cannot
survive it: INS agrees at 1,431 of 1,431 bases and at 402 of 1,431 with the
orientation flipped by hand, and the gate is what caught RLN2 being read off the
wrong strand — a minus-strand record of a minus-strand gene counts up with the
chromosome while reading the complement of it, which are two separate facts.

Offline, like everything else under `pipeline/`. Run by hand; the output is uploaded
to storage with `upload_tracks.py`.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import statistics
import subprocess
import sys
import time

BACKEND = pathlib.Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.paths import DATA, SKILLS  # noqa: E402
from pipeline.targets import TARGETS, Target  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
GENCODE_DIR = HERE / "gencode"
CHECKPOINT_DIR = HERE / "checkpoints"

# What the header claims, and what `check_assets.py` and `GeneImpact.fromJson`
# check it against. A changed scorer or annotation build is a new asset, not a
# quietly different one.
ASSEMBLY = "GRCh38"
ANNOTATION = "GENCODE v46"
SCORER = "AVI_SCORE"
SCORE_UNITS = "phred"
ALT_ORDER = "ACGT minus wildtype"

# The Atlas meters requests per minute, and a single wide interval fans out into
# many of them: a 189 kb pull died with RESOURCE_EXHAUSTED where 20 kb and 100 kb
# went through. Windows are cut to CHUNK_BP, paced, and retried on the way back.
CHUNK_BP = 20_000
# Two drawn runs closer than this are pulled as one window. The scores between
# them are thrown away, but one request for a 1,200 bp window beats two for a
# 500 bp one when the quota counts requests rather than bases.
MERGE_GAP_BP = 1_500
PACE_SECONDS = 1.5
MAX_ATTEMPTS = 7
# How far an exon may be slid to find its own sequence. Covers the UTR-boundary
# disagreements between RefSeqGene and GENCODE — TP53's is three bases — without
# being wide enough to land on a different exon.
RESOLVE_PAD = 64
# A RefSeqGene record and the primary assembly are not the same bases
# everywhere: NG_012232 carries a G where GRCh38 has an A inside dystrophin's
# exon 21, with the other 180 bases of that exon exact. So an exon is placed
# where it fits best rather than where it fits perfectly, and the whole record
# is allowed a few differences — but only a few. A mapping that is simply wrong
# disagrees at about three bases in four, nowhere near these.
EXON_MISMATCH_RATE = 0.02
RECORD_MISMATCH_RATE = 0.005

# Phred is calibrated genome-wide: 20 is the top 1% of all SNVs, 10 the top 10%.
# The app buckets on the same numbers, so they are written down once here and
# once in `gene_impact.dart`, and the test checks the boundaries.
HIGH_PHRED = 20.0
MIDDLE_PHRED = 10.0

BASES = ("A", "C", "G", "T")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

# The skill that owns the Atlas. Its GTF command is the only annotation source
# allowed here: gene models have to match the scores and the website, so Ensembl
# or UCSC would be a different set of exons.
DEFAULT_SKILL = SKILLS / "alphagenome-variant-impact-score"


class BakeError(RuntimeError):
    """A gate said no. Nothing is written."""


# ------------------------------------------------------------------ GENCODE


def gencode(gene: str, skill: pathlib.Path, refresh: bool = False) -> dict:
    """The MANE Select transcript for `gene`, cached on disk.

    Cached because the GTF command ingests a 318 MB feather to answer, and
    because a bake should be reproducible from the repo alone — the same reason
    `targets.py` pins UniProt regions rather than fetching them.
    """
    GENCODE_DIR.mkdir(parents=True, exist_ok=True)
    cached = GENCODE_DIR / f"{gene}.json"
    if cached.exists() and not refresh:
        return json.loads(cached.read_text())

    script = skill / "scripts" / "alphagenome_atlas_avi.py"
    if not script.exists():
        raise BakeError(
            f"no GENCODE cache for {gene} and no skill at {skill}. "
            f"Pass --skill, or restore {cached.relative_to(BACKEND)}."
        )
    result = subprocess.run(
        ["uv", "run", str(script), "gtf", "--gene", gene, "--exons", "--format", "json"],
        cwd=skill, capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        raise BakeError(f"gtf failed for {gene}: {result.stderr.strip()[-400:]}")
    text = result.stdout
    start = text.find("[")
    if start < 0:
        raise BakeError(f"gtf returned no JSON for {gene}")
    entries = json.loads(text[start:])
    mane = [e for e in entries if e.get("is_mane_select")]
    if len(mane) != 1:
        raise BakeError(
            f"{gene}: expected exactly one MANE Select transcript, found {len(mane)}"
        )
    cached.write_text(json.dumps(mane[0], indent=2, sort_keys=True) + "\n")
    return mane[0]


# ------------------------------------------------------- the coordinate map


def record_exons(record: dict) -> list[tuple[int, int]]:
    """The record's exons as (5' end, 3' end) in transcript order.

    Order comes from position and strand. The `number` field cannot be used:
    NG_017013 leaves four of TP53's eleven exons unnumbered and skips 2, 4 and 9.
    """
    reverse = record["location"].get("strand") == -1
    spans = sorted(record["exons"], key=lambda e: e["start"], reverse=reverse)
    return [
        (e["end"], e["start"]) if reverse else (e["start"], e["end"])
        for e in spans
    ]


def mane_exons(mane: dict) -> list[tuple[int, int]]:
    """GENCODE's exons as (5' end, 3' end) in transcript order.

    Its `start` is 0-based where everything else here is 1-based, which is why
    every exon is one base wider than the record's until it is corrected.
    """
    minus = mane["strand"] == "-"
    spans = sorted(mane["exons"], key=lambda e: e["exon_number"])
    out = []
    for e in spans:
        low, high = e["start"] + 1, e["end"]
        out.append((high, low) if minus else (low, high))
    return out


def exon_pairs(record: dict, mane: dict, gene: str):
    """The record's exons beside GENCODE's, in transcript order, with both steps."""
    rec = record_exons(record)
    gen = mane_exons(mane)
    if len(rec) != len(gen):
        raise BakeError(
            f"{gene}: record has {len(rec)} exons, MANE Select has {len(gen)}"
        )
    rec_step = -1 if record["location"].get("strand") == -1 else 1
    gen_step = -1 if mane["strand"] == "-" else 1
    return rec, gen, rec_step, gen_step


def seed_placements(rec, gen, gen_step) -> list[int]:
    """Where each exon's transcript-5' base sits, before the sequence confirms it.

    Anchored on the donor, which an intron shares and both annotations therefore
    agree on, and sized by the record's own width rather than GENCODE's. The two
    disagree at UTR boundaries — TP53's first exon is 174 bases to GENCODE's 114
    and its second is 99 to GENCODE's 102 — and it is the record that says what
    the app draws.
    """
    out: list[int] = []
    last = len(rec) - 1
    for k, ((l5, l3), (g5, g3)) in enumerate(zip(rec, gen)):
        width = abs(l3 - l5) + 1
        out.append(g5 if k == last else g3 - (width - 1) * gen_step)
    return out


def map_from(rec, rec_step, gen_step, placements) -> dict[int, int]:
    """Every drawn local position to its genomic one, given where the exons sit."""
    mapping: dict[int, int] = {}
    for k, (l5, l3) in enumerate(rec):
        width = abs(l3 - l5) + 1
        start = placements[k]
        for j in range(width):
            mapping[l5 + j * rec_step] = start + j * gen_step
        if k == len(rec) - 1:
            break
        donor_local, donor_genomic = l3, start + (width - 1) * gen_step
        acceptor_local, acceptor_genomic = rec[k + 1][0], placements[k + 1]
        drawn = abs(acceptor_local - donor_local) - 1
        if drawn <= 0:
            raise BakeError(f"exons {k + 1} and {k + 2} are not separated")
        # `compress()` keeps a shortened intron's own first half and last half,
        # so the head runs on from the donor and the tail runs back from the
        # acceptor. An intron that was never shortened has the two meet.
        head, tail = drawn // 2, drawn - drawn // 2
        for j in range(head):
            mapping[donor_local + (j + 1) * rec_step] = donor_genomic + (j + 1) * gen_step
        for j in range(tail):
            mapping[acceptor_local - (j + 1) * rec_step] = (
                acceptor_genomic - (j + 1) * gen_step
            )
    return mapping


def resolve_placements(record, rec, rec_step, gen_step, seeds, reference,
                       complemented, gene, pad) -> list[int]:
    """Move each exon to where its own sequence says it is.

    The seed is right wherever the two annotations agree on the donor, which is
    everywhere except a handful of UTR ends. Rather than trust that, each exon is
    slid within `pad` bases until the letters the app draws are the letters the
    chromosome has. A unique answer is required: no match means the locus is
    wrong, and several means the exon is too repetitive to place this way.
    """
    resolved: list[int] = []
    differences = 0
    for k, ((l5, l3), seed) in enumerate(zip(rec, seeds)):
        width = abs(l3 - l5) + 1
        drawn = [drawn_base(record, l5 + j * rec_step) for j in range(width)]
        allowed = int(width * EXON_MISMATCH_RATE)
        scored: list[tuple[int, int]] = []
        for delta in range(-pad, pad + 1):
            start = seed + delta * gen_step
            misses = 0
            for j, letter in enumerate(drawn):
                entry = reference.get(start + j * gen_step)
                if entry is None:
                    misses = width + 1
                    break
                wanted = entry.translate(COMPLEMENT) if complemented else entry
                if letter != wanted:
                    misses += 1
                    if misses > allowed:
                        break
            if misses <= allowed:
                scored.append((misses, start))
        if not scored:
            raise BakeError(
                f"{gene}: exon {k + 1} ({width} bases) matches nothing within "
                f"{pad} bases of {seed}"
            )
        fewest = min(m for m, _ in scored)
        best = [start for m, start in scored if m == fewest]
        if len(best) > 1:
            # A tie on distance as well is a genuinely ambiguous exon; a tie
            # broken by distance is the seed being right and a repeat nearby.
            nearest = min(best, key=lambda h: abs(h - seed))
            if sum(1 for h in best if abs(h - seed) == abs(nearest - seed)) > 1:
                raise BakeError(
                    f"{gene}: exon {k + 1} fits {len(best)} places equally well"
                )
            best = [nearest]
        differences += fewest
        resolved.append(best[0])
    moved = [k + 1 for k, (a, b) in enumerate(zip(seeds, resolved)) if a != b]
    if moved:
        print(f"    exon{'s' if len(moved) > 1 else ''} "
              f"{', '.join(str(m) for m in moved)} moved to match the reference")
    if differences:
        print(f"    {differences} exonic base(s) differ from the assembly")
    return resolved


def runs_of(mapping: dict[int, int]) -> list[dict]:
    """The map as runs, which is how the asset carries it.

    A run is a stretch where local and genomic both step by one, so the whole of
    a 24,000-base gene is a few hundred numbers rather than 24,000 pairs.
    """
    out: list[dict] = []
    for local in sorted(mapping):
        genomic = mapping[local]
        if out:
            last = out[-1]
            expected_local = last["local"] + last["length"]
            expected_genomic = last["genomic"] + last["step"] * last["length"]
            if local == expected_local and genomic == expected_genomic:
                last["length"] += 1
                continue
            if local == expected_local and last["length"] == 1:
                last["step"] = 1 if genomic > last["genomic"] else -1
                if genomic == last["genomic"] + last["step"]:
                    last["length"] += 1
                    continue
        out.append({"local": local, "genomic": genomic, "step": 1, "length": 1})
    for run in out:
        if run["length"] == 1:
            run["step"] = 1
    return out


# ------------------------------------------------------------- the Atlas pull


def windows(positions: set[int], chunk: int, merge_gap: int,
            pad: int = 0) -> list[tuple[int, int]]:
    """Genomic positions coalesced into intervals to ask for."""
    ordered = sorted(positions)
    merged: list[list[int]] = []
    for pos in ordered:
        if merged and pos - merged[-1][1] <= merge_gap:
            merged[-1][1] = pos
        else:
            merged.append([pos, pos])
    if pad:
        merged = [[max(1, low - pad), high + pad] for low, high in merged]
        squashed: list[list[int]] = []
        for low, high in merged:
            if squashed and low <= squashed[-1][1] + 1:
                squashed[-1][1] = max(squashed[-1][1], high)
            else:
                squashed.append([low, high])
        merged = squashed
    out: list[tuple[int, int]] = []
    for low, high in merged:
        while high - low + 1 > chunk:
            out.append((low, low + chunk - 1))
            low += chunk
        out.append((low, high))
    return out


def phred_of(quantile: float) -> float:
    """The Atlas's own calibration: Phred = -10 log10(1 - CDF)."""
    return -10.0 * math.log10(max(1e-7, 1.0 - quantile))


def pull(client, chromosome: str, low: int, high: int) -> dict[int, dict[str, object]]:
    """Every SNV the Atlas holds in a closed 1-based interval.

    Retried rather than paced optimistically: the per-minute quota is shared with
    whatever else is using the key, so the only safe assumption is that any call
    can bounce.
    """
    from alphagenome.data import genome

    interval = genome.Interval(chromosome=chromosome, start=low - 1, end=high)
    delay = PACE_SECONDS
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            scores = client.query_interval(interval, requested_scorers=[SCORER])
            break
        except Exception as error:  # noqa: BLE001 - grpc raises its own family
            text = str(error)
            # The quota is the expected one. The rest are the transport giving
            # up — a dropped route, an unavailable backend — and a bake that is
            # an hour long will meet them; neither says the request was wrong.
            retriable = any(
                marker in text
                for marker in (
                    "RESOURCE_EXHAUSTED",
                    "UNAVAILABLE",
                    "DEADLINE_EXCEEDED",
                    "failed to connect",
                    "No route to host",
                    "Broken pipe",
                )
            )
            if attempt == MAX_ATTEMPTS or not retriable:
                raise BakeError(
                    f"{chromosome}:{low}-{high} failed after {attempt} "
                    f"attempt{'s' if attempt > 1 else ''}: {text[:300]}"
                ) from error
            reason = "quota reached" if "RESOURCE_EXHAUSTED" in text else "transport failed"
            print(f"      {reason}, waiting {delay:.0f}s", flush=True)
            time.sleep(delay)
            delay = min(delay * 2, 120)
    else:  # pragma: no cover - the loop always breaks or raises
        raise BakeError("unreachable")

    if SCORER not in scores:
        raise BakeError(f"no {SCORER} for {chromosome}:{low}-{high}")
    adata = scores[SCORER]
    if "quantiles" not in adata.layers:
        raise BakeError(f"no quantiles for {chromosome}:{low}-{high}")

    import numpy as np

    quantiles = np.ravel(adata.layers["quantiles"])
    out: dict[int, dict[str, object]] = {}
    for index, variant in enumerate(adata.obs["variant"]):
        reference = str(variant.reference_bases)
        alternate = str(variant.alternate_bases)
        if len(reference) != 1 or len(alternate) != 1:
            continue
        entry = out.setdefault(int(variant.position), {"ref": reference})
        entry[alternate] = round(phred_of(float(quantiles[index])), 1)
    return out


def scores_for(client, chromosome: str, needed: set[int], slug: str,
               resume: bool, pad: int = 0) -> dict[int, dict[str, object]]:
    """Every needed position, pulled and checkpointed window by window."""
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
    checkpoint = CHECKPOINT_DIR / f"{slug}.jsonl"
    done: dict[str, dict] = {}
    if resume and checkpoint.exists():
        for line in checkpoint.read_text().splitlines():
            if not line.strip():
                continue
            record = json.loads(line)
            done[record["w"]] = record["d"]
        print(f"    resuming: {len(done)} window(s) already pulled")

    plan = windows(needed, CHUNK_BP, MERGE_GAP_BP, pad)
    bases = sum(high - low + 1 for low, high in plan)
    print(f"    {len(plan)} window(s), {bases:,} bp")

    gathered: dict[int, dict[str, object]] = {}
    with checkpoint.open("a") as log:
        for number, (low, high) in enumerate(plan, 1):
            key = f"{chromosome}:{low}-{high}"
            if key in done:
                payload = done[key]
            else:
                print(
                    f"    [{number}/{len(plan)}] {key} ({high - low + 1:,} bp)",
                    flush=True,
                )
                payload = {
                    str(pos): entry
                    for pos, entry in pull(client, chromosome, low, high).items()
                }
                log.write(json.dumps({"w": key, "d": payload}) + "\n")
                log.flush()
                time.sleep(PACE_SECONDS)
            for pos, entry in payload.items():
                gathered[int(pos)] = entry
    return gathered


# ------------------------------------------------------------------- gates


def drawn_base(record: dict, local: int) -> str:
    """The letter the app draws at a local position.

    R2.1: a minus-strand record's sequence is already reverse-complemented, so
    it is indexed from the far end and never complemented a second time.
    """
    loc = record["location"]
    if loc.get("strand") == -1:
        return record["sequence"][loc["end"] - local]
    return record["sequence"][local - loc["start"]]


def check_sequence(record: dict, mapping: dict[int, int], gathered: dict,
                   complemented: bool, gene: str) -> int:
    """The gate the whole coordinate map rests on."""
    checked = disagreed = 0
    examples: list[str] = []
    for local, genomic in mapping.items():
        entry = gathered.get(genomic)
        if entry is None:
            continue
        reference = str(entry["ref"])
        wanted = reference.translate(COMPLEMENT) if complemented else reference
        checked += 1
        if drawn_base(record, local) != wanted:
            disagreed += 1
            if len(examples) < 5:
                examples.append(
                    f"local {local} -> {genomic}: "
                    f"draws {drawn_base(record, local)}, reference {reference}"
                )
    if checked == 0:
        raise BakeError(f"{gene}: nothing to check the coordinate map against")
    rate = disagreed / checked
    if rate > RECORD_MISMATCH_RATE:
        raise BakeError(
            f"{gene}: the coordinate map is wrong — {disagreed} of {checked} "
            f"bases ({rate:.1%}) disagree with the reference:\n      "
            + "\n      ".join(examples)
        )
    if disagreed:
        # Not an error, and not silent either. These are places the record and
        # the assembly genuinely hold different bases; the Atlas scored the
        # assembly's, so the app cannot claim an exact score for the letter it
        # draws, and the position falls through to an estimate below.
        print(
            f"    sequence gate: {checked - disagreed:,} of {checked:,} bases "
            f"agree; {disagreed} differ from the assembly"
        )
        for example in examples:
            print(f"      {example}")
    else:
        print(f"    sequence gate: {checked:,} of {checked:,} bases agree")
    return disagreed


def intron_spans(record: dict) -> list[tuple[int, int]]:
    """Drawn intron interiors, as local (low, high) pairs."""
    spans = sorted(record["exons"], key=lambda e: e["start"])
    return [
        (a["end"] + 1, b["start"] - 1)
        for a, b in zip(spans, spans[1:])
        if b["start"] - a["end"] > 1
    ]


def check_biology(record: dict, positions: dict[int, list[float]], gene: str) -> dict:
    """Splice boundaries and exons must outscore intron interiors.

    Not a statistical claim — a smoke alarm. A coordinate map that passed the
    sequence gate by accident, or a gene pulled against the wrong locus, shows
    up here as an intron that scores like an exon.
    """
    exonic: list[float] = []
    for exon in record["exons"]:
        for local in range(exon["start"], exon["end"] + 1):
            if local in positions:
                exonic.append(max(positions[local]))

    junction: list[float] = []
    interior: list[float] = []
    for low, high in intron_spans(record):
        edge = min(8, (high - low + 1) // 2)
        for local in range(low, high + 1):
            if local not in positions:
                continue
            value = max(positions[local])
            near = local < low + edge or local > high - edge
            (junction if near else interior).append(value)

    if not exonic or not interior:
        raise BakeError(f"{gene}: too few scored bases to check the biology")

    summary = {
        "exon_median_phred": round(statistics.median(exonic), 2),
        "splice_median_phred": round(statistics.median(junction), 2) if junction else None,
        "intron_median_phred": round(statistics.median(interior), 2),
    }
    if summary["exon_median_phred"] <= summary["intron_median_phred"]:
        raise BakeError(
            f"{gene}: exons score no higher than intron interiors "
            f"({summary['exon_median_phred']} vs {summary['intron_median_phred']})"
        )
    if junction and summary["splice_median_phred"] <= summary["intron_median_phred"]:
        raise BakeError(
            f"{gene}: splice boundaries score no higher than intron interiors "
            f"({summary['splice_median_phred']} vs {summary['intron_median_phred']})"
        )
    return summary


# ------------------------------------------------------------------- baking


def bake(target: Target, client, skill: pathlib.Path, resume: bool,
         refresh_gencode: bool) -> None:
    record = json.loads((DATA / target.mock_asset).read_text())
    mane = gencode(target.gene, skill, refresh=refresh_gencode)
    chromosome = mane["chromosome"]
    print(f"    {chromosome} {mane['strand']} {mane['transcript_id']} "
          f"({mane['num_exons']} exons)")

    rec, gen, rec_step, gen_step = exon_pairs(record, mane, target.gene)
    # Two different facts, and conflating them is a real bug the sequence gate
    # caught: `orientation` is whether the record's coordinates count the same
    # way as the chromosome's, and `complemented` is whether its letters are the
    # other strand's. A minus-strand record of a minus-strand gene counts up
    # with the chromosome and still reads the complement of it.
    orientation = rec_step * gen_step
    complemented = mane["strand"] == "-"
    seeds = seed_placements(rec, gen, gen_step)
    provisional = map_from(rec, rec_step, gen_step, seeds)
    print(f"    {len(provisional):,} drawn bases, orientation {orientation:+d}"
          f"{', complemented' if complemented else ''}")

    gathered = scores_for(
        client, chromosome, set(provisional.values()), target.slug, resume,
        pad=RESOLVE_PAD,
    )
    reference = {pos: str(entry["ref"]) for pos, entry in gathered.items()}

    placements = resolve_placements(
        record, rec, rec_step, gen_step, seeds, reference, complemented,
        target.gene, RESOLVE_PAD,
    )
    mapping = map_from(rec, rec_step, gen_step, placements)
    loc = record["location"]
    drawable = loc["end"] - loc["start"] + 1
    if len(mapping) != drawable:
        raise BakeError(
            f"{target.gene}: mapped {len(mapping)} of {drawable} drawn bases"
        )

    outstanding = {g for g in mapping.values() if g not in gathered}
    if outstanding:
        print(f"    {len(outstanding):,} base(s) fell outside the first pull")
        gathered.update(
            scores_for(client, chromosome, outstanding, target.slug, resume)
        )

    differences = check_sequence(record, mapping, gathered, complemented, target.gene)

    positions: dict[int, list[float]] = {}
    unscored = 0
    for local, genomic in sorted(mapping.items()):
        entry = gathered.get(genomic)
        wildtype = drawn_base(record, local)
        if entry is None:
            unscored += 1
            continue
        # The app's letters read along the transcript; the Atlas filed its
        # alternates against the chromosome's plus strand. A gene on the minus
        # strand has to ask for the complement of the substitution it means.
        values = []
        for alt in (b for b in BASES if b != wildtype):
            key = alt.translate(COMPLEMENT) if complemented else alt
            value = entry.get(key)
            if value is None:
                break
            values.append(value)
        if len(values) != 3:
            unscored += 1
            continue
        positions[local] = values

    if not positions:
        raise BakeError(f"{target.gene}: nothing scored")
    coverage = len(positions) / len(mapping)
    if coverage < 0.98:
        raise BakeError(
            f"{target.gene}: only {coverage:.1%} of drawn bases scored "
            f"({unscored:,} without a full set of three)"
        )

    summary = check_biology(record, positions, target.gene)
    print(
        f"    exons {summary['exon_median_phred']} · "
        f"splice {summary['splice_median_phred']} · "
        f"introns {summary['intron_median_phred']} (median Phred)"
    )

    asset = {
        "gene": target.gene,
        "uniprot": target.uniprot,
        "accession": target.source.accession,
        "assembly": ASSEMBLY,
        "annotation": ANNOTATION,
        "chromosome": chromosome,
        "transcript": mane["transcript_id"],
        "orientation": orientation,
        "complemented": complemented,
        "scorer": SCORER,
        "score_units": SCORE_UNITS,
        "alt_order": ALT_ORDER,
        "high_phred": HIGH_PHRED,
        "middle_phred": MIDDLE_PHRED,
        "start": loc["start"],
        # The drawn letters in increasing record position, which is how every
        # reader indexes them: `sequence[position - start]`. A minus-strand
        # record stores its own sequence in transcript order, read from the far
        # end (R2.1), so copying it here gave RLN2's and GCG's pages the wrong
        # letter at three bases in four.
        "sequence": "".join(
            drawn_base(record, local) for local in range(loc["start"], loc["end"] + 1)
        ),
        "generation": {
            "drawn_bases": len(mapping),
            "scored_bases": len(positions),
            "unscored_bases": unscored,
            "assembly_differences": differences,
            **summary,
        },
        "runs": runs_of(mapping),
        "positions": {str(local): values for local, values in sorted(positions.items())},
    }

    out = DATA / target.impact_asset
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(asset, separators=(",", ":")) + "\n")
    print(f"    wrote {target.impact_asset} ({out.stat().st_size / 1024:.0f} KB)")

    checkpoint = CHECKPOINT_DIR / f"{target.slug}.jsonl"
    if checkpoint.exists():
        checkpoint.unlink()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--all", action="store_true", help="every scored row")
    group.add_argument("--target", help="one row, by slug")
    parser.add_argument(
        "--skill", type=pathlib.Path, default=DEFAULT_SKILL,
        help="the alphagenome-variant-impact-score skill, for GENCODE",
    )
    parser.add_argument(
        "--no-resume", action="store_true",
        help="ignore any checkpoint and pull every window again",
    )
    parser.add_argument(
        "--refresh-gencode", action="store_true",
        help="re-query GENCODE rather than reading the committed cache",
    )
    parser.add_argument(
        "--map-only", action="store_true",
        help="report the coordinate map without scoring; the sequence gate "
             "needs the reference, so this is the seed map, not the final one",
    )
    args = parser.parse_args()

    if args.target:
        chosen = [t for t in TARGETS if t.slug == args.target]
        if not chosen:
            print(f"no target with slug {args.target!r}", file=sys.stderr)
            return 2
    else:
        chosen = [t for t in TARGETS if t.impact_scored]

    if args.map_only:
        failed = 0
        for number, target in enumerate(chosen, 1):
            try:
                record = json.loads((DATA / target.mock_asset).read_text())
                mane = gencode(target.gene, args.skill, refresh=args.refresh_gencode)
                rec, gen, rec_step, gen_step = exon_pairs(record, mane, target.gene)
                orientation = rec_step * gen_step
                mapping = map_from(
                    rec, rec_step, gen_step, seed_placements(rec, gen, gen_step)
                )
                runs = runs_of(mapping)
                print(
                    f"[{number}/{len(chosen)}] {target.gene:6s} {mane['chromosome']:>5s}"
                    f" {mane['strand']} orientation {orientation:+d}"
                    f"  {len(mapping):6,d} bases  {len(runs):4d} runs"
                )
            except BakeError as error:
                failed += 1
                print(f"[{number}/{len(chosen)}] {target.gene:6s} FAILED: {error}",
                      file=sys.stderr)
        print(f"\n{len(chosen) - failed} of {len(chosen)} mapped.")
        return 1 if failed else 0

    key = os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        print(
            "ALPHAGENOME_API_KEY is not set. Put it in ~/.env or the environment.",
            file=sys.stderr,
        )
        return 2

    from alphagenome.atlas import atlas

    client = atlas.create(key)

    failures: list[str] = []
    for number, target in enumerate(chosen, 1):
        print(f"[{number}/{len(chosen)}] {target.gene} ({target.slug})", flush=True)
        try:
            bake(target, client, args.skill, not args.no_resume, args.refresh_gencode)
        except BakeError as error:
            print(f"    FAILED: {error}", file=sys.stderr, flush=True)
            failures.append(target.slug)

    if failures:
        print(f"\n{len(failures)} failed: {', '.join(failures)}", file=sys.stderr)
        return 1
    print(f"\n{len(chosen)} baked.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
