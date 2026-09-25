# Baking the impact tracks

One track per gene, `assets/impact/<slug>_avi.json` under `pipeline/data/`:
AlphaGenome Variant Impact for every base the gene page can reach, filed under
the coordinate the app already uses. Offline like everything else under
`pipeline/` — it runs by hand, its output is uploaded to storage with
`upload_tracks.py --kind impact`, and the app never calls the Atlas at tap time.

```sh
export ALPHAGENOME_API_KEY=...          # or put it in ~/.env
pipeline/impact/venv/bin/python -u pipeline/impact/bake_impact.py --all
```

`--target <slug>` bakes one row and leaves the rest byte for byte.
`--map-only` builds and prints the coordinate map without scoring anything,
which is the cheap way to check a new row before spending any quota.

The environment is `alphagenome` and nothing else:

```sh
uv venv pipeline/impact/venv --python 3.12
uv pip install --python pipeline/impact/venv/bin/python -r pipeline/impact/requirements.txt
```

## The coordinate problem

AVI is indexed by GRCh38 `chr:pos`. The app has no chromosome anywhere: every
position it holds is a 1-based offset inside the `NG_` or `NC_` record the gene
was lifted from, and DMD, APP and CFTR have had their introns compressed to fit
the gene page's 24,000 bp budget on top of that.

The map back is recoverable without a per-gene table, because of two things the
rest of the pipeline already guarantees. The record is clipped to its transcript
(R1.4), so every drawn base is exonic or intronic and nothing flanks. And
`compress()` keeps a shortened intron's own first `d // 2` and last `d - d // 2`
bases rather than a middle slice, so every drawn base is still a real base.

So each exon is paired with GENCODE's MANE Select exon in transcript order —
by position and strand, never by the record's `number` field, which is null four
times in TP53 and skips 2, 4 and 9 — anchored on the end it shares with an
intron, and a compressed intron's head runs forward from the exon before it
while its tail runs back from the exon after it. The first exon's 5' end and the
last exon's 3' end are left to float, because RefSeqGene and GENCODE disagree
there: CFTR's first exon is 185 bases to GENCODE's 124.

**Two facts, not one.** `orientation` is whether the record's coordinates count
the same way as the chromosome's; `complemented` is whether its letters are the
other strand's. A minus-strand record of a minus-strand gene counts *up* with
the chromosome and still reads the *complement* of it. Conflating them read RLN2
and GCG off the wrong strand, and the sequence gate is what caught it.

**The asset's `sequence` is in record order, not the record's order.** It lists
the drawn letters by increasing record position, since the app reads it as
`sequence[position - start]`. A minus-strand record stores its own letters from
the far end (R2.1), so there the two are reverses of each other. Copying the
record's string as stored left every score filed correctly and still gave
RLN2's and GCG's base sheets the wrong wildtype at about three bases in four.
`check_assets.py` and the catalog ClinVar test compare it with the page's
letters, base for base.

## Gates

Each one blocks the write.

| gate | what it proves |
|---|---|
| sequence | every mapped position's AVI `ref` is the base the app draws there, complemented where the gene is on the minus strand |
| exon correspondence | the record and MANE Select have the same exons, and the same widths for the internal ones |
| coverage | at least 98% of drawn bases carry a full set of three substitutions |
| biology | exons and splice boundaries both outscore intron interiors |

The sequence gate is the one the whole map rests on, and it is sharp: INS agrees
at 1,431 of 1,431 bases, and at 402 of 1,431 with the orientation flipped by
hand. A mapping that is simply wrong disagrees at about three bases in four.

**Assembly differences are counted, not forbidden.** A RefSeqGene record and the
primary assembly are not the same bases everywhere, so up to 0.5% of a record
may differ; each one is printed, and the count goes in the asset's header. There
are three across all twenty, all in dystrophin — and they are exactly residues
882, 2366 and 2937, the three `uniprot_variants` that `targets.py` already
declares for it, reached from the other direction. At such a position the Atlas
scored the assembly's base rather than the record's, so no exact score can be
claimed for the letter the app draws; the position is left out and the app
borrows its nearest neighbour's, saying so.

## Costs and quota

The Atlas meters **requests per minute**, and a single wide interval fans out
into many of them: a 189 kb pull dies with `RESOURCE_EXHAUSTED` where 20 kb and
100 kb go through. Windows are cut to 20 kb, paced, and retried with a backoff
that also covers the transport giving up over a run this long.

Only the drawable bases are pulled, coalesced into windows — for dystrophin that
is about 120 kb of the gene's real 2.09 Mb span, in 68 windows. The whole
catalog is roughly 300 kb of transfer; wall clock is set by the quota, not the
bytes, so expect 15–30 minutes for `--all`. Each gene checkpoints to
`pipeline/impact/checkpoints/<slug>.jsonl` as it goes and resumes from there, so a
run interrupted by the quota picks up where it stopped. The checkpoint is
deleted once that gene's asset is written.

GENCODE answers are cached in `pipeline/impact/gencode/<GENE>.json` and committed,
for the same reason `targets.py` pins UniProt regions: a bake should be
reproducible from the repo. `--refresh-gencode` re-queries.

## What a correct result looks like

```
[1/20] INS (insulin)
    chr11 - ENST00000381330.5 (3 exons)
    1,431 drawn bases, orientation -1, complemented
    1 window(s), 1,559 bp
    sequence gate: 1,431 of 1,431 bases agree
    exons 17.8 · splice 21.9 · introns 5.3 (median Phred)
    wrote assets/impact/insulin_avi.json (33 KB)
```

Across the twenty: 4.2 MB, every drawn base scored but dystrophin's three, and
splice boundaries between 11.9 and 26.1 against intron interiors between 2.5 and
9.7. `python3 pipeline/check_assets.py` re-checks the shipped files against the
records, and `test/features/gene_lookup/impact/` checks them from the Dart side.

## Research use only

AVI predicts molecular effect — splicing, expression, chromatin, coding. Nothing
the bake writes or the app draws is a clinical statement, and the panel's words
are kept molecular throughout.
