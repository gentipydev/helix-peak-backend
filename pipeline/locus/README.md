# Baking `locus`: where each gene lies, in the body and by band

The app's zoom runs from a body down to a gene, through an organ, a tissue, a
cell, its nucleus and a chromosome. Its last level before the gene is the
chromosome the gene is on, with the gene marked at its band, and its levels on
the way down go through the organ and cells the gene is read in. That needs,
per protein, the chromosome's bands and where the gene falls among them, and
where in the body its RNA is found, and nothing in the records says any of it:
a RefSeqGene record is a stretch of sequence with its own coordinates, not a
place on a chromosome or in a body.

This is a separate track kind, `locus`. Nothing is added to `targets.py`,
`curated/catalog.json` or the `protein` table (there is no cytoband column),
and no other track is touched.

```sh
.venv/Scripts/python pipeline/locus/bake_locus.py --all      # writes pipeline/data/assets/locus/
.venv/Scripts/python pipeline/locus/check_locus.py           # every check below; exits 1 on any failure
.venv/Scripts/python pipeline/locus/verify_locus.py          # each place held to HGNC's (network)
.venv/Scripts/python pipeline/upload_tracks.py --kind locus --dry-run
```

Upload only when asked: `DATABASE_URL`, `SUPABASE_URL` and `SUPABASE_SERVICE_KEY`
set, `migrations/0007_locus.sql` applied (after `0006`), then the same command
without `--dry-run`. Objects go to the `tracks` bucket as
`locus/<slug>.<sha12>.json`.

Which proteins get the track: all twenty (`LOCUS_TARGETS` in `bake_locus.py`).

## Where it comes from

1. **The gene's span**, from MANE, found the way the protein index finds it.
   The protein's UniProt entry (`pipeline/uniprot.py`) is joined to MANE's
   current summary (`pipeline/mane.py`, which the index loader reads through
   too) by `rows_for_entry` in `app/protein_index.py`, and the buildable row
   for the table's gene is the one. Its GRCh38 sequence, start, end and strand
   are MANE Select's.
2. **The chromosome's name**, from UCSC's own `chromAlias` table: the RefSeq
   accession MANE gives (`NC_000011.10`) is UCSC's `chr11`.
3. **The bands**, from UCSC's `cytoBand` table for hg38: every band of that
   chromosome, with its Giemsa stain.
4. **Where in the body**, from the Human Protein Atlas (`proteinatlas.org`,
   the gene's own summary by its Ensembl gene, which MANE gives): how
   specific the gene's RNA is to a tissue and to a single cell type, in the
   Atlas's categories, and the tissues and cell types it names, highest
   first. The Atlas's release is read from its releases page: its version,
   release date and the Ensembl version it is built on. The Atlas is licensed
   CC BY 4.0: the app has to name it wherever it shows these.

It rides in this track rather than one of its own because the app's
`TrackKind` enum is closed over the walk's tests: a new kind would mean
editing a walk-test helper, and the zoom is the only reader of either.

Both UCSC tables are read through the Genome Browser's REST API
(`api.genome.ucsc.edu/getData/track`), which serves the same tables as the
download server; the download server does not answer from every network.

### Moving the MANE download

`scripts/load_protein_index.py` downloaded MANE's summary inline. That block
moved unchanged into `pipeline/mane.py`'s `current_summary`, which both the
loader and this bake call. The join and the parse were already shared.

The move was proven on a fresh build: one download of UniProt 2026_03, MANE
v1.5 and LRG_RefSeqGene, built by the loader before and after the move, gave
the same 20,543 rows, 292,164 terms, duplicates and gene names; with the
catalog laid over them (twenty listed, no failures, no notes) the rows and all
292,202 terms equal the stored `protein_index` and `protein_index_term`, read
without writing.

## What it holds

One JSON object per protein. Every position is 1-based and inclusive on
GRCh38, as MANE and GenBank count; UCSC's 0-based band starts are made 1-based.

| field | what |
|---|---|
| `chromosome`, `sequence`, `length` | `11`, `NC_000011.10`, and its last base |
| `span` | the gene's start, end and strand, MANE Select's |
| `transcript` | that transcript, RefSeq and Ensembl |
| `bands` | every band of the chromosome in order from the end of the short arm: name, start, end and stain (`gneg`, `gpos25` to `gpos100`, `acen`, `gvar`, `stalk`) |
| `band` | the band or bands the span lies in, and where they begin and end |
| `locus` | the place as cytogenetics writes it: `11p15.5`; `Xp21.2-p21.1` across two bands |
| `expression` | the Ensembl gene, and for `tissue` and for `cell_type` the Atlas's `specificity` (`Tissue enriched`, `Group enriched`, `Tissue enhanced`, `Low tissue specificity`, `Not detected`, and the cell-type counterparts), its `distribution`, and the `specific` ones it names with their levels (nTPM for a tissue, nCPM for a cell type), highest first; none where nothing is specific |
| `sources` | `cytoband` and `chrom_alias`: the table, its genome, when UCSC last updated it (its version: UCSC tables carry no other), how many rows were read and their digest; `mane`: the release; `uniprot`: the release, its date and the entry's version; `hpa`: the Atlas version, its release date, its Ensembl version, its licence and the gene's URL |

A band is a stain pattern seen down a microscope at low resolution, millions of
bases long. The track says where the gene lies; it does not say the gene can be
seen there, and the app must not say so either.

## What the check holds

`check_locus.py`, and the uploader before it writes, hold each track to:

- its own protein: the table's slug, gene and accession;
- GRCh38 and hg38, on chromosome 1 to 22, X or Y, named by its RefSeq accession;
- the whole chromosome: bands from base 1 to its last without a gap or an
  overlap, the short arm's before the long arm's, each with a cytoBand stain;
- the gene inside it: the band names are exactly those the span overlaps, the
  band's boundaries are theirs, and the locus is written from them;
- the record, where the record says where it is: a record cut from a
  chromosome (oxytocin, relaxin, glucagon, amylase) is cut from this one,
  around this span, and a table row naming its transcript names this one;
- the Atlas as it reads itself: its categories only, levels above zero and
  highest first, and tissues or cell types named exactly when the gene is
  specific to some;
- a version for every source, and the Atlas's licence.

`verify_locus.py` holds each baked place to the one HGNC curates for the gene
symbol, read at the time: HGNC may write a coarser place, never a different one.

## What it found

All twenty resolve to a band, and all twenty agree with HGNC, written the same
way. Dystrophin's 2.1 million bases run across two bands, Xp21.2 and Xp21.1.

| chromosome | genes |
|---|---|
| 11 | INS (p15.5) and HBB (p15.4): neighbours at the tip of the short arm |
| 17 | TP53 (p13.1), UBB (p11.2), GH1 (q23.3) |
| 20 | OXT, AVP and PRNP, all in p13 |
| 21 | APP (q21.3), SOD1 (q22.11) |
| 7 | EPO (q22.1), CFTR (q31.2), LEP (q32.1) |
| X | DMD (p21.2-p21.1) |
| 1, 2, 6, 9, 12, 22 | AMY1A (1p21.1), GCG (2q24.2), TNF (6p21.33), RLN2 (9p24.1), LYZ (12q15), MB (22q12.3) |

From the Atlas, version 25.1: sixteen are specific to a tissue in some way,
from insulin and glucagon in the pancreas and hemoglobin in the bone marrow to
growth hormone in the pituitary; TP53, UBB, DMD and APP are read in every
tissue it measured. Hemoglobin's cell type is `Erythrocytes`: mature red
cells, which have no nucleus, so the zoom lands in the marrow precursor that
still has one and says so. The Atlas's top cell type is not always a cell of
its top tissue (myoglobin's is thymic myoid cells, the muscle-like cells of the
thymus, where its tissue is skeletal muscle), so the app states each as the
Atlas's reading rather than zooming one into the other.
