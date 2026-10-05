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
| `expression` | the Ensembl gene, and for `tissue` and for `cell_type` the Atlas's `specificity` (`Tissue enriched`, `Group enriched`, `Tissue enhanced`, `Low tissue specificity`, `Not detected`, and the cell-type counterparts), its `distribution`, and the `specific` ones it names with their levels (nTPM for a tissue, nCPM for a cell type), highest first; none where nothing is specific. Since schema 2 also `tissue_cell_type`, the Atlas's tissue cell type pairs as it lists them (`{tissue, cell_type}`, from `Pancreas - Beta cells`); `subcellular`, its `main` and `additional` locations of the protein in a cell (`[]` where it has none); and `secretome`, where it is secreted to, or null |
| `path` | schema 2: the zoom's one way down for this gene, chosen by `path_of` from those readings: `tissue` and how it was chosen (`tissue_from`), the `cell_type` in it, its Atlas `cell_class` and how it was chosen (`cell_from`), and `lands_in`, the cell with a nucleus the zoom enters where that cell type has none |
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
  specific to some, each by a name the Atlas lists; its tissue cell type
  pairs, each in a tissue it lists; its subcellular and secretome words;
- the path the rule takes, recomputed from those readings, with a class for
  its cell;
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
still has one and says so.

The Atlas's top cell type is not always a cell of its top tissue: myoglobin's
is thymic myoid cells, the muscle-like cells of the thymus, where its tissue
is skeletal muscle; tumour necrosis factor's is microglia, a cell of the
brain, where its tissue is bone marrow. Schema 1 left the app to take each
on its own, and seven of the twenty zooms went down into a cell of another
organ. Schema 2 chooses one path per gene (below).

## Schema 2: one path per gene

Schema 2 adds keys and changes none: `expression.tissue_cell_type`,
`expression.subcellular`, `expression.secretome` and `path`. All of them are
read from the same Atlas answer the bake already fetched, which carries 119
fields of which schema 1 kept six. Changing an existing track's format
departs from the contract on purpose, as folding's schema 2 did: the keys are
only added, a schema 1 reader ignores them, and the zoom is the track's one
reader. Every re-baked object has new bytes and a new name.

`path_of` takes the path from the readings alone, the same way for every gene:

1. The organ is the tissue the RNA is highest in, without its sample number.
2. The cell is one that lives there: the first cell type the Atlas finds the
   gene enriched in within that tissue (its tissue cell type pairs), else the
   single cell type it is highest in among those whose home is that tissue,
   else none, and the app draws the tissue's own cells and says so.
3. A gene specific to no tissue goes to the home of the single cell type it
   is highest in (the first with a home), else to its first pair, else
   nowhere in particular.
4. A cell type with no nucleus lands in the one that still has one.

The rule reads two tables of the Atlas's own vocabularies, once for every
gene: `tissues.py`, its 37 consensus tissues and the 19 tissues its pairs
name, each made one of those; and `cell_types.py`, its 154 single cell types
with their class as the Atlas gives it (one of fifteen) and the tissues each
lives in, the classes of the pairs' own cell type names, the two cell types
with no nucleus, and its subcellular and secretome words. A name the tables
lack stops the bake, so a release that renames a tissue cannot lead a zoom
astray. The homes are general biology; a cell type found along nerves or in
connective tissue everywhere, or in a tissue the consensus reading does not
sample (the lacrimal gland, the conjunctiva), has none.

The twenty paths, Protein Atlas 25.1:

| protein | organ (from) | cell (from) | protein in the cell | secreted |
|---|---|---|---|---|
| insulin | pancreas (tissue) | Beta cells (pair) | | to blood |
| glucagon | pancreas (tissue) | Alpha cells (pair) | endoplasmic reticulum; vesicles | to blood |
| cftr | pancreas (tissue) | Pancreatic duct cells (single cell) | | |
| amylase | salivary gland (tissue) | Salivary acinar cells (single cell) | Golgi apparatus, cytosol | to digestive system |
| lysozyme | salivary gland (tissue) | Minor salivary glandular cells (pair) | Golgi apparatus, actin filaments; nucleoplasm | to blood |
| hemoglobin | bone marrow (tissue) | Erythrocytes (single cell), lands in erythroblasts | | |
| tnf | bone marrow (tissue) | none in it: the app draws the marrow's cells | | to blood |
| erythropoietin | liver (tissue) | Hepatocytes (single cell) | | to blood |
| sod1 | liver (tissue) | Hepatocytes (pair) | nucleoplasm; cytosol | |
| leptin | adipose tissue (tissue) | Adipocytes (Subcutaneous) (pair) | vesicles, plasma membrane | to blood |
| myoglobin | skeletal muscle (tissue) | Skeletal myocytes (pair) | | |
| dystrophin | skeletal muscle (its cell's home) | Myonuclei (single cell) | | |
| oxytocin | brain (tissue) | Other brain neurons (single cell) | | to blood |
| vasopressin | brain (tissue) | Other brain neurons (single cell) | | to blood |
| prion | choroid plexus (tissue) | none in it: the app draws the plexus's cells | nuclear membrane, cytosol; vesicles | |
| somatotropin | pituitary gland (tissue) | Somatotropes (pair) | | to blood |
| relaxin | fallopian tube (tissue) | Fallopian tube ciliated cells (single cell) | | to blood |
| ubiquitin | testis (its cell's home) | Late primary spermatocytes (single cell) | cytosol; acrosome, equatorial segment | |
| app | lymphoid tissue (its cell's home) | Lymphatic endothelial cells (single cell) | Golgi apparatus; vesicles | to blood |
| p53 | stomach (its pair) | Mitotic cells (Stomach) (pair) | nucleoplasm; vesicles, cytosol | |

Main locations come before the semicolon, additional ones after it. A blank
means the Atlas gives none.
