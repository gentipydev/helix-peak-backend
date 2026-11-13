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
.venv/bin/python -I pipeline/locus/snapshot.py               # the pinned inputs, once per machine
.venv/Scripts/python pipeline/locus/bake_locus.py --all      # writes pipeline/data/assets/locus/
.venv/Scripts/python pipeline/locus/check_locus.py           # every check below; exits 1 on any failure
.venv/bin/python -I pipeline/locus/verify_locus.py           # each place held to NCBI Gene's (network)
.venv/bin/python -I pipeline/locus/survey.py --live --out OUT  # every protein the search offers
.venv/Scripts/python pipeline/upload_tracks.py --kind locus --dry-run
```

Upload only when asked: `DATABASE_URL`, `SUPABASE_URL` and `SUPABASE_SERVICE_KEY`
set, `migrations/0007_locus.sql` applied (after `0006`), then the same command
without `--dry-run`. Objects go to the `tracks` bucket as
`locus/<slug>.<sha12>.json`.

Which proteins get the track:
- every curated one (`LOCUS_TARGETS` in `bake_locus.py`), baked and uploaded by hand as above;
- every protein built on demand, made by the resolver as it builds the protein (`worker.make_locus`, schema 3), or by `scripts/queue_locus.py` for one built before that. Both read the same snapshot (`pipeline/resolver/README.md`).

Seven of the proteins the search offers get no zoom: their locus is `refused`, with the reader's sentence ([Schema 3](#schema-3-every-protein-the-search-offers)).

## Where it comes from

Every source but UniProt is read from one pinned snapshot, never fetched at
bake time: `snapshot.py` downloads the Atlas 25.1 (from its own archive site,
`v25.proteinatlas.org`), UCSC's tables, MANE v1.5's summary, HGNC's monthly
archive of 2026-10-06 and NCBI's GRCh38.p14 scaffold placements, and writes
`pipeline/data/locus_inputs/locus-inputs.<sha12>.jsonl.gz` (git ignores it).
`inputs.py` pins the sha256 of its text, `INPUTS_SHA256`, and refuses any
other file; the same sources on the same date (`INPUTS_CREATED`) give the
same text, so a machine without it builds it again. The Atlas's site serves
its newest release, and a bake that read it gene by gene would have moved a
curated zoom to 26.0 unasked. Built 2026-11-13: `bb475d41dcf9`, 20,162 Atlas
genes, 42 MB of text, 9.4 MB on disk. It reproduces every curated protein's
stored `expression` and `path`.

1. **The gene's span**, from MANE, found the way the protein index finds it.
   For a curated protein, its UniProt entry (`pipeline/uniprot.py`, the one
   source still fetched) is joined to the snapshot's MANE v1.5 summary by
   `rows_for_entry` in `app/protein_index.py`, and the buildable row for the
   table's gene is the one. A protein built on demand needs no fetch: its
   index row (`protein_index`, which the loader made by that same join) is
   that row. Its GRCh38 sequence, start, end and strand are MANE Select's. A
   gene MANE puts on a patch or alternate locus scaffold is carried to its
   chromosome ([Schema 3](#schema-3-every-protein-the-search-offers)).
2. **The chromosome's name**, from UCSC's own `chromAlias` table: the RefSeq
   accession MANE gives (`NC_000011.10`) is UCSC's `chr11`.
3. **The bands**, from UCSC's `cytoBand` table for hg38: every band of that
   chromosome, with its Giemsa stain.
4. **Where in the body**, from the Human Protein Atlas 25.1 as the snapshot
   holds it: the gene's entry in the Atlas's bulk download, found by its
   Ensembl gene ([the match](#schema-3-every-protein-the-search-offers)).
   From it come:
   - how specific the gene's RNA is to a tissue and to a single cell type,
     in the Atlas's categories;
   - the tissues and cell types it names, highest first;
   - since schema 3, its consensus RNA level in every tissue and its level
     in every single cell type.

   The release line (version, release date, the Ensembl version it is built
   on) is the releases page's, read as the snapshot was built. The Atlas is
   licensed CC BY 4.0: the app has to name it wherever it shows these.

It rides in this track rather than one of its own because the app's
`TrackKind` enum is closed over the walk's tests: a new kind would mean
editing a walk-test helper, and the zoom is the only reader of either.

UCSC's tables are read through the Genome Browser's REST API
(`api.genome.ucsc.edu/getData/track`) as the snapshot is built; it serves the
same tables as the download server, which does not answer from every network.

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
| `sources` | `cytoband` and `chrom_alias`: the table, its genome, when UCSC last updated it (its version: UCSC tables carry no other), how many rows were read and their digest; `mane`: the release; `uniprot`: the release, its date and the entry's version; `hpa`: the Atlas version, its release date, its Ensembl version, its licence and the gene's URL (on `v25.proteinatlas.org`); since schema 3, `locus_inputs`: the snapshot's sha256 and date |
| schema 3: `expression.atlas_match`, `atlas_gene` | how the Atlas entry was found (`mane`, `hgnc` or `uniprot`), and the Atlas's own name for the gene (AFG2A's is `SPATA5`) |
| schema 3: `expression.tissue_levels` | the consensus nTPM in every tissue, its sample tissues counted toward their group (`CONSENSUS_GROUPS`), the highest of each: what rule v3 reads |
| schema 3: `expression.cell_levels` | where rule v3 chose the cell by level: the nCPM of each single cell type living in the chosen tissue |
| schema 3: `placed_from` | only for a gene MANE puts on a patch or alternate locus scaffold: the scaffold's RefSeq accession, its UCSC name and NCBI's name for it, the gene's span on it, the placement (GRCh38.p14, orientation, the aligned stretch on the scaffold and on the chromosome) and the method (`proportional`) |

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
  its cell: rule v2 for a schema 2 track, v3 for schema 3;
- a version for every source, and the Atlas's licence;
- schema 3:
  - how the Atlas entry was matched and named, and its v25 URL;
  - a level for every one of the 37 tissues, in order;
  - `cell_levels` exactly where the rule goes by level, naming every cell type that lives in the tissue it goes to;
  - a path to a tissue: a gene the Atlas reads in no tissue is refused, never stored;
  - `sources.locus_inputs`;
  - `placed_from`, its place worked out again from the block's own numbers;
  - for a record cut from a scaffold (an `NW_` or `NT_` source), `placed_from` naming that scaffold, and its span inside the slice. The rule once skipped any source but an `NC_` one.

`verify_locus.py` holds each baked place to the `map_location` NCBI Gene
gives the gene, by the GeneID MANE gives its transcript and never by its
symbol, read in bulk. It classes every difference:
- `equal`;
- `coarser`;
- `range`: NCBI ranges over every transcript;
- `multi`;
- `par`;
- `placed`: carried from a scaffold;
- `none`;
- `wrong`, of which there must be none.

HGNC, written by hand and sometimes stale, is counted as a second opinion and never held to.

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

The pairs' cell type names have no class from the Atlas, so each takes the
class of its counterpart among the single cell types (`Spermatocytes` as
`Late primary spermatocytes`); a qualifier (`Mitotic cells (Stomach)`) or a
cluster's number (`Neutrophils_2`) is the bare name's. The table once held
twenty of those names, and run over the
Atlas's whole download (19,300 genes with an entry) the rule stopped on 860
of them, at 29 names it lacked. It now has all 96 the pairs name in 25.1,
and the rule takes a path for every one of the 19,300.

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

## Schema 3: every protein the search offers

Schema 3 (2026-11-13) takes the zoom to every protein the search offers to
build (7,311 that day), not the curated alone. It adds keys and changes
none, and a schema 2 reader ignores them.

The curated stay schema 2, their stored objects untouched; rule v3 gives
every one of them the path it has (a test holds v3 to v2 on all 54).
`bake_locus.py --all` now writes schema 3 for them, with the same
`expression` and `path` field for field. A locus re-bake was never
byte-for-byte: `retrieved` was the day it ran. Uploading the curated again
would change their format, so it waits for the user's word.

**A built protein's payload is a pure function of three inputs:**
- its index row;
- the index's release;
- the snapshot.

`retrieved` is the snapshot's date, `sources.uniprot` and `sources.mane` are
the index release's, and `sources.locus_inputs` names the snapshot. So the
census's bytes and the resolver's are the same bytes, which was checked live
on B2M, PFN1 and GFAP.

**Patch scaffolds.** MANE puts 16 genes' transcripts on GRCh38 patch or
alternate locus scaffolds; 12 of them are offered.
- Such a gene is carried to the chromosome NCBI places its scaffold on. The
  source is GRCh38.p14's `all_alt_scaffold_placement.txt`: proportionally
  within the aligned stretch, its strand turned where the scaffold lies
  reversed (KIR2DS4).
- `sequence` is the chromosome's `NC_`, and `placed_from` says where it was
  read and how it was carried.
- UCSC's `altLocations` and `fixLocations`, kept in the snapshot, are a cross-check only.
- NCBI Gene's place holds all 12: APOBEC3A_B, C4B_2, CASP8AP2, FCGBP, GBA3,
  IFNL4, KIR2DS4, POLR2A, RYBP, SAMD1, SHANK3, SLC37A4.

**The Atlas entry** is found by the first of these that gives one:
- `mane`: the Ensembl gene MANE puts the transcript on;
- `hgnc`: the Ensembl gene HGNC's archive gives the row's HGNC id (a gene MANE annotates on a scaffold carries the scaffold's gene id);
- `uniprot`: the one Atlas entry naming the UniProt accession, where exactly one does.

Of the 7,304 placed: 7,290 / 11 / 3. The Atlas names 27 of them otherwise
(`atlas_gene`).

**Rule v3** (the user's decision of 2026-11-13; `path_of_v3`) keeps rule
v2's first two ways as they are, and reads its last two by level:
- A gene specific to no tissue and to no cell type with a home, with tissue
  cell type pairs only, goes to the pair in the tissue its RNA is highest
  in. Where two tissues are level, the earlier pair in the Atlas's order
  wins. v2 took the first pair, and the Atlas lists them by tissue,
  alphabetically: 141 proteins move.
- A gene with none of these goes to the tissue its RNA is highest in
  (`tissue_from: tissue_level`); of two level, the later name. It goes into
  the cell type living there it is highest in (`cell_from:
  single_cell_level`), the later name where two are level, else the
  tissue's own cells. v2 gave these 255 no tissue at all.

A sample tissue counts toward its group (`CONSENSUS_GROUPS` in `tissues.py`:
the brain's eight regions, the intestine's four, the lymphoid tissues' five),
at the group's highest. 396 paths differ from v2, each to another organ,
and no curated one does. Six genes are settled by a tie-break (HEMK2,
L3MBTL2, METTL16, NSD1, RBM12, TAF2).

**Refusals.** Seven offered proteins get no zoom: their `locus` row is
`refused` with the sentence the reader is shown (`LocusRefused`).

| Sentence | Genes |
|---|---|
| "The Human Protein Atlas has no reading for this gene." (`NO_ATLAS_ENTRY`) | GSTT1, KIR2DS2, HLA-DRB3, LILRA3: the Atlas has no entry |
| "The Human Protein Atlas detects this gene's RNA in no tissue." (`NOT_DETECTED`) | BECN2, MRGPRX4, TAS2R16 |

## The census

`survey.py` bakes and checks every protein the search offers, as the
resolver would. It writes only under `--out`.

```sh
.venv/bin/python -I pipeline/locus/survey.py --live --out OUT \
    --gene-info Homo_sapiens.gene_info.gz --gff MANE.GRCh38.v1.5.refseq_genomic.gff.gz
cd ../helix-peek && ZOOM_CATALOGUE=OUT flutter test test/zoom_catalogue_check.dart
```

- `--live` reads the index read-only; `--save-dump` / `--dump` keep a copy to run from again.
- `--gff` writes, for the app's check, a record per protein from MANE's exons with made-up bases. The check stages and paints every stop of every protein at 390×844 and 375×667.
- Re-run it after any index reload (`scripts/load_protein_index.py` reads UniProt's current release).

On 2026-11-13 (UniProt 2026_04, MANE v1.5, snapshot `bb475d41dcf9`):

| | |
|---|---|
| Offered | 7,311 |
| Tracks / refused | 7,304 / 7 |
| Failures / check problems | 0 / 0 |
| NCBI: equal, placed, multi, range, PAR, none (THPO) | 7,247 / 12 / 22 / 13 / 9 / 1 |
| **NCBI wrong** | **0** |
| Paths changed from v2 / curated among them | 396 / 0 |
| No cell (the tissue's own cells drawn) | 858 |
| Organs / organ and cell pairs | 37 / 270 |
| App check: proteins / failures / plates cut short | 7,299 / 0 / 0 |

Five proteins are offered but have no record: their exons alone are over the
gene page's 24,000 bp, so the resolver refuses to build them (GRIN2B, LRRC7,
RPS6KA5, SMAD2, ZBTB20).

## Atlas 26.0 is its own project

The Atlas's site serves 26.0 since this snapshot was built. Moving to it is
one new snapshot, pinned, plus a review. Measured on 2026-11-13:
- 3,480 paths change, 8 of them curated zooms;
- the rule stops on 185 proteins, at names its tables lack;
- four new cell names (cDCs, Basophils, Bone marrow MSCs, Enteroendocrine
  cells) need a class and a home.

Never read the live Atlas gene by gene.
