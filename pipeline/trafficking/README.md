# Baking `trafficking`: where in the cell each protein goes

The app's cell scene lights up the compartments a protein passes through,
and works that route out from the protein's sequence features
(`lib/features/lab/trafficking/domain/` in the app). Four of the features are
in the records already, as regions of the precursor: the signal peptide, the
disulfides, the cuts, and the prion's GPI-anchor signal. The fifth is not. A
region table says what a precursor is cut into, not where any of it sits, so no
record says whether a stretch crosses a membrane. Without that, nineteen of the
twenty routes stop at "unknown": secreted or held in the membrane, after a
signal peptide; left in the cytosol or anchored in the ER, without one.

This is a separate track kind, `trafficking`, and it carries what the records
cannot. Nothing is added to `targets.py`, `curated/catalog.json` or the
`protein` table, and no other track is touched.

```sh
.venv/Scripts/python pipeline/fetch_tracks.py --kind record              # the records, digest-checked
.venv/Scripts/python pipeline/trafficking/bake_trafficking.py --all      # writes pipeline/data/assets/trafficking/
.venv/Scripts/python pipeline/trafficking/check_trafficking.py           # every check below; exits 1 on any failure
.venv/Scripts/python pipeline/upload_tracks.py --kind trafficking --dry-run
```

Upload only when asked: `DATABASE_URL`, `SUPABASE_URL` and `SUPABASE_SERVICE_KEY`
set, `migrations/0004_trafficking.sql` applied (after `0003`), then the same
command without `--dry-run`. Objects go to the `tracks` bucket as
`trafficking/<slug>.<sha12>.json`.

Which proteins get the track: all twenty (`TRAFFICKING_TARGETS` in
`bake_trafficking.py`).

## What it holds

One JSON object per protein, read from its UniProtKB entry:

| field | from UniProt | what |
|---|---|---|
| `transmembrane` | `Transmembrane` features | each span, in precursor numbering, with UniProt's description and evidence codes; `[]` where UniProt annotates none |
| `gpi_anchor` | the `Lipidation` whose description starts `GPI-anchor` | the anchored residue, and the signal after it that is cut off in its place (`site + 1` to the last residue); `null` where there is none |
| `location` | the `SUBCELLULAR LOCATION` comments | each place, with its topology and orientation, and the molecule a comment is about where it names one (an isoform, or a chain cut from the precursor, as TNF's soluble form is) |

Evidence codes are kept because they say how a span is known: CFTR's twelve
helices were observed (`ECO:0000269`), TNF's signal anchor was predicted from
sequence (`ECO:0000255`), and the prion protein's anchor was placed by
similarity (`ECO:0000250`).

Each payload also names its protein (`slug`, `gene`, `uniprot`), its length
(`residues`), the entry's `entry_version` and `sequence_version`, and where it
came from: `release` (`2026_03`), `release_date` (from the service's
`X-UniProt-Release-Date` header) and `retrieved`, the day the bake fetched it.
The uploader copies those into the row's provenance. Two bakes within one
release differ in `retrieved` alone.

## Where the entries come from

`pipeline/uniprot.py`, which the record builder shares: `fetch_entry` reads
`rest.uniprot.org/uniprotkb/<accession>.json` with the release it came from,
and `mock/build_gene_record.py`'s `canonical_sequence` is that same fetch. The
fetch moved out of the record builder unchanged first, and the twenty records
it bakes were re-baked after each step and still match the sha256 on their
stored rows.

## What the check holds

`check_trafficking.py`, and the uploader before it writes, hold each track to:

- its own protein: the table's slug, gene and accession;
- the record's length: `residues` is the table's `aa` and the length of the
  protein the stored record translates, so a span numbered on UniProt's
  sequence lands on the residue the walk draws;
- the chain: every span, and the anchor's site and signal, within
  1..residues, the spans in order and apart;
- the region table: where `targets.py` names a `GPI-anchor signal`, the
  anchor's signal is that region, and where it names none there is no anchor;
- its provenance: a release, a release date and a retrieval date.

## What it found

| | proteins |
|---|---|
| crosses a membrane | CFTR (12 helices), APP (one, 702-722, type I), TNF (one, 36-56, a type II signal anchor) |
| GPI-anchored | the prion protein, at Ser230; 231-253 is cut off |
| neither | the other sixteen |

Hemoglobin's entry names no subcellular location, and ubiquitin's names one
only for the ubiquitin chain it is cut into, not for the precursor.
