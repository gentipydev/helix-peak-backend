# Baking `assets/mock/gene_*.json`

Twenty records, one per protein, each of them exactly what the backend's
`GET /gene/{id}/{gene}` answers — because the parsing *is* the backend's:
`build_gene_record.py` imports `app.genbank_parser.extract_gene` from this
repo rather than reimplementing it. A record that parses through a different
parser is a record that can drift from the contract it stands in for.

They are the walk's `record` track: written under `pipeline/data/assets/mock/`
and uploaded with `upload_tracks.py --kind record`.

```sh
NCBI_EMAIL=you@example.com \
  .venv/bin/python pipeline/mock/build_gene_record.py --target leptin
```

`--all` rebuilds every record from whatever NCBI serves today, so a new
protein is baked with `--target`, which leaves the others byte for byte.

The backend's venv, because that is where biopython is. Fetched flat files are
cached under `$HELIXPEEK_GB_CACHE` (a temp directory by default), which is not
politeness alone: dystrophin's RefSeqGene is 2.2 Mb of GenBank.

## The gates

Nothing is written until all of these pass, and each one has caught something.

**The frame.** The CDS is translated straight out of the payload's own sequence
and coordinates, and has to come back equal to the record's `/translation`. It
is the only check that proves a coordinate remap did not shift a reading frame,
and it caught exactly that: minus-strand compression was complementing the
bases a second time, which the relaxin record reported as residue 1 being Y
instead of M.

**UniProt.** The record's protein is compared with `rest.uniprot.org`, and the
differences have to be *exactly* the ones the row declares. Not a tolerance:
NG_012232 and P11532 disagree at three of dystrophin's 3,685 residues, those
three are named in `targets.py`, and a fourth would stop the bake.

**Overlap.** Mature peptides have to be disjoint. The stage that draws them
cuts the precursor into pieces and cannot draw one base in two of them. A record
whose peptides overlap as alternatives — glucagon's pancreatic and intestinal
products, APP's secretase fragments — declares `mature_peptides=False`, and its
walk stops at the precursor until there is a page for alternatives.

## Six things the live service does not do yet

All six are marked at their call sites. Until the backend learns them, these
are the places the fixture and the service would differ.

1. **Exons from the transcript.** A chromosome slice — which is how OXT and
   RLN2 are reached, neither having a RefSeqGene — annotates `mRNA` and `CDS`
   but no `exon` features, and the gene page's caption counts exons. A
   transcript's segments are its exons by definition, so they are named rather
   than invented.
2. **`/product` qualifiers.** NCBI writes a UniProt feature's whole note into
   `/product` on some records: oxytocin's arrives as `Oxytocin.
   /evidence=ECO:0000269|PubMed:13591312. /id=PRO_0000020495`. The app draws
   `/product` as a label on the molecule, so the trailing qualifiers are cut.
3. **Unannotated connecting peptides.** RLN2's record annotates relaxin's B
   and A chains but not the C-peptide between them, so the mature-peptide page
   counted all 108 residues from B to A as one cut site. A removed region in
   `targets.py` that sits between two annotated peptides and overlaps none of
   them is added as a peptide, and its translation is checked against the
   protein's.
4. **Unannotated proproteins.** The protein page names what the signal peptide
   leaves from the record's `proprotein` feature, which only insulin's
   RefSeqGene has. `proprotein` in `targets.py` names it for oxytocin
   (oxytocin-neurophysin 1) and relaxin (prorelaxin), and the residues after
   the signal peptide are added under that name, translated and checked.
   Vasopressin's RefSeqGene names none either, and APP's precursor is named
   this way because its walk stops there.
5. **A lone chain behind a signal peptide.** NG_021471 annotates
   erythropoietin's signal peptide and no chain, which left its protein page
   calling the chain 'proprotein' and its walk without a mature page. Where a
   record has a signal peptide and no peptide at all, the table's kept regions
   after the signal peptide are added as chains. They have to tile the rest of
   the precursor exactly, and each is translated and checked.
6. **The transcript, where two share one CDS.** AMY1A's slice holds two mRNAs
   around identical CDS coordinates, differing only in an untranslated first
   exon. The CDS no longer identifies the transcript, so `transcript_id` in
   `Source` names the MANE Select one.

## Shortened introns

Only where the gene cannot otherwise be drawn: dystrophin, whose transcript
spans 2.1 Mb, APP (290 kb) and CFTR (189 kb), against a 24,000-base budget that
comes from `AnatomyLayout.fit` sizing the gene page to one screen with a
two-point floor on a cell. Every other record keeps every base of the span it is
clipped to.

What is shortened is introns and only introns. Exons keep every base, so the
transcript page, the protein page and the constraint track are the real, whole
molecule; a shortened intron keeps its own first and last bases rather than a
middle slice, so the donor `GT` and acceptor `AG` that make it an intron are
still the bases on screen. The scale goes into the payload as `intron_scale`,
the real gene as `real_span_bp`, and each intron's real length, in transcript
order, as `real_intron_bp`: an intron shortened to the 60-base floor was
shortened by less than the scale, so the scale cannot give its length back. The
gene page's caption says "Introns are 99.3% of the gene, drawn shortened; exons
are to scale", its badge reads the real span, and a tapped intron reads its real
length and how much of it is drawn.

## What a correct bake produces

| slug | gene | record | span | exons | protein |
|---|---|---|---:|---:|---:|
| insulin | INS | NG_007114 | 1,431 | 3 | 110 |
| hemoglobin | HBB | NG_059281 | 1,608 | 3 | 147 |
| myoglobin | MB | NG_007075 | 10,566 | 3 | 154 |
| p53 | TP53 | NG_017013 | 19,149 | 11 | 393 |
| lysozyme | LYZ | NG_008195 | 5,880 | 4 | 148 |
| relaxin | RLN2 | NC_000009.12 slice | 4,853 | 2 | 185 |
| oxytocin | OXT | NC_000020.11 slice | 898 | 3 | 125 |
| somatotropin | GH1 | NG_011676 | 1,637 | 5 | 217 |
| ubiquitin | UBB | NG_023320 | 1,650 | 2 | 229 |
| dystrophin | DMD | NG_012232 | 24,000 | 79 | 3,685 |
| vasopressin | AVP | NG_008663 | 2,169 | 3 | 164 |
| glucagon | GCG | NC_000002.12 slice | 9,366 | 6 | 180 |
| app | APP | NG_007376 | 24,000 | 18 | 770 |
| cftr | CFTR | NG_016465 | 24,000 | 27 | 1,480 |
| erythropoietin | EPO | NG_021471 | 3,233 | 5 | 193 |
| leptin | LEP | NG_007450 | 16,352 | 3 | 167 |
| tnf | TNF | NG_007462 | 2,772 | 4 | 233 |
| sod1 | SOD1 | NG_008689 | 9,310 | 5 | 154 |
| amylase | AMY1A | NC_000001.11 slice | 8,795 | 11 | 511 |
| prion | PRNP | NG_009087 | 15,133 | 2 | 253 |

Three records are shortened: dystrophin's real span is 2,092,329, APP's 290,221
and CFTR's 188,703. Every span
is the transcript's own: `clip_to_transcript` narrows a gene feature wider than
its transcript, which is what brings relaxin's 39 kb locus down to 4,853 bases.

A record with more than one transcript needs the row to name which: RLN2's
chromosome slice carries six CDS features, four of them predicted `XP_` models,
and `extract_gene` takes whichever comes first in the file. `protein_id` in
`Source` picks one; the mRNA is found rather than named, being the one whose
segments contain every segment of the chosen CDS — unless two do, as AMY1A's
do, and then `transcript_id` names it.
