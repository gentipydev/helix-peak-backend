# ClinVar observed evidence

Every catalog gene ships a snapshot: INS since September 21, 2026, and the
other nineteen since September 22. A row added without one says "not yet
included" in its About sheet. Nothing calls NCBI at app runtime: the snapshots
are stored tracks, like the AVI and ESM ones -- baked here into `pipeline/data/`,
sent with `upload_tracks.py --kind clinvar`, and fetched by the app from storage.
The bake reads the record and the AVI map from `pipeline/data/`; run
`fetch_tracks.py` first when they are not there.

```sh
NCBI_EMAIL=you@example.com \
  .venv/bin/python -u pipeline/clinvar/bake_clinvar.py --target <slug>
# Every row with a snapshot, as the other bakers' --all:
NCBI_EMAIL=you@example.com \
  .venv/bin/python -u pipeline/clinvar/bake_clinvar.py --all
# Replay exactly the downloaded source without another network request:
.venv/bin/python pipeline/clinvar/bake_clinvar.py --all --replay
.venv/bin/python -m pytest pipeline/clinvar
python3 pipeline/check_assets.py
```

The gene's AVI track has to be baked first: records are placed through its
coordinate map. `NCBI_API_KEY` is read if set; the bake is sequential either
way. The nineteen baked on 2026-09-22 took about seventeen minutes, dystrophin's
12,041 records about seven of them, and left about 1 GB of XML in the cache.

The baker uses Biopython Entrez: a paginated `GENE[gene]` search followed by
`efetch(rettype="vcv", is_variationid="true")` in batches of 100. The current VCV schema separates
germline, somatic clinical impact and oncogenicity classifications. This overlay
quotes **germline** classifications only. It retains the aggregate classification
and each condition's RCV classification, review status, evaluation date, VCV
accession/version, update date and collection methods. It neither requests nor
infers individual patient identities or patient counts.

Conditions are named, not described. A top-level `traits` map gives each
condition name the identifiers ClinVar attaches to it: the RCV's MedGen
concept, and from the record's own ConditionList that concept's OMIM entry (a
MIM number, or a `PS` phenotypic series), MONDO class and symbol. The symbol is
the one OMIM cross-references for the condition's own entry, else ClinVar's
preferred one; records list symbols in no fixed order, so ties go by rank
(preferred, fewest other sources, shortest), never by position. ClinVar's
placeholders and names it has not mapped (INS: "INS-related disorder") carry
none, and one name with two concepts aborts the bake. Definitions are not
taken: MedGen gives type 2 diabetes the WFS1 GeneReviews summary.

Raw XML and a retrieval manifest live in the ignored `cache/<GENE>/` directory.
The shipped JSON includes the UTC retrieval time, query, total search count,
exclusion counts, and each XML batch's SHA-256. Search pagination, duplicate IDs,
failed fetches and incomplete responses fail before replacing the app asset.
Biopython retries a failed request but not a failed read of its body, so each
batch is fetched, read and parsed as one attempt and retried whole, with a
backoff that doubles from ten seconds, until it holds exactly the records asked
for; CFTR's bake needed one such retry. `--replay` verifies the raw XML
checksums and the full set of returned IDs.

Mapping uses the existing verified GRCh38 coordinate runs, independently of AVI
scores. A direct `SimpleAllele/Location` must identify exactly one substitution
on the correct chromosome, and the genomic reference must match the actual
drawn base after strand conversion. Gene spans, complex records, other variant
types, omitted intron middles, other assemblies and reference mismatches cannot
be silently attached to a nearby base. A record inside the gene's span that the
map leaves out sits in the undrawn middle of an intron drawn shortened (DMD,
APP, CFTR), and is counted as `intron_not_drawn` rather than as outside the
gene. Coding position and precursor residue
are derived through the selected transcript's spliced CDS, never by subtracting
genomic endpoints or by trusting an arbitrary protein isoform's HGVS string.
The app rechecks the current gene sequence and residue mapping before display.

What the twenty hold. These are **record counts**, not patient counts, and the
classes are the app's colour groups (R9.1), not ClinVar's wording, which the
app always quotes. Absence means no mapped SNV in this snapshot; it does not
imply that the position has no ClinVar record or a benign effect.

| gene | records | mapped | P/LP | Conflicting | VUS | B/LB | Other |
|---|---:|---:|---:|---:|---:|---:|---:|
| INS | 242 | 166 | 37 | 24 | 54 | 42 | 9 |
| HBB | 1,963 | 1,437 | 187 | 89 | 232 | 772 | 157 |
| MB | 43 | 24 | 1 | 0 | 20 | 3 | 0 |
| TP53 | 4,023 | 2,708 | 421 | 515 | 711 | 1,060 | 1 |
| LYZ | 118 | 93 | 4 | 6 | 55 | 28 | 0 |
| RLN2 | 243 | 42 | 0 | 0 | 35 | 7 | 0 |
| OXT | 56 | 11 | 0 | 0 | 10 | 1 | 0 |
| GH1 | 227 | 183 | 25 | 15 | 90 | 53 | 0 |
| UBB | 66 | 14 | 0 | 0 | 5 | 9 | 0 |
| DMD | 12,041 | 8,051 | 1,160 | 706 | 2,950 | 3,235 | 0 |
| AVP | 196 | 141 | 34 | 5 | 74 | 28 | 0 |
| GCG | 27 | 4 | 0 | 0 | 2 | 2 | 0 |
| APP | 720 | 554 | 25 | 30 | 249 | 246 | 4 |
| CFTR | 6,466 | 4,930 | 783 | 363 | 2,268 | 1,440 | 76 |
| EPO | 91 | 50 | 1 | 1 | 29 | 19 | 0 |
| LEP | 150 | 115 | 7 | 4 | 60 | 44 | 0 |
| TNF | 36 | 16 | 1 | 0 | 9 | 6 | 0 |
| SOD1 | 391 | 290 | 127 | 19 | 73 | 71 | 0 |
| AMY1A | 49 | 12 | 0 | 0 | 10 | 2 | 0 |
| PRNP | 234 | 180 | 26 | 5 | 84 | 64 | 1 |

Of the 27,382 records searched, 19,021 map to a single drawn base. The rest
are other variant types (4,669), records with no GRCh38 location on the gene's
chromosome at all (2,248; almost all deletions, duplications and copy-number
changes whose breakpoints ClinVar does not place),
records inside an intron drawn shortened (772), records with no germline
classification (312; somatic-only records), records outside the drawn gene
(290), complex records or ones filed under another gene (66), and dystrophin's
four reference mismatches, at the assembly differences R2.5 already names.
Small genes lose most of their records to the first two: GCG keeps 4 of 27,
AMY1A 12 of 49, RLN2 42 of 243.

In the app every record is read once, against its own allele: the exact AVI
alternative at its base and, for a missense record, the ESM score of its own
amino acid (`VariantEvidence`, `lib/.../domain/entities/variant_evidence.dart`).
A residue or base sheet lists only the records at that position, each as one
row carrying the model its bars do not show; an opened row is the one place
both numbers and the one-line reading meet. The overview draws every record as
one mark at its real position — height its own AVI, standing on the residue's
ESM constraint, coloured by class — with records off the protein on a drawing
of the gene, and lists them grouped by region. Class colour is for grouping
only (rules in `docs/protein-pipeline-rules.md` §9, R9.1); the verbatim
classification is always what is shown. ClinVar puts terms off the Mendelian
axis after a semicolon — CFTR's `Pathogenic; drug response`, HBB's `Pathogenic;
other` — and those are read term by term like the slash-joined ones, which
moved 17 CFTR and 28 HBB records out of Other. The list builds its rows as they
come on screen: dystrophin's overview opens with 139 of its 8,051 built. No
regression, enrichment statistic,
diagnostic verdict or count of agreement with ClinVar is generated: ClinVar
ascertainment is selective and classifications can incorporate computational
evidence.

Sources:

- [NCBI programmatic access and VCV examples](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/)
- [NCBI review status definitions](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
- [What ClinVar collects](https://www.ncbi.nlm.nih.gov/clinvar/intro/)
- [Evidence and computational predictions in submissions](https://www.ncbi.nlm.nih.gov/clinvar/docs/faq_submitters/)
