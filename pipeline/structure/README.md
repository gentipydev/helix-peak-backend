# Baking `assets/models/*.glb`

Everything here exists to regenerate the curated rows' models: the folds the
structure page renders. It is offline, run by hand, and run once per protein. The `.glb`s it
writes land under `pipeline/data/assets/models/`; storage keeps each one beside
the `.fsceneb` compiled from it, which is what the app fetches.

Two things here are not run by hand: `experimental.py` makes the model of a
protein built on demand from an experimental entry chosen by rule, with its
fold animation ([below](#experimental-entries-for-proteins-built-on-demand)),
and `alphafold.py` from AlphaFold DB where no entry passes
([below](#alphafold-db-models-for-proteins-built-on-demand)); the resolver's
worker runs both. They use this bake's own export, and change none of the
twenty.

It is kept because a stored binary that cannot be regenerated is a
liability. If a model ever needs to change — a different colour split, a
different span, a lighter mesh — this is the only record of how it was made.

```sh
brew install pymol                 # 3.1.0; build-time only, nothing ships
python3 -m venv pipeline/structure/venv
pipeline/structure/venv/bin/pip install -r pipeline/structure/requirements.txt

pipeline/structure/venv/bin/python pipeline/structure/bake.py --all   # from the repo root
```

Every curated row's model is cut from a file in the repository, never fetched
by a bake (`pdb.source_of`): the twenty from their entries' own legacy files,
`structures/<ENTRY>.pdb`, and a row added since from a working copy of its
entry, `structures/working/<slug>.pdb` ([below](#curated-rows-cut-from-a-working-copy)).
Each bake writes its intermediates to `output/<slug>/` — the generated `.pml`,
one `.obj` per chain — and copies the finished model into
`pipeline/data/assets/models/`. `output/` and `venv/` are gitignored.

The app ships no model. It fetches a compiled `.fsceneb` from the `models`
bucket, and the `.glb` it was compiled from is stored beside it, named in the
structure row's provenance. `upload_tracks.py --kind structure` compiles each
scene on the way up (`compiled_scene`), as the worker compiles a model made on
demand: with flutter_scene's own importer, from the app's checkout, which is
what this runs:

```sh
cd ../helix-peek && dart run flutter_scene:import -i model.glb -o model.fsceneb
```

Run on the twenty's stored `.glb`s it gives their stored `.fsceneb`s byte for
byte (2026-10-17: `upload_tracks.py --kind structure --dry-run --against` the
service, twenty of twenty the same). The container is versioned, so the
importer has to be the flutter_scene the stored scenes were compiled by
(0.23.0, `alphafold.FLUTTER_SCENE`).

## Which entry, and why

All experimental, for every curated row: the twenty chosen by hand, the rows added since by policy 1 and cut from their working copies. AlphaFold is not used for any of them:
`AF-P01308` predicts insulin's 110-residue *preproprotein* — the third page's
molecule, not this one — at a mean pLDDT of 52.9 with 51% of residues below 50,
which would draw a coil of low-confidence loops exactly where the reader has
been promised a hormone. (Over insulin's mature span the mean is 48.0, and
`alphafold.py` refuses it at its gate.)

| slug | entry | Å | chains → nodes | notes |
|---|---|---:|---|---|
| insulin | `3I40` | 1.85 | A→`chainA`, B→`chainB` | human; `4INS` is porcine |
| hemoglobin | `2DN1` | 1.25 | β→`chainA` | the gene's chain only; its alpha partners are HBA1/HBA2's |
| myoglobin | `3RGK` | 1.65 | A→`chainA` | |
| p53 | `2OCJ` | 2.05 | A 96-289→`chainA` | most of the rest is disordered |
| lysozyme | `1REX` | 1.50 | A→`chainA` | |
| relaxin | `6RLX` | 1.50 | A→`chainA`, B→`chainB` | insulin's architecture |
| oxytocin | `7RYC` | cryo-EM | L→`chainA` | the ligand of a receptor complex |
| somatotropin | `1HGU` | 2.50 | A→`chainA` | |
| ubiquitin | `1UBQ` | 1.80 | A→`chainA` | one of the precursor's three |
| dystrophin | `1DXX` | 2.60 | A 9-246→`chainA` | the other 3,447 residues are solved only in fragments |
| vasopressin | `7KH0` | cryo-EM | L→`chainA` | the ligand of a receptor complex, as oxytocin's is; drawn as a tube |
| glucagon | `6LMK` | cryo-EM | E→`chainA` | human glucagon on its receptor; no bridges, so no `bonds` node |
| app | `4PWQ` | 1.40 | A→`chainA` | the E1 domain, 28-189; no structure covers the whole protein |
| cftr | `5UAK` | cryo-EM | A→`chainA` | wild type; the entries sharper than it carry E1371Q. Sampling 1 holds 1,139 residues near half a megabyte |
| erythropoietin | `1EER` | 1.90 | A→`chainA` | the glycosylation-site mutant every EPO entry is |
| leptin | `1AX8` | 2.40 | A→`chainA` | carries W121E |
| tnf | `7JRA` | 2.10 | A→`chainA` | one chain of the soluble trimer; apo `1TNF` has Leu where UniProt has Asp219 |
| sod1 | `2C9V` | 1.07 | A→`chainA` | one chain of the dimer; numbered from Ala2 |
| amylase | `1SMD` | 1.60 | A→`chainA` | wild type; its neighbours in the PDB are point mutants |
| prion | `4KML` | 1.50 | A→`chainA` | the folded domain, 117-225; its nanobody (chain B) is not exported |
| phospholamban | `1ZLL` | NMR | A→`chainA` | working copy; one chain of the unphosphorylated pentamer |
| s100b | `5D7F` | 1.30 | A→`chainA` | working copy; one chain of the calcium-bound pair; its RAGE peptide is not exported |
| hcgbeta | `7FIG` | cryo-EM | Y→`chainA` | working copy; hCG's beta chain on its receptor; the alpha chain is CGA's |
| ppbp | `1F9P` | 1.93 | A→`chainA` | working copy; CTAP-III, 44-128, the one chain of its entry |
| pappa | `8HGG` | cryo-EM | C→`chainA` | working copy; one PAPP-A chain of the 2:2 complex with proMBP; sampling 1 |
| cystatinb | `4N6V` | 1.80 | 2→`chainA` | working copy; one of the crystal's ten chains |
| tyrp1 | `9EY8` | 2.20 | AA→`chainA` | working copy; the part inside the melanosome, 25-471; sampling 2 |
| recoverin | `2D8N` | 2.20 | A→`chainA` | working copy; a selenomethionine crystal, drawn as Met; sampling 5 |
| thyroglobulin | `7B75` | cryo-EM | A→`chainA` | working copy; one chain of the dimer; sampling 1, the largest scene (1.17 MB) |
| uromodulin | `9NU1` | cryo-EM | A 329-580→`chainA` | working copy; the ZP domain; 6TQL (292-584), ranked first, fails the app's check |

Every deviation from UniProt in these entries, and why each was chosen over the
alternatives, is in `docs/protein-verification.md`.

The node names are the contract with `structure_view.dart`, which looks them up
to assign a material each. Renaming them here takes their colours off there.

## Curated rows cut from a working copy

A row added to the curated list since the twenty is drafted from its gene by
`scripts/draft_curated.py`: its entry chosen by the on-demand path's rule
(policy 1, [below](#experimental-entries-for-proteins-built-on-demand)), baked
and passed by the app's own check of a fold, as a protein built on demand's
is. What that bake was made from is stored as the row's file:

- `structures/working/<slug>.pdb`, the working copy (`normalize.working_copy`):
  one PDB-format text in the record's numbering, one chain a node (`A` to `C`,
  which the row's `Chain`s name), holding only what is drawn, its first line
  saying so. The twenty's own readers read it as they read an entry.
- `structures/working/<slug>.json`, where it came from: the entry, its model,
  the policy, the entry's own chains and authors' numbering, and the digests of
  the record, model, scene and fold the app's check passed.

`bake.py` and `../folding/bake_folding.py` cut such a row from its working
copy (`pdb.export_source`), with no dashed loop across a gap (R5.4,
`experimental.GAP_DASHES`); its fold is numbered as the record, its ribbon
pruned as the on-demand path prunes one, and its shape named as the working
copy's. `check_bundle.py` then holds the row to what the draft passed: the
`.glb` baked here has the digest the sidecar records, and the app's own check
(`test/fold_bundle_check.dart`) passes it with the curated row, the fold and
the scene compiled from it. `../bake_curated.py` runs all of it in order.

The twenty are not re-cut from working copies: theirs would not be the same
bytes (insulin's and relaxin's chains come out in another order, bridges in
another, the depositors' helices and strands give way to PDBe's, alternate
locations and modified residues are read otherwise, and the fidelity gate
refuses EPO and somatotropin). Their legacy files stay, and they bake to the
digests their rows carry, the proof re-run on 2026-10-17.

## Five things PyMOL does that the script works around

Do not "simplify" these back. Each was found by checking the output.

1. **The OBJ exporter ignores its selection argument.** `save chainA.obj,
   chainA` writes the whole visible scene; exporting three objects that way
   produced three byte-identical files. Each object is exported by hiding
   everything else first.
2. **It cannot export sticks, and exports spheres as single degenerate
   triangles.** Sticks come out as zero vertices and zero faces, so the
   disulfides — the entire point of the page, where there are any — would
   silently be missing from a plain `show cartoon` export. `bake.py` builds
   them as `CA->CB->SG->SG->CB->CA` rods straight from the crystallographic
   coordinates, with the pairs taken from the file's own `SSBOND` records and
   filtered to the chains and span actually exported. They start at CA, not
   CB, because the ribbon runs through CA: a bridge from CB stops about 1.5 Å
   short of it and floats beside the chain.
3. **It exports in camera space, not model space.** Coordinates came out near
   (-5.5, 7.2, -1.5) against a PDB frame near (-20.8, -0.6, -12.9), which would
   have put the generated rods in a different frame from the ribbons. The
   script pins an identity view and puts the rotation origin at the molecule
   centre, so the export is a pure translation. `verify_frame.py` proves it:
   every CA lands within 0.25 Å of the exported ribbon (mean 0.18 Å).
4. **It echoes the script's source before running it.** The export origin is
   printed with a format string, so the format string itself comes past on
   stdout first. The reader matches a line of numbers at the start of a line,
   which the echo never is.
5. **Nothing normalises the model but us.** `bake.py` centres it and scales the
   longest axis to exactly 1.0, so the Dart side carries no scale constant and
   the camera can be framed once. See `_framingMargin` in `structure_view.dart`.

## Shared with the other bakers

Three modules here hold what every baker that starts from these twenty
entries needs, so that each reads them one way:

- `pdb.py`: the entry as the bake reads it (`pdb_path`, `atoms`, `ssbonds`),
  and what it exported (`ca_atoms`, `bridge_atoms`).
- `glb.py`: a stored model's meshes, node by node, with numpy alone.
- `frame.py`: where a point of the entry lands in a stored model,
  `(p - centre) / L`. Exact from the bridges' joints where the model has
  bridges, fitted to the ribbon's CA atoms where it has none.

They need numpy and nothing else, so a baker that reads the stored model runs
on the backend's `.venv` without PyMOL. The three moved out of `bake.py`,
`verify_frame.py` and the AR bake (deleted with the app's lab on 2026-10-09)
unchanged.

`folding/` reads the entry through `pdb.py` too, but takes its frame from
`bake.py` itself: `export` and `normalisation`, re-run and held to the stored
model vertex for vertex. That is exact on all twenty, where `frame.py`'s ribbon
fit is not (see `folding/README.md`).

For an export outside the twenty, `write_pml` and `export` take a `source` (an
assembly built from an entry, or a working copy of one, loaded in place of the
entry itself), an `origin` (one export origin for two exports, so that they
share a frame) and `dashes` (false leaves out PyMOL's dashed loop across a
short gap), and `verify_frame.py` takes `--pdb`, `--out`, `--origin` and
`--chain`. Unset, each is what the twenty are baked and audited with:
re-baking all twenty after they were added gave the same `.glb` files byte
for byte, the same `.pml` scripts and the same `verify_frame.py` report, and
so did re-baking them after `dashes` was added (2026-10-12: every `.glb` the
digest its row names).

### Re-baking on Windows

PyMOL 3.1.0 is conda-forge's `pymol-open-source=3.1.0`, installed with
micromamba beside numpy 2.5.3, scipy 1.18.1 and trimesh 5.1.0 as pinned; put
its `Scripts` directory on `PATH` and run `bake.py` with its Python. It
reproduces the stored models' structure exactly (the same nodes, vertex counts
and triangle counts) but not their bytes: none of the twenty matches storage's
sha256. Every vertex lands within 2.8e-6 model units of the stored one, about
0.00005 Å, and the two tube-drawn peptides wind three and nine zero-area
triangles the other way. That is an x86-64 build's floating point against the
arm64 Mac the models were baked on.

So on Windows a refactor is proved by baking before and after the change on
the same machine and comparing those two byte for byte. Moving the shared
modules out was proved that way: twenty of twenty `.glb`s identical, the same
PyMOL `.obj` intermediates, the same log, and `verify_frame.py` printing the
same distances. Comparing with storage's sha256 needs the Mac.

## What a correct bake produces

Insulin is the reference. Its ribbons reproduce the original hand-made model
exactly; its bridges are heavier, since they now run down to CA:

| node | vertices | faces |
|---|---:|---:|
| `chainA` | 8,256 | 2,752 |
| `chainB` | 11,364 | 3,788 |
| `bonds` | 894 | 1,680 |

- Disulfide SG-SG distances 2.03, 2.04, 2.04 Å (textbook bond lengths); across
  the first ten, every generated bridge falls between 2.00 and 2.06 Å. The
  second ten's run 1.98 to 2.18 Å, the long ends being SOD1's one bridge and
  one of amylase's five, both straight from their entries' coordinates.
- Bounding box before normalising: 24.31 × 18.66 × 20.33 Å.
- After: 1.0000 × 0.7678 × 0.8362, centred at the origin.
- `insulin.glb` is about 594 KB. All twenty together are 20 MB of `.glb` and
  6.0 MB of compiled scene, none of it over the 900 KB a scene is allowed.

Mesh size is `cartoon_sampling`, set per row. Insulin's 8 is PyMOL's own
default and gives about 394 vertices per residue; the longer chains drop to 6,
5 or 4 to hold each compiled scene near half a megabyte. A peptide with no
secondary structure to draw — oxytocin's nine residues — takes
`representation="tube"` instead, which is `cartoon tube` and not
`cartoon_trace_atoms`: the latter threads the tube through the side chains.

## Experimental entries, for proteins built on demand

A protein resolved on demand (`pipeline/resolver/`) has nobody to choose an
entry for it. A rule does, over every PDB entry of the protein, and AlphaFold
DB's model ([below](#alphafold-db-models-for-proteins-built-on-demand)) is for
the protein a whole, successful search finds nothing suitable for. The Mac's
worker runs it for every protein built on demand (`pipeline/resolver/`, its
README has the publication and the backfill); `scripts/dry_fold.py` runs it
for any protein and writes a report and nothing else.

| | |
|---|---|
| `pdb_search.py` | what there is: every polymer entity RCSB's search maps to the UniProt accession (paged, experimental only), and of each, from RCSB's data API, its source organism, its SIFTS alignment, each chain's unplaced residues and its residue-by-residue validation; UniProt's domains, and Pfam's |
| `pdb_policy.py` | which of it may be drawn, each rejection with its reason, and the order of what may: **policy 1** |
| `normalize.py` | the entry's own file (PDBe's SIFTS-enriched mmCIF, read by `mmcif.py`) held to the choice residue by residue, and from it the one model, one conformation and one numbering a bake reads |
| `pdb_words.py` | the fold page's words, from the entry's own numbers |
| `fold_job.py` | all of it for one protein, ending `chosen`, `none` or `incomplete` |
| `experimental.py` | the entry chosen, baked into the model and its fold, and held to the app's own check ([below](#baking-the-entry-chosen)) |
| `fold_report.py` | what a search found, as the dry run's report and as an audit's record |

An outage is none of those three. A search that could not be finished, a page
that did not arrive, a file that cannot be read: each is raised, and asked
again. "Nothing suitable" is said only of a search that ran to its end.

**Policy 1**, frozen on 2026-10-12. Every number is a constant of
`pdb_policy.py` (and `normalize.py`, for what only the file shows); a change
to one is policy 2, and what policy 1 decided is not reused.

| may be drawn | |
|---|---|
| method | X-ray, neutron or electron crystallography, cryo-EM, solution or solid-state NMR. No computed or integrative model |
| organism | the whole entity's source is human. The expression host is not read; a fusion with another organism's protein is out |
| mapping | SIFTS's, to the canonical sequence, in one unbroken stretch. An insertion, a deletion or a chimera splits it, and stretches are never joined. An entity SIFTS has not mapped is found by sequence (95% identity) and laid on the record only where it fits once, without a gap |
| coverage | 70% of a mature product placed. Otherwise 90% of an annotated domain of 30 residues or more (UniProt's Domain, Repeat and Zinc finger; Pfam's where UniProt names none), drawn with at most 10 residues either side |
| fidelity | at most 2% of the residues drawn are not the record's, rounded down, so a peptide under 50 residues may differ nowhere |
| crystal | 3.0 Å or better, R-free 0.30 or under; of the residues drawn, at most 10% RSRZ outliers and 20% with RSCC under 0.8 |
| cryo-EM | 4.0 Å or better; at most half of the residues drawn fitting the map poorly (a Q-score more than 0.10 under what the resolution predicts) |
| NMR | at most 5% Ramachandran outliers |
| any method | at most 30% of the residues drawn in a clash |
| its own file | one model (an NMR entry's representative conformer, never an average), one conformation a residue, no break in the chain, chains drawn together touching, a bridge only where its sulfurs are within 2.5 Å |

Coverage counts residues the entry places, each once: not the construct's
length, and not its copies. Chains are drawn together (at most three, one
node each) only where they are different products of the one precursor, in
one entry, and touch.

Ranked, among what may be drawn, in this order:

1. A mature product before a domain.
2. More of a precursor's products drawn together before fewer.
3. Coverage, in tenths below the best: a few flexible residues never beat a
   better experiment.
4. A mature construct before one that places five or more residues the
   precursor loses.
5. Fewer residues that are not the gene's.
6. A method with a resolution before one without.
7. A region measured residue by residue before one that is not.
8. The finer resolution, then the better fit.
9. Entry id, then chain id, so the same inputs choose the same chains.

Three of these were the user's to settle, on the report's numbers:

- **5 before 6 to 8: wild type first.** R5.4's "prefer wild-type". CFTR draws
  9DW5 (3.8 Å, the gene's own residues) and not 9MXL (2.1 Å, one substitution).
  With the experiment first, 11 of the 52 picks change.
- **6 before 7: a resolution first.** An old crystal with no structure factors
  cannot be measured residue by residue, and an NMR model can. Relaxin keeps
  6RLX (1.5 Å), the entry chosen for it by hand; the other order draws 2MV1.
- **Cryo-EM: half, not 30%.** A hormone on its receptor and a loop in a 3 to
  4 Å map sit under the map's average. At 30% oxytocin has no entry left, and
  mTOR loses its fullest (8ERA, 2,177 residues placed).

**The words** (`pdb_words.py`). The header names the entry and the span, so
the sentence carries the method, the resolution and how much is placed, in
two lines at most:

| drawn | label | sentence |
|---|---|---|
| one chain | the fold | X-ray, 1.65 Å. 149 of these 153 residues placed. |
| one chain of a symmetric assembly of this protein alone | the subunit | X-ray, 1.07 Å. One chain of two; every residue placed. |
| a precursor's chains together | the mature chains | X-ray, 0.92 Å. Two chains, joined by two disulfide bridges. |
| one of a precursor's products | one of its products | X-ray, 1.10 Å. Glucagon, one of six products; 28 of 29 placed. |
| a domain | the P53 DNA-binding domain | X-ray, 1.32 Å. Only this part is drawn; 192 of these 200 residues placed. |

A residue the entry does not place is never called disordered: an experiment
that did not see it has not shown that. A partner, a receptor or an antibody
fragment beside the chain is named to a screen reader only, as not drawn.

**The dry run.** Against live RCSB, PDBe and UniProt; the database is read in
read-only transactions; nothing is queued, stored or uploaded:

```sh
.venv/bin/python scripts/dry_fold.py --out /tmp/folds TTR INSL3     # resolved as the worker would
.venv/bin/python scripts/dry_fold.py --out /tmp/folds --curated --built
```

`report.md` and `report.json` give, for each protein, the winner, the
runners-up and why the winner beat the next, every rejection by its reason,
each candidate's fate against its own file, and the words. On 2026-10-12, over
the twenty (as biological cases: their entries stay the ones chosen by hand),
the nine proteins built by then, and 23 chosen to be hard:

| | |
|---|---|
| 43 | an entry passes |
| 7 | none does: no entry at all (CHRDL2, OCA2); fragments under 70% of the product (APLN, GAL, PTH); chains that fit their maps poorly (PIEZO1); and erythropoietin, whose 1EER differs at 5 residues where 3 are allowed and whose other two entries fail on quality |
| 2 | the resolver refuses the protein before any search (CALCA, KCNQ1) |

`fixtures/pdb/` stands in for the network in the tests, which hold each rule
on recorded answers and on small files written where one thing is wrong at a
time.

### Baking the entry chosen

A candidate whose file passes is baked before it wins, and the next is tried
where the bake refuses it by a documented rule (`fold_job.choose(draw=...)`).
`experimental.py`'s docstring has each step; in short:

- **One working copy.** `normalize.working_copy` writes the chosen chains as
  one PDB-format file in the record's numbering, chains A to C, with the
  entry's own helices and strands (PDBe's DSSP; PyMOL's `dss` where the file
  assigns none, written into the copy so both bakes read the same), its
  `SSBOND`s, `REMARK 465` and `SEQADV`. The twenty's readers read it as they
  read an entry.
- **The model, by the twenty's bake.** `bake.export(source=...)` and
  `bake.assemble`, unchanged; the scene compiled by the app's importer;
  sampling lowered until it fits 1.3 MB. A gap the entry does not place is
  left open: PyMOL's dashed loop across one of up to ten residues is not
  drawn (`GAP_DASHES`, `bake.write_pml(dashes=False)`). It is a curve made up
  where the entry has no atoms, and the app cannot lay it on the chain where
  it curves.
- **The fold, from the same export.** `folding/bake_folding.payload(entry=,
  offset=0)` in the export's own frame, so the track names the very `.glb` it
  was made with. Then each chain is named as the entry names it (`pdb_chain`),
  with the author's numbering offset where one number does for every residue,
  and `null` where it does not (an insertion code); the whole map is in the
  audit.
- **Held to the twenty's rule.** `folding/check_folding.problems_of` with the
  record, the model and the copy, its bridges by R5.6's reach (sulfurs within
  2.5 Å). A ribbon measure that check would turn down (a strand's first
  residue far from where PyMOL flattens its sheet) is dropped first, as one
  `ribbon_at` cannot make is: the app smooths the strand there.
- **Held to the app's rules.** `test/fold_bundle_check.dart`, run from the
  app's checkout with `FOLD_BUNDLE`, reads the bundle as the fold page will:
  the scene equals the `.glb` as the compiler bakes it, the track parses for
  the row's own chains, every ring binds to its stretch of chain (in order to a
  quarter of a residue), every rod to its bridge, the fold is the stored mesh
  float for float once the bridges close, and no vertex jumps between two of
  480 frames. It writes `check.json`: `ok`, `unsuitable` (the entry breaks a
  rule of the animation: the next candidate) or `fault` (no verdict: tried
  again).

By hand, writing nothing anywhere but `--out`:

```sh
pipeline/structure/venv/bin/python pipeline/structure/experimental.py \
    --service https://helix-peak-backend.onrender.com --slug ttr --out /tmp/ttr
.venv/bin/python scripts/dry_fold.py --out /tmp/folds --curated --built --bake
```

Measured on 2026-10-12 over the dry run's 52 proteins, on this Mac: every
one of the 43 the rule chooses for bakes and passes both checks, as its first
candidate. Two keep fewer ribbon measures (HERC2 two, HGF one).

| | min | median | max |
|---|---:|---:|---:|
| scene | 27 kB | 332 kB | 973 kB (mTOR, 2,549 residues, sampling 1) |
| `.glb` | 81 kB | 1.1 MB | 3.5 MB |
| folding track | 2.8 kB | 26 kB | 345 kB |
| bake | 0.7 s | 1.0 s | 8.4 s |
| the app's check | 1.8 s | 2.1 s | 3.4 s |
| binding, in the check | 27 ms | 74 ms | 124 ms |
| 480 frames, in the check | 19 ms | 114 ms | 498 ms |

The check's times are a Mac's CPU under `flutter test`, not a phone's. 42 are
cartoons and one a tube (oxytocin); 42 take PDBe's DSSP and one PyMOL's.

## AlphaFold DB models, for proteins built on demand

A protein resolved on demand (`pipeline/resolver/`) has no entry chosen by
hand. Where a whole search finds no experimental entry that passes
([above](#experimental-entries-for-proteins-built-on-demand)), its fold page
draws AlphaFold DB's model of its UniProt sequence, still, under rule R5.4 as
the user settled it on 2026-10-08 (for every protein built on demand, until
2026-10-12). `alphafold.py`'s docstring has each rule; in short:

| | |
|---|---|
| entry | the canonical one, `AF-<accession>-F1`: the API also answers with one for each isoform |
| sequence | the record's protein, the same length, and no more residues different than the resolver allows the record |
| span | the mature chain or chains: first kept region to the last |
| gate | no model where the mean pLDDT over the span is under 50 |
| colour | pLDDT's four bands (90, 70, 50), one node each: `plddtVeryHigh`, `plddtConfident`, `plddtLow`, `plddtVeryLow` |
| bridges | UniProt's pairs, each drawn only where the model's sulfurs are within 2.5 Å, as the `bonds` node |
| shape | a cartoon; a tube where PyMOL finds no helix or strand |
| size | sampling by length, lowered until the compiled scene fits `alphafold.BUDGET_BYTES` (1.3 MB), refused if none does |

The worker bakes it (the resolver's README). By hand, from what the service
serves of a protein, writing nothing anywhere but `--out`:

```sh
pipeline/structure/venv/bin/python pipeline/structure/alphafold.py \
    --service https://helix-peak-backend.onrender.com --slug oca2 --out /tmp/oca2
```

That leaves the model file, the bake's copy of it (`span.pdb`), PyMOL's export,
`oca2.glb`, `oca2.fsceneb` and `described.json`: the fold page's seven words
(`chrome`), its `chains` and the provenance. It needs PyMOL, this directory's
`venv`, and the app's checkout beside this repository for the importer (`--app`,
`--dart`). A bake is about two seconds, and gives the same bytes each time.

How it is made, where that is not what the twenty's bake does:

- **The bands are cut from the mesh, not asked of PyMOL.** The chain is
  exported whole by `bake.export`, exactly as the twenty's are, and each
  triangle then goes to the band of the residue it lies on. Shown a part at a
  time, PyMOL draws another cartoon, cut short at every break: on B2M the parts
  have 180 fewer triangles than the whole's 8,298.
- **A triangle's residue is the nearer end of the nearest stretch between two
  consecutive CA atoms**, which is where a cartoon changes colour. Nearest CA
  alone starves a strand's residues: flat arrows are smoothed past their atoms.
- **The bake is handed a copy of the model**: the span's atoms, and `SSBOND`
  records for the bridges that pass. An AlphaFold file names none, and this
  way `build_bonds` reads them as it reads an entry's. Nothing in `bake.py`
  knows the difference.
- **The frame's audit is `verify_frame.py`'s measure, not its numbers.** "Every
  CA within 0.25 Å" holds only without strands (B2M: mean 0.44, max 2.82 Å). It
  is held to `folding/check_folding.py`'s tolerances instead: no CA further
  than 3 Å, and a quarter at least within 0.25 Å (0.6 Å on a tube). A turned
  or shifted frame fails that by tens of ångströms. And it is taken with a k-d
  tree: `verify_frame.py`'s CA-by-vertex matrix is 30 GB on a long model.
- **The mean pLDDT is over the span, from the file's own B-factors**, so on a
  whole chain it can differ from the API's in the second decimal.

Measured on 2026-10-08, by hand:

| protein | span | mean pLDDT | sampling | scene |
|---|---|---:|---:|---:|
| sarcolipin | 1–31 | 91.6 | 8 | 120 kB |
| B2M | 21–119 | 97.0 | 8 | 210 kB, bridge 45–100 at 2.02 Å |
| OCA2 | 1–838 | 73.8 | 1 | 384 kB |
| TSC2 | 1–1,807 | 67.9 | 1 | 767 kB |
| filamin A | 1–2,647 | 76.5 | 1 | 914 kB |
| mTOR | 1–2,549 | 78.0 | 1 | 1.22 MB |

At sampling 1 a helical protein costs about 480 B a residue and a strand-rich
one about 350 B, so the 900 kB the twenty's scenes are held to would stop
between 1,900 and 2,600 residues, short of AlphaFold DB's own limit of 2,700,
and refuse some 240 buildable proteins as too large. A model made on demand
is therefore held to 1.3 MB (`alphafold.BUDGET_BYTES`, the user's decision of
2026-10-09), which the longest chain AlphaFold DB models fits; the twenty keep
their 900 kB. Insulin (48.0) and
BRCA1 (41.6) are refused at the gate. Of 150 reviewed human proteins sampled
that day, AlphaFold DB's sequence was UniProt's current one for 147; the other
three differ in length, and are refused as another sequence.

`test_alphafold.py` holds the rules with numpy alone (`fixtures/alphafold/`
stands in for the network), and the mesh and the whole bake where trimesh,
PyMOL and the importer are.

## Licensing

PDB coordinate files are distributed by RCSB without restriction. If a model
appears anywhere citable, cite the entry.

AlphaFold DB's models are CC BY 4.0 (EMBL-EBI and Google DeepMind), so a page
that draws one credits it: the app's About sheet does, where it names the
model. Cite Jumper et al. 2021 (Nature 596:583) and Varadi et al. 2024
(Nucleic Acids Res. 52:D368). Each stored model's provenance carries the entry,
its version, the licence and both citations.
