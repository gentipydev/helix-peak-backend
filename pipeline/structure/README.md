# Baking `assets/models/*.glb`

Everything here exists to regenerate twenty files: the folds the structure page
renders. It is offline, run by hand, and run once per protein. The `.glb`s it
writes land under `pipeline/data/assets/models/`; storage keeps each one beside
the `.fsceneb` compiled from it, which is what the app fetches.

One thing here is not run by hand: `alphafold.py` makes the model of a protein
built on demand, from AlphaFold DB, and the resolver's worker runs it
([below](#alphafold-db-models-for-proteins-built-on-demand)). It uses this
bake's own export, and changes none of the twenty.

It is kept because a stored binary that cannot be regenerated is a
liability. If a model ever needs to change — a different colour split, a
different span, a lighter mesh — this is the only record of how it was made.

```sh
brew install pymol                 # 3.1.0; build-time only, nothing ships
python3 -m venv pipeline/structure/venv
pipeline/structure/venv/bin/pip install -r pipeline/structure/requirements.txt

pipeline/structure/venv/bin/python pipeline/structure/bake.py --all   # from the repo root
```

Entries are downloaded from RCSB into `structures/` on first use and kept.
Each bake writes its intermediates to `output/<slug>/` — the generated `.pml`,
one `.obj` per chain — and copies the finished model into
`pipeline/data/assets/models/`. `output/` and `venv/` are gitignored.

The app ships no model. It fetches a compiled `.fsceneb` from the `models`
bucket, and the `.glb` it was compiled from is stored beside it, named in the
structure row's provenance. Nothing compiles the twenty's by hand, so
`upload_tracks.py --kind structure` refuses every target: their scenes are the
ones the app's build hook compiled while it still shipped them. A scene is
compiled now by flutter_scene's own importer, from the app's checkout:

```sh
cd ../helix-peek && dart run flutter_scene:import -i model.glb -o model.fsceneb
```

Run on insulin's stored `.glb` it gives the stored `.fsceneb` byte for byte.
The container is versioned, so the importer has to be the flutter_scene the
stored scenes were compiled by (0.23.0, `alphafold.FLUTTER_SCENE`).

## Which entry, and why

All experimental, for the twenty. AlphaFold is not used for any of them:
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

Every deviation from UniProt in these entries, and why each was chosen over the
alternatives, is in `docs/protein-verification.md`.

The node names are the contract with `structure_view.dart`, which looks them up
to assign a material each. Renaming them here takes their colours off there.

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

They need numpy and nothing else, so a baker that reads the stored model, as
`structure_ar/` does, runs on the backend's `.venv` without PyMOL. The three
moved out of `bake.py`, `structure_ar/bake_ar.py` and `verify_frame.py`
unchanged.

`folding/` reads the entry through `pdb.py` too, but takes its frame from
`bake.py` itself: `export` and `normalisation`, re-run and held to the stored
model vertex for vertex. That is exact on all twenty, where `frame.py`'s ribbon
fit is not (see `folding/README.md`).

For an export outside the twenty, `write_pml` and `export` take a `source` (an
assembly built from an entry, loaded in place of the entry itself) and an
`origin` (one export origin for two exports, so that they share a frame), and
`verify_frame.py` takes `--pdb`, `--out`, `--origin` and `--chain`. Unset, each
is what the twenty are baked and audited with: re-baking all twenty after they
were added gave the same `.glb` files byte for byte, the same `.pml` scripts
and the same `verify_frame.py` report.

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

## AlphaFold DB models, for proteins built on demand

A protein resolved on demand (`pipeline/resolver/`) has no entry chosen by
hand. Its fold page draws AlphaFold DB's model of its UniProt sequence, under
rule R5.4 as the user settled it on 2026-10-08. `alphafold.py`'s docstring has
each rule; in short:

| | |
|---|---|
| entry | the canonical one, `AF-<accession>-F1`: the API also answers with one for each isoform |
| sequence | the record's protein, the same length, and no more residues different than the resolver allows the record |
| span | the mature chain or chains: first kept region to the last |
| gate | no model where the mean pLDDT over the span is under 50 |
| colour | pLDDT's four bands (90, 70, 50), one node each: `plddtVeryHigh`, `plddtConfident`, `plddtLow`, `plddtVeryLow` |
| bridges | UniProt's pairs, each drawn only where the model's sulfurs are within 2.5 Å, as the `bonds` node |
| shape | a cartoon; a tube where PyMOL finds no helix or strand |
| size | sampling by length, lowered until the compiled scene fits `FSCENEB_BUDGET_BYTES`, refused if none does |

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
| filamin A | 1–2,647 | 76.5 | 1 | 914 kB: over the budget, refused |
| mTOR | 1–2,549 | 78.0 | 1 | 1.22 MB: over the budget, refused |

At sampling 1 a helical protein costs about 480 B a residue and a strand-rich
one about 350 B, so the 900 kB every scene is held to stops between 1,900 and
2,600 residues, short of AlphaFold DB's own limit of 2,700. Insulin (48.0) and
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
