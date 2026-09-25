# Baking `assets/models/*.glb`

Everything here exists to regenerate twenty files: the folds the structure page
renders. It is offline, run by hand, and run once per protein. The `.glb`s it
writes land under `pipeline/data/assets/models/`; storage keeps each one beside
the `.fsceneb` compiled from it, which is what the app fetches.

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
structure row's provenance. Nothing here compiles a `.fsceneb` yet, so
`upload_tracks.py --kind structure` refuses every target until Phase 9 of
HANDOFF-ONDEMAND.md compiles them in a worker.

## Which entry, and why

All experimental. AlphaFold is not used anywhere: `AF-P01308` predicts insulin's
110-residue *preproprotein* — the third page's molecule, not this one — at a
mean pLDDT of 52.9 with 51% of residues below 50, which would draw a coil of
low-confidence loops exactly where the reader has been promised a hormone.

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

## Licensing

PDB coordinate files are distributed by RCSB without restriction. If a model
appears anywhere citable, cite the entry.
