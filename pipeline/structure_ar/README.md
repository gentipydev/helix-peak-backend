# Baking `structure_ar`: each fold as a USDZ, at its real size

iOS AR Quick Look reads USDZ and nothing else, and the walk's `structure` track
is a `.glb` normalised to a longest axis of 1.0, which throws the molecule's
size away. This is a separate track kind, `structure_ar`: the same meshes, back
at their real size, packaged as USDZ, with that size in angstroms in the row's
provenance. The `structure` rows, objects and provenance are not touched.

```sh
.venv/Scripts/python -m pip install -r pipeline/structure_ar/requirements.txt
.venv/Scripts/python pipeline/fetch_tracks.py --kind structure   # the stored .glbs, digest-checked
.venv/Scripts/python pipeline/structure_ar/bake_ar.py --all      # writes pipeline/data/assets/models_ar/
.venv/Scripts/python pipeline/structure_ar/check_ar.py           # every check below; exits 1 on any failure
.venv/Scripts/python pipeline/upload_tracks.py --kind structure_ar --dry-run
```

Upload only when asked: `DATABASE_URL`, `SUPABASE_URL` and `SUPABASE_SERVICE_KEY`
set, `migrations/0003_structure_ar.sql` applied, then the same command without
`--dry-run`. Objects go to the `models` bucket as
`structure_ar/<slug>.<sha12>.usdz` (`model/vnd.usdz+zip`), with its `.glb`
beside it.

Which proteins get the track: every target with a structure (`AR_TARGETS` in
`bake_ar.py`), which is all twenty.

## Tooling

Nothing here needs PyMOL or trimesh. The meshes are read from the stored `.glb`
with numpy (`glb.py`), and the USDZ is written by Pixar's own `usd-core`
(26.8, `pip install usd-core`, wheels for Windows, macOS and Linux). usd-core
also brings the validators `check_ar.py` runs, `UsdzPackageValidator` and
`RootPackageValidator` among them. Nothing else was needed: no Reality
Converter, no Xcode, no `usdzconvert`.

## Where the meshes come from

The stored `.glb` of the `structure` track, never a new PyMOL run: the walk's
fold and the AR fold are the same triangles. Its sha256 is written into the
USDZ's metadata and checked against the file beside it.

No code is shared with `structure/bake.py`. Sharing would have to be proved
by re-baking its twenty models byte for byte, which needs PyMOL, and PyMOL is
not installed here. The two small readers the size needs (SSBOND pairs, atom
coordinates) are written again in `bake_ar.py`, reading the file exactly as the
structure bake does.

## How big it really is

`structure/bake.py` exports the cartoon as a pure translation of the PDB frame,
then centres it and divides by its longest axis, L. The real size is the
stored model's extent times L, and L is found two ways:

- **Exactly, from the bridges.** Where the model has a `bonds` node, each bridge
  was built as five rods and four joint spheres centred on its CB, SG, SG and
  CB atoms. The centroid of each sphere's vertices is that atom, in the model;
  the PDB gives it in angstroms. Twelve or more points fix the scale and shift
  exactly (residual under 0.001 A). Thirteen of the twenty.
- **By fitting the ribbon to its CA atoms**, for the seven with no bridges
  (hemoglobin, myoglobin, p53, ubiquitin, dystrophin, glucagon, CFTR): the
  chains and span the bake exported, each CA moved to its nearest model
  vertex, the least-squares scale and shift solved, and repeated until it
  settles. Measured against the exact answer on the thirteen that have both,
  it lands within 1.2% for a cartoon. It is 8-9% out for the two tube-drawn
  peptides, oxytocin and vasopressin, which is why a model with bridges never
  uses it.

Insulin is the check: its bake logged "24.31 x 18.66 x 20.33 A" before
normalising, and this bake gives 24.31 x 18.66 x 20.33 A.

| slug | Å (x × y × z) | size from |
|---|---|---|
| insulin | 24.31 × 18.66 × 20.33 | bridges |
| hemoglobin | 21.18 × 39.98 × 39.39 | CA fit, rms 0.19 Å |
| myoglobin | 41.01 × 32.66 × 37.99 | CA fit, rms 0.20 Å |
| p53 | 41.34 × 36.80 × 42.60 | CA fit, rms 0.51 Å |
| lysozyme | 32.01 × 38.28 × 30.70 | bridges |
| relaxin | 27.82 × 23.87 × 21.03 | bridges |
| oxytocin | 16.59 × 5.75 × 12.24 | bridges |
| somatotropin | 44.25 × 51.19 × 38.29 | bridges |
| ubiquitin | 25.17 × 29.07 × 32.57 | CA fit, rms 0.62 Å |
| dystrophin | 47.28 × 66.38 × 69.05 | CA fit, rms 0.19 Å |
| vasopressin | 15.34 × 4.76 × 11.61 | bridges |
| glucagon | 29.58 × 31.86 × 6.00 | CA fit, rms 0.27 Å |
| app | 33.15 × 43.68 × 49.59 | bridges |
| cftr | 89.18 × 103.28 × 99.55 | CA fit, rms 0.35 Å |
| erythropoietin | 49.96 × 36.39 × 34.09 | bridges |
| leptin | 29.28 × 34.54 × 42.46 | bridges |
| tnf | 55.08 × 26.47 × 47.30 | bridges |
| sod1 | 38.87 × 42.94 × 31.12 | bridges |
| amylase | 51.77 × 69.99 × 51.52 | bridges |
| prion | 30.04 × 26.80 × 41.89 | bridges |

## What the USDZ is

- One unit is one centimetre (`metersPerUnit` 0.01) and the points are in
  angstroms, so AR Quick Look shows 1 Å as 1 cm. Y up.
- `/Fold` is the default prim. Under it, one mesh per node of the `.glb`
  (`chainA`, `chainB`, `bonds`), with the same names, vertices, normals and
  triangles, at one scale.
- It rests on its lowest point, centred over the anchor, because AR Quick
  Look puts the model's origin on the floor it finds.
- No colour is baked. Which colour a chain takes is the app's decision
  (`targets.Chain`: "a colour written down twice is a colour that will
  disagree"). Every mesh has the walk's matte material (UsdPreviewSurface,
  roughness 0.65, metallic 0) at USD's default neutral base colour.
- The root layer's `customLayerData["helixpeek"]` holds the provenance the
  upload writes to the row: the entry, the nodes, `bbox_angstrom`, the scale,
  how the size was found and its residual, the source `.glb`'s sha256, and
  the usd-core version.
- Two bakes of the same `.glb` give the same bytes.

## Beside it, a `.glb` for Android

Android's Scene Viewer takes a `.glb`, in metres. The same placed meshes go
beside each USDZ as `<slug>.glb`, every point times 0.01, so 1 Å is 1 cm there
too. The upload puts it next to the USDZ and names it in the row's
`provenance.glb` (`path`, `sha256`, `bytes`), as a `structure` row names the
`.glb` beside its `.fsceneb`. `check_ar.py` holds it to the USDZ's own points.
