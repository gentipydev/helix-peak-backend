# Baking `folding`: each fold's CA trace, for the fold animation

The walk's fold page plays the animation, drawn in the same scene as the
model it ends on (the Lab's own fold screen is gone).

The `structure` track is a baked mesh. It draws the finished fold, but it
cannot morph and it carries no residues, so an animation of a chain folding
has nothing to move. This track is what it moves: per residue of each chain
the structure exports, where its CA sits in the finished fold, whether that is
in a helix, a strand or a coil, and whether it has a place at all.

It is a separate track kind, `folding`. Nothing is added to `targets.py`,
`curated/catalog.json` or the `protein` table, and no other track is touched.

```sh
.venv/Scripts/python pipeline/fetch_tracks.py --kind record --kind structure   # digest-checked
<structure bake python> pipeline/folding/bake_folding.py --all                  # PyMOL on PATH
.venv/Scripts/python pipeline/folding/check_folding.py                          # exits 1 on any failure
.venv/Scripts/python pipeline/upload_tracks.py --kind folding --dry-run
```

The bake re-runs the structure bake's own export, so it runs in that bake's
environment: `pipeline/structure/venv` on the Mac, or on Windows conda-forge's
PyMOL 3.1.0 (see `structure/README.md`). The check, and the uploader that runs
it before it writes, need only numpy, on the backend's `.venv`.

Upload only when asked: `DATABASE_URL`, `SUPABASE_URL` and `SUPABASE_SERVICE_KEY`
set, `migrations/0005_folding.sql` applied (after `0004`), then the same command
without `--dry-run`. Objects go to the `tracks` bucket as
`folding/<slug>.<sha12>.json`.

Which proteins get the track: every target with a structure, all twenty
(`FOLDING_TARGETS` in `bake_folding.py`).

## What it holds

One JSON object per protein: its `slug`, `gene`, `uniprot` and `pdb`, the
`frame`, and one entry in `chains` per node of the structure model
(`chainA`, `chainB`), each with the entry's chain, the `offset` from the
entry's numbering to the precursor's, the `first` and `last` residue, the
positions where the entry's residue is not the gene's (`entry_differs`), and
`residues`, one line each:

```json
{"n": 24, "aa": "K", "state": "disordered"},
{"n": 118, "aa": "A", "state": "ordered", "ss": "strand", "ca": [-0.02149, 0.24273, -0.31324]},
```

- **Which residues.** Each chain covers the mature chain it is cut from, in the
  precursor's numbering, the walk's: the kept region of a cleaved precursor
  (insulin's A chain, 90-110), or an uncut chain less what the table removes
  from its ends (myoglobin, 2-154, not its globin-fold domain), clipped to the
  span the structure exports (p53's core, 96-289). Every residue of it is
  listed once.
- **`state`.** `ordered`: the entry locates its CA. `disordered`: it was in the
  experiment and has no place in it -- listed in REMARK 465, or modelled
  without its CA -- between residues that do. `absent`: beyond the ends of
  what the entry holds at all. A fold animation should leave `disordered`
  residues loose and never draw `absent` ones.
- **`ca`**, for ordered residues: the CA in the stored structure model's frame,
  centred with its longest axis 1.0, to five places.
- **`ss`**, for ordered residues: `helix`, `strand` or `coil`, from the entry's
  own HELIX and SHEET records. Every helix class is a helix (3-10 included:
  SOD1 has no other), every strand of every sheet a strand.
- **`aa`**: the record's letter, what the gene makes. Where the entry carries
  another (an engineered mutation, a sequence conflict, each declared in its
  SEQADV records) the coordinates are the entry's and the letter the gene's.

Schema 2 adds two fields, so that an animation drawn from the track can end on
the model rather than near it:

- **`bridges`**: each disulfide the model draws, as the atoms its rods run
  through, in the same frame:

  ```json
  {"a": 31, "a_node": "chainB", "b": 96, "b_node": "chainA", "path": [[ca], [cb], [sg], [sg], [cb], [ca]]}
  ```

  From the entry's SSBOND records and `pdb.atoms`, the structure bake's own
  reading, and only where the model has a `bonds` node: none for the eight
  that draw none, whatever UniProt lists. Numbered in the precursor, the lower
  first.
- **`cartoon`**: what the chains are drawn as, in angstroms: `representation`
  (`cartoon` or `tube`), `helix` and `strand` as `half_width` and
  `half_thickness`, `loop_radius`, `tube_radius` and the bridges'
  `rod_radius`. PyMOL 3.1.0's settings (`cartoon_oval_length` 1.35 and
  `_width` 0.25, `cartoon_rect_length` 1.4 and `_width` 0.4,
  `cartoon_loop_radius` 0.2) and the structure bake's radii. Measured on the
  stored models they are half-extents: insulin's helices are 2.7 A wide and
  0.5 A thick, its loops 0.4 A across. The bake refuses to run against a PyMOL
  whose settings differ.

Schema 2 changes an existing track's format, which the backend's contract
(`CLAUDE.md`) otherwise reserves for a new kind. It was done in place on
purpose (2026-09-28): the keys are only added, a schema-1 reader ignores
them, and the fold page is the track's one reader.

## The frame

The structure bake exports PyMOL's cartoon as a pure translation of the PDB
frame, then centres it and divides by its longest axis (`bake.py`'s `export`
and `normalisation`). This re-runs both and holds the re-export to the stored
`.glb` vertex for vertex: on Windows every vertex lands within 3e-6 model units
of the stored one, and a frame that does not match within 1e-5 is refused. So
each CA sits where the fold page's model says it does, and the last frame of a
fold animation lands on exactly the fold the walk draws. The frame's
`centre_angstrom` and `length_angstrom` put a PDB point `p` at
`(p - centre) / length`; `glb_sha256` names the stored model, and `bounds` is
its bounding box in model units, which is what a viewer frames it by. The
app's fold page frames the loaded model with `PerspectiveCamera.framing` on
exactly that box, so a fold drawn from this track alone is framed the same way
without loading the model.

A fit of the CA atoms to the stored ribbon, which is how `structure_ar` sizes a
model with no bridges, was measured against this frame first. It is exact where
the model has bridges and within 0.25 A for five of the other seven, but
ubiquitin's comes out 2.1% large and glucagon's single helix 3.9% small, which
moves its ends by 1.5 A. A lone helix barely fixes the ribbon's scale along its
axis. That fit is not used here.

## What the check holds

`check_folding.py`, and the uploader before it writes, hold each track to:

- its own protein and entry, and the structure's chains, node for node;
- the mature chain, residue for residue: its own reading of the region table,
  clipped to the exported span, so the residue count is the chain's length,
  and the record's letters, its mature peptide's where it names one;
- the entry: each ordered CA is the entry's CA for that residue, placed by the
  frame; a letter where the entry differs is one the entry declares;
- a chain: every CA 2.8 to 4.4 A from the next (3.8 A for a trans peptide bond,
  2.9 A for amylase's two cis ones, and 1HGU, at 2.5 A resolution, stretches one
  to 4.33 A), and across residues with no place, no more than 3.8 A a residue;
- the stored model: every helix and coil CA within 0.25 A of its own chain's
  ribbon (0.6 A on the two tube-drawn peptides), which is `verify_frame.py`'s
  criterion for insulin, now held on all twenty. A strand's CA, which PyMOL's
  flattened arrows smooth past, within 3 A;
- states that say only what they can, with nothing absent between two
  residues that are there;
- the cartoon the model is drawn with, and exactly the model's bridges: one
  per SSBOND pair where it draws them, each a pair of ordered cysteines whose
  path starts and ends on their CA atoms, with a cysteine's bond lengths, and
  whose CB and SG atoms are vertices of the stored `bonds` node, where its
  rods end.

## What it found

| slug | residues | ordered | helix | strand | disordered | absent | bridges | bytes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| insulin | 51 | 51 | 30 | 0 | 0 | 0 | 3 | 8,512 |
| hemoglobin | 146 | 145 | 127 | 0 | 1 | 0 | 0 | 16,082 |
| myoglobin | 153 | 149 | 132 | 0 | 4 | 0 | 0 | 16,810 |
| p53 | 194 | 194 | 24 | 64 | 0 | 0 | 0 | 21,080 |
| lysozyme | 130 | 130 | 55 | 8 | 0 | 0 | 4 | 16,797 |
| relaxin | 53 | 51 | 42 | 4 | 1 | 1 | 3 | 8,542 |
| oxytocin | 9 | 9 | 0 | 0 | 0 | 0 | 1 | 2,807 |
| somatotropin | 191 | 186 | 104 | 0 | 5 | 0 | 2 | 22,678 |
| ubiquitin | 76 | 76 | 16 | 33 | 0 | 0 | 0 | 9,026 |
| dystrophin | 238 | 238 | 167 | 3 | 0 | 0 | 0 | 25,674 |
| vasopressin | 9 | 9 | 0 | 0 | 0 | 0 | 1 | 2,818 |
| glucagon | 29 | 29 | 28 | 0 | 0 | 0 | 0 | 4,256 |
| app | 753 | 162 | 27 | 58 | 11 | 580 | 6 | 50,920 |
| cftr | 1,480 | 1,139 | 807 | 91 | 341 | 0 | 0 | 136,212 |
| erythropoietin | 166 | 166 | 105 | 6 | 0 | 0 | 2 | 19,810 |
| leptin | 146 | 130 | 95 | 0 | 16 | 0 | 1 | 16,040 |
| tnf | 157 | 153 | 5 | 78 | 4 | 0 | 1 | 17,805 |
| sod1 | 153 | 153 | 17 | 58 | 0 | 0 | 1 | 17,356 |
| amylase | 496 | 496 | 134 | 85 | 0 | 0 | 5 | 54,662 |
| prion | 208 | 109 | 71 | 13 | 98 | 1 | 1 | 18,323 |

486,210 bytes for all twenty, 3,775 CA atoms. The playbook's 9,991 is every CA
in the twenty files, every chain; the folds export 3,775 of them.

- **The prion protein**: 23-230, of which only 117-225 has a shape. 24-116 and
  226-230 were in the crystal and never located, so they are `disordered`;
  23 is `absent`, before the construct 4KML crystallised.
- **CFTR**: 341 residues with no place, most of them the regulatory region,
  which the fold page already says is too mobile to resolve.
- **APP**: 580 `absent`, everything after the E1 domain, which is all the entry
  holds; its 11 `disordered` are the ends of that construct.
- **Relaxin**: 6RLX's B chain starts at Ser26, so Asp25 is `absent`.
