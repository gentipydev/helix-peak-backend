# Assemblies: a molecule of more than one gene, as a morph pair

The catalog is twenty proteins, one gene each, and the walk draws only what a
gene makes (R5.1): the `hemoglobin` row is HBB's beta chain alone. The lab's
oxygen feature needs the whole tetramer, two alpha chains from HBA1 and two
beta chains from HBB, moving between its two states. That molecule cannot be
a catalog row, so it is a row of `ASSEMBLIES` in `assemblies.py` instead,
separate from `TARGETS`, and nothing that seeds or serves the catalog reads
it.

```sh
.venv/Scripts/python pipeline/fetch_tracks.py --kind record         # the HBB record, to check the numbering
<structure bake python> pipeline/assemblies/bake_assembly.py --all  # PyMOL on PATH
.venv/Scripts/python pipeline/assemblies/check_assembly.py
.venv/Scripts/python pipeline/assemblies/seed_assemblies.py --dry-run
.venv/Scripts/python pipeline/assemblies/upload_assemblies.py --dry-run
```

## Its own tables, not a namespaced slug

`migrations/0006_assemblies.sql` adds `assembly` and `assembly_track`, the
second in `protein_track`'s shape. A namespaced slug in `protein` was the other
way, and it does not work here: `protein_track.slug` is a foreign key into
`protein`, so a track for the tetramer would need a `protein` row, and every
`protein` row is one `/catalog`, `/catalog/search` and `/proteins/suggest`
serve. Keeping it out would mean a filter in each of them, and a filter
forgotten once puts a twenty-first protein in the catalog. Tables of their own
cannot leak: no catalog query names them (pinned by
`tests/test_assemblies_api.py`), and the service reaches an assembly only by
its own slug, at `/assembly/{slug}/tracks`.

`seed_assemblies.py` writes those two tables and no others. After a write, and
in `--check`, it asserts that the catalog did not move: `seed_catalog.py
--check` finds no difference in the twenty's live rows, and `GET /catalog`
returns exactly twenty. Both held on 2026-09-27, with 0006 not yet applied.

## The pair, and whether its frames agree

Hemoglobin A at 1.25 A in both states, from one study (Park, Yokoyama,
Shibayama, Shiro and Tame, 2006):

| state | entry | what the entry holds | how the tetramer is built |
|---|---|---|---|
| tense | `2DN2`, deoxy | four chains, alpha A and C, beta B and D | the entry itself |
| relaxed | `2DN1`, oxy | one alpha (A), one beta (B), oxygen on each haem | biological assembly 1: A and B, and their copy under the entry's second BIOMT operator (a two-fold: x and y swapped, z negated), lettered C and D. Four chains, checked. |

**As deposited, they do not share a frame.** The two are different crystals
(`2DN2` is P2<sub>1</sub>, `2DN1` P4<sub>1</sub>2<sub>1</sub>2), and their
coordinates agree about nothing: matched residue for residue, the 570 CA atoms
both place are 58.6 A apart (RMSD), and the best superposition of the whole
tetramer turns one onto the other by 126.7 degrees. A morph between those
coordinates would be the crystal lattice turning, not the molecule. And no
other pair would be different: the T to R change does not happen inside a
crystal, so a deoxy and an oxy tetramer always come from two lattices.

So the pair is brought into one frame, as T and R hemoglobin always are to be
compared: the relaxed tetramer is superposed onto the tense one, holding one
alpha-beta dimer still (Baldwin and Chothia's convention). The tense entry's
frame is the pair's.

| | CA RMSD |
|---|---:|
| as deposited | 58.61 A |
| alpha1-beta1, superposed on itself | 0.93 A |
| alpha2-beta2, after that, which turns 14.1 degrees | 5.19 A |
| the tetramer, after that | 3.73 A |
| the tetramer, best fit on all of it | 2.41 A |

The 0.93 A is the change inside a dimer; the 14.1 degree turn of the other is
the quaternary switch the textbooks put at about 15.

**Both exports share that frame, proved by `verify_frame.py`.** Each state is
exported by the structure bake's own `export`, with one origin between them
(`write_pml`'s `origin`), and `verify_frame.py` audits each export against that
one origin: every CA of every chain, in both states, lies within 0.25 A of its
own exported ribbon, as insulin's does. Two exports that are each a pure
translation of their coordinates, by the same translation, of coordinates in
one frame, are one frame.

## What the bake writes

Under `pipeline/data/assets/assemblies/hemoglobin-a/`:

- `tense.glb` and `relaxed.glb` (3.6 MB each): the two states' cartoons, cut
  to one box together (the structure bake's `normalisation` over both, longest
  axis 1.0), nodes `alpha1`, `beta1`, `alpha2`, `beta2`.
- `hemoglobin-a_morph.json` (171 kB), what a viewer animates: for each
  subunit, each residue both states place (a few chain ends are unplaced in
  one; `dropped` lists them), numbered as the precursor numbers it, with its
  CA in each state in that one frame and each state's own secondary
  structure; each subunit's haem iron in each state, and the oxygen bound to it
  in the relaxed one. The frame carries its centre, its length in angstroms,
  the box both models fill, the superposition and its RMSDs, and both audits.

The entries number each chain from Val1, after the initiator methionine,
which the precursor numbers Val2; the bake checks the beta chains' letters
against the catalog's own HBB record at that offset.

`check_assembly.py`, and the uploader before it writes, rebuild both molecules
from the entries and hold the morph to them: four chains in each state, every
CA the entry's in the stated frame, a chain of 2.8 to 4.4 A steps, one haem
iron per subunit per state, oxygen only in the oxy state, both audits within
0.25 A, and the models beside it the ones it names.

Nothing is uploaded and 0006 is not applied.
