# The protein pipeline

The bake tools. They were `helix-peek/tool/` until Phase 3 of
`HANDOFF-ONDEMAND.md`, when the app stopped carrying protein data and storage
became the place the tracks live. Everything here runs by hand, offline, and the
web service never imports it (yet: Phase 6 turns the record builder into a
library the resolver calls).

One table drives the twenty curated proteins: [`targets.py`](targets.py). A row
says where the gene comes from, which regions and disulfides the precursor has,
which chains of which PDB entry the fold is cut from, and whether it is scored.
The rows' hand-written prose -- display name, summary, chain name, facts, chain
tints and the fold page's sentences -- is [`curated/catalog.json`](curated/catalog.json),
which is where it is edited. `check_assets.py` holds the two to each other.

Each row produces these tracks, written under `pipeline/data/` in the layout the
app's bundle once had, then uploaded:

| track | baker | environment | what it is |
|---|---|---|---|
| `assets/mock/gene_<gene>.json` | [`mock/build_gene_record.py`](mock/build_gene_record.py) | `.venv` | the gene record the walk opens on (`record`) |
| `assets/constraint/<slug>_esm_constraint.json` | [`constraint/score_protein.py`](constraint/score_protein.py) | `pipeline/.esm-venv` | ESM-2 masked marginals, per residue |
| `assets/impact/<slug>_avi.json` | [`impact/bake_impact.py`](impact/bake_impact.py) | `pipeline/impact/venv` | AlphaGenome Variant Impact, per base |
| `assets/impact_explanations/<slug>.json` | [`impact/bake_explanations.py`](impact/bake_explanations.py) | `pipeline/impact/venv` | AVI feature attributions (insulin, hemoglobin, CFTR) |
| `assets/clinvar/<slug>_clinvar.json` | [`clinvar/bake_clinvar.py`](clinvar/bake_clinvar.py) | `.venv` | ClinVar single-base variants on the drawn gene |
| `assets/models/<slug>.glb` | [`structure/bake.py`](structure/bake.py) | `pipeline/structure/venv` | the fold, as named meshes |

Each directory's README has its detail, how to make its environment, and what a
correct result looks like.

## Storage is the source

The tracks the app reads are the objects in Supabase Storage, named by the
`protein_track` rows. To have files to work on, bring them back:

```sh
python3 pipeline/fetch_tracks.py            # all twenty, every family, sha256-verified
python3 pipeline/check_assets.py            # the tracks against each other and the tables
python3 pipeline/check_assets.py --against https://helix-peak-backend.onrender.com
```

Against fetched files the offline check proves the stored tracks agree with each
other and with the tables. It compares storage with something else only when
`pipeline/data/` holds a fresh bake, before it is uploaded.

## Rebaking one protein

From the repository root. The record comes first: the constraint scorer and the
impact bake read it, and the ClinVar bake places records through the impact
track's coordinate map. Structures are independent.

```sh
NCBI_EMAIL=you@example.com .venv/bin/python pipeline/mock/build_gene_record.py --target <slug>
pipeline/.esm-venv/bin/python -u pipeline/constraint/score_protein.py --target <slug>
ALPHAGENOME_API_KEY=... pipeline/impact/venv/bin/python -u pipeline/impact/bake_impact.py --target <slug>
NCBI_EMAIL=you@example.com .venv/bin/python -u pipeline/clinvar/bake_clinvar.py --target <slug>
pipeline/structure/venv/bin/python pipeline/structure/bake.py --target <slug>
python3 pipeline/check_assets.py
set -a && . ./.env && set +a
.venv/bin/python pipeline/upload_tracks.py --kind <family> --target <slug> --dry-run
.venv/bin/python pipeline/upload_tracks.py --kind <family> --target <slug>
```

`--all` re-fetches every record from NCBI as it is today; `--target` leaves the
others byte for byte. The uploader validates before it writes, stores plain
bytes at `<kind>/<slug>.<sha12>.<ext>` with immutable caching, and deletes an
object only once no row points at it. `--kind structure` refuses until Phase 9
compiles `.fsceneb` containers in a worker.

Rough costs on an M-series laptop: records under a minute; structures about a
minute; scores about eighty-five minutes for all twenty, forty-five of them
dystrophin's and twenty-one CFTR's; impact tracks fifteen to thirty minutes, set
by the Atlas's per-minute quota; ClinVar about seventeen minutes, most of it
dystrophin's 12,000 records and CFTR's 6,500.

## The catalog rows

```sh
python3 pipeline/seed_catalog.py --dry-run                     # print the rows
set -a && . ./.env && set +a
.venv/bin/python pipeline/seed_catalog.py --check              # diff against the live rows
.venv/bin/python pipeline/seed_catalog.py                      # upsert them
```

A real run upserts `protein` rows and replaces their aliases; it inserts only
the track rows that are missing, so it never undoes an upload.

## Tests

`pytest` at the repository root runs the service's suite and the pipeline's own
tests (`pytest.ini`). The two that hold stored tracks to their gates skip, and
say so, until `fetch_tracks.py` has filled `pipeline/data/`.

## Where things are

`paths.py` names every location outside the code, and each can be moved:

| variable | default | what |
|---|---|---|
| `HELIXPEEK_DATA` | `pipeline/data` | baked and fetched tracks |
| `HELIXPEEK_CLIENT` | `../helix-peek` | the app checkout, for `fetch_tracks.py --verify-client` |
| `ALPHAGENOME_SKILLS` | `../.claude/skills` | the Atlas skills the impact bakes call |

## Adding a protein

Read [the pipeline rules](../../helix-peek/docs/protein-pipeline-rules.md)
first: they are what the first ten had to be corrected to agree on. Check the
gene against its record before writing anything --
[protein-verification.md](../../helix-peek/docs/protein-verification.md) is how
the second ten were checked. Then add a row to `targets.py` and one to
`curated/catalog.json`, bake with `--target <slug>`, run `check_assets.py`,
upload, and seed. Phase 6 replaces the hand-written row with a resolver.

Two things to check before committing to a protein:

- **The gene has to fit one screen.** The gene page sizes to a single screen
  with a two-point floor on a cell, so past about 24,000 bases there is no
  picture left. `build_gene_record.py` shortens introns past that and records
  the scale, but a gene whose *exons* alone exceed the budget has no honest page
  and the script says so rather than guessing.
- **The structure has to exist.** Where the full length has never been solved
  in one piece -- p53's disordered stretches, the other 3,447 residues of
  dystrophin, all but APP's E1 domain -- the row names the span that has, and
  the fold page's sentence says what is missing.
