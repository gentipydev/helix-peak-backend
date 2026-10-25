# HANDOFF: ten more curated proteins, through a pipeline that adds more

Prepared 2026-10-17. Starts from app `c0a704c`, backend `d55644a` (both local,
not pushed). Read both CLAUDE.md files and `helix-peek/docs/protein-pipeline-rules.md`
first. Under each phase, write a short **What actually happened** when it is done.

## Decided by the user (2026-10-17)

1. **The ten below**, each landing in a cell type none of the twenty reach.
2. **The curated list grows, often.** The "exactly twenty" rule is gone (app
   `c0a704c`, backend `d55644a`). A session may add a curated protein to the
   app's `test/fixtures/` (its files, its row after the others), checked by
   `python3 tool/fixtures_added_only.py <start>`. The walk's tests then run on it.
   Leave the pipeline so the next batch is a list of genes plus hand-written words.
3. **Curated proteins get hand-written** `summary`, structure `sentence` and
   `semantics`, held to R4.5. Built-on-demand proteins stay templated.
4. **This session executes, live writes included**: `--dry-run` every upload and
   seed first, write only with every gate green. Never push. Ask before a
   migration, a deploy, money or secrets.

## The ten

Organ and cell come from the locus bake's `path_of` (Protein Atlas 25.1). The
entry is PDB policy 1's choice. Every one passed `scripts/dry_fold.py`, which
also builds the record under the 24,000-base page budget. None was built on
demand on 2026-10-17 (check `/proteins/built` again before seeding).

| order | slug (suggested) | gene | UniProt | organ → cell (drawn as) | entry | notes |
|---|---|---|---|---|---|---|
| 20 | phospholamban | PLN | P26678 | heart → cardiomyocytes (cardiomyocyte) | 1ZLL NMR, 52/52 | one chain of five |
| 21 | s100b | S100B | P04271 | brain → Bergmann glia (glial) | 5D7F X-ray 1.3 Å, 90/91 | |
| 22 | hcgbeta | CGB3 | P0DN86 | placenta → syncytiotrophoblasts (trophoblast) | 7FIG cryo-EM 3.9 Å, 110/145 | 6 bridges; CGB3/5/8 are copies, check the record |
| 23 | pf4 | PF4 | P02776 | bone marrow → platelets, lands in a megakaryocyte | 1F9Q X-ray 2.0 Å, 62/70 | |
| 24 | pappa | PAPPA | Q13219 | placenta → decidual stromal cells (spindle) | 8HGG cryo-EM 3.6 Å, 1,484/1,546 | 36 bridges; 1,627 aa; introns shortened |
| 25 | cystatinb | CSTB | P04080 | esophagus → esophageal apical cells (squamous class) | 4N6V X-ray 1.8 Å, 91/98 | |
| 26 | tyrp1 | TYRP1 | P17643 | skin → melanocytes (pigment class) | 9EY8 X-ray 2.2 Å, 447/513 | 7 bridges; **no legacy .pdb file**, mmCIF only |
| 27 | recoverin | RCVRN | P35243 | retina → rod photoreceptors (neuron) | 2D8N X-ray 2.2 Å, 196/199 | |
| 28 | thyroglobulin | TG | P01266 | thyroid → thyroid glandular cells (acinar) | 7B75 cryo-EM 3.2 Å, 2,483/2,749 | one chain of two; 57 bridges; 2,768 aa; introns shortened |
| 29 | uromodulin | UMOD | P07911 | kidney → ascending loop of Henle (epithelial) | 6TQL cryo-EM 4.0 Å, 292–584 | domain only, like p53 and APP |

Re-run the check with
`NCBI_EMAIL=gentiannurce@gmail.com .venv/bin/python scripts/dry_fold.py --out <dir> PLN S100B CGB3 PF4 PAPPA CSTB TYRP1 RCVRN TG UMOD`.
If one fails its verification, the alternates are PPBP (for PF4), A2M (lung
fibroblasts) or CNN1 (for PAPPA), and PRSS1 (for UMOD). Swap only with the user.

## What the pipeline lacks today (Phase 1 fixes it)

- **Rows are written by hand.** The pieces exist already:
  - `pipeline/resolver/resolve.py` `resolve()` derives a `Target` from UniProt
    and MANE (source, regions, disulfides, declared variants).
  - `structure/fold_job.choose` and `experimental._bake_target` derive its
    `Structure`.
  - `structure/pdb_words.describe` drafts the fold words.
- **Curated folds are baked from hand-kept legacy `.pdb` files** in
  `pipeline/structure/structures/` (author numbering). 9EY8 has none, so a new
  row needs the on-demand path's working copy (`normalize.working_copy`,
  precursor numbering), stored there.
- **`upload_tracks.py --kind structure` refuses every curated row.** The app's
  build hook no longer makes `.fsceneb` files. The worker compiles them with
  `alphafold.compile_scene(glb, out, app, dart)`.
- **No tool writes the app's test copies.** `test/fixtures/catalog.json` is a
  `/catalog` response. Per protein, the walk tests read `mock/gene_<gene>.json`,
  `constraint/`, `impact/`, `clinvar/` and `models/<slug>.glb`. The zoom tests
  read `test/features/zoom/fixtures/<slug>_locus.json`.

## Phases

1. **Pipeline** (offline, with tests, committed):
   - `scripts/draft_curated.py GENE…`: resolve and choose as `dry_fold.py` does,
     store the working copy, print the `targets.py` row and the
     `curated/catalog.json` entry (templated words marked for hand edit).
   - The twenty's `structure/bake.py` and `folding/bake_folding.py` bake from a
     working copy. Prove the twenty unchanged with backend CLAUDE.md's sha256
     proof for every baker you touch.
   - The structure upload compiles the scene for a curated row.
   - One runner, `pipeline/bake_curated.py --target <slug>`: record, structure
     and folding, constraint, impact, ClinVar, locus, then `check_assets.py`,
     each in its own environment.
   - A fixture sync (backend script) that writes one curated protein's stored
     tracks and its `/catalog` row into the app checkout.
   - `pipeline/README.md` "Adding a protein" becomes that sequence.

   **What actually happened (2026-10-17).** Backend `3fe60bf..9d6c8dc`, eight commits,
   not pushed; the worker restarted on them.
   - The twenty keep their legacy `.pdb` files: re-cut from working copies they would
     not be the same bytes (node and SSBOND order, DSSP vs depositors' records, altlocs,
     modified residues; the fidelity gate refuses EPO and somatotropin). A new row is cut
     from `structures/working/<slug>.pdb` (+ `<slug>.json` sidecar): `pdb.source_of`,
     `export_source`; `bake.py`, `bake_folding.py` and `check_folding.reading_of` read it
     (no dashes, offset 0, `experimental._unmeasured`, reach 2.5 Å). The on-demand path's
     code is untouched.
   - `upload_tracks.py --kind structure` compiles scenes (`compiled_scene`), exits 1 on
     a refusal, `--against` compares a dry run with the service. `fetch_tracks.py
     --check`; a not-yet-seeded slug no longer fails fetch or `--verify-client`.
   - Guards: `check_assets.py` facts rule (record-derived; holds for the twenty), chain
     rule, display rule, `DRAFT:` refusal, `--target`; `seed_catalog.py` refuses drafted
     words and a row over a protein built on demand (s100b and pf4's slugs are their
     genes'); `check_locus.py --target`; `check_built.py` no longer pins twenty.
   - New: `scripts/draft_curated.py`, `pipeline/bake_curated.py`,
     `pipeline/structure/check_bundle.py`, `pipeline/sync_fixtures.py` (in `pipeline/`,
     not `scripts/`, to be importable by its tests).
   - Proof, from a scratch copy of the stored tracks: `bake_folding.py --all` then
     `bake.py --all`, `fetch_tracks.py --check` 143 of 143 the same; `check_folding` 0
     problems; `check_assets` 20 check out; the twenty's scenes compiled from the re-baked
     models twenty of twenty the same as served. `sync_fixtures --check` from local and
     from live: the twenty's 144 fixtures equal the app's. pytest 850 passed offline; 539
     passed with PyMOL in a throwaway 3.12 venv.
   - Order corrected: the seed must come before the uploads (`protein_track` references
     `protein`).
2. **Verify and write the ten.** Check each gene against its record as the
   second ten were, recorded as a new section of
   `helix-peek/docs/protein-verification.md`: span, exons, neighbouring genes,
   processing (which regions are kept, `mature_peptides`), structure. Add rows at
   orders 20–29, then the prose per R4.5. PAPPA and TG join the compressed-intron
   genes: update the app CLAUDE.md section that names DMD, APP and CFTR.

   **What actually happened (2026-10-17).** `draft_curated.py` drafted nine on their
   planned entries. UMOD's 6TQL fails the app's fold check (a step of 1.16, where 1.1 is
   allowed); the user chose policy 1's next entry, 9NU1 (the ZP domain, 329–580). PAPPA
   draws 8HGG's chain C (chain D narrowly fails the ring check). Verified against the
   records in `protein-verification.md` ("The third ten"). By the rules: PF4's and UMOD's
   records name nested products, so `mature_peptides=False`; PAPPA's unannotated
   propeptide (23–81) is the N-terminal extension (R3.6), named "Propeptide" in the table,
   with proprotein "pro-PAPP-A"; UMOD's 615–640 is the GPI-anchor signal. RefSeqGenes
   unversioned without transcript ids, as the twenty's; slices keep MANE's span. Words
   hand-written per R4.5 and clear of the search tests' pinned words.
3. **Bake** all ten. ESM-2 takes longest: TG about 35 min, PAPPA about 15.
   Then: `check_assets.py` over every curated row, pytest green, the twenty's
   stored sha256 untouched.

   **What actually happened (2026-10-17/18).** `bake_curated.py` for each. Two red gates,
   each the user's call: PLN fails ESM-2's alignment gate (40.0% / 46.7%) and ships
   unscored (`scored=False`); PF4 fails AVI's biology check (exons 7.4 = introns 7.4) and
   was replaced by PPBP (1F9P, CTAP-III), which passed everything. ESM-2: TG 35 min,
   PAPPA 20. Every model equals its draft's and passed the app's check again. Gates:
   `check_assets.py` 30 check out (29 scored), `check_folding`/`check_locus` 30 clean,
   `fetch_tracks.py --check` the twenty's 143 files as stored, pytest 866 offline + 542
   with PyMOL, the ESM freeze unchanged. Backend `237c86f`. Then the app's tests ran on
   the ten from local fixtures before anything was live (the stop for the user's review
   of the words, approved as written, with the further test and zoom changes below).
4. **Live.** Upload each kind with `--target` for the new slugs only, dry run
   first. Then `seed_catalog.py` (dry run, real, `--check` clean). Then
   `LIVE_BACKEND=https://helix-peak-backend.onrender.com flutter test test/live_backend_check.dart`
   (curl `/health` first; Render sleeps).

   **What actually happened (2026-10-18).** The order is seed first: `protein_track`
   references `protein`. Seed `--target` the ten (65 aliases, 60 track rows), then every
   kind uploaded (record, structure with scenes compiled, folding, constraint 9, impact,
   ClinVar, locus), nothing refused. `seed_catalog.py --check` 30 rows, 0 differences;
   `check_assets.py --against` 30 rows agree; `/catalog` 30 in order;
   `/proteins/resolve/PLN` ready at `phospholamban`; `fetch_tracks.py --check` 212 the
   same; `live_backend_check.dart` passes (after it, too, reads ESM-2 only where scored).
5. **App.** Sync the ten's fixtures. `flutter analyze` and `flutter test` must
   pass: the walk tests now run on the ten, and a failure is a row to fix,
   never a test. Then the three guards in the app CLAUDE.md. Golden pages are
   made only on x86, so list the new PNGs for the user to generate.

   **What actually happened (2026-10-18).** App `372cb54..a94426a`, not pushed. Fixtures
   from storage (`sync_fixtures.py`: the local pre-flight's 70 files the same). The tests
   met four things the twenty never had, each approved by the user: an unscored curated
   protein (five walk tests read ESM-2 only where scored), a mature count over a thousand
   ("(1,546)", one walk-test regex), a protein in the nucleoli (CSTB; three zoom tests),
   and the skin and eye reached for the first time (their drawings had no outline;
   `organ_art.dart` gives them theirs). Plus the two walk edits approved up front and the
   zoom tables. `flutter test` 1,748 passed; the three guards hold (the walk diff is the
   approved edits). 44 new golden PNGs owed to the x86 host.
6. **Report** to the user, and update the memory notes on unpushed commits and
   this proposal (`ten-more-curated-proposal`).

   **What actually happened (2026-10-18).** Reported. The thirty reach all 18 of the
   zoom's cell kinds (the twenty reached 13). Memory notes updated:
   `ten-more-curated-proposal`, `local-main-ahead-of-origin` (app 56 ahead, backend 11,
   GitHub unmoved), `live-db-migration-state`, `zoom-cell-coverage`,
   `goldens-need-x86-host`, `protein-pipeline-goal`. Owed to the user: the 44 golden PNGs
   on the x86 host, and pushing when they choose (the backend's commits change no service
   code; the ten are live already).

## Environment

- Every bake environment and command is in `pipeline/README.md`.
- `NCBI_EMAIL=gentiannurce@gmail.com` goes in the command's environment only.
- `ALPHAGENOME_API_KEY` is in `~/.env`; the Atlas meters per minute, so chunk.
- Never install into `pipeline/.esm-venv`.
- The resolver worker runs this backend tree against production. Read its
  queue first, and `scripts/resolver_worker.sh restart` after changing `pipeline/`.

**Out of scope:** the zoom draws the new organs (heart, placenta, esophagus,
skin, retina, thyroid, kidney) with its plain organ surface, since only the
twenty's organs have studies. All five new cell drawings (cardiomyocyte, glial,
trophoblast, megakaryocyte, spindle) are the older plain ones; only endocrine
and neuron cells are drawn part by part. Studies and cell recipes for them are zoom
work the user asks for by name.
