# Session plan: everything the Mac worker does, on Modal

**Goal.** A protein built on demand gets everything it gets today (record, fold, AVI, ClinVar and ESM-2 scores) from Modal, with the Mac turned off. The Mac stays the machine for curated data and re-bake proofs, not a server.

**Written against** backend `main` at `e415000` and app `main` at `f8ab13c` (2026-10-23). Re-check both before starting: if `pipeline/resolver/` or `test/fold_bundle_check.dart` moved since, re-read §1.

**Scope.** Backend only: `pipeline/resolver/modal_app.py`, one read-only helper in `store.py`, tests, docs. No baker, no `worker.py` bake function, no `local_worker.py` class and no app code changes. The app's checkout is only read, as the Mac's worker reads it.

**Before the session:**
- **A payment method and a spend limit** in the Modal workspace (Usage & Billing). Without them Modal refuses the GPU function, which is what stopped the deploy on 2026-10-07.
- **Run the session on the Mac.** The deploy has to carry two things that exist only there: the app's checkout and the AlphaGenome skill.

---

## 1. What the analysis found

### What Modal does today

`modal_app.py` is unchanged since Phase 6 and has never been deployed.
- `sweep` calls `worker.sweep(conn, storage)` with `structures` and `evidence` off. Modal-built proteins would get no fold, AVI or ClinVar.
- `score` runs ESM-2 in-process on an L4 (`score_in_process`, CUDA, the reader's Stop read from the database).
- `wake` is the endpoint the service calls.

### What the Mac's worker does (`local_worker.py`), one protein at a time

| Step | How the Mac runs it | What it needs |
|---|---|---|
| Resolve | `worker.sweep(structures=True, evidence=True)` | Biopython, psycopg. Already in `resolver_image` |
| Fold | `Modeller` → `model_local.py` as a subprocess in `pipeline/structure/venv` | PyMOL **3.1.0** (`experimental.PYMOL`, version-checked); numpy 2.5.3, scipy 1.18.1, trimesh 5.1.0; Dart and Flutter; **the app's checkout**, whose `pubspec.lock` must pin flutter_scene **0.23.0** (`alphafold.FLUTTER_SCENE`, checked by `importer_unready`); `flutter test --no-pub test/fold_bundle_check.dart` (headless: reads the compiled scene's bytes, no GPU); RCSB, PDBe, UniProt and AlphaFold DB over the network |
| AVI | `Evidencer` → `impact_local.py` as a subprocess in `pipeline/impact/venv` | the `alphagenome` client (`requirements.txt` says only `>=0.1.0`); `uv`; the AlphaGenome skill (`ALPHAGENOME_SKILLS`, default `<workspace>/.claude/skills`, **outside both repos**); `ALPHAGENOME_API_KEY`; about 2.5 GB of memory for the GENCODE lookup; `Settings.impact_state` kept between tries (checkpoints, GENCODE answers) |
| ClinVar | `ClinVarBaker` → `worker.clinvar_in_process` | Biopython only |
| ESM-2 | `Scorer` → `score_local.py` in `.esm-venv` (`mps`) | On Modal: `score_in_process` in `scorer_image` (`cuda`). Already done |

**Order** is enforced by the claim, not the caller (`store.WAITS_FOR`): ClinVar waits for AVI, and ESM-2 waits for both. A score container started too early simply finds nothing it can claim.

### The key finding: the Mac's classes can be reused unchanged

`Modeller`, `Evidencer` and `ClinVarBaker` are plain subprocess wrappers configured entirely by `Settings`:
- paths: `structure_python`, `dart`, `flutter`, `app`, `impact_python`, `uv`;
- `alphagenome_key`;
- `state` and the directories derived from it.

Nothing in them is macOS-specific. The macOS parts are `Awake` (`caffeinate`), launchd and the home-directory defaults in `load_settings`, and none of those are needed. Each class also has an `unready()` check that runs its driver's `--check`, so no bake is claimed while one could not be made.

So Modal builds a `Settings` that points at its image's paths and calls the same classes. **Same code path as the Mac, so the same behaviour by construction.**

### Other facts that shape the design

- **A Linux Flutter recipe already exists:** the app's `tool/goldens/Dockerfile` (Ubuntu 24.04, Flutter 3.47.5, `flutter precache --linux`).
- **The scene compiler works on Linux.** In the cloud sandbox, `dart run flutter_scene:import` ran headless on Linux and gave byte-identical output on repeat runs.
- **PyMOL on x86-64 does not give the Mac's bytes** (`pipeline/structure/README.md`, "Re-baking on Windows"). That matters for nothing stored: no built fold is re-checked against a digest.
- **Today's Modal images upload only `.py` files** (`not_code`). The worker path has never run on Modal, so any non-`.py` file it opens at runtime would be missing there. The obvious readers of `curated/catalog.json` (`seed_catalog.curated_rows`, `check_assets.curated_catalog`) don't look like they're on the worker path. Prove that by running it (§3, steps 5–6), not by reading.
- **Per-protein timeouts already sit under the reaper's limits:**

  | step | per-protein timeout | reaper's stale limit |
  |---|---|---|
  | fold | 15 min | 20 min |
  | AVI | 50 min | 60 min |
  | scoring | — | 90 min |

  Modal's function timeouts must allow several proteins per call, and must never hold one past its stale limit.

---

## 2. Design

### 2.1 Two new images in `modal_app.py`

**`fold_image`:**
1. Ubuntu x86-64 with Flutter 3.47.5, following `tool/goldens/Dockerfile` (`Image.from_dockerfile` on a copy, or the same `run_commands`).
2. micromamba from conda-forge with `python=3.12` and `pymol-open-source=3.1.0`, then pip-installed numpy 2.5.3, scipy 1.18.1, trimesh 5.1.0, plus the worker's packages (`biopython`, `psycopg[binary]`, `python-dotenv`).
3. **The app's checkout** from `paths.CLIENT` (`HELIXPEEK_CLIENT`, default `<workspace>/helix-peek`):
   - copied with `copy=True` into `/opt/helix-peek`;
   - leaving out `.git`, `build/`, `.dart_tool/`, the platform folders and `test/goldens/**/*.png`;
   - then `flutter pub get` there at image build.
4. **Deploy-time guards** in `modal_app.py`, so a mismatched app fails the deploy rather than every protein's fold:
   - the checkout exists;
   - its `pubspec.lock` pins `alphafold.FLUTTER_SCENE`;
   - `test/fold_bundle_check.dart` exists.

**`evidence_image`:**
1. Debian slim, Python 3.12.
2. `alphagenome==<the version in the Mac's pipeline/impact/venv>`, written into a new lock file. Also `python-dotenv`, numpy and the worker's packages.
3. **`uv`**, pinned to the Mac's version.
4. **The AlphaGenome skill** copied from `paths.SKILLS` into `/opt/skills`, with `ALPHAGENOME_SKILLS=/opt/skills`.
5. **Warmed at image build:** run one GENCODE lookup (`uv run … gtf --gene INS --exons --format json`) so the skill's own dependencies and its GENCODE data are inside the image, not downloaded at every cold start. Find on the Mac where the skill keeps that data, and copy it in if a build-time run can't fetch it.
6. Deploy-time guard: the skill folder exists.

Both images keep `add_local_python_source("app", "pipeline", ignore=not_code)`.

### 2.2 `Settings` for Modal: `modal_settings(state: Path) -> local_worker.Settings`

A pure function, tested:
- credentials from the environment (the secret);
- `structure_python` / `impact_python`: the image's interpreter;
- `dart`: the SDK's own binary, `/opt/flutter/bin/cache/dart-sdk/bin/dart` (not the `flutter/bin/dart` wrapper, for the reason `default_dart` gives);
- `flutter`: `/opt/flutter/bin/flutter`;
- `app`: `/opt/helix-peek`;
- `uv`: its absolute path;
- `structures` and `evidence` true;
- `state=/tmp/worker`, with `/tmp/worker/impact` a symlink into a **Modal Volume**, `helix-peak-resolver-state`. AVI checkpoints and GENCODE answers then survive between containers, as `impact_state` does on the Mac. Commit the volume after each call.
- `log_file`: unused, under `/tmp`.

### 2.3 Functions

| Function | Image | Resources | Does |
|---|---|---|---|
| `sweep` (every 5 min, and via `wake`) | `resolver_image` | CPU | `worker.sweep(conn, storage, structures=S, evidence=E)`, then spawns per §2.4 |
| `fold` (new) | `fold_image` | CPU, 4 GB, 60 min, `max_containers=1` | `Modeller(settings, stop).unready()`: if not ready, print why and claim nothing; else `worker.structure_all(conn, storage, model=modeller)`. Then spawn `score` if a constraint bake is claimable |
| `evidence` (new) | `evidence_image` | CPU, 6 GB, 60 min, `max_containers=1` | `Evidencer.unready()` as above; `worker.impact_all(avi=Evidencer)`, then `worker.clinvar_all(clinvar=ClinVarBaker)`. Then spawn `score` if a constraint bake is claimable |
| `score` | `scorer_image` | L4 | Unchanged |
| `wake` | `resolver_image` | | Unchanged |
| `preflight` (new, `modal run` only) | each image | | Runs each class's `unready()` and the scorer's import in its own image, and prints ready, or why not |

- **`stop`** is a `threading.Event` set on SIGTERM, so a preempted container hands its bake back as the Mac's does. A reader's Stop of a scoring is already read from the database.
- **Secrets.** `helix-peak-resolver` stays as it is, on every function. A new **`helix-peak-alphagenome`** secret holds only `ALPHAGENOME_API_KEY` and is attached to `evidence` alone. The Mac likewise hands the key to the AVI bake alone.
- **Switches.** `S` and `E` come from `RESOLVER_STRUCTURES` / `RESOLVER_EVIDENCE` in the `helix-peak-resolver` secret: the names the Mac's `.env` uses, read at run time, so they flip without a redeploy. **On Modal both default to off** until preflight passes, because a switch on with nothing able to bake leaves tracks `pending`.

### 2.4 What `sweep` starts

A pure function, `spawns(counts) -> set[str]`, tested offline:
- `fold` if any structure bake is queued;
- `evidence` if any `impact` or `clinvar` bake is queued;
- `score` if a constraint bake is **claimable now**.

That needs **`store.claimable_bakes(conn, kind) -> int`**: the count of what `claim_bake` would take, with the same `WAITS_FOR` and `resolver_version` conditions, read-only. Without it, a GPU container would start for a scoring still waiting on its evidence. With it, the GPU starts only when it has work, and `fold`/`evidence` start it when they finish.

### 2.5 What does not change

- **The bakers**, `worker.py`'s bake functions, the store's claim semantics and the service.
- **The app.**
- **`local_worker.py`:** its classes are imported, not edited. If one needs a change to run on Linux, that is a finding: stop, say what, and make it the smallest default-preserving change, with a test.
- **The Mac's worker** still runs as it does. Both can run at once (every claim is `FOR UPDATE SKIP LOCKED`); day to day, run one.

---

## 3. Steps, in order

0. **On the Mac, record versions for the locks:**
   - `pipeline/impact/venv/bin/pip freeze | grep -i alphagenome`
   - `uv --version`
   - the skill folder's version or date, and where it keeps GENCODE data
   - `pipeline/structure/venv/bin/pip freeze`
   - the PyMOL the Mac reports
1. **Images and deploy-time guards** in `modal_app.py`, plus `pipeline/resolver/requirements-modal-*.txt` locks where pip pins are needed.
2. **Code:** `modal_settings`, the `fold`/`evidence`/`preflight` functions, `spawns()` and the `sweep` change, and `store.claimable_bakes`.
3. **Tests (offline):**
   - `test_modal_app.py`: the six functions exist; the AlphaGenome secret is on `evidence` only (keep each function's options in a module-level dict so this is testable); the guards raise on a wrong pin or a missing skill; `not_code` still sends only code.
   - `spawns()` truth table; `modal_settings` paths and the volume symlink.
   - `store.claimable_bakes` in `test_worker_pg.py`: waiting, claimable, running, a curated row ignored.
   - `.venv/bin/pytest` green. Postgres tests green with `RESOLVER_TEST_DATABASE_URL`.
4. **Deploy with both switches off.** `modal secret create helix-peak-alphagenome ALPHAGENOME_API_KEY=…`, then `modal deploy pipeline/resolver/modal_app.py`. Note the first build time and image sizes.
5. **Preflight:** `modal run pipeline/resolver/modal_app.py::preflight` must print ready for fold, evidence and the scorer. Fix the images until it does.
6. **Switch over, one protein at a time:**
   1. `scripts/resolver_worker.sh stop` on the Mac.
   2. Set `RESOLVER_STRUCTURES=1` and `RESOLVER_EVIDENCE=1` in the secret.
   3. Build, each followed by `.venv/bin/python scripts/check_built.py <GENE> cuda`:
      - a small protein with an experimental entry (a fold from the PDB, with its folding track);
      - one that falls to AlphaFold DB's model;
      - one whose fold is refused (over 2,700 residues);
      - a minus-strand gene, for AVI and ClinVar;
      - a long protein, for GPU time and the evidence-then-score order.
7. **Parity against the Mac**, on three proteins the Mac built (one PDB fold, one AlphaFold fold, one with ClinVar records):
   - `scripts/queue_structures.py --refresh <GENE>` and `scripts/queue_evidence.py <GENE>`, made by Modal.
   - **Must be identical:** the entry chosen and the search's `selection`; the row's `chrome` words and `chains`; the scene's node names and vertex and face counts.
   - **Small differences are expected:** coordinates may differ by about 1e-5 (PyMOL on x86-64), which is a finding to note, not a failure.
   - **The evidence must match the Mac's**, apart from the retrieval date: the AVI track (its `runs` and scores) and the ClinVar track's variants.
   - Anything else is a parity bug: roll back with `scripts/fold_rollback.py`, and fix before going on.
8. **Backfill:** `queue_structures.py --dry-run --all` and `queue_evidence.py --dry-run`, then for real, for any protein the Mac never gave a fold or evidence.
9. **Retire the Mac's agent:** `scripts/resolver_worker.sh uninstall`. Set the wake on Render if not yet: `MODAL_WAKE_URL`, `MODAL_KEY`, `MODAL_SECRET`.
10. **Docs:**
    - resolver README: rewrite "Moving to Modal" as done; the Mac section becomes "the Mac as a fallback worker";
    - `modal_app.py` docstring;
    - backend `CLAUDE.md` ("is not": the worker runs on Modal);
    - a line in the app's `CLAUDE.md` that Modal's fold image carries the app's checkout, so **a change to the flutter_scene pin or to `test/fold_bundle_check.dart` means redeploying Modal**.

## 4. Done means

- pytest is green offline, and with Postgres.
- Preflight is ready on Modal.
- The five builds in step 6 pass `check_built.py … cuda`.
- Parity in step 7 holds.
- The Mac's agent is uninstalled, and a build made with it uninstalled has its fold, AVI, ClinVar and scores.
- **No change** under `pipeline/structure`, `pipeline/impact`, `pipeline/clinvar`, `pipeline/constraint`, or to the bake functions in `pipeline/resolver/worker.py` (`git diff --stat`).
- `seed_catalog.py --check` is clean, and `/catalog` is unchanged.

## 5. Risks to watch

- **Images.** Flutter, PyMOL and the GENCODE data make big images with slow first builds and slower cold starts for `fold` and `evidence`. Keep `max_containers=1`, and don't add `min_containers` unless waits hurt.
- **Deploy coupling.** The fold image carries the app's checkout and the evidence image the skill. A deploy from a machine with a different app commit changes what folds are checked against. The guards catch the flutter_scene pin, not the check's content, so deploy from the commit you released.
- **The AlphaGenome skill** may download at run time what the warm-up didn't cache. Watch the first `evidence` call's log, and keep its data on the Volume if needed.
- **`.py`-only uploads** may miss a data file the worker reads. Preflight and the first builds would show it as a broken bake. Fix with an explicit `add_local_file` for that file, never by dropping `not_code`.
- **Cost.** L4 at about $0.80 an hour, plus CPU minutes for folds and evidence. `RESOLVES_PER_DAY` and the spend limit bound it.
- **AlphaGenome quota** is shared with anything else using the key. A spent quota shows as retried AVI bakes.

## 6. Not in this session

- Making the curated proteins' re-bake proofs run on Modal. They compare against bytes your Mac made, so they stay on the Mac.
- Any change to how folds, AVI or ClinVar are made.
