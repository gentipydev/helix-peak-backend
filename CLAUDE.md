# helix-peak-backend

FastAPI service behind the Helix Peek Flutter app. It lifts one gene out of an NCBI
GenBank record, serves the protein catalog and each protein's track index from Supabase,
and suggests from every reviewed human protein. `pipeline/` holds the offline bake tools.

## THE CONTRACT — read before changing anything

Baker CODE may be shared and refactored. The DATA SHAPE may not change:

- No new column on the `protein` table.
- No new field in `pipeline/targets.py` or `pipeline/curated/catalog.json`.
- No change to an existing track's bytes, format or provenance.
- No edit to an applied migration (`0001`, `0002`). A change is the next numbered file.
- New data is a new track kind: its own `protein_track` rows, objects and provenance.

The data shape holds while all of these hold: re-baked tracks match the stored sha256;
`pipeline/seed_catalog.py --check` diffs clean; `/catalog` returns exactly twenty.

### The proof that makes a baker refactor safe: re-bake, compare sha256

A stored track's digest is the `sha256` on its row (`GET /protein/{slug}/tracks`).

1. `python pipeline/fetch_tracks.py`: storage's copy into `pipeline/data/`, digest-checked.
   Bakes read each other's output, so their inputs must be the stored bytes.
2. Re-bake with `--target <slug>` (commands in `pipeline/README.md`).
3. `python pipeline/check_assets.py`, then `... --against https://helix-peak-backend.onrender.com`.
   `check_records` holds each served record to its row's sha256 and then to the re-baked
   file byte for byte. `check_scenes` holds each `.fsceneb` to its row's sha256.
4. check_assets compares bytes only for `record`. For other kinds, compare the re-baked
   file's sha256 with that kind's `sha256` on `/protein/{slug}/tracks`. For `structure`,
   the row's digest is the compiled `.fsceneb`, so compare the `.glb` with `provenance.glb.sha256`.

Don't run `fetch_tracks.py` between re-bake and compare: it overwrites any file whose digest
differs from storage. A mismatch means the refactor changed data. Stop; never upload to fix it.

## What the service is, and is not

- Is: `/gene/{id}/{gene}` (read-through cache of whole records), `/catalog`, `/catalog/search`,
  `/protein/{slug}`, `/protein/{slug}/tracks`, `/proteins/suggest`, `/proteins/built`, `/proteins/resolve` (POST,
  GET `/{gene}`, and POST `/{gene}/stop`), `.../impact-explanations`, `/health`, `/health/db`.
  It writes `genbank_record`, and one `resolve_request` row per new ask. It names bytes and never carries
  them: a ready track resolves to a public Supabase storage URL the client fetches directly.
- Is not: authenticated (there is no auth), rate limited (Biopython's ~0.37 s spacing between
  Entrez calls is the only throttle; resolve requests alone are capped per day,
  `RESOLVES_PER_DAY`) or CORS-restricted (`allow_origins=["*"]`). It holds no
  service_role key; rows are reached through `DATABASE_URL` only. It resolves nothing itself:
  a request is a row the resolver (`pipeline/resolver/`) works, on Modal or, until it is deployed
  there, on the dev Mac. The service only wakes Modal's (`MODAL_WAKE_URL`, proxy auth); the Mac's
  asks the queue every minute. `/catalog/search`'s `candidates` stays empty.
- Only `NCBI_EMAIL` is required; with no `DATABASE_URL`, `/gene` still works and the catalog is 503.
- `app/` stays Python 3.9-compatible (`Optional`/`List`, never `X | Y`); the image runs 3.12.

## app/ module map

- `main.py`: wiring. Entrez email and tool, socket timeout, CORS, lifespan (pool plus `genbank_record`), health endpoints.
- `config.py`: pydantic-settings. `NCBI_EMAIL` has no default, so importing without it fails by design.
- `db.py`: optional sync psycopg pool, opened in the lifespan. `None` when `DATABASE_URL` is unset.
- `entrez_client.py`: the one `Entrez.efetch` call. Fetch-as-text and parse stay apart so the cache can store text.
- `record_cache.py`: read-through cache of whole GenBank records by accession. Best-effort and self-creating.
- `genbank_parser.py`: `extract_gene(record, gene)`. Pure, exact `/gene` match, 1-based inclusive.
- `catalog.py`: reads `protein` and `protein_alias`. Raises `CatalogUnavailable` (503) and never answers empty.
- `tracks.py`: `protein_track` rows become `{kind: Track}` with public storage URLs. Any kind with a row is served; `KINDS` is the floor.
- `impact_explanations.py`: AVI explanations. Storage redirect first, then the local directory.
- `protein_index.py`: pure. `normalize()`, UniProt/MANE parsing, index rows and terms. Shared with `scripts/load_protein_index.py`.
- `suggest.py`: `/proteins/suggest`. Ranked prefix tiers over `protein_index_term`, near misses last;
  `building` marks a protein whose build is under way, which the app watches rather than opens.
  And `/proteins/built`: every install's builds as the same suggestions, newest first, paged by
  (`resolved_at`, slug); one under way is left out, one stopped is kept.
- `resolves.py`: `/proteins/resolve`. What an ask means (ready, pending, refused, failed, unavailable,
  buildable), the `resolve_request` row, the daily cap and the best-effort wake, and how far a
  build has got (`build`: step, residues scored, seconds left), read from the request and the
  constraint bake's progress (`0010`), which the worker writes as the scorer prints it; and
  stopping a build for the install that asked (`asker`, `0011`).
- `schemas.py`: pydantic response models that mirror the Dart entities field for field.
- `router.py`: every non-health endpoint, with the NCBI 404/502 and catalog 503 mappings.

The app's lab was deleted on 2026-10-09, and with it every baker and endpoint only the lab
read: `pipeline/audio`, `structure_ar`, `trafficking`, `assemblies`, and `/assembly/{slug}/tracks`.
Their tables and kinds stay in the schema (`0006`'s `assembly` and `assembly_track`, and
`protein_track_kind_known`), empty: a migration here only ever adds.

## Migrations: why one table creates itself

`genbank_record` is created by `record_cache.create_schema()` in the lifespan
(`create table if not exists`): one table, no history to migrate, a deploy that cannot
half-apply. A failure is logged and swallowed, and every read then misses.

Every other change is a `migrations/000N_*.sql`, applied by hand once and in order:
`.venv/bin/python scripts/apply_migration.py migrations/000N_x.sql` (or `psql -f`). Run from a
web worker on every boot, multi-object DDL could half-apply, run concurrently with itself
across workers, and fail where the only honest response is to keep serving.
Each file is one transaction, additive, with RLS on and no policies.

## The web service must not import pipeline/

The Dockerfile copies only `requirements.txt` and `app/`, so `import pipeline` in `app/`
breaks the deployed service at import. Pipeline code also needs bake environments (ESM,
AlphaGenome, PyMOL) and write credentials. `config.impact_explanations_dir` defaults to a
path under `pipeline/data/`; that is the only link, and it is a path, not an import.

`pipeline/resolver/` (Phase 6) is the one part of `pipeline/` that runs unattended: on Modal,
with the uploader's credentials in a Modal secret, calling the record builder and the ESM-2
scorer unchanged. A protein it resolves is a `protein` row with a null `catalog_order` and
`resolver_version` 1, so it never joins `/catalog`. It may never write over a curated row
(the upsert is guarded), and a change to a baker it calls needs the same sha256 proof:
`structure/bake.py` is one of them now (its `export` and `assemble` make every model).

Until Modal has a payment method, the worker is a launchd agent on the dev Mac:
`local_worker.py`, the same `worker.sweep` and `worker.score_all`, with `.env` for the secret
(`scripts/resolver_worker.sh status|stop|start|restart`; the runbook is its README). Three
things follow for anyone working in this checkout:

- It runs this working tree against production. `restart` it after a change under `pipeline/`,
  and `stop` it before leaving the tree somewhere it cannot run from (a long rebase, an old
  branch). While its scorer cannot start it claims no bake, so nothing is refused meanwhile.
- Nothing is installed into `pipeline/.esm-venv`, which baked the twenty. The worker lives in
  `.worker-venv` and runs the scorer there as a subprocess (`score_local.py`), so what that
  imports stays `scoring.py`, `score_protein.py`, `targets.py` and `paths.py`.
- `local_worker.py` imports nothing under `pipeline` at module level: `paths.DATA` is fixed at
  first import, and the worker's is its own directory, never `pipeline/data/`.
- It also makes each protein's fold, from AlphaFold DB's model (`pipeline/structure/alphafold.py`),
  between resolving and scoring: `model_local.py` in `pipeline/structure/venv`, with PyMOL and the
  scene importer of the app's checkout (`../helix-peek`, read and run from, never written). The
  importer is the flutter_scene that checkout's `pubspec.lock` pins, which must be the one every
  stored scene was compiled by (`alphafold.FLUTTER_SCENE`); otherwise no model is claimed.
  `RESOLVER_STRUCTURES=0` in `.env` turns it off. It writes the row's `structure` column, the one
  column a bake writes on `protein`, and only where `catalog_order is null`.

## Checklist: adding a track kind

1. A baker in `pipeline/<kind>/`, with README, requirements and tests. Which proteins get the
   kind is recorded there, never as a field in `targets.py` or `catalog.json`. If it
   shares code with an existing baker, re-bake that baker's tracks and prove the sha256 unchanged.
2. `migrations/000N_<kind>.sql` drops and re-adds the `protein_track_kind_known` check with
   the new value. Without it, the insert is rejected. Never edit `0001`.
3. No `app/` change. `/protein/{slug}/tracks` and the catalog's `tracks` maps serve any kind
   that has a row, as the row says (pinned by `test_catalog.py`). `KINDS` is only the walk's
   floor, reported `absent` when there is no row. Leave it alone.
4. `pipeline/upload_tracks.py`: `--kind` choices, `asset_of`, `validate` (which returns the
   provenance), and the bucket, suffix and content type if it is not JSON. The upload
   upserts the row as `ready`. A protein with no row reads `absent`.
5. `pipeline/fetch_tracks.py`: `KINDS` and `asset_path`.
6. A verification script for the new payload. Checks on existing kinds do not change.
7. Done: pytest green, `seed_catalog.py --check` clean, `/catalog` still twenty, existing
   sha256 unchanged. Upload only when asked, with `--dry-run` first.
8. The client skips kinds it cannot classify (`TrackKind.fromWire`), so no installed app breaks.

## Tests run offline

`.venv/bin/pytest` (on Windows, `.venv\Scripts\python -m pytest`) runs `tests/` and `pipeline/`
(`pytest.ini`). pytest must never hit the network: not NCBI, not Supabase, not Render.

- `conftest.py` sets `NCBI_EMAIL` before importing the app.
- Every Entrez call is patched: `efetch_returning(text)` / `efetch_raising(exc)` replace
  `entrez_client.Entrez.efetch`, fed `tests/fixtures/ng_007114.gb` (a real NCBI response) or
  `tests/ncbi_errors.py`. `efetch_raising(AssertionError(...))` proves a path never calls NCBI.
- The database is a `FakePool` set on `db.pool` that answers by SQL shape (`test_catalog.py`,
  `test_record_cache.py`, `test_suggest.py`, `test_built.py`). `no_pool` sets it to `None`.
- `client` is `TestClient(app)` outside a `with`, so the lifespan never runs and the real
  `DATABASE_URL` in `.env` is never dialled. Never write `with TestClient(app)`.
- Tests that need stored tracks skip until `fetch_tracks.py`, run by hand, fills `pipeline/data/`.
- `pipeline/resolver/test_worker_pg.py` runs the resolver's, `/proteins/resolve`'s and
  `/proteins/built`'s SQL against a real Postgres, and skips unless `RESOLVER_TEST_DATABASE_URL`
  names one it may create a scratch database on (a local server, never Supabase). `test_modal_app.py` skips without the Modal client.

## NCBI error mapping

| Upstream | We return |
|---|---|
| 200 with a GenBank record | 200 (404 if no feature carries that `/gene`) |
| 200 with error text (`F a i l e d  t o  u n d e r s t a n d`) | 404, because `SeqIO.read` raises `ValueError` |
| HTTP 400 | 404, because `id` is the only input that can be invalid |
| Any other HTTP status | 502, with the upstream body in `detail` |
| DNS failure, refused connection, timeout (`URLError` / `OSError`) | 502 |

503 means our own source is unreadable (catalog, protein index, explanations directory), never
answered as empty or 404: that would say a protein does not exist when nobody could look.

A valid response is recognised positively: it parses as GenBank, or it is not a record.
NCBI answers a bad id with plain text and HTTP 200, so an empty-body check would pass that
text through as a 200. The cache parses before it writes, and evicts a cached row that
fails to parse. `_upstream_errors` keeps its clauses most-specific-first, because
`HTTPError` < `URLError` < `OSError`.

Git: commit locally and never push, in any form. The user does all pushing.
