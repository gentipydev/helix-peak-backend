# The resolver: any protein, on demand (Phase 6)

A reader searches for a protein the catalog does not list. `/proteins/suggest`
already says whether it is `buildable`; asking for it (`POST /proteins/resolve`)
queues a request, and this directory turns the request into a protein the walk
can open:

```
app (search field) ──POST /proteins/resolve──▶ service ──row──▶ resolve_request
                                                  │
                                                  └──wake──▶ Modal: sweep (CPU, every 5 min anyway)
                                                               │  UniProt entry + GenBank record
                                                               │  → protein row, aliases, record track (ready)
                                                               │  → constraint track (pending) + bake_job
                                                               └──spawn──▶ Modal: score (one L4)
                                                                             ESM-2 650M → constraint track (ready)
app ──GET /proteins/resolve/{gene}──▶ ready ──▶ /gene/<slug>, the walk, as for the twenty
```

The service never imports this directory and holds no storage key: it writes
one `resolve_request` row and calls a URL. Everything else runs on Modal, with
the uploader's credentials.

## What a protein resolved on demand is (resolver version 1)

`resolve.py` turns a buildable `protein_index` entry and its UniProt entry into
the `Target` row the twenty were written by hand as, and hands it to the bakers
that already exist, unchanged: `mock/build_gene_record.build` for the record and
`constraint/score_protein.score_protein` for ESM-2. Its docstring has every
rule; in short:

- **Where the gene is**: the entry's RefSeqGene, or the span MANE puts its
  transcript on, with the MANE Select transcript and protein named so a record
  holding several transcripts draws MANE's.
- **Regions**: UniProt's processing features. Initiator methionine, signal and
  transit peptides and propeptides are removed; chains and peptides are kept.
  Domains are not read (they overlap; regions tile). Overlapping processing
  features drop the peptides inside chains, then fall back to leader plus one
  chain. Short names follow `targets.py`: `S`, `T`, `Met`, the letter a piece's
  name gives it ("Insulin B chain" is `B`), else numbers.
- **Disulfides**: UniProt's pairs whose residues are cysteines in the record.
- **No mature-peptide page** (`mature_peptides` False): the walk stops at the
  precursor.
- **Alleles**: the record's protein may differ from UniProt's at up to 1% of
  residues (at least three); more, or a different length, is refused.
- **Prose** is templated from the index ("Made by BRCA1 on chromosome 17. Built
  on demand from UniProt P38398 and MANE Select NM_007294.4, not curated by
  hand."), and `provenance.prose` says `templated`.
- **Tracks**: `record` ready, `constraint` pending until the GPU has scored it,
  and nothing else (`impact`, `clinvar`, `structure` read `absent`).
- **ESM-2** is the twenty's model and revision, `facebook/esm2_t33_650M_UR50D`,
  run through the same scorer and its alignment gate. The app accepts only that
  model's tracks, so no app change is needed for them. A protein that fails the
  gate has its constraint track `refused`, with the reason.

On insulin, which is one of the twenty, the automatic row reproduces the
hand-written one: the same cut, short names, cleavage sites, numbering and
disulfides (`test_resolve.py`).

## Deploying, once

From the repository root, in order.

1. **The table.** Apply the migration, as every other:

   ```sh
   .venv/bin/python scripts/apply_migration.py migrations/0009_resolve_request.sql
   ```

2. **The secret.** In the Modal workspace (modal.com → Secrets → Custom), create
   `helix-peak-resolver` with four keys: `DATABASE_URL`, `SUPABASE_URL`,
   `SUPABASE_SERVICE_KEY` (the uploader's three) and `NCBI_EMAIL`.

3. **The app.**

   ```sh
   python3 -m venv .modal-venv && .modal-venv/bin/pip install -r pipeline/resolver/requirements.txt
   .modal-venv/bin/modal token new          # once, in a browser
   .modal-venv/bin/modal deploy pipeline/resolver/modal_app.py
   ```

   The first deploy builds the scorer's image, which downloads the 650M weights
   into it (a few minutes). The output prints the `wake` URL,
   `https://<workspace>--helix-peak-resolver-wake.modal.run`.

4. **The wake.** In the Modal workspace (Settings → Proxy Auth Tokens), create a
   token. On Render, set `MODAL_WAKE_URL` to the wake URL and `MODAL_KEY` /
   `MODAL_SECRET` to the token's key and secret. Optional: without them the
   five-minute schedule still takes every request, just later.

5. **The cap.** `RESOLVES_PER_DAY` (default 50) on Render bounds a day's
   requests, which is what bounds the GPU bill.

### The first build

```sh
curl -s -X POST https://helix-peak-backend.onrender.com/proteins/resolve \
     -H 'content-type: application/json' -d '{"gene": "B2M"}'     # 202, pending
curl -s https://helix-peak-backend.onrender.com/proteins/resolve/B2M  # pending … ready
curl -s https://helix-peak-backend.onrender.com/protein/b2m/tracks   # record ready, constraint pending … ready
```

`modal app logs helix-peak-resolver` shows each sweep's summary and each bake.

### The scorer's one change, and its proof

`constraint/score_protein.py` now picks CUDA where there is one. The machine
that baked the twenty has none, so its path is unchanged; to prove it, on that
Mac:

```sh
pipeline/.esm-venv/bin/python -u pipeline/constraint/score_protein.py --target insulin
shasum -a 256 pipeline/data/assets/constraint/insulin_esm_constraint.json
```

and compare with the `constraint` row's `sha256` on `/protein/insulin/tracks`
(after `fetch_tracks.py`, and without running it again in between).

## Operating it

```sql
-- what is waiting, running, refused or failed
select gene, state, reason, attempts, requested_at from resolve_request
order by requested_at desc limit 20;

-- the constraint bakes
select slug, state, attempts, error from bake_job order by requested_at desc limit 20;

-- let a newer resolver try a refused gene again
delete from resolve_request where lower(gene) = lower('TTN') and state = 'refused';
```

A worker that dies holding a request (or a bake) is noticed by the next sweep
after 30 minutes (90 for a bake) and put back, or failed once it has been tried
three times.

## Tests

```sh
.venv/bin/pytest pipeline/resolver          # offline: the resolver, the worker's
                                            # decisions, and the Modal app's shape
RESOLVER_TEST_DATABASE_URL=postgresql://postgres@127.0.0.1:5432/postgres \
    .venv/bin/pytest pipeline/resolver/test_worker_pg.py
```

The second runs every statement `store.py` and `app/resolves.py` make against a
scratch database with all migrations applied, from request to ready track, and
the service reading what was written. It needs a Postgres it may create a
database on, and skips without one. `test_modal_app.py` skips where the Modal
client is not installed.
