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
one `resolve_request` row and calls a URL. Everything else is the worker's,
with the uploader's credentials: on Modal, as drawn, or until it is deployed
there on a Mac, which makes the same two calls from a launchd agent
([below](#running-the-worker-on-a-mac-temporary)).

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

  The gate was written to catch scores filed a position off, and each of the
  twenty clears it, the lowest at 64.9%. It also falls short where ESM-2 knows
  little about a protein, and it cannot tell the two apart. Of 34 buildable
  proteins of 150 residues or fewer, resolved live and scored on the dev Mac on
  2026-10-07, nine fell short: six of the eight of 40 residues or fewer, but
  also AKAIN1 (69 residues, 38.7%) and FAM24B (94, 43.8%), while RPL41, 25
  residues, cleared it at 95.2%. So the reason says what was measured, how
  often the model prefers the residue that is there to its neighbour's, and
  does not say the scores were misfiled (`worker._refusal`).

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

## Running the worker on a Mac (temporary)

Modal creates no GPU function for a workspace without a payment method: on
2026-10-07 step 3 above was refused at its last line, with the table and the
secret already in place. Until it is deployed, a Mac is the worker.
`local_worker.py` makes the two calls `modal_app.py` makes, `worker.sweep` and
`worker.score_all`, from a launchd agent. It asks the queue every 60 seconds
and scores on the Mac's GPU (`mps`) in `pipeline/.esm-venv`, the environment
that baked the twenty.

Once, from the repository root:

```sh
/opt/homebrew/bin/python3.12 -m venv .worker-venv
.worker-venv/bin/pip install -r requirements.txt
scripts/resolver_worker.sh install
```

It needs `.env` with the secret's four keys, `pipeline/.esm-venv`, and the 650M
weights in the Hugging Face cache, where a bake of the twenty leaves them.
Nothing is installed into `.esm-venv`: the worker has an environment of its
own and runs the scorer there as a subprocess (`score_local.py`).

| `scripts/resolver_worker.sh` | |
|---|---|
| `status` | whether it runs, what waits in the queue, its last ten lines |
| `logs` | follow its log |
| `stop` | stop it, and keep it stopped across logins |
| `start` | let it run again |
| `restart` | it runs the working tree as it was when it started, so: after any change under `pipeline/`. Within a minute of its last start this takes up to a minute (launchd's `ThrottleInterval`) |
| `install`, `uninstall` | write the agent and start it; stop it and remove the agent |

- **The agent** is `~/Library/LaunchAgents/com.helixpeek.resolver-worker.plist`.
  It starts at login and is started again if it crashes. It holds no
  credential; the worker reads `.env`.
- **The log** is `~/Library/Logs/HelixPeek/resolver-worker.log`: a line a
  protein, and a line a minute while there is nothing to do. 1 MB, and five
  older files. `resolver-worker.launchd.log` beside it has only what was said
  before the log was open.
- **Its files** are in `~/Library/Application Support/HelixPeek/resolver-worker/`
  and are emptied after each build, as a container's `/tmp` is. It never
  writes under `pipeline/data/`.
- **Overrides** go in `.env`: `RESOLVER_POLL_SECONDS` (60),
  `RESOLVER_SCORE_TIMEOUT` (3600; at most 5100, because a sweep takes a bake
  that has run 90 minutes for a dead worker's) and `RESOLVER_ESM_PYTHON`.
- **A protein it built** is checked by `.venv/bin/python scripts/check_built.py
  <GENE> mps`: the service, the stored bytes against their sha256, the rows and
  UniProt's features, ending `ALL OK`. `scripts/live_resolved_check.dart` reads
  it through the app; its header says how to run it.

**It works while the Mac is awake and its user is logged in.** A request made
while it sleeps waits on the queue: the app says "Building…" and, after 15
minutes, that the build is taking longer than it should. The worker keeps the
Mac from idle sleep while it has work (`caffeinate -i`) and at no other time,
so whether a build *starts* unattended is the Mac's own setting. `pmset -g
custom` on the dev Mac, 2026-10-08: never on the power adapter, after 1 minute
on battery (System Settings → Battery → Options, "Prevent automatic sleeping
on power adapter when the display is off"). Closing the lid sleeps a MacBook
whatever that says, unless an external display is attached.

A laptop is stopped, loses its network and has its files moved under it, and a
container is not. What the worker does about each is in `local_worker.py`. What
an operator sees of it:

- A stop or a restart lands between two proteins. A scorer stopped part-way
  has its bake back on the queue at once, with one of its three tries used.
- A protein that fails is tried again a poll later, not at once.
- No bake is claimed while the scorer cannot start (a checkout without
  `score_local.py`, the weights cleared from the cache). The log says "bake(s)
  left on the queue", and why, each minute until it can.
- macOS guards `~/Desktop`, where the repositories are. The agent's Python
  was not refused them on 2026-10-08 (macOS 26.6). Should it ever be,
  `resolver-worker.launchd.log` says `Operation not permitted` each minute:
  give the Python the environments link to (`readlink -f
  .worker-venv/bin/python`) Full Disk Access in System Settings → Privacy &
  Security, or move the repositories.

Measured on the dev Mac (M5, 24 GB), start of the scorer's process to its
track: sarcolipin, 31 residues, 5 s; beta-2-microglobulin, 119, 12 s; ANTXR1,
564, 203 s, the first protein the agent built unattended (2026-10-08). Scored
again through the worker, B2M's stored record gives its stored track byte for
byte.

Those three grow with the square of the length. Carried on, which is an
estimate and not a measurement, 1,000 residues is about 11 minutes, each
residue past ESM-2's 1,022-residue window adds about 0.6 s, and the 60-minute
default is reached near 5,000 residues. A protein that does not finish in time
is tried three times, an hour each, and its track then reads "Scoring failed 3
times: ESM-2 was still scoring after 60 minutes and was stopped."

## Moving to Modal

Nothing in the code changes: Modal's `sweep` and `score` call what the Mac's
worker calls.

1. In the Modal workspace (Usage & Billing), add the payment method and set a
   spend limit, a monthly cap past which Modal stops workloads. An L4 is $0.80
   an hour.
2. If `.env` has changed since the secret was made (2026-10-07):
   `.modal-venv/bin/modal secret create helix-peak-resolver --from-dotenv .env --force`
3. `.modal-venv/bin/modal deploy pipeline/resolver/modal_app.py`
4. `scripts/resolver_worker.sh uninstall`
5. Optionally the wake, step 4 of "Deploying, once".

The two may run at once during the switch: every claim is `FOR UPDATE SKIP
LOCKED`, so they never take the same request or the same bake. Day to day, run
one.

A track scored after the switch names `cuda` as its `generation.device` where
the Mac's name `mps`. The app does not read it, and nothing is scored again.
`.venv/bin/python scripts/check_built.py <GENE> cuda` checks the first protein
Modal builds.

Back again: `.modal-venv/bin/modal app stop helix-peak-resolver`, then
`scripts/resolver_worker.sh install`.

## Operating it

```sql
-- what is waiting, running, refused or failed
select gene, state, reason, attempts, requested_at from resolve_request
order by requested_at desc limit 20;

-- the constraint bakes
select slug, state, attempts, error from bake_job order by requested_at desc limit 20;

-- let a newer resolver try a refused gene again
delete from resolve_request where lower(gene) = lower('TTN') and state = 'refused';

-- score again a protein whose scoring broke three times ("Scoring failed 3
-- times: ..."), once what broke it is mended
update protein_track set state = 'pending', reason = null, updated_at = now()
where slug = 'ttn' and kind = 'constraint' and state = 'refused';
insert into bake_job (slug, kind) values ('ttn', 'constraint');
```

On a Mac, `scripts/resolver_worker.sh status` prints the counts by state, the
last five requests and every bake still waiting, with its last error.

`.venv/bin/python scripts/dry_resolve.py GENE ...` shows what the resolver would
build for a gene, and whether the upload gate takes its record, writing
nothing. It does not run ESM-2.

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

`test_local_worker.py` is the Mac's worker with no torch and no database: a
script stands in for the scorer's environment and two lists for the queue.
`test_worker_pg.py` runs its cycle on the scratch database as well: a request
taken to a scored protein, a stop mid-score, a scorer that cannot start.
