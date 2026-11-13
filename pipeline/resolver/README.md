# The resolver: any protein, on demand (Phase 6)

A reader searches for a protein the catalog does not list. `/proteins/suggest`
already says whether it is `buildable`; asking for it (`POST /proteins/resolve`)
queues a request, and this directory turns the request into a protein the walk
can open. Since 2026-11-13 the index offers only proteins expected to build
whole, their fold from the PDB (`app/protein_index.py`, `withhold`); every
other is `unavailable`, with the sentence why:

```
app (search field) ──POST /proteins/resolve──▶ service ──row──▶ resolve_request
                                                  │
                                                  └──wake──▶ Modal: sweep (CPU, every 5 min anyway)
                                                               │  UniProt entry + GenBank record
                                                               │  → protein row, aliases, record track (ready)
                                                               │  → structure, impact, clinvar, constraint
                                                               │    tracks (pending), a bake_job each
                                                               ├──spawn──▶ Modal: fold (CPU)
                                                               │             PDB entry, else AlphaFold DB
                                                               │             → structure (+ folding) track (ready)
                                                               ├──spawn──▶ Modal: evidence (CPU)
                                                               │             AlphaGenome → impact, NCBI → clinvar
                                                               │             └──spawn──▶ score, once both are over
                                                               └──spawn──▶ Modal: score (one GPU: T4, else L4, else A10)
                                                                             ESM-2 650M → constraint track (ready)
app ──GET /proteins/resolve/{gene}──▶ state + build ──▶ /gene/<slug>, the walk, as for the twenty
```

The service never imports this directory and holds no storage key: it writes
one `resolve_request` row and calls a URL. Everything else is the worker's,
with the uploader's credentials: on Modal, as drawn, or on a Mac, which makes
the same calls from a launchd agent
([below](#running-the-worker-on-a-mac-temporary)). The two may run side by side.

Each makes a protein's fold as well as its scores, from an experimental entry
chosen by rule, with the fold animation the twenty play, and from AlphaFold
DB's model only where a whole search finds none
([below](#the-fold-an-experimental-entry-else-alphafold-dbs-model)), and its
variant evidence, AlphaGenome's AVI and ClinVar's records
([below](#the-variant-evidence-avi-and-clinvar)). Until 2026-11-13 only the
Mac's could; Modal's images now hold the same tools
([Moving to Modal](#moving-to-modal)).

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
- **Tracks**: `record` ready; `locus`, the zoom's way in, ready in the same
  transaction, made from the index row and the pinned snapshot of the locus
  bake's inputs (`worker.make_locus`, `pipeline/locus/README.md`) or refused
  with the reader's sentence for the seven genes the Atlas cannot place, and
  absent only where the worker could not read the snapshot (its outcome says
  why); `constraint` pending until the GPU has scored it,
  where the worker makes models `structure` pending until its model is made
  (seconds), and where it makes the variant evidence `impact` and `clinvar`
  pending until AlphaGenome's scores and ClinVar's records are in (minutes).
  Nothing else; where a worker is set to make none of the three
  (`RESOLVER_STRUCTURES=0`, `RESOLVER_EVIDENCE=0`), they read `absent`.
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

## What a reader waiting sees

ESM-2 takes minutes (838 residues, 9 on the dev Mac), so `GET
/proteins/resolve/{gene}` says how far a build has got as well as what it is
(`build`, `BuildReport` in `app/schemas.py`): the step (`queued`, `record`,
`evidence`, `scoring`, `check`, `structure`, `done`), the proteins the worker
takes first, the residues scored and of how many, the seconds the scorer
estimates are left, the seconds since it was asked for, the transcript the
record was built from, and its fold's state (`fold`: `pending`, `ready`,
`refused`, or none where no model was asked for). A step that ended without
its result carries the track's reason.

`structure` holds a build whose scores are over while its first fold is still
being made: a fold is made before the scoring, in seconds to a minute, so this
is rare, and it lasts `resolves.FOLD_HOLD` (20 minutes) at most. A bake still
queued after that is one no worker is making (models switched off, or the
modeller unable to start), and the build opens without it, its fold page saying
a model is on its way. A protein that draws a model already is never held
while another is made for it (a backfill), and a stopped build is never held.
An app older than the step reads the report as none, and opens the protein once
its ESM-2 track is settled, as it did before.

The scorer prints a line every 25 residues ("Scored 450/838 (283.3s, 4.1 min
left)"). `worker.score_next` hands the scorer a `Report`, and the line reaches
it either way: on the Mac `local_worker.Scorer` reads the scorer's output every
5 s, and on Modal `worker.score_in_process` reads what it prints in-process.
`store.note_progress` writes it on the bake (`progress_*`, `0010`), where it
counts only for the run that wrote it. A write that fails costs the bake
nothing. The scorer itself is unchanged.

`/proteins/suggest` marks a protein `building` while a request is queued or
running, its constraint bake is not over, or its build is held for its first
fold. It is `ready` from the moment its
row is written, and a walk opened before ESM-2 is done would keep the pending
track until the app restarts, so the app watches it instead.

The reader who asked may stop a build (`POST /proteins/resolve/{gene}/stop`,
`0011`). There are no accounts: the app sends its install's own random id with
the ask (`asker`), the request keeps it, the bake its resolution queues is
given it, and only that id may stop either. Stopped while queued, the request
is `stopped` and nothing was built. Stopped while ESM-2 scores, the bake is
`stopped` and the track `absent`; the worker hears it the next time it writes
its progress, or, on the Mac, within 5 s while the model loads, and ends the
scorer (`worker.Stopped`). The protein opens without scores, suggest marks it
`stopped`, and asking for it again queues its scoring afresh, for whoever
asked. A record being written (seconds) cannot be stopped, and a stop that
lands as the track is made is too late for it.

Once stopped, the same reader decides (`0014`; the app's card offers Resume,
Keep and Discard, and comes back after a restart until one is chosen):

- **Resume** is asking again, as above.
- **Keep** (`POST /proteins/resolve/{gene}/keep`) keeps the protein as it is,
  without its ESM-2 track: `kept_at` on the stopped bake. It still opens
  unscored, and suggest still marks it `stopped`.
- **Discard** (`POST /proteins/resolve/{gene}/discard`) deletes everything
  built for it: the service deletes its rows at once (the protein, with its
  aliases, tracks, bakes and folds, and every request for its gene), refusing
  while a part of it is still being made, and writes a `discard_request` row.
  The service holds no storage key, so the next sweep deletes its objects
  (`forget.purge_discarded`, first thing in `worker.sweep`, on Modal and on the
  Mac): those named as its own (`<folder>/<slug>.<sha12>.<ext>`) that no row
  names and that were made before the discard, so a build of the same gene
  after it keeps every object. A build stopped before anything was built has
  only its request, which goes.
- `POST /proteins/stopped` (`asker` in the body, never in a URL) lists one
  install's stopped builds with no decision yet: what its app pins again.

Only the install that stopped a build may keep or discard it, and only a
stopped build is either: a curated protein, a protein built whole, or one kept
already answers 409 with the reason.

## The fold: an experimental entry, else AlphaFold DB's model

A protein's fold page draws a model, and the twenty's are cut from entries
chosen by hand, each with the fold animation that ends on it. A protein
resolved on demand has its entry chosen by rule (policy 1,
`pipeline/structure/pdb_policy.py`; the structure README has the rule, what it
chose over 52 proteins, and the bake), and folds on its page as the twenty do.
Only where a whole, successful search finds no entry that passes is it drawn
from AlphaFold DB's model (`pipeline/structure/alphafold.py`), still, coloured
by pLDDT, as before. A search cut short, or a service not answering, is never
read as none.

- **Queued with the scoring.** `store.write_resolution(..., structures=True)`
  puts the `structure` track `pending` behind a `bake_job` of its own, in the
  transaction that queues the `constraint` one, for the same `asker`. Only a
  worker that can make a model asks (`worker.sweep(structures=True)`): the
  Mac's and Modal's both do unless switched off, so a track is never left
  pending where nothing will bake it.
- **Made before the scoring.** `worker.structure_next` claims it: the Mac's
  cycle runs it between resolving and scoring, and on Modal the `fold`
  function runs it, beside the evidence, as soon as the sweep has resolved
  the protein. `model_local.py` (in the
  structure bake's environment) searches the PDB, ranks every way of drawing
  the protein, holds the best to its own file, and bakes it into the model and
  its `folding` track from one working copy, by the twenty's own bakers; then
  it holds the pair to the app's own check (`test/fold_bundle_check.dart`, run
  from the app's checkout), which binds the fold to the model as the walk's
  fold page will. A candidate that fails a documented rule gives way to the
  next; a whole search with none left falls to AlphaFold DB's model. Seconds
  to a minute (the dry run's 43: 5.6 s to 70 s, mTOR's 2,549 residues the
  longest), so the fold page is ready long before a long protein's ESM-2
  track is.
- **Remembered between attempts** (`fold_selection`, 0012). A search that found
  nothing for the same inputs under the same rule is not made again for 28
  days (four PDB releases): AlphaFold DB's model is made straight away. A
  search cut short at its bounds (eight entry files, three bakes an attempt)
  keeps how far it got, goes back on the queue without spending a try, and
  the next attempt goes on from there where the search finds the same
  entries. An experimental fold that stands for the protein's record, regions
  and rule is kept until a search is asked to start afresh
  (`queue_structures.py --refresh`).
- **Stored as the twenty's are.** The compiled `.fsceneb` in the `models`
  bucket as `structure/<slug>.<sha12>.fsceneb`, the `.glb` beside it and named
  in the provenance (`upload_tracks.structure_row`, the uploader's own
  branch); an experimental fold's `folding` track in `tracks` as
  `folding/<slug>.<sha12>.json`, whose frame names that `.glb` by its digest;
  and the bake's audit beside it, `fold_audit/<slug>.<sha12>.json`, named in
  the structure track's provenance (the entry and its file's digest, the
  record's, every candidate and why it was passed over, the residue map, the
  tools). No object is ever deleted.
- **Published whole, or not at all** (`store.publish_fold`). The row's
  `structure` column (`{"chrome": ..., "chains": ...}`, the fold page's seven
  words templated from the model's own numbers, and a tint for each node), the
  structure track, the folding track (ready for an experimental fold, absent
  for AlphaFold's), the generation that records them (`fold_generation`), the
  search that chose them, and the bake done, in one transaction. It writes
  only while the bake is still this worker's at the attempt it claimed (one
  the reaper put back, or another worker took, writes nothing:
  `store.Superseded`), only while the protein's stored record, regions and
  disulfides are the ones the fold was made for, and never on one of the
  twenty.
- **Refused, with the sentence the fold page shows** (`worker.Unmodelled`):
  where no experimental entry passes and AlphaFold DB has no model either, say
  so, the PDB's part first ("The PDB holds no experimental structure of it.
  AlphaFold DB has no model of proteins over 2,700 residues; this one has
  4,834."). A protein that draws a model already keeps it. Anything else (a
  service down, PyMOL, the importer or the app's check breaking) is tried
  three times, a cycle apart; then a track with no model says "The fold could
  not be made this time: its bake broke each of the 3 times it was tried.",
  the job's `error` keeps what broke, and a model already drawn stays.
- **Not stopped.** A reader's Stop ends the ESM-2 scoring. The fold is made all
  the same: a protein stopped while scoring opens with its fold and without
  its scores.
- **On the build card** (`BuildReport.fold`, and the `structure` step above).

`scripts/check_built.py` holds a built protein's fold too: the scene and its
`.glb` to their digests, the nodes to the row's chains; for an experimental
fold its `folding` track, the `.glb` its frame names, its chains and bridges,
the audit beside it, and the generation the database records as drawn; for
AlphaFold's, the span, the gate and the bridges; or, refused, the sentence why.

## The variant evidence: AVI and ClinVar

The twenty's walk reads AlphaGenome's per-base Variant Impact scores (`impact`)
and ClinVar's records (`clinvar`); since 2026-10-09 a protein built on demand
gets both too (HANDOFF-AVI-CLINVAR.md has the user's decisions). They are made
by the bakers that made the twenty's, unchanged (`impact/bake_impact.py`,
`clinvar/bake_clinvar.py`), so nothing of the twenty's bytes moved.

- **Queued with the scoring.** `store.write_resolution(..., evidence=True)`
  puts both tracks `pending` behind a bake each (`store.queue_evidence`), for
  the same `asker`. Only a worker that can make them asks: the Mac's, and
  Modal's (its `evidence` function), each unless switched off.
- **In a build's order, held by the claim.** ClinVar's records are placed on
  the record's bases by the map AVI's track carries (`runs`, `chromosome`,
  `complemented`), so a `clinvar` bake is claimed only once its protein's
  `impact` bake is over, and a `constraint` bake only once both are
  (`store.WAITS_FOR`). The long step stays last, and a build that opens once
  its scores are in has its evidence too.
- **AVI**: `local_worker.Evidencer` runs `impact_local.py` in the AVI bake's
  own environment (`pipeline/impact/venv` on the Mac; on Modal the same
  packages, `pipeline/impact/requirements-lock.txt`), with
  `ALPHAGENOME_API_KEY` (from `.env`, or the Modal secret: handed to that bake
  alone) and `uv` on its PATH for the AlphaGenome skill's GENCODE lookup.
  GENCODE's answers and the pull's checkpoints are kept in the worker's state
  between tries, never in `pipeline/impact/` (on Modal, for as long as the
  container lives). One of the bake's own gates saying no (the exons do not
  pair with GENCODE's MANE Select ones, the record's bases are not GRCh38's
  where GENCODE puts the gene, too few were scored) refuses the track with the
  gate's words (`worker.Unplaced`); the Atlas or the lookup not answering is
  tried again, three times a cycle apart. A gene with no intron skips only the
  bake's biology alarm (exons must outscore intron interiors, which it cannot
  have), the user's choice of 2026-10-09; the sequence gate still holds every
  base, and the track's (unread) `generation` says null for the two intron
  medians.
- **ClinVar**, in the worker's own process (Biopython is all it needs):
  `worker.clinvar_in_process` runs the baker with its raw responses in the
  worker's folder, deleted once the bake is over; the track names the query,
  the day and each batch's sha256. NCBI not answering in full is tried again
  (`worker.Unfetched`), and the baker's own `ValueError` is a refusal. Where
  AVI was refused, so is ClinVar, its sentence ending in AVI's.
- **Through the upload gate.** `upload_tracks.validate` holds each to
  `check_assets.py`'s own checks against the record (ClinVar's against AVI's
  map) before a row says ready. The twenty's 40 stored payloads pass them.
- **On the build card**: the step is `evidence` while either bake is queued or
  running, and never carries a reason. No Stop is offered there: Stop ends
  ESM-2.
- **Checked** by `scripts/check_built.py`, as it stands: ready, each to its
  digest and to `check_assets.py` again; refused, with the sentence why.

Measured on the dev Mac (2026-10-09, `HANDOFF-AVI-CLINVAR.md` Phase 1): GENCODE's
lookup, 1.8 s a gene at about 2.5 GB; TTR's AVI, 6.9 s for one Atlas window;
its ClinVar, 17 s for 503 records in 6 batches, 28.6 MB of raw responses. Of
twenty genes tried (the eight built then and twelve more, both strands,
chromosome slices, compressed introns, chrX, chrY and a PAR gene) all twenty
mapped every drawn base.

## Deploying, once

From the repository root, in order.

1. **The tables.** Apply the migrations, as every other:

   ```sh
   .venv/bin/python scripts/apply_migration.py migrations/0009_resolve_request.sql
   .venv/bin/python scripts/apply_migration.py migrations/0010_bake_progress.sql
   .venv/bin/python scripts/apply_migration.py migrations/0011_build_stop.sql
   ```

2. **The secret.** In the Modal workspace, `helix-peak-resolver` with five
   keys: `DATABASE_URL`, `SUPABASE_URL`, `SUPABASE_SERVICE_KEY` (the
   uploader's three), `NCBI_EMAIL` and `ALPHAGENOME_API_KEY` (the Atlas's,
   which only the AVI bake is handed). From `.env`, those five lines alone,
   never the whole file and never printed:

   ```sh
   five=$(mktemp)    # readable by you alone
   grep -E '^(DATABASE_URL|SUPABASE_URL|SUPABASE_SERVICE_KEY|NCBI_EMAIL|ALPHAGENOME_API_KEY)=' .env > "$five"
   .modal-venv/bin/modal secret create helix-peak-resolver --from-dotenv "$five" --force
   rm "$five"
   ```

3. **The app.**

   ```sh
   python3 -m venv .modal-venv && .modal-venv/bin/pip install -r pipeline/resolver/requirements.txt
   .modal-venv/bin/modal token new          # once, in a browser
   .venv/bin/python -I pipeline/locus/snapshot.py   # the locus inputs, once
   .modal-venv/bin/modal deploy pipeline/resolver/modal_app.py
   ```

   The first deploy builds the four images: the scorer's downloads the 650M
   weights into it, and the fold's Flutter (1.6 GB) and PyMOL (a few minutes
   each). Run it from a checkout with the app's beside it (`../helix-peek`,
   clean) and the AlphaGenome skill (`../.claude/skills`), which the fold's and
   the evidence's images copy. The output prints the `wake` URL,
   `https://<workspace>--helix-peak-resolver-wake.<region>.modal.run`.

4. **The wake.** In the Modal workspace (Settings → Proxy Auth Tokens), create a
   token. On Render, set `MODAL_WAKE_URL` to the wake URL and `MODAL_KEY` /
   `MODAL_SECRET` to the token's key and secret. Optional: without them the
   five-minute schedule still takes every request, just later.

5. **The caps and the budget.** `RESOLVES_PER_DAY` (default 20) and
   `RESOLVES_PER_INSTALL_PER_DAY` (default 5) on Render bound a day's asks, a
   stopped protein's scoring asked again included, counted in `spend_ledger`
   (0015), which no discard deletes. What bounds the bill is the budget in
   `app/budget.py` ($2 a rolling day, $10 a calendar month): each of `score`,
   `fold` and `evidence` reserves its worst case (its timeout at Modal's rates,
   on the dearest GPU) before it claims anything and settles what it cost when
   it ends, and the sweep starts none whose worst case would not fit; its bake
   waits, queued. `SECURITY.md` at the repository's root has the numbers.

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

A Mac can be the worker too, beside Modal or instead of it. It was the only
one from 2026-10-08 to 2026-11-13, while the Modal workspace had no payment
method. `local_worker.py` makes the calls `modal_app.py` makes,
`worker.sweep`, the fold's, the evidence's and `worker.score_all`, from a
launchd agent. It asks the queue every 60 seconds
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

To make AVI tracks it needs the AVI bake's environment (`pipeline/impact/venv`,
its README), `ALPHAGENOME_API_KEY` in `.env`, the AlphaGenome skill under
`../.claude/skills/` and `uv` (`~/.local/bin/uv`, or `RESOLVER_UV`). It asks
first (`impact_local.py --check`: the client, the key, the skill under `uv`,
and one 1-base Atlas request), and claims no AVI bake while it could not.

To make models it also needs the structure bake's environment
(`pipeline/structure/venv`, and PyMOL 3.1.0: its README), and the app's
checkout beside this one with a Flutter SDK, for the scene importer and for the
app's own check of a fold. It runs the bake there as a subprocess too
(`model_local.py`), having first asked whether one could be made (`--check`:
the packages, PyMOL and its version, the importer on a box, Flutter, and the
app's `test/fold_bundle_check.dart`), and only once the database has 0012's
tables. While it could not, no model is claimed, and the log says why each
minute. The app's checkout is only read and run from: `dart run
flutter_scene:import`, with the flutter_scene its `pubspec.lock` pins, which
has to be the one every stored scene was compiled by, and `flutter test --no-pub
test/fold_bundle_check.dart`, which writes the test's own build files under its
`.dart_tool/` and `build/`, both ignored by git. Where someone's `flutter run` holds Flutter's start-up lock, the
check waits for it, up to ten minutes; a check that does not finish is no
verdict, and the bake is tried again.

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
  that has run 90 minutes for a dead worker's) and `RESOLVER_ESM_PYTHON`. For
  the models: `RESOLVER_STRUCTURES=0` makes none (no structure bake is queued
  or claimed, as on Modal), and `RESOLVER_STRUCTURE_PYTHON`, `RESOLVER_DART`
  (default: the SDK inside `~/flutter`, not `flutter/bin/dart`, which takes
  Flutter's start-up lock), `RESOLVER_FLUTTER` (`~/flutter/bin/flutter`) and
  `RESOLVER_APP_DIR` (`../helix-peek`) say where its four tools are. For the evidence: `RESOLVER_EVIDENCE=0` makes none (none
  queued, none claimed; any already queued waits, and so does its protein's
  scoring), and `RESOLVER_IMPACT_PYTHON` and `RESOLVER_UV` say where its two
  tools are.
- **A protein it built** is checked by `.venv/bin/python scripts/check_built.py
  <GENE> mps`: the service, the stored bytes against their sha256, the rows and
  UniProt's features, ending `ALL OK`. One whose ESM-2 track is not ready (the
  scorer refused it, or its scoring was stopped) is checked for its model and
  fold alone, `scripts/check_built.py <GENE> fold`. The app's
  `test/live_resolved_check.dart` reads it as the app does; its header says how
  to run it.

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
564, 203 s, the first protein the agent built unattended (2026-10-08); OCA2,
838, 534 s, asked for from the app the same day. Scored again through the
worker, B2M's stored record gives its stored track byte for byte.

A model, from the modeller's start to its stored objects: AlphaFold DB's,
about two seconds at any length (OCA2's 838 residues, 2.0 s; a refusal, under
one). An experimental one, search and all, measured over the dry run's 43 on
2026-10-12: 5.6 s to 70 s, nine seconds for most; the bake 0.7 s to 8.4 s and
the app's check about two seconds of that, after a first check that compiles
the app's test (a minute or so, once after each change to the app).

Those three grow with the square of the length. Carried on, which is an
estimate and not a measurement, 1,000 residues is about 11 minutes, each
residue past ESM-2's 1,022-residue window adds about 0.6 s, and the 60-minute
default is reached near 5,000 residues. A protein that does not finish in time
is tried three times, an hour each, and its track then reads "Scoring failed 3
times: ESM-2 was still scoring after 60 minutes and was stopped."

## Moving to Modal

Nothing in the code changes: Modal's functions call what the Mac's worker
calls, through the same drivers (`local_worker.Modeller`, `Evidencer`,
`ClinVarBaker`) and the same bakers. Since 2026-11-13 Modal makes everything
the Mac's worker makes, so a protein is built with the Mac off.

- **`fold`** runs in an image that holds the structure bake's environment
  (`pipeline/structure/requirements.txt`), PyMOL 3.1.0 from conda-forge in a
  conda environment of its own (on the path as `pymol`, so its numpy never
  meets the bake's), Flutter 3.47.3 (held to its framework and engine
  revisions as the image is built) and the part of the app's checkout the
  scene importer and the app's check read (`modal_app.APP_FILES`,
  `APP_DIRS`). Never the app's `.env`: an empty file stands in for the asset
  `pubspec.yaml` declares. The image's build fetches the app's packages as
  its lockfile pins them, compiles the check, and runs the build hooks and
  the importer once, so a failure shows there and not on a protein.
- **`evidence`** runs in an image that holds the AVI bake's environment,
  pinned whole from the one that baked the twenty's
  (`pipeline/impact/requirements-lock.txt`), `uv`, the AlphaGenome skill
  synced from its lock with its own Python 3.14, and GENCODE v46's table,
  fetched once (`ALPHAGENOME_GTF_PATH`) rather than for each protein. It
  makes one AVI track and then one ClinVar track a call, and starts `score`
  once a scoring no longer waits for them.
- Each claims nothing while its check says no bake could be made in its
  image (`model_local --check`, `impact_local --check`: one Atlas request),
  exactly as the Mac's worker. A broken image then costs one check a sweep
  while a bake waits, and the function's log says why.
- The sweep queues a fold and the evidence for each protein it resolves, and
  starts each function only where a bake of its kind waits that it can take
  (`store.claimable_bakes`): `fold` and `evidence` side by side, then `score`
  from `evidence`. `RESOLVER_STRUCTURES=0` or `RESOLVER_EVIDENCE=0` in the
  shell that runs `modal deploy` turns either off, read as the Mac's `.env`
  reads them; a plain deploy makes both.
- The fold image is the app's checkout as it was at the deploy, so deploy
  again after a change to `lib/shared/folding/`, `test/fold_bundle_check.dart`
  or the flutter_scene the app pins, from a clean app tree.
- **`sweep`** makes each protein's locus, the zoom's way in, as it resolves
  it (`worker.make_locus`), from the pinned snapshot of the locus bake's
  inputs (`pipeline/locus/inputs.py`), which its image carries
  (`modal_app.LOCUS_INPUTS`). The snapshot is data, which the code's upload
  leaves behind, so it goes up on its own, and a deploy needs it built first
  on the deploying machine (`.venv/bin/python -I pipeline/locus/snapshot.py`,
  which says whether what it built is the pinned one). Without it the deploy
  stops and names the file. Its digest is checked as a sweep loads it; a
  sweep that cannot read it resolves the protein all the same, its locus
  `absent`, which `scripts/queue_locus.py` makes afterwards. The load takes
  the sweep's process to about 164 MB, so it asks for 512 MiB
  (`modal_app.SWEEP_MEMORY`).
- A model made there is x86-64 PyMOL's, whose bytes are not this Mac's
  (`pipeline/structure/README.md`, "Re-baking on Windows"). That matters for
  nothing stored: no model made on demand is ever re-baked against a digest,
  and the app's check holds each model to its own scene.

1. In the Modal workspace (Usage & Billing), a payment method and a spend
   limit, a monthly cap past which Modal stops workloads. `score` takes a T4
   ($0.59 an hour), and when Modal has none to give an L4 ($0.80), then an
   A10 ($1.10) (`modal_app.SCORER_GPU`, in that order). The fold and the
   evidence run on CPUs, cents a month.
2. The secret's five keys, step 2 of "Deploying, once".
3. The locus inputs' snapshot, once on the deploying machine:
   `.venv/bin/python -I pipeline/locus/snapshot.py`.
4. `.modal-venv/bin/modal deploy pipeline/resolver/modal_app.py`
5. Optionally the wake, step 4 of "Deploying, once", and `RESOLVES_PER_DAY`.
6. The Mac's worker may stay beside it, or go: `scripts/resolver_worker.sh uninstall`.

The two may run at once: every claim is `FOR UPDATE SKIP LOCKED`, so they never
take the same request or the same bake, and either may make any part of a
protein the other resolved. Modal starts its GPU only for a constraint bake it
can take (`store.claimable_bakes`), never for one still waiting on its AVI and
ClinVar or for one of the twenty's, and works one protein a call, as the Mac
does: a request or bake put back on the queue waits for the next sweep rather
than spending its three tries at once. Two workers baking AVI at once share one
Atlas quota; the bake paces itself.

A track scored after the switch names `cuda` as its `generation.device` where
the Mac's name `mps`. The app does not read it, and nothing is scored again.
`.venv/bin/python scripts/check_built.py <GENE> cuda` checks a protein Modal
builds.

Back again: `.modal-venv/bin/modal app stop helix-peak-resolver`, then
`scripts/resolver_worker.sh install`.

## Operating it

```sql
-- what is waiting, running, refused or failed
select gene, state, reason, attempts, requested_at from resolve_request
order by requested_at desc limit 20;

-- the constraint bakes, and how far a running one has got
select slug, state, attempts, error, progress_done, progress_total, progress_at
from bake_job order by requested_at desc limit 20;

-- let a newer resolver try a refused gene again
delete from resolve_request where lower(gene) = lower('TTN') and state = 'refused';

-- score again a protein whose scoring broke three times ("Scoring failed 3
-- times: ..."), once what broke it is mended
update protein_track set state = 'pending', reason = null, updated_at = now()
where slug = 'ttn' and kind = 'constraint' and state = 'refused';
insert into bake_job (slug, kind) values ('ttn', 'constraint');

-- what each protein built on demand draws now, and how its last search ended
select g.slug, g.source, g.state, g.created_at, s.outcome, s.candidate, s.expires_at
from fold_generation g
left join lateral (select outcome, candidate, expires_at from fold_selection f
                   where f.slug = g.slug order by decided_at desc, id desc limit 1) s on true
where g.active order by g.created_at desc;
```

Make the fold of proteins built without one (before the worker made models, or
by Modal), or again (a backfill), without SQL by hand:

```sh
.venv/bin/python scripts/queue_structures.py --dry-run --all     # what it would queue; reads only
.venv/bin/python scripts/queue_structures.py TTR PRL             # or --all
.venv/bin/python scripts/queue_structures.py --refresh TTR       # search afresh, past what stands
```

One transaction a protein (`store.queue_structure`), one of the twenty refused.
A model it draws stays drawn until the one that replaces it is published, the
tracks and words together, and where nothing better can be drawn it keeps it;
its build is never held. An experimental fold that stands for the same record,
regions and rule is left as it is unless `--refresh`.

Every fold published is a generation (`fold_generation`), and so is the one a
protein drew before its first since 0012: an AlphaFold model, kept when the
next is published, or the sentence why it had none, kept when its fold is
queued again, before its track says a model is on its way. Put one back, whole:

```sh
.venv/bin/python scripts/fold_rollback.py TTR --list       # its folds, newest first
.venv/bin/python scripts/fold_rollback.py TTR              # the one before the fold drawn now
.venv/bin/python scripts/fold_rollback.py TTR --to 7       # that one, back or forward
```

One transaction (`store.roll_back_fold`): the row's words, the structure track
and the folding track as that generation wrote them (AlphaFold's has none, so
it reads absent again). Refused while a structure bake for it is queued or
running. An app reads it the next time the protein is opened.

AVI and ClinVar for a protein built without them (before 2026-10-09, or by
Modal), or again after a refusal, without SQL by hand:

```sh
.venv/bin/python scripts/queue_evidence.py --dry-run TTR PRL   # what it would queue; reads only
.venv/bin/python scripts/queue_evidence.py TTR PRL             # or --all
```

One transaction a protein, a ready track left as it is, and one of the twenty
refused. A protein opened while its evidence is pending keeps it missing until
the app restarts, so queue it while nobody walks it.

The locus for a protein built before the resolver made one (2026-11-13), or
again:

```sh
.venv/bin/python scripts/queue_locus.py --dry-run --all   # what it would write; reads only
.venv/bin/python scripts/queue_locus.py B2M PFN1          # or --all
```

- Nothing is queued: a locus is a pure function of its index row and the
  pinned snapshot, so each is made at once (`worker.make_locus`) and written
  in a transaction of its own, ready or refused.
- Only a protein built on demand is touched (`catalog_order` null).
- It needs the snapshot on this machine (`pipeline/locus/snapshot.py`).
- On 2026-11-13 it gave B2M and PFN1 theirs.
- `check_built.py` holds a protein's locus to the bytes the snapshot makes
  from its index row today, and calls an absent one wrong.

On a Mac, `scripts/resolver_worker.sh status` prints the counts by state, what
the proteins built on demand draw now and how their last searches ended
(`folds  4 pdb ready, 2 alphafold ready, 1 none refused; searches 4 chosen,
3 none`), the last five requests and every bake still waiting, with its last
error.

`.venv/bin/python scripts/dry_resolve.py GENE ...` shows what the resolver would
build for a gene, and whether the upload gate takes its record, writing
nothing. It does not run ESM-2.

Delete proteins built on demand, whole: their rows (the `protein` row, and by
cascade its aliases, tracks, bakes, fold selections and generations), every
request for their genes, and every object of theirs in storage, those their
rows name and any other a re-bake left under their names:

```sh
.venv/bin/python scripts/forget_built.py --all --dry-run    # what would go; reads only
.venv/bin/python scripts/forget_built.py TTR PRL            # or --all
```

A curated protein is never deleted: naming one refuses the lot. Refused while a
request or a bake for one of them is queued or running. The rows go in one
transaction, then the objects. A gene deleted can be asked for again, if the
index offers it. On 2026-11-13 it deleted all eleven built so far, at the
user's word.

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
client is not installed (`.modal-venv/bin/python -m pytest pipeline/resolver`
runs it). It runs each Modal function's body on fakes: what the sweep starts
for each kind of bake it can take, and only where switched on; the fold and
the evidence claiming nothing while their check says no; each one's own
restart, and the GPU started once a scoring no longer waits; every call ending
before the reaper's patience with its claim; and what goes up of the app's
checkout, never its `.env`.

`test_local_worker.py` is the Mac's worker with no torch and no database: a
script stands in for the scorer's environment and two lists for the queue.
`test_worker_pg.py` runs its cycle on the scratch database as well: a request
taken to a scored protein, a stop mid-score, a scorer that cannot start.

The fold is tested the same three ways: `test_worker.py` holds the worker's
gate on a model and on an experimental fold's track and audit,
`test_worker_pg.py` every statement from a queued model to a ready one, a
refused one and a reaped one, and the service serving it, and
`test_local_worker.py` the modeller with a script in the place of the structure
bake's interpreter, and the model driver with the network stood in for.
`test_fold_pg.py` holds publication: a bake taken from its worker and a
protein changed mid-bake write nothing, a publication interrupted leaves the
fold drawn before it whole, a curated row is untouched, a fold is put back and
forward, a search cut short or finding nothing is remembered, a fold that
stands is kept, a backfill keeps its model, the hold and its bound, and the two
operator commands. The bake itself is `pipeline/structure/test_alphafold.py`
and `test_experimental.py`.

So is the evidence, with `test_impact_local.py` for the AVI bake's driver: its
verdicts told from its failures, and the bake pointed at the worker's folders,
with no AlphaGenome client installed. `test_worker_pg.py` holds the order the
claim keeps, the cascade, the reaper's patience and the build card's step.
