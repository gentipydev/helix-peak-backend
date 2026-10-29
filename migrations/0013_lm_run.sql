-- A recorded training run of the de novo peptide model, which the app plays
-- stage by stage: where its bytes are, as a track's row says where a track is.
--
-- Run by hand, once, against the Supabase project, after 0001 to 0012:
--
--     .venv/bin/python scripts/apply_migration.py migrations/0013_lm_run.sql
--
-- A run is not one protein's: it is trained on the peptides of several
-- (insulin's, relaxin's and vasopressin's, `pipeline/denovo/curated_peptides.py`),
-- so it cannot be a `protein_track` row, whose slug references one protein.
-- This table is that row's shape for a run. `pipeline/denovo/publish_trace.py`
-- writes it, and `GET /peptide-lm/runs/{name}` serves it; the bytes stay in
-- storage and the client fetches them there.
--
-- No `protein` column and no existing row, object or provenance is touched.
-- Everything here is additive: one table, RLS on and no policies.

begin;

create table lm_run (
    -- Which run: `curated` is the curated proteins' peptides.
    name         text primary key,

    -- Where the recorded run is, as `protein_track` says where a track is:
    -- the path carries the digest, so new bytes are a new object.
    bucket       text not null,
    object_path  text not null,
    bytes        integer not null,
    sha256       text not null,
    format       text not null default 'json',

    -- What made it: the model, the steps and checkpoints, the engine and the
    -- versions, the records the peptides were read from, the commit and date.
    provenance   jsonb not null default '{}'::jsonb,

    updated_at   timestamptz not null default now()
);

alter table lm_run enable row level security;

commit;
