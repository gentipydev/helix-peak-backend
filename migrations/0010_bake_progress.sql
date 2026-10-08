-- How far a running bake has got, so a reader waiting on a protein built on
-- demand sees the residues scored and the time left rather than a spinner.
--
-- Run by hand, once, against the Supabase project, after 0001 to 0009:
--
--     psql "$DATABASE_URL" -f migrations/0010_bake_progress.sql
--
-- The scorer already says it as it goes ("Scored 450/838 (283.3s, 4.1 min
-- left)"); the worker reads those lines and writes the last one here
-- (`store.note_progress`). `GET /proteins/resolve/{gene}` serves it.
--
-- Progress belongs to the run that wrote it: it counts only while
-- `progress_at` is not before the job's `started_at`. A job claimed again, or
-- put back on the queue, needs nothing cleared, and a worker that predates
-- this file runs as it did, writing nothing here.
--
-- No `protein` column and no existing row, object or provenance is touched.
-- Everything here is additive: four nullable columns on `bake_job`.

begin;

alter table bake_job
    -- Units done and in all: residues scored, for a constraint bake.
    add column progress_done  integer,
    add column progress_total integer,
    -- Seconds the baker estimated were left when it wrote this.
    add column progress_left  real,
    add column progress_at    timestamptz;

commit;
