-- A reader may stop a build they asked for: one that has run long enough
-- (HERC2, 4,834 residues, an hour of ESM-2 on the dev Mac) to be worth
-- giving up.
--
-- Run by hand, once, against the Supabase project, after 0001 to 0010:
--
--     psql "$DATABASE_URL" -f migrations/0011_build_stop.sql
--
-- There are no accounts, so who asked is an install's own random id, sent with
-- the ask (`asker`). Only that id may stop it. A request keeps the id; the
-- constraint bake its resolution queues is given the same one, as is a bake
-- queued again when the reader asks for a stopped protein to be scored.
-- Requests and bakes from before this file have none, and no reader can stop
-- them.
--
-- `stopped` is a new state for both. A request stopped while queued built
-- nothing, and asking again starts over. A bake stopped while queued or
-- running leaves the protein built without its ESM-2 track (`absent`, nothing
-- is coming), and asking again queues the scoring afresh. The worker notices a
-- stop the next time it writes its progress, and ends the scorer.
--
-- No `protein` column and no existing row, object or provenance is touched.
-- Everything here is additive: two nullable columns, and each check dropped
-- and re-added with one more value, in one transaction, so no moment exists
-- without it.

begin;

alter table resolve_request add column asker text;
alter table bake_job add column asker text;

alter table resolve_request drop constraint resolve_request_state_known;
alter table resolve_request add constraint resolve_request_state_known check (state in (
    'queued', 'running', 'done', 'refused', 'failed', 'stopped'
));

alter table bake_job drop constraint bake_job_state_known;
alter table bake_job add constraint bake_job_state_known check (state in (
    'queued', 'running', 'done', 'failed', 'stopped'
));

commit;
