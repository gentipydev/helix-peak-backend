-- What the reader decides once a build they stopped is stopped: Resume, Keep
-- or Discard (the user's design of 2026-11-13). The app shows the decision on
-- the stopped build's card, and the card comes back after the app restarts
-- until one is made, so the decision is the service's to remember.
--
-- Run by hand, once, against the Supabase project, after 0001 to 0013:
--
--     .venv/bin/python scripts/apply_migration.py migrations/0014_build_decisions.sql
--
-- Resume needs nothing new: asking again queues the scoring afresh (0011), and
-- a stopped bake is then no longer the protein's latest.
--
-- Keep is said on the stopped constraint bake: `kept_at`, when the reader
-- chose to keep the protein as it is, without its ESM-2 track. A stopped bake
-- with no `kept_at` is a decision still to make (`POST /proteins/stopped`).
--
-- Discard deletes the protein's rows at once (the service holds the database
-- URL alone), and leaves its objects in storage to the resolver, which holds
-- the storage key: each `discard_request` row is one protein whose objects are
-- to go. The resolver deletes only objects named as that protein's that no row
-- names and that were made before `requested_at`, so a protein built again
-- after its discard keeps every object of its new build. `purged_at` is when
-- that was done, and `objects` how many went. A build stopped before anything
-- was built has nothing in storage: its row has no slug and is purged at once.
--
-- No `protein` column and no existing row, object or provenance is touched.
-- Everything here is additive: one nullable column, and one table with RLS on
-- and no policies.

begin;

alter table bake_job add column kept_at timestamptz;

create table discard_request (
    id            bigserial primary key,

    -- The gene as the app asked for it, and the slug its protein had: null
    -- where nothing was built.
    gene          text not null,
    slug          text,

    -- The install that discarded it: the one that asked for the build.
    asker         text not null,

    requested_at  timestamptz not null default now(),

    -- When the resolver deleted its objects, and how many.
    purged_at     timestamptz,
    objects       integer
);

create index discard_request_unpurged on discard_request (requested_at)
    where purged_at is null;

alter table discard_request enable row level security;

commit;
