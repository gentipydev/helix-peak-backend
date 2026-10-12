-- A protein built on demand draws an experimental structure where one passes
-- the rule (policy 1, `pipeline/structure/pdb_policy.py`), with the fold
-- animation the twenty play, and AlphaFold DB's model only where a whole,
-- successful search finds none (HANDOFF-PDB-FOLDING.md). What each search
-- decided, and each fold that was published, is kept here: so that the next
-- search can go on from where one stopped, a search that found nothing is not
-- made again for a while, every outcome can be counted, and a fold can be put
-- back as it was.
--
-- Run by hand, once, against the Supabase project, after 0001 to 0011:
--
--     .venv/bin/python scripts/apply_migration.py migrations/0012_fold_generation.sql
--
-- No `protein` column and no existing row, object or provenance is touched.
-- A fold is still the `structure` track, the row's `structure` column and a
-- `folding` track, as the twenty's are; these tables record them, and no
-- reader reads these tables. Everything here is additive: one column on
-- `bake_job` with a default, two tables, their indexes, RLS on and no
-- policies, the posture every table in 0001 takes.

begin;

-- A structure bake asked to search again, past what an earlier one decided
-- and past the fold the protein draws (`scripts/queue_structures.py --refresh`).
alter table bake_job add column refresh boolean not null default false;

create table fold_selection (
    id            bigserial primary key,
    slug          text not null references protein (slug) on delete cascade,
    job_id        bigint references bake_job (id) on delete set null,

    -- The rule it was decided by, and what it was decided for: a digest of the
    -- protein's letters, products, domains and disulfides
    -- (`fold_report.inputs_digest`), and the stored record they were read from.
    policy        integer not null,
    inputs        text not null,
    record_sha256 text not null,

    -- chosen      an experimental entry passed, baked and checked: it is drawn
    -- none        the whole search ran, and nothing passed
    -- incomplete  one attempt reached its bounds; `refused` is how far it got,
    --             and the next attempt goes on from there
    outcome       text not null,

    -- The search: its query, RCSB's weekly release, how many entities it
    -- found, and the digest of their ids.
    search        jsonb,
    -- The candidate drawn, as `fold_job.candidate_id` names it, and its entry.
    candidate     text,
    entry         text,
    -- Candidates turned down so far, by id: what a resumed search skips.
    refused       jsonb not null default '[]'::jsonb,
    -- The whole record of the search (`fold_report.selection`): every
    -- rejection by its reason, each candidate examined, and why the winner won.
    detail        jsonb not null default '{}'::jsonb,
    decided_at    timestamptz not null default now(),
    -- When a search that found nothing is to be made again: the PDB grows
    -- every week.
    expires_at    timestamptz,

    constraint fold_selection_outcome_known check (outcome in ('chosen', 'none', 'incomplete')),
    constraint fold_selection_chosen_names_one check (
        outcome <> 'chosen' or (candidate is not null and entry is not null)
    ),
    constraint fold_selection_none_expires check (outcome <> 'none' or expires_at is not null)
);

create index fold_selection_slug_idx on fold_selection (slug, decided_at desc, id desc);

create table fold_generation (
    id              bigserial primary key,
    slug            text not null references protein (slug) on delete cascade,

    -- pdb        an experimental entry, with its folding track
    -- alphafold  AlphaFold DB's model, drawn still
    -- none       no model: the structure track says why
    source          text not null,
    -- The structure track as published: ready, or refused.
    state           text not null,
    job_id          bigint references bake_job (id) on delete set null,
    selection_id    bigint references fold_selection (id) on delete set null,

    -- Exactly what was written, so that it can be written again: the protein
    -- row's `structure` column, the structure track's row and the folding
    -- track's (null: none, read as absent). No object is ever deleted, so
    -- every generation's objects stay where its rows point.
    structure       jsonb,
    structure_track jsonb not null,
    folding_track   jsonb,
    -- Where the bake's audit is: its bucket, path, sha256 and size.
    audit           jsonb,

    -- What it was made from: the stored record, the protein's regions and
    -- disulfides (an md5 the database takes of them), and the rule.
    record_sha256   text,
    regions         text,
    policy          integer,
    created_at      timestamptz not null default now(),
    -- The one the protein draws now.
    active          boolean not null default false,

    constraint fold_generation_source_known check (source in ('pdb', 'alphafold', 'none')),
    constraint fold_generation_state_known check (state in ('ready', 'refused')),
    constraint fold_generation_pdb_folds check (source <> 'pdb' or folding_track is not null)
);

create unique index fold_generation_active_idx on fold_generation (slug) where active;
create index fold_generation_slug_idx on fold_generation (slug, id desc);

alter table fold_selection enable row level security;
alter table fold_generation enable row level security;

commit;
