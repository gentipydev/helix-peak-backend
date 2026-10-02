-- Phase 6: a reader asks for a protein the catalog does not list, and it is
-- resolved on demand.
--
-- Run by hand, once, against the Supabase project, after 0001 to 0008:
--
--     psql "$DATABASE_URL" -f migrations/0009_resolve_request.sql
--
-- `bake_job` cannot hold the request: its `slug` references `protein`, and a
-- protein that has not been resolved has no row yet. Resolving is what writes
-- that row, from the protein's `protein_index` entry and its UniProt features,
-- with the record track beside it; the bakes that follow (ESM-2 constraint)
-- are `bake_job` rows as 0001 planned. So this table is the step before: one
-- row per request, written by the service (`POST /resolve`) and worked by
-- `pipeline/resolver/` on Modal.
--
-- No `protein` column and no existing row, object or provenance is touched.
-- Everything here is additive: one table, two partial indexes, RLS on and no
-- policies, the posture every table in 0001 takes.

begin;

create table resolve_request (
    id               bigserial primary key,

    -- What was asked for, as `protein_index` names it. A buildable entry is
    -- unique by lower(gene) (`protein_index_buildable_gene`), so the gene is
    -- the request's key and the UniProt accession records which entry it was.
    gene             text not null,
    uniprot          text not null,

    -- The slug the resolved protein takes: lower(gene), the rule 0001 set for
    -- every protein resolved after the twenty.
    slug             text not null,

    -- queued    asked for, and not started
    -- running   a worker holds it
    -- done      the protein row and its record track are written
    -- refused   the resolver declined, and `reason` is the sentence why
    -- failed    it broke more often than it is retried; `reason` says how
    state            text not null default 'queued',
    reason           text,
    attempts         integer not null default 0,

    -- Which resolver finished it. A refusal stands for the version that made
    -- it; a later version may decide differently.
    resolver_version integer,

    requested_at     timestamptz not null default now(),
    started_at       timestamptz,
    finished_at      timestamptz,

    constraint resolve_request_state_known check (state in (
        'queued', 'running', 'done', 'refused', 'failed'
    )),
    constraint resolve_request_said_why check (
        state not in ('refused', 'failed') or reason is not null
    )
);

-- Two readers asking for BRCA1 in the same minute make one request, not two.
-- Partial, so the history of finished requests stays.
create unique index resolve_request_inflight_idx on resolve_request (lower(gene))
    where state in ('queued', 'running');

-- The worker's work list, oldest first.
create index resolve_request_queue_idx on resolve_request (requested_at)
    where state = 'queued';

-- The daily cap counts requests by when they were made.
create index resolve_request_requested_idx on resolve_request (requested_at);

alter table resolve_request enable row level security;

commit;
