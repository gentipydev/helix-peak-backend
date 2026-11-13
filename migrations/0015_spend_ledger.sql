-- What builds on demand have cost, kept where nothing that deletes a build can
-- reach it (2026-11-13, the security review before the app's release).
--
-- Run by hand, once, against the Supabase project, after 0001 to 0014:
--
--     .venv/bin/python scripts/apply_migration.py migrations/0015_spend_ledger.sql
--
-- The service counts a day's asks to cap them (`RESOLVES_PER_DAY`, and a cap
-- for each install). It used to count `resolve_request` rows, which a discard
-- deletes, so ask, stop and discard in a loop never reached the cap; and
-- asking again for a stopped protein queued its GPU scoring without counting
-- at all. Each is now a row here, which no discard, forget or cascade touches.
--
-- The resolver on Modal writes a row for each call that costs money (`score`
-- on a GPU, `fold` and `evidence` on CPUs) before it claims anything: `usd` is
-- what the call could cost at most, its timeout at Modal's rates, and is
-- replaced by what it did cost once it ends (`finished_at`). A call that dies
-- leaves its worst case standing, which is what it may have cost. The budget
-- (`app/budget.py`) is the sum of `usd` in the last day and this month.
--
-- Nothing else is touched. One table with RLS on and no policies.

begin;

create table spend_ledger (
    id           bigserial primary key,
    at           timestamptz not null default now(),

    -- `ask`: a request queued by the service; `resume`: a stopped protein's
    -- scoring queued again. `score`, `fold`, `evidence`: a call on Modal.
    kind         text not null,

    -- The gene asked for, and the install that asked (as `resolve_request`).
    gene         text,
    asker        text,

    -- A call on Modal: what it may cost, then what it cost, and how long it
    -- ran once it has ended.
    usd          numeric(10, 6),
    seconds      double precision,
    finished_at  timestamptz,

    constraint spend_ledger_kind_known
        check (kind in ('ask', 'resume', 'score', 'fold', 'evidence')),
    constraint spend_ledger_usd_known
        check ((kind in ('ask', 'resume')) = (usd is null))
);

create index spend_ledger_kind_at on spend_ledger (kind, at);
create index spend_ledger_asker_at on spend_ledger (asker, at) where asker is not null;

alter table spend_ledger enable row level security;

commit;
