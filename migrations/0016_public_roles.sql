-- Supabase's public roles get nothing in `public` (2026-11-13, the security
-- review before the app's release).
--
-- Run by hand, once, against the Supabase project, after 0015:
--
--     .venv/bin/python scripts/apply_migration.py migrations/0016_public_roles.sql
--
-- Nothing reads these tables through Supabase's Data API: the service and the
-- resolver reach them over DATABASE_URL, and the app reads storage's public
-- URLs alone. But Supabase grants `anon` and `authenticated` every privilege
-- on every table `postgres` makes there, so row level security with no
-- policies was the only thing between the anon key and the rows, and one
-- policy, `genbank_record`'s "Public read access", let them read that table.
-- Here those grants go, for the tables there are and for those to come, and
-- that policy with them. `service_role` and `postgres` keep theirs.
--
-- On a Postgres that is not Supabase's (the opt-in SQL tests), the roles do
-- not exist and nothing is done for them. `genbank_record` is made by the
-- service at startup (`app/record_cache.py`), so it may not be there yet.

begin;

do $$
declare
    role_name text;
begin
    foreach role_name in array array['anon', 'authenticated'] loop
        if exists (select 1 from pg_roles where rolname = role_name) then
            execute format('revoke all on all tables in schema public from %I', role_name);
            execute format('revoke all on all sequences in schema public from %I', role_name);
            execute format('alter default privileges for role postgres in schema public '
                           'revoke all on tables from %I', role_name);
            execute format('alter default privileges for role postgres in schema public '
                           'revoke all on sequences from %I', role_name);
        end if;
    end loop;

    if to_regclass('public.genbank_record') is not null then
        alter table genbank_record enable row level security;
        drop policy if exists "Public read access" on genbank_record;
    end if;
end
$$;

commit;
