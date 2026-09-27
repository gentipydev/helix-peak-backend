-- Assemblies: molecules made of more than one gene's chains, kept apart from
-- the catalog.
--
-- Run by hand, once, against the Supabase project, after 0001 to 0005:
--
--     psql "$DATABASE_URL" -f migrations/0006_assemblies.sql
--
-- The catalog is twenty proteins, one gene each: `protein` is one gene's
-- product and every `protein_track` row belongs to one. The hemoglobin
-- tetramer is two genes' chains, so it cannot be a `protein` row, and a track
-- of it cannot be a `protein_track` row: that table's slug is a foreign key
-- into `protein`, and a row there is one `/catalog`, `/catalog/search` and
-- `/proteins/suggest` would all serve. So an assembly has tables of its own,
-- which none of them reads. The catalog stays twenty, and no `protein` column
-- and no existing row, object or provenance is touched.
--
-- `pipeline/assemblies/seed_assemblies.py` writes `assembly`, and
-- `upload_assemblies.py` the `assembly_track` rows, the way `seed_catalog.py`
-- and `upload_tracks.py` write theirs. `assembly_track` is `protein_track`'s
-- shape, state for state.
--
-- Everything here is additive, in one transaction, RLS on and no policies.

begin;

create table assembly (
    slug        text primary key,
    display     text not null,
    -- Each chain: its node, its gene, its UniProt accession and its letter.
    subunits    jsonb not null,
    -- Each state: its name, its PDB entry, what is bound, how it is built.
    states      jsonb not null,
    provenance  jsonb not null,
    updated_at  timestamptz not null default now()
);

create table assembly_track (
    slug             text not null references assembly (slug) on delete cascade,
    kind             text not null,
    state            text not null,
    reason           text,

    bucket           text,
    object_path      text,
    bytes            bigint,
    sha256           text,
    content_encoding text,
    format           text not null,

    provenance       jsonb not null,
    updated_at       timestamptz not null default now(),

    primary key (slug, kind),
    constraint assembly_track_kind_known check (kind in ('morph')),
    constraint assembly_track_state_known check (state in (
        'ready', 'pending', 'absent', 'refused'
    )),
    constraint assembly_track_ready_has_object check (
        state <> 'ready' or (bucket is not null and object_path is not null)
    ),
    constraint assembly_track_refused_has_reason check (
        state <> 'refused' or reason is not null
    )
);

alter table assembly enable row level security;
alter table assembly_track enable row level security;

commit;
