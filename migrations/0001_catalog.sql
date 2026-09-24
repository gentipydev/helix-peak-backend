-- The catalog, the track index and the bake queue.
--
-- Run by hand, once, against the Supabase project:
--
--     psql "$DATABASE_URL" -f migrations/0001_catalog.sql
--
-- Not startup DDL. `genbank_record` creates itself in `app/record_cache.py`
-- because it is one table with no history to migrate and a deploy that cannot
-- half-apply. This file is five tables, an extension and two partial indexes;
-- applying that from a web worker on every boot would be a deploy that can
-- half-apply, run concurrently with itself, and fail in a place where the only
-- honest response is to keep serving.
--
-- Everything here is additive. `genbank_record` is untouched.

begin;


-- One protein: what `tool/targets.py` holds on the bake side and
-- `protein_catalog.dart` holds on the Dart side, in the one place both will
-- read it from instead.
--
-- `slug` is the key for every storage object, every URL and every legacy asset
-- path, so it is written down once and never derived from a display name --
-- HBA1 and HBA2 would both want "hemoglobin alpha". New proteins take
-- `lower(gene)`; the twenty that shipped keep the slug they shipped with, and
-- `protein_alias` carries the tie.
create table protein (
    slug             text primary key,
    gene             text not null,
    uniprot          text not null,

    -- A RefSeqGene accession where one exists, otherwise the chromosome
    -- accession the slice is cut from. The three slice columns are null
    -- together or set together; set means `accession` alone does not identify
    -- the record and the fetch carries seq_start/seq_stop/strand.
    accession        text not null,
    slice_start      integer,
    slice_end        integer,
    slice_strand     smallint,

    -- MANE Select. Written down because the CDS does not always identify the
    -- transcript: AMY1A's slice carries two mRNAs around identical CDS
    -- coordinates and only one of them is MANE Select.
    transcript_id    text,
    protein_id       text,

    taxon_id         integer not null default 9606,
    display          text not null,
    summary          text not null,

    -- What an uncut protein calls its one chain. Null for a precursor that is
    -- cut, which names its pieces instead.
    chain_name       text,

    -- Whether the record's mature peptides are a cleavage series worth a page,
    -- or fragments and alternatives that one page of disjoint chains cannot
    -- draw. False stops the walk at the precursor.
    mature_peptides  boolean not null default true,

    residues         integer not null,
    exons            integer not null,
    chains           integer not null,
    bridges          integer not null,

    -- [{label, short, start, end, origin, kept}], a gapless tiling of 1..residues.
    regions          jsonb not null,
    -- [[31, 96], ...] in precursor numbering.
    disulfides       jsonb not null,
    -- StructureChrome plus the chain->node->tint mapping the fold page paints by.
    -- Null where no entry passed the picker; the structure track says why.
    structure        jsonb,

    -- Where each field came from, which positions the record and UniProt
    -- disagree on, and whether the prose was written or templated. A resolved
    -- protein must be able to say how it was resolved.
    provenance       jsonb not null,
    resolver_version integer not null,
    resolved_at      timestamptz not null default now(),

    -- The order an empty search shows. The twenty that shipped run insulin
    -- first and then small to large, which is a reading order and not a
    -- sorting of any field; it is written down because it cannot be derived.
    -- A protein resolved on demand has no place in that sequence and sorts
    -- after all of them, by name.
    catalog_order    integer,

    -- One row per gene per organism. HBA1 and HBA2 are two rows and are never
    -- merged.
    unique (gene, taxon_id)
);

create index protein_gene_idx on protein (gene);


-- Everything a search can match against, flattened. `matching()` ranks in four
-- tiers -- exact, prefix, contains, summary -- and this is the table that
-- serves the first three without a scan.
create table protein_alias (
    slug  text not null references protein (slug) on delete cascade,
    alias text not null,
    -- One row per string `matching()` can rank against: the five names it
    -- compares whole (display, gene, slug, uniprot, accession), each word of
    -- the display name, UniProt's synonyms, and the slug a protein shipped
    -- with before its gene symbol became the rule.
    kind  text not null,   -- display | gene | slug | uniprot | accession | word | synonym | legacy_slug
    primary key (slug, alias, kind)
);

create extension if not exists pg_trgm;

-- text_pattern_ops serves the prefix tier (`like 'ins%'`), which the default
-- opclass will not use on a non-C collation. The trigram index serves the
-- contains tier and near misses.
create index protein_alias_prefix_idx on protein_alias (lower(alias) text_pattern_ops);
create index protein_alias_trgm_idx   on protein_alias using gin (lower(alias) gin_trgm_ops);


-- Where one track for one protein is, and what state it is in.
--
-- This replaces the four booleans in `protein_catalog.dart`. A boolean can say
-- "there is no ClinVar snapshot for this gene" but it cannot tell an unbaked
-- gene from one whose bake is four minutes out, and it cannot say that the
-- gene's exons exceed the page budget. Those are different sentences on
-- screen, so they are different states here.
create table protein_track (
    slug             text not null references protein (slug) on delete cascade,
    kind             text not null,
    -- ready    the object is there and `object_path` points at it
    -- pending  a bake_job is queued or running; arriving, not absent
    -- absent   nothing was asked for, and nothing is coming
    -- refused  the pipeline declined, and `reason` is the sentence why
    state            text not null,
    reason           text,

    bucket           text,
    object_path      text,
    bytes            bigint,
    -- The client's cache key. Content is immutable, so a changed sha is a new
    -- file and there is no invalidation to get wrong.
    sha256           text,
    content_encoding text,
    format           text not null,   -- json | packed-v1 | glb

    -- Model, revision, method and vocabulary for a constraint track; scorer and
    -- units for impact; retrieval date, query and exclusions for ClinVar; entry,
    -- resolution, SEQADV and pLDDT for a structure. The About sheet reads this
    -- rather than a string constant, because two constraint tracks can now come
    -- from two different models.
    provenance       jsonb not null,
    updated_at       timestamptz not null default now(),

    primary key (slug, kind),
    constraint protein_track_kind_known check (kind in (
        'record', 'constraint', 'impact', 'clinvar', 'structure', 'impact_explanations'
    )),
    constraint protein_track_state_known check (state in (
        'ready', 'pending', 'absent', 'refused'
    )),
    -- A ready track has somewhere to read it from; a refused one has a sentence.
    constraint protein_track_ready_has_object check (
        state <> 'ready' or (bucket is not null and object_path is not null)
    ),
    constraint protein_track_refused_has_reason check (
        state <> 'refused' or reason is not null
    )
);

-- The worker's work list.
create index protein_track_pending_idx on protein_track (kind) where state = 'pending';


create table bake_job (
    id           bigserial primary key,
    slug         text not null references protein (slug) on delete cascade,
    kind         text not null,
    state        text not null default 'queued',   -- queued | running | done | failed
    attempts     integer not null default 0,
    error        text,
    requested_at timestamptz not null default now(),
    started_at   timestamptz,
    finished_at  timestamptz,

    constraint bake_job_state_known check (state in ('queued', 'running', 'done', 'failed'))
);

-- Two people searching CFTR in the same minute enqueue one job, not two. The
-- index is partial so the history of finished jobs stays, and a protein can be
-- rebaked later without deleting what happened last time.
create unique index bake_job_inflight_idx on bake_job (slug, kind)
    where state in ('queued', 'running');

create index bake_job_queue_idx on bake_job (requested_at) where state = 'queued';


-- The client never holds a Supabase key; it reads blobs from public buckets over
-- plain HTTPS and everything else through the backend, which holds service_role.
-- So: row level security on, and no policies at all. That denies anon and
-- authenticated outright and leaves service_role unaffected, which is the
-- posture that stays correct if an anon key ever does end up in a build.
alter table protein enable row level security;
alter table protein_alias enable row level security;
alter table protein_track enable row level security;
alter table bake_job enable row level security;

commit;
