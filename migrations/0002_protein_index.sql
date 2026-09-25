-- The protein index: every reviewed human protein, whether or not the app can
-- walk it yet. Search suggests from here.
--
-- Run by hand, once, against the Supabase project:
--
--     psql "$DATABASE_URL" -f migrations/0002_protein_index.sql
--
-- and filled by `scripts/load_protein_index.py`, which replaces every row in
-- one transaction whenever UniProt or MANE publishes a release.
--
-- Everything here is additive. The catalog tables are untouched: a protein the
-- app can walk is a `protein` row, and this index only says what could become
-- one.

begin;


-- One row per (UniProt entry, gene). Usually that is one row per entry, but
-- not always: P69905 is made by both HBA1 and HBA2, and each gene is its own
-- walk, so the entry is two rows. An entry with no MANE Select transcript is
-- one row under UniProt's own gene name, or '' where UniProt names none.
--
-- `gene` is MANE's symbol wherever MANE has the transcript, because that is the
-- gene a build would fetch. `buildable` says whether the pipeline can start
-- from this row alone, and `unavailable_reason` says why not where it cannot.
create table protein_index (
    uniprot            text not null,
    gene               text not null default '',

    name               text not null,          -- UniProt's recommended name
    length             integer not null,       -- residues in the canonical sequence
    annotation_score   smallint not null,      -- UniProt's 1-5
    existence          smallint not null,      -- 1 protein level ... 5 uncertain

    gene_id            integer,                -- NCBI Gene
    hgnc_id            text,

    -- The MANE Select transcript and protein, versioned as MANE writes them.
    refseq_nuc         text,                   -- NM_000207.3
    refseq_prot        text,                   -- NP_000198.1
    ensembl_nuc        text,                   -- ENST00000381330.5
    ensembl_prot       text,
    -- The UniProt isoform MANE's protein is, where it is not the canonical one.
    mane_isoform       text,                   -- Q13625-3

    -- Where the MANE transcript sits on GRCh38, 1-based inclusive.
    chrom_acc          text,                   -- NC_000011.10
    chrom_start        integer,
    chrom_end          integer,
    chrom_strand       smallint,

    -- The RefSeqGene record whose LRG_RefSeqGene row lists `refseq_prot`, in
    -- any of that file's categories. Null where there is none, and a build
    -- then reads a chromosome slice instead (R1.2).
    refseqgene         text,                   -- NG_007114.1

    buildable          boolean not null,
    unavailable_reason text,

    primary key (uniprot, gene),
    constraint protein_index_reason check (buildable or unavailable_reason is not null),
    constraint protein_index_strand check (chrom_strand is null or chrom_strand in (1, -1))
);

-- A build is keyed by its gene, so two buildable rows may not share one.
create unique index protein_index_buildable_gene on protein_index (lower(gene)) where buildable;


-- Every string a search can match, normalised the way `app/protein_index.py`
-- normalises a query: lower-case ASCII letters and digits, single spaces,
-- Greek letters spelled out. "Insulin-degrading enzyme" is stored as
-- "insulin degrading enzyme", and so are each of its words.
--
-- The column is collated "C" so a prefix is a plain range, `term >= 'ins' and
-- term < 'int'`, which the btree serves under a generic plan too. psycopg
-- prepares a statement it has run five times, and a prepared `like $1` cannot
-- use a btree at all, because the planner no longer knows the pattern is a
-- prefix. The trigram index serves near misses ("insuln").
create table protein_index_term (
    uniprot text not null,
    gene    text not null,
    term    text collate "C" not null,
    -- 0 an accession, 1 the gene symbol, 2 a gene synonym, 3 the protein's
    -- name (UniProt's recommended one, MANE's gene name, the app's own), 4
    -- any other name (alternative, short, CD, INN, what it is cut into), 5 one
    -- word of a name. The search ranks by it after exact matches, and only a
    -- whole term is ever an exact match: "ins" is a word of an inositol
    -- kinase's name, and must not rank it beside insulin.
    kind    smallint not null,

    primary key (uniprot, gene, term, kind),
    foreign key (uniprot, gene) references protein_index (uniprot, gene) on delete cascade,
    constraint protein_index_term_kind check (kind between 0 and 5)
);

create index protein_index_term_term on protein_index_term (term);
create index protein_index_term_trgm on protein_index_term using gin (term gin_trgm_ops);


-- Which releases the index was built from, for the answer to say so. One row.
create table protein_index_release (
    only_row   boolean primary key default true check (only_row),
    uniprot    text not null,                  -- 2026_03
    mane       text not null,                  -- v1.5
    entries    integer not null,
    buildable  integer not null,
    loaded_at  timestamptz not null default now()
);


-- Same posture as 0001: row level security on and no policies, so only the
-- backend's own role reaches these rows.
alter table protein_index enable row level security;
alter table protein_index_term enable row level security;
alter table protein_index_release enable row level security;

commit;
