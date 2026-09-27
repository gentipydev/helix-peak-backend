-- The locus track kind: where on its chromosome each protein's gene lies, by
-- cytogenetic band, with every band of that chromosome, for the lab's zoom.
--
-- Run by hand, once, against the Supabase project, after 0001 to 0006:
--
--     psql "$DATABASE_URL" -f migrations/0007_locus.sql
--
-- A new kind is a new value in the check `protein_track` holds `kind` to, and
-- nothing else: the rows come from `pipeline/upload_tracks.py --kind locus`,
-- and the service serves any kind a row names. No `protein` column and no
-- existing row, object or provenance is touched. In particular there is no
-- cytoband column: the band is a track, like every other new datum.
--
-- 0006 created the assembly tables and left this check as 0005 left it, so
-- the list below is 0005's with one more value.
--
-- Everything here is additive: the check is dropped and re-added with one more
-- value, in one transaction, so no moment exists without it.

begin;

alter table protein_track drop constraint protein_track_kind_known;
alter table protein_track add constraint protein_track_kind_known check (kind in (
    'record', 'constraint', 'impact', 'clinvar', 'structure', 'impact_explanations',
    'structure_ar', 'trafficking', 'folding', 'locus'
));

commit;
