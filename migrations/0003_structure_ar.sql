-- The structure_ar track kind: each fold as a USDZ at its real size, for AR.
--
-- Run by hand, once, against the Supabase project, after 0001 and 0002:
--
--     psql "$DATABASE_URL" -f migrations/0003_structure_ar.sql
--
-- A new kind is a new value in the check `protein_track` holds `kind` to, and
-- nothing else: the rows come from `pipeline/upload_tracks.py --kind
-- structure_ar`, and the service serves any kind a row names. The `structure`
-- rows the walk reads, their objects and their provenance are untouched.
--
-- Everything here is additive: the check is dropped and re-added with one more
-- value, in one transaction, so no moment exists without it.

begin;

alter table protein_track drop constraint protein_track_kind_known;
alter table protein_track add constraint protein_track_kind_known check (kind in (
    'record', 'constraint', 'impact', 'clinvar', 'structure', 'impact_explanations',
    'structure_ar'
));

commit;
