-- The folding track kind: each fold's CA trace, residue by residue, with its
-- secondary structure, in the stored structure model's frame, for the lab's
-- fold animation.
--
-- Run by hand, once, against the Supabase project, after 0001 to 0004:
--
--     psql "$DATABASE_URL" -f migrations/0005_folding.sql
--
-- A new kind is a new value in the check `protein_track` holds `kind` to, and
-- nothing else: the rows come from `pipeline/upload_tracks.py --kind
-- folding`, and the service serves any kind a row names. No `protein` column
-- and no existing row, object or provenance is touched.
--
-- Everything here is additive: the check is dropped and re-added with one more
-- value, in one transaction, so no moment exists without it.

begin;

alter table protein_track drop constraint protein_track_kind_known;
alter table protein_track add constraint protein_track_kind_known check (kind in (
    'record', 'constraint', 'impact', 'clinvar', 'structure', 'impact_explanations',
    'structure_ar', 'trafficking', 'folding'
));

commit;
