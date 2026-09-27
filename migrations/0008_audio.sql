-- The audio track kind: each protein as sound, one note a residue, for the
-- lab's Listen, with the map from each residue to the millisecond its note
-- starts carried inside the file.
--
-- Run by hand, once, against the Supabase project, after 0001 to 0007:
--
--     psql "$DATABASE_URL" -f migrations/0008_audio.sql
--
-- A new kind is a new value in the check `protein_track` holds `kind` to, and
-- nothing else: the rows come from `pipeline/upload_tracks.py --kind audio`,
-- and the service serves any kind a row names. No `protein` column and no
-- existing row, object or provenance is touched.
--
-- The playbook named this file 0007_audio; 0007 is the locus track's, so this
-- is the next number, and its list is 0007's with one more value.
--
-- Everything here is additive: the check is dropped and re-added with one more
-- value, in one transaction, so no moment exists without it.

begin;

alter table protein_track drop constraint protein_track_kind_known;
alter table protein_track add constraint protein_track_kind_known check (kind in (
    'record', 'constraint', 'impact', 'clinvar', 'structure', 'impact_explanations',
    'structure_ar', 'trafficking', 'folding', 'locus', 'audio'
));

commit;
