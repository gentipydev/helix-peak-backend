# PDB fixtures

What the tests of the experimental-structure choice read in place of the
network, recorded on 2026-10-12.

- `6rlx_cut.cif`, `2mv1_cut.cif`, `2k6t_cut.cif`, `1ppb_cut.cif`: PDBe's
  SIFTS-enriched mmCIF (`<id>_updated.cif`) of four entries, cut down to the
  categories `normalize.examine` reads (`entry`, `exptl`, `struct_conf_type`, `entity_poly`,
  `pdbx_poly_seq_scheme`, `pdbx_sifts_xref_db` with one row a residue,
  `atom_site`, `struct_conf`, `struct_sheet_range`, the disulfides of
  `struct_conn`, `pdbx_struct_mod_residue`, `pdbx_nmr_representative`) and
  written out again; values are as served.
  - 6RLX: human relaxin-2, X-ray, both chains; a pyroglutamate.
  - 2MV1: human relaxin-2, NMR, its first two conformers.
  - 2K6T: human INSL3, NMR, its first three conformers.
  - 1PPB: human thrombin, X-ray, in chymotrypsinogen numbering (insertion
    codes), its two protein chains.
- `graphql_*.json`: RCSB data API answers to `pdb_search._QUERY`, for the
  entities each name says (relaxin: 2MV1 and 6RLX; insl3: 2K6T; thrombin: 1PPB;
  ttr: 1BM7, 8C86, 1TTA; b2m: 1K5N's; cftr: 9DW5's; sod1: 2C9V's).
- `uniprot_*.json`: UniProt entries (release 2026_03) cut down to their
  sequence and the features the choice reads.
- `pfam_*.json`: InterPro's Pfam entries on p53 and the prion protein, as served.

PDB data are free of copyright (CC0 1.0); UniProt and InterPro are CC BY 4.0.
