# AlphaFold DB fixtures

What `test_alphafold.py` reads in place of the network, fetched on 2026-10-08.

- `AF-O00631-F1-model_v6.pdb`: AlphaFold DB's model of human sarcolipin (SLN,
  UniProt O00631), 31 residues, as served.
- `O00631.json`, `Q9H6X2.json`: the prediction API's answers for sarcolipin and
  for ANTXR1, cut down to the keys the bake reads and a few it is checked
  against. ANTXR1's keeps two of its five isoform entries, put first, and the
  canonical one.

AlphaFold DB data are CC BY 4.0 (EMBL-EBI and Google DeepMind). Jumper et al.
2021, Nature 596:583; Varadi et al. 2024, Nucleic Acids Res. 52:D368.
