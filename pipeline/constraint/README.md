# Offline constraint scores

Twenty tracks, one for each protein in the catalog, scored with
`facebook/esm2_t33_650M_UR50D` (MIT licensed). A row added before its track is
baked says `scored=False` in `pipeline/targets.py`: `--all` passes it by, `--target`
refuses it, and the app draws its protein page without the conservation toolbar.
The checkpoint revision is pinned in `score_protein.py`. Weights are downloaded
once to Hugging Face's normal cache, outside the app. This directory is never a
runtime dependency.

From the Flutter project root, with Python 3.12:

```sh
python3.12 -m venv pipeline/.esm-venv
pipeline/.esm-venv/bin/python -m pip install -r pipeline/constraint/requirements-lock.txt
pipeline/.esm-venv/bin/python -m unittest discover -s pipeline/constraint -t . -v
pipeline/.esm-venv/bin/python -u pipeline/constraint/score_protein.py --all
```

`requirements.txt` records direct dependencies; `requirements-lock.txt` records
the environment used for generation. Apple Silicon uses MPS, with CPU fallback
on machines without MPS. Forward passes use float32 and evaluation mode.

About fifty minutes for the first ten on an M-series laptop, of which dystrophin
is forty-five, and thirty-four for the ten after them, of which CFTR is
twenty-one. `--target <slug>` bakes one, and can be repeated.

## The sequence is not an argument

`score_protein` reads the protein out of `assets/mock/gene_<gene>.json` — the
record the app draws on the page this track colours. So the two cannot describe
different molecules, which is the failure this arrangement exists to prevent:
`AnatomyScreen` only applies a track where the stage's own letters are the
track's sequence, and a mismatch is not an error, it is a page that quietly
loses its colour. Build the gene records first.

Region boundaries and disulfide pairs come from `pipeline/targets.py` (they are
UniProt features, pinned) and are written into the asset, which is what lets
the Dart side name a residue's domain and its bonding partner without knowing
which protein it is looking at. Changing a boundary does not need a re-score:

```sh
pipeline/.esm-venv/bin/python pipeline/constraint/score_protein.py --metadata-only --all
```

## Method

One residue masked per pass. The scorer verifies the complete token sequence,
the single mask, BOS/EOS, and every unchanged context token on every pass. Each
substitution is the natural log probability ratio against the wildtype, rounded
to three decimals. The wildtype entry is included and is exactly zero.

Entropy is computed from all 33 softmax probabilities, including special and
noncanonical tokens. Float64 CPU reductions are used for entropy and log
ratios. Conservation is `(max_entropy - entropy) / (max_entropy - min_entropy)`
across the entire sequence. There is no per-domain normalization or manual
reweighting — so conservation is a rank within one protein and means nothing
across two. Raw entropy is included in the asset for auditing.

**Windows, for dystrophin and CFTR.** ESM-2's positions are rotary, so a
sequence past the trained 1,022-residue context runs rather than raising; it is
just being asked about distances it never saw in training. Dystrophin's 3,685
residues and CFTR's 1,480 are therefore scored in windows of 1,022 centred on
the mask, so every residue gets the largest context the model was trained to use
and the ones away from the termini get it on both sides. Every other protein
takes one window that is the whole of it, down the identical code path. The
asset records which under `context`.

## The gate

Publication is gated, and the gate runs before any file is created.

Universally: the track has to be the sequence, position for position;
conservation has to be finite in [0, 1]; the wildtype has to score exactly
zero. These are what the Dart loader checks anyway, checked earlier and with a
better message.

For every protein, also **alignment**: wherever a residue differs from the one
beside it, the model has to score the residue that is there above its
neighbour's more often than not, on the side before and on the side after. What
that is for is offset bugs — a mask or context offset files each score under
the wrong residue, and a track shifted by one position scores each residue's
neighbour as its own, which inverts the measure. The token checks on every pass
rule most of those out; this catches what gets past them, from the scores
themselves. Runs of one amino acid are not counted, since a Q beside a Q says
nothing about which one a score belongs to. Every shipped track clears it with
room: the lowest is relaxin, at 64.9% on its worse side.

Where a protein has disulfides, the bonded residues have to be cysteines, and
their ranks by lowest entropy are printed as findings but not gated. They used
to be the gate, on the argument that they are the most constrained residues a
protein has, and the rule was rewritten three times, each time because the data
was right and the rule was not:

- *all within the top decile* cannot be satisfied at all by oxytocin, whose
  precursor has sixteen bonded cysteines in 125 residues;
- *all within two places each* then rejected growth hormone over C215, which
  ranks 28/217 while its other three rank 1, 2 and 4 — the small C-terminal
  loop being genuinely more tolerant than the bridge holding the helix bundle;
- *the median in the top decile, none past halfway* then rejected three of the
  second ten: SOD1's bridge ranks 46 and 51 of 154, under the histidines that
  hold its copper and zinc; leptin's ranks 24 and 33 of 167, under its core
  leucines; amylase's C85 and C130 rank 449 and 265 of 511, while its other
  eight bonded cysteines rank 1–41.

A proxy made of biology can be beaten by biology, and it checked nothing at all
on the seven proteins in the catalog that have no disulfides. Alignment is the thing the
cysteines stood in for, measured directly.

The scorer also prints regional means and the most constrained positions, so
the proposed pattern can be checked against the actual result before UI work.
See [`verification.md`](verification.md) for the measured results.

## Independent audit

After generation, a second implementation re-derives every number on CPU:

```sh
pipeline/.esm-venv/bin/python -u pipeline/constraint/verify_cpu.py --target insulin
```

It takes residue offsets from the tokenizer's special-token mask rather than
from `i + 1`, masks token IDs directly, checks log ratios as raw logit
differences, uses NumPy entropy, and re-derives the alignment measure from its
own logits. It reads cached weights and never modifies the asset. Full-length
tracks only — auditing dystrophin's or CFTR's would mean sharing the windowing,
which is the one thing an independent check must not do.

## Sources

- [UniProt](https://rest.uniprot.org/) for sequences, regions and disulfides
- [ESM-2 checkpoint](https://huggingface.co/facebook/esm2_t33_650M_UR50D)
- [Meta's masked-marginal implementation](https://github.com/facebookresearch/esm/blob/main/examples/variant-prediction/predict.py)
- [Meier et al. 2021](https://proceedings.neurips.cc/paper/2021/hash/f51338d736f95dd42427296047067694-Abstract.html)
