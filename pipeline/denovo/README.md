# A de novo protein model, from scratch

A learning project: a small language model that writes peptide sequences one
residue at a time, built up a file at a time. Offline only. Nothing here is a
track, nothing is uploaded, and the service never imports it.

## Not ESM-2's kind of model

ESM-2 is masked: trained to fill a hidden residue from both sides. That is what
makes it a scorer (the constraint track is its masked marginals), but it has no
first residue to start from and no order to write in. This model is
autoregressive, like ProtGPT2 and ProGen: each residue is predicted from the
ones before it, so sampling one at a time writes a sequence from `<bos>` to
`<eos>`.

## Steps

1. `alphabet.py`: the token set, 23 ids, and the lengths modelled. Done.
2. `tokenizer.py`: a sequence to ids and back, framed by `<bos>`/`<eos>`. Done.
3. `fetch_peptides.py` + `sequences.py` + `dedup.py`: the corpus. Done.
4. `bigram.py` + `bigram_torch.py` + `baselines.py`: the loss to beat. Done.
5. `helix_peek_LM.py` + `train.py` + `sweep.py`: the transformer. Done.
6. Sampling: temperature and top-k, from `<bos>` until `<eos>`. **Not written.**

## The corpus: 10,416 peptides, and why not 20,000

Reviewed UniProt does not hold 20,000 peptides of 10-50 residues. It holds
12,423; non-fragment, 9,460; leaving out the toxin keyword, 8,273; of those,
8,245 use only the twenty, and **a third of them are exact copies of each
other**. Whole reviewed entries alone come to about 5,500 distinct peptides.

The rest come from inside longer precursors. An entry's PEPTIDE and PROPEPTIDE
features name the mature products cut out of it -- the insulin chains, each
neuropeptide of a prohormone -- and those are peptides in every sense that
matters here. Together:

| source | count |
|---|---|
| reviewed whole entries, 10-50 aa, non-fragment, not toxin | 8,273 |
| mature products (PEPTIDE 9,745 + PROPEPTIDE 3,230) | 12,975 |
| of those, 10-50 canonical residues | 9,072 |
| **corpus, after range, alphabet and exact duplicates** | **10,416** |
| clusters, at 3-mer Jaccard >= 0.5 | 6,750 |

Rebuild it with `python fetch_peptides.py`; `data/provenance.json` records the
queries, the UniProt release and the corpus sha256. Reaching 20,000 means
unreviewed sequences -- 1,076,603 exist at this length, nearly all predicted
small ORFs rather than peptides anyone has seen.

Toxin-keyword entries are left out. At this length they are the most redundant
families in the database, so keeping them would skew the residue distribution
and flatter the model through near-duplicates.

### The split is by cluster, never by sequence

A third of the corpus starts out as exact duplicates and much of the rest is a
substitution or two apart. Split that at random and training holds a copy of
half the validation set. `dedup.py` clusters on k-mer Jaccard -- one matrix
multiply on the GPU, seconds for the whole corpus -- and holds out whole
clusters. `python dedup.py --report` prints the sizes and asserts that no
cluster lands on both sides.

## The loss to beat

Nats per token on the held-out clusters, scored over the 21 tokens that can
actually follow (`python baselines.py`):

| model | nats/token | perplexity |
|---|---|---|
| uniform | 3.0445 | 21.00 |
| unigram | 2.9493 | 19.09 |
| **bigram** | **2.8918** | **18.03** |

A bigram barely beats counting residues, which is the usual finding for protein
sequence: adjacent residues say little. The transformer is what tests whether a
longer context says more.

## What this Mac does

MacBook Pro 14", Apple M5, 10-core CPU, 24 GB unified memory, torch 2.14.0 on
MPS. One training step, forward through AdamW, at 51 positions:

| device | model | batch | ms/step | tokens/s | s/epoch @20k |
|---|---|---|---|---|---|
| mps | 128/4/4/512 | 16 | 5.8 | 144k | 7.2 |
| **mps** | **128/4/4/512** | **64** | **17.7** | **188k** | **5.5** |
| mps | 128/4/4/512 | 256 | 78.9 | 169k | 6.2 |
| mps | 128/4/4/512 | 1024 | 371 | 143k | 7.3 |
| cpu | 128/4/4/512 | 64 | 64.9 | 51k | 20.3 |
| mps | 256/8/8/1024 | 16 | 30.0 | 28k | 37.5 |

Three things follow. The GPU is never saturated by a model this small, so
throughput hardly moves with batch size and **64 is the fastest per epoch** --
256 and 1024 are slower. **MPS is about four times CPU**, which is the whole
difference between a loop that moves tensors to a device and one that does not.
And memory is not a constraint and will not become one: the tokenised corpus is
about 5 MB and a batch of 1024 stays under 2 GB of the 24.

The real run is 9,367 training peptides, 147 steps an epoch, **2.3 s an epoch**.
A hundred epochs is four minutes. The machine was never the limit here; how many
real peptides exist is.

## The hyperparameters, and why

The model is 456,887 parameters against roughly 290k training tokens -- about
1.6 parameters per token, where compute-optimal is nearer 0.05. It memorises, so
the choices that matter are the ones that decide when to stop and how hard to
regularise, not the architecture.

| property | was | is | why |
|---|---|---|---|
| positional capacity | 100 | **51** | `<bos>` + 50 residues is the longest input there can be. Positions 52-99 never saw a gradient. |
| d_model / heads / layers / ff | 128 / 4 / 4 / 512 | **96 / 4 / 4 / 384** | swept; head width stays >= 16 at every size swept, so 4 heads throughout |
| dropout | 0.1 | **0.2** | swept; at two parameters per token this matters more than depth |
| learning rate | 3e-4 | **1e-3** | swept; 3e-4 is a GPT-2-scale value, and at this width the usual choice is nearer 1e-3 |
| weight decay | 0.01 | **0.1, matrices only** | it was decaying LayerNorm gains, biases and the 23-row embedding |
| scheduler | cosine to 1e-5 | **+ 100-step warmup** | `T_max` was already right, in steps; warmup matters once the rate rises |
| batch size | 16 | **64** | measured fastest per epoch, and a quieter gradient |
| epochs | 100 fixed | **200 cap, stop after 10 idle** | best validation lands near epoch 100-170; a 60-epoch cap cut every run short while it was still improving |
| gradient clipping | 1.0 | **1.0** | sound |
| seed | 0 | **0, 1, 2** | differences here are hundredths of a nat |

### What the sweep measured

27 trainings, 3 seeds each, 200-epoch cap, 99 minutes on this Mac
(`sweep.py`; `runs/sweep/runs.tsv` has every row).

Capacity, at dropout 0.2 and lr 1e-3, mean validation loss over 3 seeds:

| d_model / layers / ff | parameters | val loss |
|---|---|---|
| 64 / 4 / 256 | 206,295 | 2.2200 +/- 0.0402 |
| **96 / 4 / 384** | **456,887** | **2.1114 +/- 0.0085** |
| 128 / 2 / 512 | 409,239 | 2.1641 +/- 0.0087 |
| 128 / 3 / 512 | 607,511 | 2.1277 +/- 0.0158 |
| 128 / 4 / 512 | 805,783 | 2.0939 +/- 0.0225 |

Depth beats width: 96 with 4 layers beats 128 with 2 layers on fewer parameters.
The largest shape has the lowest mean, but only by 0.0175 nats, inside the pooled
seed spread of 0.0192 -- so the rule in `choose` takes the smaller model, and the
choice between those two is genuinely unresolved at 3 seeds.

Dropout 0.1 / 0.2 / 0.3 gave 2.1617 / 2.1114 / 2.1494, and learning rate
3e-4 / 1e-3 / 3e-3 gave 2.2025 / 2.1114 / 2.1195: 1e-3 and 3e-3 are a tie inside
seed noise, 3e-4 is clearly worse.

The 60-epoch cap tried first was binding -- every run was still improving when it
hit it. At 200 the best epochs land between 57 and 185 and early stopping, not
the cap, ends each run.

The bigram baseline on the same cluster split is 2.8918, so the transformer buys
**0.78 nats** (perplexity 18.0 down to 8.3).

`train.py` also scores residues and `<eos>` apart. Predicting where a peptide
ends is easier than predicting the next residue -- about 1.5 nats against 2.2 --
and mixing them hides which one a change improved.

### Why the loss is masked

`<pad>` and `<bos>` can never be the next token, so `train.py` sets their logits
to `-inf` before scoring, exactly as `bigram_torch.masked` does. That makes the
transformer's loss the same quantity the bigram baseline reports; otherwise the
model would spend probability mass the baseline never spends and the two numbers
could not be compared. The gather on a padded position is then an infinite loss,
so `ignore_index` has to go with it.

A checkpoint's `lm_head` rows for `<pad>` and `<bos>` are therefore untrained.
Anything that samples from this model must apply the same mask.

## Running it

```
python fetch_peptides.py          # corpus + provenance into data/
python baselines.py               # uniform, unigram, bigram on the same split
python dedup.py --report          # cluster sizes, and the split holds
python train.py --tag default     # one run, best weights into runs/default/
python sweep.py                   # 33 runs, three seeds, runs/sweep/
```

Tests are offline and never touch the network: `pytest pipeline/denovo` in this
directory's venv. The pure-Python parts (`alphabet`, `tokenizer`, `bigram`) also
run in the backend's own `.venv`, which has no torch; everything that imports
torch skips there.
