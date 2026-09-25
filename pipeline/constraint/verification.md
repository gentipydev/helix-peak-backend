# Scoring verification

Three dates. The insulin section below is the original run of 2026-09-14 and is
unchanged; the nine that follow were scored on 2026-09-16, when the catalog
grew from one protein to ten; the ten after them on 2026-09-17, the day they
were added, which is also when the gate changed.

**Insulin was re-scored on 2026-09-16 and reproduced exactly.** Same model, same
revision, same method, and the same regional means to six decimals — 0.411320,
0.851111, 0.116584, 0.825591 — against a generator that had been rewritten
around a target table in between. The asset grew from 34,596 to 35,237 bytes,
all of it the region and disulfide tables the Dart side now reads.

---

# P01308 scoring verification — 2026-09-14

The complete precursor was scored on Apple Silicon MPS in float32 in 8.5 seconds
after model loading. The asset has **110 positions and 34,596 bytes**, with
20 canonical scores per position and exactly zero for every wildtype score.
Every masked input passed the full context/BOS/EOS/token-offset check.

## Cysteine gate: passed

All six cysteines rank in the **top six**, stronger than the predeclared top-11
threshold. Ranking uses unrounded full-vocabulary entropy.

| Precursor residue | Chain residue | Conservation | Rank / 110 |
| --- | --- | ---: | ---: |
| C31 | B7 | 0.999669 | 2 |
| C43 | B19 | 1.000000 | 1 |
| C95 | A6 | 0.997829 | 6 |
| C96 | A7 | 0.999595 | 3 |
| C100 | A11 | 0.999266 | 4 |
| C109 | A20 | 0.998725 | 5 |

## Independent CPU audit: passed

`verify_cpu.py` scored all 110 masked positions in 26.7 seconds. It derives
offsets from the tokenizer's special-token map instead of assuming `index + 1`,
masks IDs instead of inserting a text token, uses NumPy/SciPy instead of Torch
for the reduction, and compares raw logit differences to the saved log ratios.

- Largest entropy difference: **0.000020663** nats.
- Largest conservation difference: **0.000008557**.
- Largest difference between rounded substitution scores: **0.001**.
- All six cysteine ranks are identical on CPU and MPS.

The independent run supports the saved data and rules out a token offset or
MPS-specific explanation for the regional pattern below.

## Regional result differs from the proposed expectation

| Region | Precursor residues | Mean conservation |
| --- | --- | ---: |
| Signal peptide | 1–24 | 0.411320 |
| B chain | 25–54 | 0.851111 |
| C-peptide | 57–87 | 0.116584 |
| A chain | 90–110 | 0.825591 |

The mature chains are much more constrained and the cysteines are the highest
spikes. **C-peptide is less constrained than the signal peptide**, not
intermediate between the signal peptide and mature chains. There is no broad
rise across the signal peptide's end: positions 21–24 are 0.196, 0.500, 0.109,
and 0.659. B and A also contain individual tolerant positions.

These are model measurements, not an observed alignment or proof that a
residue never varies. No scores were adjusted to create the expected pattern.


---

# The other nine — 2026-09-16

All nine scored with the same model, revision, method and normalization.
Dystrophin is the only one not scored in a single pass; see the context column.

| Protein | Residues | Context | Bytes | Disulfides | Bonded ranks |
| --- | ---: | --- | ---: | ---: | --- |
| Oxytocin | 125 | full | 39,917 | 8 | 1–14, 18, 20 of 125 |
| Hemoglobin β | 147 | full | 46,022 | 0 | — |
| Lysozyme | 148 | full | 47,302 | 4 | 1, 3, 4, 5, 7, 10, 11, 12 of 148 |
| Myoglobin | 154 | full | 48,535 | 0 | — |
| Relaxin | 185 | full | 57,556 | 3 | 1, 2, 3, 4, 5, 6 of 185 |
| Somatotropin | 217 | full | 67,145 | 2 | 1, 2, 4, **28** of 217 |
| Ubiquitin | 229 | full | 74,350 | 0 | — |
| p53 | 393 | full | 121,632 | 0 | — |
| Dystrophin | 3,685 | windows of 1,022 | 1,144,355 | 0 | — |

## The gate had to be rewritten twice, and both times the data was right

The original gate was insulin's: every disulfide cysteine in the top decile.
Neither failure it produced was a scoring fault.

**Oxytocin cannot satisfy it arithmetically.** Its precursor has sixteen bonded
cysteines in 125 residues, and a decile has twelve places. What the scores
actually say is stronger than the gate was asking: fourteen of the sixteen are
the fourteen most constrained residues in the whole precursor, in order, and the
other two rank 18 and 20.

**Growth hormone's C215 ranks 28/217** while C79, C191 and C208 rank 4, 1 and 2.
That is the small C-terminal loop, two residues from the end of the chain, being
genuinely more tolerant than the bridge that pins the helix bundle — conservation
0.732 against 0.977, 1.000 and 0.981. An offset bug scrambles every position at
once, so three cysteines at ranks 1, 2 and 4 rule one out by themselves.

The gate now asks that the **median** bonded cysteine rank in the top decile and
that none fall past the halfway mark. Insulin's six still answer to its original
eleven; nothing was relaxed to make a particular protein pass. (A day later the
second ten broke this form too, and the gate stopped reading offsets off the
cysteines at all; see the last section.)

## Regional results

Conservation is min-max normalised within one protein, so a number here is a
rank inside that molecule and is not comparable across the table.

**Relaxin repeats insulin's pattern.** Signal peptide 0.164, B chain 0.443,
C-peptide **0.085**, A chain 0.462. The connecting peptide is again the least
constrained region of the precursor, below the signal peptide — the same
ordering insulin measured and the same one that contradicted the expectation
recorded above. Two proteins built the same way and cut the same way now
measure the same way.

**Ubiquitin is the most constrained protein in the set, and scores itself three
times.** The precursor is three tandem copies of the same 76 residues, scored
independently in three different contexts: 0.992, 0.993, 0.988. Three
near-identical answers to the same question asked three times is a check on the
method that no single protein can give. The 229th residue, a lone C-terminal
cysteine that is cut off and discarded, scores 0.000 — the least constrained
position in the protein, and the only one that is not part of a tag.

**p53's order is its disorder.** DNA-binding domain 0.508, proline-rich region
0.350, tetramerisation 0.381, transactivation domain **0.092**, and the
unannotated stretches between them 0.157 and 0.057. The part that grips DNA is
the part the model is sure about; the ends that have no structure are the ends
it is not.

**Dystrophin's domains stand out of its linkers.** Calponin-homology 1 0.590 and
CH2 0.537, the twenty-four spectrin repeats around 0.435, and the sequence
between named domains 0.086–0.145. Averaged over 3,685 residues that separation
is the clearest single result in the table, and it is one nobody chose: the
region boundaries are UniProt's and the scores never saw them.

**Both globins put the haem pocket at the top, without being told there is
one.** Haemoglobin's five most constrained residues are H64, F43, H93, Y146 and
D100; myoglobin's are H94, H65, F44, H98 and F34. UniProt annotates H64 and H93
as haemoglobin's distal and proximal haem-binding residues, and H65 and H94 as
myoglobin's — the top two in each case, and nothing in the scoring saw those
annotations. Two proteins that share a fold and 25% of their sequence agree on
which residues the fold is for.

Both also make their initiator methionine unremarkable — rank 25/147 and 6/154 —
which is the right answer for a residue whose identity is fixed by the genetic
code rather than by the molecule, and which is cut off it.

These are model measurements, not observed alignments or proof that a residue
never varies. No scores were adjusted in any protein.


---

# The second ten — 2026-09-17

All ten scored the day they were added, with the same model, revision, method
and normalization as the first ten. CFTR is the second protein, after
dystrophin, too long for one pass, and is scored in windows the same way.

| Protein | Residues | Context | Bytes | Disulfides | Bonded ranks | Alignment, before / after |
| --- | ---: | --- | ---: | ---: | --- | --- |
| SOD1 | 154 | full | 48,469 | 1 | **46, 51** of 154 | 94.5% / 91.7% |
| Vasopressin | 164 | full | 52,126 | 8 | 1–15, 54 of 164 | 85.4% / 80.1% |
| Leptin | 167 | full | 51,527 | 1 | **24, 33** of 167 | 71.2% / 71.2% |
| Glucagon | 180 | full | 56,801 | 0 | — | 75.1% / 75.1% |
| Erythropoietin | 193 | full | 59,943 | 2 | 1, 4, 5, 91 of 193 | 76.8% / 80.2% |
| TNF-alpha | 233 | full | 72,601 | 1 | 1, 2 of 233 | 83.3% / 80.5% |
| Prion protein | 253 | full | 79,572 | 1 | 1, 5 of 253 | 82.7% / 87.4% |
| Amylase | 511 | full | 161,352 | 5 | 1, 4, 5, 7, 12, 14, 33, 41, **265, 449** of 511 | 95.7% / 95.9% |
| APP | 770 | full | 240,047 | 9 | 1–9, 12, 13, 15, 17, 19, 28, 33, 38, 55 of 770 | 78.7% / 78.1% |
| CFTR | 1,480 | windows of 1,022 | 455,757 | 0 | — | 80.3% / 79.8% |

*Alignment* here is whether each score sits on the residue it was measured for,
not a sequence alignment: wherever a residue differs from the one beside it, how
often the model scores the residue that is there above its neighbour's.

## The gate had to be rewritten a third time, and the data was right again

Under the gate as it stood — the median bonded cysteine in the top decile, none
past halfway — three of the ten would not have been written. A dry run of the
identical loop, writing nothing, found them before any of this was baked.

**SOD1's bridge is not the most constrained thing about it.** C58 and C147 rank
51 and 46 of 154. Four of its ten most constrained residues are histidines
UniProt annotates as its metal ligands — H47 and H49 hold the copper, H72 and
H81 the zinc — and H64, which holds both, ranks 12. The site the enzyme works
with outranks the bridge that stabilises it.

**Leptin's one bridge sits under its helix core.** C117 and C167 rank 33 and 24
of 167, below the leucines and isoleucines that fill its top ten (L154, L34,
I38, L150, L101, L79, L147). Its scale is also skewed: D100 alone sits at 1.0
and M1 next at 0.94, so min-max leaves the chain's mean at 0.195.

**Amylase's bridges are not alike.** Eight of its ten bonded cysteines rank 1–41
of 511, and C399 is the most constrained residue in the protein; C85 and C130,
bonded to each other, rank 449 and 265. Its catalytic proton donor, E248, ranks
second.

None of the three is an offset bug, and the evidence for that was never going to
be the cysteines. A track filed one residue out scores each residue's neighbour
as its own, so the model would prefer the neighbour's residue at most positions;
in these three it prefers the residue that is there at 71–96% of them. The gate
now measures that, on every protein, and refuses a track that does not clear 50%
on both sides. The first ten's shipped tracks, re-checked from their stored
scores rather than re-scored:

| Protein | Alignment, before / after |
| --- | --- |
| Insulin | 85.9% / 80.8% |
| Oxytocin | 92.9% / 94.7% |
| Ubiquitin | 99.5% / 100.0% |
| Lysozyme | 97.2% / 97.2% |
| Hemoglobin β | 92.6% / 91.1% |
| Myoglobin | 95.9% / 93.2% |
| Relaxin | 68.5% / 64.9% |
| Somatotropin | 79.7% / 79.7% |
| p53 | 73.0% / 73.6% |
| Dystrophin | 81.8% / 82.2% |

A run of one amino acid is not counted, since a Q beside a Q says nothing about
which one a score belongs to. The unit tests build a track shifted one residue
each way and the gate refuses both; they also run the gate over every track on
disk. Bonded cysteine ranks are still printed on every bake, and recorded above.

## TNF's bridge was missing from the table

`targets.py` had no disulfide for TNF, though its fold page names one and draws
it from 7JRA's `SSBOND` record. UniProt P01375 has Cys145–Cys177, and it was
added before scoring. The two rank 1 and 2 of 233. TNF's other two cysteines,
C30 in the cytoplasmic tail and C49 in the membrane anchor, rank 166 and 109.

## Regional results

Conservation is min-max normalised within one protein, so a number here is a
rank inside that molecule and is not comparable across the table.

**Vasopressin repeats oxytocin.** Fifteen of its sixteen bonded cysteines are
the fifteen most constrained residues in the precursor; oxytocin's were fourteen
of sixteen. The hormone is 0.913 and neurophysin 2 0.722. Copeptin, which
oxytocin's precursor does not have, is 0.221, and the lone arginine it is cut
at, 0.057.

**Glucagon's hormones stand out of its spacers.** Glucagon 0.491, GLP-1 0.475
and GLP-2 0.397, against the glicentin-related polypeptide at 0.071 and
intervening peptide 1 at 0.030. Intervening peptide 2 is 0.258.

**Erythropoietin's bridges differ at one end.** C34 is its most constrained
residue, and C188 and C56 rank 4 and 5; C60, C56's partner, ranks 91.

**The prion protein's anchor signal is its most tolerant stretch.** Residues
231–253, cut off and replaced by the GPI anchor, average 0.132 against 0.674
for the chain. C214 is the most constrained residue and C179 ranks 5; four
glycines spaced four apart, G119, G123, G127 and G131, fill the rest of its top
six. The five octapeptide repeats, 51–91, do not crowd the top: the best of them
ranks 13.

**APP's domains stand out of its linkers, as dystrophin's do.** Averaged over
UniProt's domain boundaries, which the scores never saw: E1 0.899, E2 0.691 and
the Kunitz inhibitor domain 0.552, against 0.268, 0.113 and 0.057 for the
stretches between them. Seventeen of its eighteen bonded cysteines rank in the
top 38 of 770; the eighteenth, C341, ranks 55. Amyloid-beta 42 itself is 0.298,
and the 48-residue cytoplasmic tail 0.843. The track's own region table is the
row's signal peptide and chain; these means were computed for this page.

**CFTR puts its two ATP-binding loops at the top.** Six of its eight most
constrained residues sit in the two stretches UniProt annotates as binding ATP,
458–465 and 1244–1251: G458 is first of 1,480, G1249 second and K1250 fourth.
Averaged over UniProt's domains, the two nucleotide-binding domains are 0.533
and 0.509 and the two membrane domains 0.413 and 0.436; the disordered R region
between the halves is 0.177, and the stretches between domains 0.054–0.265.
G551, where the gating mutation G551D sits, ranks 14. F508, whose deletion
UniProt calls the most common cystic fibrosis mutation on Caucasian chromosomes,
ranks 369 at 0.610 — but masked marginals score substitutions, and a deletion is
not one. The track's region table is the row's one chain; these means were
computed for this page.

## Independent CPU audit

`verify_cpu.py` re-derived the nine full-length tracks on CPU, alignment
included. CFTR, like dystrophin, is scored in windows and not audited.

| Protein | Entropy | Conservation | Rounded score | Cysteine ranks and alignment |
| --- | ---: | ---: | ---: | --- |
| SOD1 | 0.000009880 | 0.000003779 | 0.001 | identical |
| Vasopressin | 0.000017352 | 0.000007013 | 0.001 | identical |
| Leptin | 0.000024283 | 0.000008539 | 0.001 | identical |
| Glucagon | 0.000007719 | 0.000002859 | 0.001 | identical (no bridges) |
| Erythropoietin | 0.000025969 | 0.000010372 | 0.001 | identical |
| TNF-alpha | 0.000020996 | 0.000007200 | 0.001 | identical |
| Prion protein | 0.000022368 | 0.000008556 | 0.001 | identical |
| Amylase | 0.000031154 | 0.000012111 | 0.001 | identical |
| APP | 0.000031587 | 0.000010966 | 0.001 | identical |

The columns are the largest differences between the CPU run and the MPS asset.

These are model measurements, not observed alignments or proof that a residue
never varies. No scores were adjusted in any protein.
