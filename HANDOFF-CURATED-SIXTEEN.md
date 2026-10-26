# HANDOFF: sixteen more curated proteins, one for each organ the zoom has no specimen of

Prepared 2026-10-25. Starts from app `bde24dd` (level with GitHub) and the backend commit that
last changed this file (`git log -1 -- HANDOFF-CURATED-SIXTEEN.md`; the backend is ahead of
GitHub only by this file's commits). Read both CLAUDE.md files,
`helix-peek/docs/protein-pipeline-rules.md` and `pipeline/README.md` ("Adding a protein")
first: that runbook is the whole method, and this file adds only the batch. Paths are the
backend's, except `helix-peek/…`, the app's checkout beside it. Under each phase, write a short
**What actually happened** when it is done.

Until this batch this file was `HANDOFF-CURATED-TEN.md`, the third ten's handoff; that record is
in git at `3907aed`. The checks of the third ten and the fourth batch are in
`helix-peek/docs/protein-verification.md`.

## Decided by the user (2026-10-25)

1. **The sixteen below**, one for each `OrganKind` with no specimen under
   `helix-peek/lib/features/zoom/presentation/scenes/organs/` (20 of the 36 have one), so each
   of those organs has a protein to show once its specimen is drawn.
2. **A red gate stops the session.** When a protein fails a gate (the app's fold check, AVI's
   biology check, ESM-2's alignment gate, a page limit, the resolver or the record build), or
   its locus bake lands it outside its organ, ask. The user chooses each time: a swap, shipping
   it unscored (as phospholamban), or dropping the organ.
3. **This session executes, live writes included**, after one stop: the user reviews the
   sixteen's words before anything is live (Phase 4). `--dry-run` every seed and upload
   first, and write only with every gate green. Never push. Ask before a migration, a deploy,
   money or secrets.
4. Standing since 2026-10-17: the curated list grows. A curated protein's files, and its
   `catalog.json` row after the others, may be added to the app's `test/fixtures/`, checked by
   `python3 tool/fixtures_added_only.py <start>`. Curated proteins get a hand-written
   `summary`, structure `sentence` and `semantics`, held to R4.5.

Not the user's decision: orders 38–53 in the table's order, after `lgi1` (37). That is what
`draft_curated.py` gives the genes in this order. The list can be reordered later, as
recoverin was.

## The sixteen

Organ and cell come from the locus bake's `path_of`, run over the Atlas bulk file
(proteinatlas.json.gz) fetched 2026-10-25. The bake reads the Atlas itself, so its `path` is
the one that counts. "Drawn as" is the zoom's `shapeOfCell`, or `shapeOfTissue` where the path
names no cell. The entry is PDB policy 1's choice in a `scripts/dry_fold.py` run of 2026-10-25
(UniProt 2026_03, RCSB search of 2026-10-21), its method and resolution as the fold page will
print them. Residues: placed / in the span drawn.

| order | slug (suggested) | gene | UniProt | organ → cell (drawn as) | RNA specificity | entry | notes |
|---|---|---|---|---|---|---|---|
| 38 | s100p | S100P | P25815 | urinary bladder → urothelial cells (epithelial) | group enriched | 1J55 X-ray 2.00 Å, 88/95 | one chain of two |
| 39 | lactalbumin | LALBA | P00709 | breast → breast lactating cells (acinar) | tissue enriched | 3B0O X-ray 1.61 Å, 119/123 | 4 bridges |
| 40 | slpi | SLPI | P03973 | cervix → no cell, the tissue's own (epithelial) | group enriched | 4DOQ X-ray 2.00 Å, 47/47 | domain only: WAP 2, 85–131; 4 bridges |
| 41 | npc2 | NPC2 | P61916 | epididymis → epididymal principal cells (epithelial) | tissue enriched | 5KWY X-ray 2.40 Å, 132/132 | 3 bridges |
| 42 | fgf19 | FGF19 | O95750 | gallbladder → cholangiocytes (epithelial) | tissue enriched | 2P23 X-ray 1.80 Å, 136/192 | 1 bridge, not UniProt's |
| 43 | fabp6 | FABP6 | P51161 | intestine → colon enteroendocrine cells (endocrine) | tissue enriched | 5L8O X-ray 2.39 Å, 125/127 | initiator methionine removed |
| 44 | surfactantd | SFTPD | P35247 | lung → alveolar cells type 2 (epithelial) | tissue enriched | 1PWB X-ray 1.40 Å, 126/126 | domain only: C-type lectin, 250–375; 2 bridges |
| 45 | padi6 | PADI6 | Q6TGC4 | ovary → oocytes (germ) | tissue enriched | 25YL cryo-EM 2.8 Å, 685/694 | one chain of two |
| 46 | casr | CASR | P41180 | parathyroid gland → no cell, the tissue's own (endocrine) | tissue enriched | 7M3G cryo-EM 2.5 Å, 804/1,059 | 9 bridges; 1,078 aa |
| 47 | acp3 | ACP3 | P15309 | prostate → prostate glandular cells (acinar) | tissue enriched | 8XJ4 cryo-EM 3.2 Å, 345/354 | one chain of two; 2 bridges |
| 48 | mpges1 | PTGES | O14684 | seminal vesicle → no cell, the tissue's own (acinar) | tissue enhanced | 5TL9 X-ray 1.20 Å, 148/152 | one chain of three |
| 49 | loxl2 | LOXL2 | Q9Y4K0 | smooth muscle → smooth muscle cells (spindle) | tissue enhanced | 5ZE3 X-ray 2.40 Å, 233/233 | domain only: SRCR 3–4, 322–554; 5 bridges |
| 50 | ckmt2 | CKMT2 | P17540 | tongue → myonuclei (myofibre) | group enriched | 4Z9M X-ray 2.10 Å, 369/380 | transit peptide 1–39 |
| 51 | glycodelin | PAEP | P09466 | endometrium (uterus) → endometrial glandular cells (acinar) | tissue enriched | 4R0B X-ray 2.45 Å, 153/162 | one chain of two; 2 bridges |
| 52 | alox12 | ALOX12 | P18054 | vagina → no cell, the tissue's own (epithelial) | tissue enhanced | 8GHE cryo-EM 2.0 Å, 662/663 | one chain of four |
| 53 | gas6 | GAS6 | Q14393 | blood vessel → no cell, the tissue's own (endothelial) | tissue enriched | 5VXZ X-ray 2.30 Å, 369/388 | domain only: laminin G-like 1–2, 288–675; 2 bridges |

Re-run the structure check with
`NCBI_EMAIL=gentiannurce@gmail.com .venv/bin/python scripts/dry_fold.py --out <dir> S100P LALBA SLPI NPC2 FGF19 FABP6 SFTPD PADI6 CASR ACP3 PTGES LOXL2 CKMT2 PAEP ALOX12 GAS6`.
The 2026-10-25 run's files, while this Mac keeps them, are in
`/private/tmp/claude-501/-Users-gentiannurce-Desktop-HelixPeek/8b7eb047-13f2-4c76-8992-68e49befd8f9/scratchpad/`:
`folds/report.md` and `folds/report.json`, and the Atlas scan, `organ_scan.json`. Draft with:

```sh
NCBI_EMAIL=gentiannurce@gmail.com .venv/bin/python scripts/draft_curated.py --out <dir> \
    S100P=s100p@1J55 LALBA=lactalbumin@3B0O SLPI=slpi@4DOQ NPC2=npc2@5KWY \
    FGF19=fgf19@2P23 FABP6=fabp6@5L8O SFTPD=surfactantd@1PWB PADI6=padi6@25YL \
    CASR=casr@7M3G ACP3=acp3@8XJ4 PTGES=mpges1@5TL9 LOXL2=loxl2@5ZE3 \
    CKMT2=ckmt2@4Z9M PAEP=glycodelin@4R0B ALOX12=alox12@8GHE GAS6=gas6@5VXZ
```

Alternates that passed the same dry run (swap only with the user):

| organ | alternate | entry | lands in |
|---|---|---|---|
| urinary bladder | MMP13 | 4FVL X-ray 2.44 Å, 368/368 | urothelial cells |
| urinary bladder | IL23A | 5MJ3 X-ray 1.74 Å, 155/170 | urothelial cells |
| cervix | S100A7 | 1PSR X-ray 1.05 Å, 100/100 | no cell |
| lung | AGER | 4P2Y X-ray 2.30 Å, 299/382 | alveolar cells type 1 |
| prostate | MSMB | 3IX0 X-ray 2.30 Å, 94/94 | prostate glandular cells |
| endometrium | PGR | 1SQN X-ray 1.45 Å, ligand-binding domain 683–923 | endometrial stromal cells |
| endometrium | ESR1 | 1XP9 X-ray 1.80 Å, ligand-binding domain 307–551 | endometrial stromal cells |
| tongue | S100A1 | 5K89 X-ray 2.25 Å, 85/94 | no cell |
| vagina | EDN3 | 6IGK X-ray 2.00 Å, endothelin-3 alone, 21 of 238 | no cell |
| smooth muscle | MYLK | 6C6M X-ray 2.50 Å, one Ig-like domain, 95 of 1,914 | smooth muscle cells |

None passed for breast, epididymis, gallbladder, intestine, ovary, parathyroid gland, seminal
vesicle or blood vessel. Refused, for any organ:

- PTH (no entry passes, and its propeptide is 6 residues) and KLK3 (a 7-residue propeptide):
  autotaxin's limit.
- ACTN1: six cell locations, and no buildable index entry.
- IDH2: the resolver refuses it (its MANE mRNA does not hold its protein).
- SPINK2, ANXA4 and SFTPC: no buildable index entry.
- No entry passes for EDNRA, ANO1, SERPINB4, MFGE8, COMP, BGN, APOA4, DEFA5, FSHR, INHA, VDR,
  CCKAR, UPK1A, WFDC2, STC2 or CSN2.
- SEMG1 and SEMG2 were not tried.

None of the sixteen or the alternates was built on demand on 2026-10-25 (`/proteins/built`
lists eleven proteins), and `/catalog` served 38 rows. Check `/proteins/built` again before
seeding: the seed refuses a row over a protein built on demand.

## What the ten and the eight met

- **The seed comes before the uploads:** `protein_track` references `protein`.
- **Three gates refused once each,** and the user chose:
  - UMOD's planned 6TQL failed the app's fold check (a step of 1.16 where 1.1 is allowed); it
    draws policy 1's next entry, 9NU1. `draft_curated.py` runs that check, so it shows at the
    draft.
  - PF4 failed AVI's `check_biology` (exons 7.4 = introns 7.4); PPBP replaced it.
  - PLN failed ESM-2's alignment gate (40.0% / 46.7%, against the more than 50% it needs) and
    ships unscored (`scored=False`).

  `bake_curated.py` runs impact after constraint, so a `check_biology` refusal arrives after
  ESM-2's time.
- **Two page limits, still open.** A region above `regionFloor` but under about 30 bases that
  appears once cannot be named (autotaxin's 8-residue propeptide), and the names round the cell
  fit three or four Atlas locations, not six (cerebellin-1). The sixteen were picked inside both:
  every region at least 10 residues, at most two locations.
- **The search tests pin words.** `search_surface_test` counts the rows the phone shows for
  "ins", so a summary carrying it ("joins", "against") pushes a suggestion off the screen. Keep
  the search tests' pinned words out of summaries.
- **A row lists the bridges its fold draws.** LGI1's row carries 9KZT's five, which UniProt does
  not annotate (`targets.py`). The fourth batch's Flag 4 in `protein-verification.md` says the
  row lists none, which is wrong: correct it in this batch's section. The stored model must be
  the one the app's check passed: taking LGI1's bridges out changed its model, and the check
  refused it.
- **A failing walk test is a row to fix.** The ten needed six walk-test edits, each approved by
  the user (among them an unscored protein and a grouped mature count over a thousand); the
  eight needed none.

## To watch in this batch

Each goes into the verification section, with what was done about it.

1. **Sixteen organs, each reached for the first time.** None has a specimen, so the zoom draws
   each as the body figure's outline grown to size, with its plain surface (`organ_art.dart`,
   `_lookOf`); the blood vessel and the parathyroid have drawings of their own. The organ test
   already passes such an organ (the ten's approved zoom-test edit). Cervix, endometrium, ovary
   and vagina are the female figure's; epididymis, prostate and seminal vesicle the male's. Six
   tissue recipes meet a curated protein for the first time: urothelial, intestinal, alveolar,
   ovarian, smooth muscle and vascular (`tissue_slide_test` draws every recipe already).
2. **Five land in no cell:** SLPI, CASR, PTGES, ALOX12 and GAS6. `path_of` finds no cell type of
   theirs that lives in the organ, so the zoom draws the tissue's own cells and says so, as it
   does for prion and TNF.
3. **FGF19's bridges.** 2P23 draws one bridge, 102–120. UniProt's 58–70 is 2.69 Å apart in it
   and is not drawn. The row's disulfides, its facts and the fold must say the same.
4. **LOXL2's bridges.** UniProt's 351–414 lies in the span drawn, but 5ZE3 does not place it.
   5ZE3 draws 364–424, 395–405, 464–530, 477–543 and 511–521.
5. **Residues the entries differ on**, each within what the fidelity gate allows: 5ZE3 has Q at
   455, 4R0B has E at 46 and 8GHE has S at 322, where each record has N. 5TL9's Cys59 is
   modified (CSO, S-hydroxycysteine).
6. **Leaders.** CKMT2's first 39 residues are a mitochondrial transit peptide, the first leader
   in a curated row that is neither a signal peptide nor a methionine. The resolver names it
   `Transit peptide (…)`, short `T` (`resolve._label`). R3.1 and R3.2 are written for a signal
   peptide, so check what the protein and mature pages call it, and record it. FABP6's chain
   starts at 2: its initiator methionine is a region of its own, as SOD1's is (R3.7).
7. **Four folds draw a domain**, as p53's, APP's and uromodulin's do: SLPI's second WAP domain,
   SFTPD's carbohydrate-recognition domain, LOXL2's SRCR 3–4 and GAS6's two laminin G-like
   domains. The structure sentence says what is not drawn.
8. **Introns shortened.** The budget holds the MANE transcript's span, exons and introns
   (`build_gene_record.intron_scale`). By Ensembl's canonical (MANE Select) spans, six are over
   24,000 bases: CASR (107,962), LOXL2 (107,224), GAS6 (43,528), ACP3 (41,435), CKMT2 (33,077)
   and PADI6 (29,504). The next longest are ALOX12 (14,706) and PTGES (14,694). The record
   decides (`GeneRecord.isIntronCompressed`). The app CLAUDE.md's compressed-intron section names
   each gene that joins: fifteen genes, if all six do.
9. **Tissue enhanced, not enriched.** That is what the Atlas says for PTGES, LOXL2 and ALOX12
   (the list given on 2026-10-25 said enriched). `path_of` reads the organ from the Atlas's
   tissue specificity either way.
10. **ESM-2's longest** is CASR, at 1,078 aa. TG's 2,768 took 35 min and PAPPA's 1,627 took 20.

## Phases

1. **Draft and verify.** Run `draft_curated.py` as above. It resolves each gene, chooses and
   bakes its fold, holds the fold to the app's check, and stores the working copy. A gene that
   does not draft still takes its order, so a swap is drafted with `--order <n>`.
   - Check each gene against its record in a new section of `protein-verification.md`, "The
     fifth batch: sixteen organs": span, exons, neighbouring genes, processing (which regions
     are kept, `mature_peptides`), structure, and the list above.
   - Settle each `# HAND:` (R3.3, R3.4).
   - Add the rows to `targets.py` and `curated/catalog.json` at orders 38–53.
   - Write every `DRAFT:` word by hand (R4.5), in the register of the 38 rows already written:
     data first, no field basics (memory `domain-expert-audience`).
   - Then `scripts/resolver_worker.sh status`, and `restart`: the worker runs this tree, and
     its scorer imports `targets.py`.

   **What actually happened (2026-10-25).** The user chose, before the session began, that a
   red gate holds its protein and every refusal comes to the one stop with the words, rather
   than stopping the session at each.
   - **All sixteen drafted on their planned entries** (16 of 16, about 15 minutes), and every
     fold passed the app's check. For CASR, 7M3G's chain B, ranked first, narrowly fails it
     (a ring 1.518 from its own stretch of chain, where 1.5 is allowed; PAPP-A's chain D
     case), and the same entry's chain A, which places 821 residues to B's 804, passes.
   - **Every translation equals its UniProt sequence**, so no row declares
     `uniprot_variants`. Each RefSeqGene holds one coding sequence for its gene, at MANE's
     version, so those rows name it unversioned, as the eight's do.
   - **Two records end a leader where UniProt does not.** FGF19's names its signal peptide 1–22
     and its chain 23–216 (UniProt: 1–24); SP-D's names a signal peptide 1–17 and no chain
     (UniProt: 1–20). In each of the 38 rows the table's signal peptide is the record's, so both
     tables follow their records (R3.2's source order), and SP-D's chain, 18–375, is filled
     from the table as erythropoietin's is (R3.9). GAS6's record also names its signal peptide
     alone, and its chain, 31–678, is filled the same way.
   - **Glycodelin stops at the precursor.** Its slice carries isoform 2's chain inside the
     MANE transcript's exons, over this one's, and the build refuses two chains on one
     stretch: `mature_peptides=False` (R3.4), as for uromodulin, its catalog `chain` named
     `glycodelin`. Every other cut record names one chain after its leader, so the walk cuts
     it (`mature_peptides=True`). α-lactalbumin's slice names its signal peptide and chain
     twice, CKMT2's its transit peptide and chain three times; the record keeps one of each.
   - **Disulfides: UniProt's, with any the fold draws that UniProt lacks.** The handoff's
     premise for FGF19 was wrong: UniProt annotates 102–120 as well as 58–70, and 2P23 forms
     102–120 (58–70's sulfurs are 2.69 Å apart). LOXL2's 351–414, which the handoff says 5ZE3
     does not place, is placed with its sulfurs 1.76 Å apart; the entry only does not record
     the bond. CASR's chain A forms 677–765, which UniProt does not annotate, and leaves 437–449
     open (3.59 Å), which its chain B forms. Dropping a UniProt pair for one chain's geometry
     would make a card contradict the knowledge base, so each row keeps UniProt's pairs and
     adds the entry's (CASR 9, as LGI1's five were added), and the fold's words name a pair
     the fold leaves out. A pair can still be dropped later without scoring again
     (`score_protein.py --metadata-only`).
   - **The words** were written by hand and held to the search tests' pinned queries by a
     port of `rank` (`protein_catalog_test`, a walk test, pins `protein`, `hormone`,
     `precursor`, `TP53` and `titin`; the search tests pin `ins`, `lysozyme`, `kinase` and
     others). So SFTPD is displayed `SP-D`, and no summary says "protein", "hormone", "ins"
     (chains, domains, proteins), "lysozyme" or "kinase".
   - `restart`ed the worker on the new tree (queue empty).
2. **Bake.** For each:
   `NCBI_EMAIL=gentiannurce@gmail.com .venv/bin/python pipeline/bake_curated.py --target <slug>`
   (`--dry-run` prints the commands, `--from <step>` resumes). Its last step checks that row
   alone. Then the gates, by hand:
   - `check_assets.py`, `folding/check_folding.py` and `locus/check_locus.py` with no
     `--target`: all 54 rows.
   - Each locus lands in its planned organ (decision 2).
   - `.venv/bin/pytest` green offline.
   - `fetch_tracks.py --check`: every other row's files as stored.

   At a red gate, stop and ask (decision 2). Commit the backend.

   **What actually happened (2026-10-25).** One protein at a time, every step but `check`
   (the words were written meanwhile, and only `check` reads them), then the checks over every
   row. 52 minutes for the sixteen; CASR's ESM-2 (1,078 residues, in windows) took 13 of them.
   - **Fifteen passed every gate.** ESM-2's alignment ran 76.8% (PADI6) to 94.4% (NPC2)
     against the 50% it asks; every AVI biology check but one passed; each fold passed the
     app's check a second time (`bundle`); each locus landed in its planned organ and cell.
   - **PTGES failed AVI's biology check** (its exons score 5.1, its intron interiors 5.3),
     though every one of its 14,694 bases agrees with GRCh38: PF4's case. ClinVar falls with it.
     Held. The handoff had no alternate for the seminal vesicle, so the Atlas scan's other
     candidates there were tried: a structure passes policy 1 for ANGPT1 (9IVT), CHEK1 (2YEX)
     and TES (2XQN), and AVI's own gates pass all three (`pipeline.resolver.impact_local` on the
     dry run's resolution, in a scratch directory; exons 13.0, 10.9 and 12.6 against introns
     5.2, 4.7 and 6.3). ANGPT1 was drafted at order 48, checked against its record, worded
     and baked whole (alignment 79.7% / 79.0%; 386 ClinVar records), in the tree in mPGES-1's
     place for the user to accept or not; mPGES-1's row and words are kept aside.
   - Glycodelin's protein page called what its signal peptide leaves "proprotein" (R3.3's
     last resort); the row now names it (`proprotein="glycodelin"`), as PAPP-A's, and the row
     was baked again.
   - Gates: `check_assets.py` 54 rows (53 scored, 54 with AVI), `check_folding.py` 54, 0
     problems, `check_locus.py` 54, 0 problems; pytest 906 passed, 96 skipped;
     `fetch_tracks.py --check` 266 files as stored, 0 differ. The numpy warnings
     `check_folding.py` prints are there for the 38 alone too.
3. **Try the app on them.** `pipeline/sync_fixtures.py --from local --target <slug>…` writes
   the rows and their tracks into the app's fixtures, as they will be served.
   - Add the sixteen to the zoom's own tables, as the eight did. In `zoom_domain_test`: the baked
     paths (a protein in no cell is tabled as prion's is) and the cell shapes (the table's
     "drawn as"). In `cell_scene_test`: the cells drawn part by part, which are the endocrine
     cells and the neurons, so FABP6 and CASR join. Any other zoom-test change goes to the user,
     as the ten's did.
   - Run `flutter analyze`, `flutter test`, and the three guards in the app CLAUDE.md against
     the start commit. A walk-test edit needs the user's word.

   **What actually happened (2026-10-25).** The sixteen's fixtures were written from the bakes
   (`sync_fixtures.py --from local`, 112 files and their `/catalog` rows) and the zoom's tables
   gained them: their baked paths and cell shapes in `zoom_domain_test`, and CASR's and FABP6's
   endocrine cells in `cell_scene_test`'s part-by-part list. `flutter analyze` clean; the three
   guards clean against `bde24dd` (no walk-test line changed, 16 proteins added to the
   fixtures and nothing else, no golden PNG moved). `flutter test`: 2,470 passed, 66 skipped,
   and every walk test passed with the sixteen in it. Five failed:
   - **Four zoom checks**, each zoom code meeting a case none of the 38 reaches:
     - Chromosomes 13 and 14 (GAS6, NPC2): their short arms are 15–16% of the chromosome
       (21's is 26%), too short for the waist and end cap, whose sizes are fixed. The stalk at
       p12 is drawn at 0.75 of the arm's width, where 0.9 is asked, and the cleft is not yet
       open at the arm's top (`chromosome_scene_test`, `chromosome_look_test`).
     - The vessel's and the parathyroid's own drawings (GAS6, CASR) have no outline, so the
       place their tissue is sampled lies outside it (`zoom_anatomy_test`).
     - LOXL2: in a smooth muscle cell's flat nucleus, chromosome 8's glowing territory lies over
       the nucleolus (`nucleolus_test`).
   - `test/features/zoom/_scratch_shots_test.dart`, a git-excluded screenshot tool of an earlier
     session that only loads with `SHOTS` set. Not this batch's; left as it was.

   Walking the sixteen's pages (a throwaway test over `AnatomyModel.derive`) found glycodelin's
   "proprotein" (Phase 2), and CKMT2's transit peptide drawn as the N-terminal extension on the
   gene page, with a protein page of one block and a mature page of its 380 residues. The
   sixteen add 75 catalog golden pages.
4. **Stop for the user.** Show them the sixteen's words (summary, structure sentence and
   semantics) and whatever Phases 1–3 left to them. Nothing goes live before their yes.

   **What actually happened (2026-10-25).** One stop, with every word, the choices made by rule
   (the bridges, the records' leaders, glycodelin, CKMT2's transit peptide, `SP-D`) and the two
   open questions. The user's answers:
   - **Angiopoietin-1 for the seminal vesicle**, in mPGES-1's place at order 48.
   - **The four zoom checks fixed generically in the zoom's code**, each change acting only
     where these cases occur, so nothing the 38 draw changes, and shown to the user before the
     app's commit.
   - **The words approved as written.**
5. **Live.** Curl `/health` first: Render sleeps.
   - `seed_catalog.py --target …`, a dry run, then for real.
   - Each kind with `upload_tracks.py --kind <kind> --target …`, a dry run first: record;
     structure (scenes compiled; `--against` the service on the dry run); folding; constraint
     (where scored); impact; clinvar; locus.
   - Then check:
     - `seed_catalog.py --check`: 54 rows, 0 differences.
     - `check_assets.py --against https://helix-peak-backend.onrender.com`: 54 rows agree.
     - `/catalog`: 54 rows in order.
     - `/proteins/resolve/<GENE>` reads `ready` at its slug.
     - `fetch_tracks.py --check`.
     - `LIVE_BACKEND=https://helix-peak-backend.onrender.com flutter test test/live_backend_check.dart`.

   **What actually happened (2026-10-25).** After the backend's commit `0ff87ca`, with no
   protein of the sixteen built on demand and the worker's queue empty:
   - **The seed.** Its dry run, then for real: 16 proteins, 92 aliases, 96 track rows.
   - **The uploads.** Every kind, a dry run first: record; structure, whose dry run against
     the service compiled every scene, each model and scene by its digest the one the draft's
     app check passed (32 of 32); folding; constraint; impact; clinvar; locus. 16 objects of
     each, and 32 for the structures.
   - **The checks.**
     - `seed_catalog.py --check`: 54 rows, 0 differences.
     - `check_assets.py --against` the service: 54 rows agree.
     - `/catalog`: 54 rows, orders 0 to 53.
     - `/proteins/resolve/<GENE>`: `ready` at its slug for all sixteen.
     - `fetch_tracks.py --check`: 378 files as stored, 0 differ.
     - `live_backend_check.dart` passes.
6. **App.** `pipeline/sync_fixtures.py --check --target …` from the service must find the local
   run's files unchanged. Commit the app: the fixtures, the zoom tables and the docs (this
   batch's verification section, the compressed-intron section). List the new golden PNGs for
   the x86 host; the ten owed 44 and the eight 34.

   **What actually happened (2026-10-25/26).**
   - `sync_fixtures.py --check` from the service: 113 fixtures, 0 differ.
   - **The four zoom checks** were fixed in the zoom's code, the user's choice:
     - the vessel's and the parathyroid's drawings given their outlines;
     - on a short arm under 3.2 half-widths, the waist and that end's rounding drawn down
       in proportion;
     - a followed territory that would lie over the nucleolus turned a quarter of the way
       round, away from it.

     A dump of every chromatid width and cleft, and every territory and nucleolus, is the
     same to the float for the thirty-eight before and after. The user saw the diff before
     the commit.
   - **The app's checks.** `flutter test` 2,472 passed; its one failure is the git-excluded
     screenshot tool that needs `SHOTS`. `flutter analyze` is clean, and the three guards are
     clean against `bde24dd` (16 proteins and 80 files added to `test/fixtures`, nothing else).
   - **The app's commits**, local:
     - `93a4587`: the zoom's three cases;
     - `7c259ce`: the fixtures and the zoom's tables;
     - `7de0d9e`: the docs, with this batch's verification section and the compressed-intron
       list, now sixteen genes.
   - **Every one of the zoom's 36 organs** now has a curated protein that reaches it.
   - **Golden pages.** 75 new catalog golden pages are owed to the x86 host
     (`test/goldens/catalog/<slug>-1-gene.png` and `-2.png` to `-5.png`). The ten's 44 and the
     eight's 34 are still owed too.
7. **Report** to the user, and update the memory notes: `unbuilt-organ-proteins`,
   `live-db-migration-state`, `zoom-cell-coverage`, `goldens-need-x86-host`,
   `protein-pipeline-goal`, `local-main-ahead-of-origin`.

## Environment

- Every bake environment and command is in `pipeline/README.md`.
- `NCBI_EMAIL=gentiannurce@gmail.com` goes in the command's environment only.
- `ALPHAGENOME_API_KEY` is in `~/.env`. AlphaGenome meters requests per minute: the impact
  bake cuts its windows to 20 kb, paces them and retries (`bake_impact.CHUNK_BP`), so bake one
  protein at a time.
- Never install into `pipeline/.esm-venv`.
- The resolver worker runs this backend tree against production. Read its queue first
  (`scripts/resolver_worker.sh status`), and `restart` it after changing `pipeline/`.
- zsh does not split a joined string, so pass slugs as an array: `args+=(--target "$s")`.

**Out of scope:**

- The sixteen organs' specimens, drawn from anatomy as the other twenty's were. That is zoom
  work the user asks for by name, organ by organ; until then each organ is the figure's outline
  with its plain surface.
- The list's order.
- The two open page limits.
- Pushing.
