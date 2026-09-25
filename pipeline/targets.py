"""The twenty curated proteins, and everything the bakers need to know.

One row per protein, read by every tool under `pipeline/`: the gene record
builder, the ESM-2 constraint scorer, the impact and ClinVar bakes and the
structure bake. The rows' hand-written prose -- what `protein_catalog.dart` held
on the Dart side until Phase 3 of HANDOFF-ONDEMAND.md -- is `curated/catalog.json`;
the two are checked against each other by `pipeline/check_assets.py`.

Provenance. Gene sources are RefSeqGene accessions where one exists, and a
`NC_*` chromosome slice where none does (OXT and RLN2 have no RefSeqGene).
Region boundaries and disulfide pairs are UniProt features in precursor
numbering, read from `rest.uniprot.org` and pinned here so a bake is
reproducible without the network; `--verify-uniprot` re-checks them against the
live entry. The two exceptions are noted where they sit.
"""

from __future__ import annotations

from dataclasses import dataclass, field


# The gene page is sized to one screen by `AnatomyLayout.fit`, which shrinks
# cells to a 2pt floor and stops being readable past about this. A gene over
# budget has its introns scaled; see `pipeline/mock/build_gene_record.py`.
GENE_PAGE_BUDGET_BP = 24_000

# ESM-2's trained context. A longer protein is scored in windows of this many
# residues centred on the mask rather than in one pass, so no score is taken
# from a position the model never saw in training. DMD and CFTR need it.
ESM_CONTEXT_RESIDUES = 1022


@dataclass(frozen=True)
class Region:
    """A named stretch of the precursor, 1-based inclusive, UniProt numbering.

    `short` is what a disulfide partner is named by: insulin's A chain is `A`,
    so its Cys96 reads `A7`. A mature protein that is all one chain leaves it
    empty and its partners are named by position alone, which is how the
    lysozyme and somatotropin literature numbers them.
    """

    label: str
    short: str
    start: int
    end: int
    # Whether this stretch is still there once the precursor has been cut.
    # Read by the panel, which has something different to say about a residue
    # in a piece that gets thrown away: tolerance there is not the same fact
    # about the molecule as tolerance in the part that goes on to work.
    kept: bool = True


@dataclass(frozen=True)
class Source:
    """Where the GenBank record comes from. A slice re-bases coordinates to 1."""

    accession: str
    seq_start: int | None = None
    seq_stop: int | None = None
    # Which transcript, where the record holds more than one. A RefSeqGene
    # usually holds exactly one and leaves this None; a chromosome slice holds
    # every model RefSeq has, predicted `XP_` ones included.
    protein_id: str | None = None
    # Which mRNA, for the case `protein_id` cannot settle: two transcripts that
    # differ only outside the coding sequence, and so share one CDS. AMY1A's
    # slice carries NM_001008221.1 and NM_004038.4 around identical CDS
    # coordinates, and only the second is MANE Select. None leaves the choice
    # to `select_isoform`, which takes the shortest mRNA holding the CDS.
    transcript_id: str | None = None


@dataclass(frozen=True)
class Chain:
    """One node of the baked `.glb`, cut from one chain of the PDB entry.

    The node names are the contract with `structure_view.dart`, which looks
    them up to give each one a material. What colour that material takes is a
    Dart decision and lives in `protein_catalog.dart`; nothing here needs to
    know, and a colour written down twice is a colour that will disagree.
    """

    node: str
    pdb_chain: str


@dataclass(frozen=True)
class Structure:
    """What the last page of the walk is cut from."""

    pdb: str
    chains: tuple[Chain, ...]
    bonds: bool
    # Restrict the export to one span, in the PDB's own residue numbering, for
    # a protein whose full length has no structure. `None` means whole chains.
    residues: tuple[int, int] | None = None
    # PyMOL's cartoon sampling, and what it costs: insulin's 51 residues come
    # out at 8 as about 394 vertices each. A longer chain drops it to hold the
    # `.fsceneb` inside the asset budget.
    sampling: int = 8
    # `cartoon` everywhere but a peptide too short to have any secondary
    # structure to draw, where a plain tube is the honest shape.
    representation: str = "cartoon"


@dataclass(frozen=True)
class Target:
    slug: str
    gene: str
    uniprot: str
    display: str
    source: Source
    structure: Structure
    # How the space between two named regions is read. A precursor is cut, so
    # the gap is a cleavage site; a single folded chain is not, so the gap is
    # just more of that chain.
    cleaved: bool
    chain_label: str
    regions: tuple[Region, ...] = ()
    disulfides: tuple[tuple[int, int], ...] = ()
    aa: int = 0
    # Positions where the genomic record's own translation differs from the
    # UniProt entry, as (position, record, uniprot). Declared rather than
    # tolerated: the bake fails unless the differences it finds are exactly
    # these, so an NCBI revision cannot slip past unnoticed. The app scores and
    # draws the record's alleles, because that is what its gene page shows.
    uniprot_variants: tuple[tuple[int, str, str], ...] = ()
    # Whether the record's `mat_peptide` features are a cleavage series worth a
    # page of their own. False where they are fragments cut back out of the
    # finished chain rather than pieces the precursor is divided into —
    # haemoglobin's hemorphins, SOD1's antimicrobial peptide — and where they
    # overlap each other as alternatives, glucagon's and APP's, which one page
    # of disjoint chains cannot draw. Such a walk stops at the precursor.
    mature_peptides: bool = True
    # What the signal peptide leaves, for a record that does not annotate it.
    # NG_007114 names proinsulin, so insulin's protein page is drawn as the
    # leader and the proprotein it becomes; oxytocin's and relaxin's chromosome
    # slices name no proprotein, and without this their pages could only call it
    # 'proprotein'. Filled into the record by `build_gene_record.py`.
    proprotein: str | None = None
    # Whether an ESM-2 constraint track is baked for this protein. All twenty
    # are; a row added without one sets this False, which the app draws as a
    # state — a protein page with no conservation toolbar — rather than meeting
    # as a missing file. `score_protein.py --all` passes an unscored row by, and
    # `check_assets.py` fails one that has a track. `scored` in
    # `protein_catalog.dart` says the same thing on the Dart side.
    scored: bool = True
    # Whether an AlphaGenome Variant Impact track is baked for this gene. It is
    # the same kind of flag as `scored` one level down: `scored` is a per-residue
    # ESM-2 track over the protein, this is a per-base AVI track over the gene
    # record. A row without one draws its nucleotide pages exactly as before —
    # a tap moves the tracer and no sheet opens — rather than meeting a missing
    # file. `impact_scored` in `protein_catalog.dart` says the same thing.
    impact_scored: bool = True
    # Whether a ClinVar snapshot is bundled for this gene, the same kind of flag
    # again. All twenty have one; a row added without one sets this False, which
    # the app draws as "not yet included" rather than as a missing file.
    # `clinvarAvailable` in `protein_catalog.dart` says the same thing.
    clinvar_available: bool = True

    @property
    def mock_asset(self) -> str:
        return f"assets/mock/gene_{self.gene.lower()}.json"

    @property
    def constraint_asset(self) -> str:
        return f"assets/constraint/{self.slug}_esm_constraint.json"

    @property
    def structure_asset(self) -> str:
        return f"assets/models/{self.slug}.glb"

    @property
    def impact_asset(self) -> str:
        return f"assets/impact/{self.slug}_avi.json"


def _spectrin() -> tuple[Region, ...]:
    """Dystrophin's 24 spectrin repeats, UniProt `Repeat` features."""
    spans = [
        (339, 447), (448, 556), (559, 667), (719, 828), (830, 934), (943, 1045),
        (1048, 1154), (1157, 1263), (1266, 1367), (1368, 1463), (1468, 1568),
        (1571, 1676), (1679, 1778), (1779, 1874), (1877, 1979), (1992, 2101),
        (2104, 2208), (2211, 2318), (2319, 2423), (2475, 2577), (2580, 2686),
        (2689, 2802), (2808, 2930), (2935, 3040),
    ]
    return tuple(
        Region(f"Spectrin repeat {n}", f"SR{n}", start, end)
        for n, (start, end) in enumerate(spans, 1)
    )


TARGETS: tuple[Target, ...] = (
    Target(
        slug="insulin", gene="INS", uniprot="P01308", aa=110,
        display="Insulin",
        source=Source("NG_007114"),
        cleaved=True, chain_label="Proinsulin",
        regions=(
            Region("Signal peptide", "S", 1, 24, kept=False),
            Region("B chain", "B", 25, 54),
            Region("C-peptide", "C", 57, 87, kept=False),
            Region("A chain", "A", 90, 110),
        ),
        disulfides=((31, 96), (43, 109), (95, 100)),
        structure=Structure(
            pdb="3I40",
            chains=(Chain("chainA", "A"), Chain("chainB", "B")),
            bonds=True,
        ),
    ),
    Target(
        slug="hemoglobin", gene="HBB", uniprot="P68871", aa=147,
        display="Hemoglobin (beta chain)",
        source=Source("NG_059281"),
        cleaved=False, chain_label="Beta chain",
        mature_peptides=False,
        regions=(
            Region("Initiator methionine", "Met", 1, 1, kept=False),
            Region("Beta chain", "", 2, 147),
        ),
        structure=Structure(
            pdb="2DN1", sampling=4,
            # Only the beta chain: the gene makes that and nothing else. The
            # alpha chain beside it in the entry is HBA1/HBA2's, and drawn in a
            # second colour it read as a second product of this gene.
            chains=(Chain("chainA", "B"),),
            bonds=False,
        ),
    ),
    Target(
        slug="myoglobin", gene="MB", uniprot="P02144", aa=154,
        display="Myoglobin",
        source=Source("NG_007075"),
        cleaved=False, chain_label="Myoglobin",
        regions=(
            Region("Initiator methionine", "Met", 1, 1, kept=False),
            Region("Globin fold", "", 2, 148),
        ),
        structure=Structure(
            pdb="3RGK", sampling=6,
            chains=(Chain("chainA", "A"),),
            bonds=False,
        ),
    ),
    Target(
        slug="p53", gene="TP53", uniprot="P04637", aa=393,
        display="p53",
        source=Source("NG_017013", protein_id="NP_000537.3"),
        cleaved=False, chain_label="p53",
        # Only `DNA-binding domain` is a UniProt feature (DNA binding 102..292).
        # The other four are the standard literature boundaries; p53's ends are
        # largely disordered (the tetramerisation domain is the folded
        # exception) and UniProt annotates them as interaction regions rather
        # than domains.
        regions=(
            Region("Transactivation domain", "TAD", 1, 61),
            Region("Proline-rich region", "PRR", 62, 101),
            Region("DNA-binding domain", "DBD", 102, 292),
            Region("Tetramerisation domain", "TET", 323, 356),
            Region("Regulatory domain", "REG", 364, 393),
        ),
        structure=Structure(
            pdb="2OCJ", residues=(96, 289), sampling=5,
            chains=(Chain("chainA", "A"),),
            bonds=False,
        ),
    ),
    Target(
        slug="lysozyme", gene="LYZ", uniprot="P61626", aa=148,
        display="Lysozyme",
        source=Source("NG_008195"),
        cleaved=True, chain_label="Lysozyme C",
        regions=(
            Region("Signal peptide", "S", 1, 18, kept=False),
            Region("Lysozyme C", "", 19, 148),
        ),
        disulfides=((24, 146), (48, 134), (83, 99), (95, 113)),
        structure=Structure(
            pdb="1REX",
            chains=(Chain("chainA", "A"),),
            bonds=True,
        ),
    ),
    Target(
        slug="relaxin", gene="RLN2", uniprot="P04090", aa=185,
        display="Relaxin",
        source=Source("NC_000009.12", 5_299_363, 5_339_825, protein_id="NP_604390.1"),
        cleaved=True, chain_label="Prorelaxin H2", proprotein="prorelaxin",
        regions=(
            Region("Signal peptide", "S", 1, 24, kept=False),
            Region("B chain", "B", 25, 53),
            Region("C-peptide", "C", 56, 157, kept=False),
            Region("A chain", "A", 162, 185),
        ),
        disulfides=((35, 172), (47, 185), (171, 176)),
        structure=Structure(
            pdb="6RLX",
            chains=(Chain("chainA", "A"), Chain("chainB", "B")),
            bonds=True,
        ),
    ),
    Target(
        slug="oxytocin", gene="OXT", uniprot="P01178", aa=125,
        display="Oxytocin",
        source=Source("NC_000020.11", 3_071_119, 3_073_016),
        cleaved=True, chain_label="Oxytocin-neurophysin 1",
        proprotein="oxytocin-neurophysin 1",
        regions=(
            Region("Signal peptide", "S", 1, 19, kept=False),
            Region("Oxytocin", "OT", 20, 28),
            Region("Neurophysin 1", "NP", 32, 125),
        ),
        disulfides=(
            (20, 25), (41, 85), (44, 58), (52, 75),
            (59, 65), (92, 104), (98, 116), (105, 110),
        ),
        structure=Structure(
            pdb="7RYC", representation="tube",
            chains=(Chain("chainA", "L"),),
            bonds=True,
        ),
    ),
    Target(
        slug="somatotropin", gene="GH1", uniprot="P01241", aa=217,
        display="Growth hormone",
        source=Source("NG_011676"),
        cleaved=True, chain_label="Somatotropin",
        regions=(
            Region("Signal peptide", "S", 1, 26, kept=False),
            Region("Somatotropin", "", 27, 217),
        ),
        disulfides=((79, 191), (208, 215)),
        structure=Structure(
            pdb="1HGU", sampling=6,
            chains=(Chain("chainA", "A"),),
            bonds=True,
        ),
    ),
    Target(
        slug="ubiquitin", gene="UBB", uniprot="P0CG47", aa=229,
        display="Ubiquitin",
        source=Source("NG_023320"),
        cleaved=True, chain_label="Polyubiquitin-B",
        regions=(
            Region("Ubiquitin 1", "U1", 1, 76),
            Region("Ubiquitin 2", "U2", 77, 152),
            Region("Ubiquitin 3", "U3", 153, 228),
            Region("C-terminal residue", "", 229, 229, kept=False),
        ),
        structure=Structure(
            pdb="1UBQ",
            chains=(Chain("chainA", "A"),),
            bonds=False,
        ),
    ),
    Target(
        slug="dystrophin", gene="DMD", uniprot="P11532", aa=3685,
        display="Dystrophin",
        source=Source("NG_012232"),
        cleaved=False, chain_label="Dystrophin",
        # NG_012232.1 carries three `misc_difference` features against the
        # UniProt entry. They are common alleles, not errors in either record.
        uniprot_variants=((882, "G", "D"), (2366, "Q", "K"), (2937, "Q", "R")),
        regions=(
            Region("Calponin-homology 1", "CH1", 15, 119),
            Region("Calponin-homology 2", "CH2", 134, 240),
            *_spectrin(),
            Region("WW domain", "WW", 3055, 3088),
            Region("ZZ zinc finger", "ZZ", 3308, 3364),
        ),
        structure=Structure(
            pdb="1DXX", residues=(9, 246), sampling=4,
            chains=(Chain("chainA", "A"),),
            bonds=False,
        ),
    ),
    # The ten after the first ten. Each was checked against its record before it
    # was added — span, exons, neighbouring genes, processing — and what that
    # found is in docs/protein-verification.md. They were scored after that, and
    # what their tracks measured is in pipeline/constraint/verification.md.
    Target(
        slug="vasopressin", gene="AVP", uniprot="P01185", aa=164,
        display="Vasopressin",
        source=Source("NG_008663"),
        cleaved=True, chain_label="Vasopressin-neurophysin 2-copeptin",
        proprotein="vasopressin-neurophysin 2-copeptin",
        # Cut as oxytocin's precursor is, and once more: copeptin comes off at a
        # lone arginine, residue 125, rather than at a basic pair.
        regions=(
            Region("Signal peptide", "S", 1, 19, kept=False),
            Region("Vasopressin", "VP", 20, 28),
            Region("Neurophysin 2", "NP", 32, 124),
            Region("Copeptin", "CP", 126, 164),
        ),
        disulfides=(
            (20, 25), (41, 85), (44, 58), (52, 75),
            (59, 65), (92, 104), (98, 116), (105, 110),
        ),
        structure=Structure(
            pdb="7KH0", representation="tube",
            chains=(Chain("chainA", "L"),),
            bonds=True,
        ),
    ),
    Target(
        slug="glucagon", gene="GCG", uniprot="P01275", aa=180,
        display="Glucagon",
        # No RefSeqGene: NCBI Gene's coordinates and 500 bases either side, as
        # for oxytocin and relaxin. The gene is on the minus strand.
        source=Source("NC_000002.12", 162_142_381, 162_152_746),
        cleaved=True, chain_label="Proglucagon",
        # The record's eight peptides overlap: glicentin holds glucagon, and
        # GLP-1 holds its two shorter forms. They are what the pancreas and the
        # gut each cut the precursor into — alternatives, where the mature page
        # draws one series — so the walk stops at the precursor, and the regions
        # below are the precursor's parts rather than chains to fill.
        mature_peptides=False,
        regions=(
            Region("Signal peptide", "S", 1, 20, kept=False),
            Region("Glicentin-related polypeptide", "GRPP", 21, 50),
            Region("Glucagon", "GCG", 53, 81),
            Region("Intervening peptide 1", "IP1", 84, 89),
            Region("Glucagon-like peptide 1", "GLP1", 92, 128),
            Region("Intervening peptide 2", "IP2", 131, 145),
            Region("Glucagon-like peptide 2", "GLP2", 146, 178),
        ),
        structure=Structure(
            pdb="6LMK",
            chains=(Chain("chainA", "E"),),
            bonds=False,
        ),
    ),
    Target(
        slug="app", gene="APP", uniprot="P05067", aa=770,
        display="Amyloid precursor protein",
        source=Source("NG_007376"),
        cleaved=True, chain_label="Amyloid-beta precursor protein",
        # Fifteen peptides, every one inside the chain the signal peptide leaves
        # and most of them inside each other: what alpha- and beta-secretase cut
        # off, and what gamma-secretase makes of the rest. Alternatives, as
        # glucagon's are, so the walk stops at the precursor too.
        mature_peptides=False,
        proprotein="amyloid-beta precursor protein",
        regions=(
            Region("Signal peptide", "S", 1, 17, kept=False),
            Region("Amyloid-beta precursor protein", "", 18, 770),
        ),
        disulfides=(
            (38, 62), (73, 117), (98, 105), (133, 187), (144, 174),
            (158, 186), (291, 341), (300, 324), (316, 337),
        ),
        structure=Structure(
            # The E1 domain, 28-189: no structure covers the whole protein.
            pdb="4PWQ", sampling=6,
            chains=(Chain("chainA", "A"),),
            bonds=True,
        ),
    ),
    Target(
        slug="cftr", gene="CFTR", uniprot="P13569", aa=1480,
        display="CFTR",
        source=Source("NG_016465"),
        cleaved=False, chain_label="CFTR",
        structure=Structure(
            # Wild type. The regulatory region, 646-843, is not resolved; the
            # better-resolved entries all carry E1371Q. 1,139 residues at the
            # lowest sampling is what holds the compiled scene under budget.
            pdb="5UAK", sampling=1,
            chains=(Chain("chainA", "A"),),
            bonds=False,
        ),
    ),
    Target(
        slug="erythropoietin", gene="EPO", uniprot="P01588", aa=193,
        display="Erythropoietin",
        source=Source("NG_021471"),
        cleaved=True, chain_label="Erythropoietin",
        # NG_021471 annotates the signal peptide and no chain after it, so
        # `fill_chain` draws the chain from here.
        regions=(
            Region("Signal peptide", "S", 1, 27, kept=False),
            Region("Erythropoietin", "", 28, 193),
        ),
        disulfides=((34, 188), (56, 60)),
        structure=Structure(
            # Every EPO entry is the glycosylation-site mutant; this one also
            # carries P148N and P149S. Listed in docs/protein-verification.md.
            pdb="1EER", sampling=6,
            chains=(Chain("chainA", "A"),),
            bonds=True,
        ),
    ),
    Target(
        slug="leptin", gene="LEP", uniprot="P41159", aa=167,
        display="Leptin",
        source=Source("NG_007450"),
        cleaved=True, chain_label="Leptin",
        regions=(
            Region("Signal peptide", "S", 1, 21, kept=False),
            Region("Leptin", "", 22, 167),
        ),
        disulfides=((117, 167),),
        structure=Structure(
            pdb="1AX8", sampling=6,
            chains=(Chain("chainA", "A"),),
            bonds=True,
        ),
    ),
    Target(
        slug="tnf", gene="TNF", uniprot="P01375", aa=233,
        display="TNF-alpha",
        source=Source("NG_007462"),
        cleaved=True, chain_label="Tumor necrosis factor",
        # TNF is a type II membrane protein: no signal peptide to cleave, and
        # ADAM17 sheds the soluble form off the part that stays in the membrane
        # (P01375, `Site 76-77 Cleavage; by ADAM17`). NG_007462 annotates
        # neither piece — its neighbour in the record, LTA, has peptides of its
        # own, which the exact `/gene` match leaves out — so both come from the
        # table, and `fill_chain` draws them.
        #
        # Written as the cut divides the precursor rather than as UniProt writes
        # it. UniProt's chains overlap — the membrane form is 1-233, the soluble
        # form 77-233, and SPPL2A and SPPL2B cut the rest into three more — and
        # a page that cuts a precursor into pieces cannot draw one residue in
        # two of them. The two pieces ADAM17 actually leaves do divide it.
        #
        # 1-76 is kept, not thrown away: it is the cytoplasmic domain and the
        # signal-anchor helix, it stays in the membrane, and membrane TNF
        # signals in both directions. Named for what it is rather than for
        # UniProt's overlapping forms.
        regions=(
            Region("Membrane anchor", "", 1, 76),
            Region("Tumor necrosis factor", "", 77, 233),
        ),
        # UniProt's one bridge (P01375), inside the soluble part: the one the
        # fold page draws. Held in precursor numbering, as every bridge is, and
        # shown as the soluble chain's own Cys69-Cys101.
        disulfides=((145, 177),),
        structure=Structure(
            # One chain of the three: TNF works as a trimer of it, which the
            # fold page says rather than draws, as haemoglobin's does.
            pdb="7JRA", sampling=6,
            chains=(Chain("chainA", "A"),),
            bonds=True,
        ),
    ),
    Target(
        slug="sod1", gene="SOD1", uniprot="P00441", aa=154,
        display="SOD1",
        source=Source("NG_008689"),
        cleaved=False, chain_label="Superoxide dismutase",
        # The record's one peptide is residues 2-21, an antimicrobial fragment
        # of the finished enzyme rather than a piece it is cut into: the
        # hemorphins' case.
        mature_peptides=False,
        regions=(
            Region("Initiator methionine", "Met", 1, 1, kept=False),
            Region("Superoxide dismutase", "", 2, 154),
        ),
        disulfides=((58, 147),),
        structure=Structure(
            # A dimer of this chain; one is drawn. The entry numbers from Ala2,
            # so its C57-C146 bridge is 58-147 here.
            pdb="2C9V", sampling=6,
            chains=(Chain("chainA", "A"),),
            bonds=True,
        ),
    ),
    Target(
        slug="amylase", gene="AMY1A", uniprot="P0DUB6", aa=511,
        display="Amylase",
        # No RefSeqGene, and two transcripts around one CDS: NM_004038.4 is MANE
        # Select, NM_001008221.1 differs only in its untranslated first exon.
        source=Source(
            "NC_000001.11", 103_655_018, 103_665_053,
            protein_id="NP_004029.2", transcript_id="NM_004038.4",
        ),
        cleaved=True, chain_label="Alpha-amylase 1A",
        regions=(
            Region("Signal peptide", "S", 1, 15, kept=False),
            Region("Alpha-amylase 1A", "", 16, 511),
        ),
        disulfides=((43, 101), (85, 130), (156, 175), (393, 399), (465, 477)),
        structure=Structure(
            pdb="1SMD", sampling=3,
            chains=(Chain("chainA", "A"),),
            bonds=True,
        ),
    ),
    Target(
        slug="prion", gene="PRNP", uniprot="P04156", aa=253,
        display="Prion protein",
        source=Source("NG_009087"),
        cleaved=True, chain_label="Major prion protein",
        regions=(
            Region("Signal peptide", "S", 1, 22, kept=False),
            Region("Major prion protein", "", 23, 230),
            Region("GPI-anchor signal", "GPI", 231, 253, kept=False),
        ),
        disulfides=((179, 214),),
        structure=Structure(
            # The folded C-terminal domain, 117-225; the first half of the chain
            # has no fixed shape to solve.
            pdb="4KML",
            chains=(Chain("chainA", "A"),),
            bonds=True,
        ),
    ),
)

BY_SLUG = {t.slug: t for t in TARGETS}
BY_GENE = {t.gene: t for t in TARGETS}


def partition(target: Target) -> list[dict]:
    """The precursor as a gapless, sorted list of regions covering 1..aa.

    Every position belongs to exactly one region, so the Dart side is a lookup
    rather than a chain of index comparisons. `origin` is what a position's
    second number is counted from; the first is always the precursor's own,
    which is UniProt's and HGVS's. The rule follows the field's conventions
    rather than the table's layout:

    - A precursor cut into two or more kept pieces numbers each piece from its
      own start: insulin's Cys96 is `A7`, polyubiquitin's Lys124 is `U2 48`.
    - Anything else is one chain. It is numbered from its first residue after
      a removed N-terminal leader, where it has one, so the prion protein's
      Cys179 is mature 157 and SOD1's Ala5 is mature 4. A domain inside that
      chain is not a numbering of its own: p53's Arg175 is 175, not the 74th
      residue of its DNA-binding domain.
    - A piece that is removed (a leader, a cut site, a trailer) counts from 1,
      because nothing is numbered from it.
    """
    named = sorted(target.regions, key=lambda r: r.start)
    for before, after in zip(named, named[1:]):
        if before.end >= after.start:
            raise ValueError(f"{target.slug}: regions {before.label}/{after.label} overlap")
    if named and named[-1].end > target.aa:
        raise ValueError(f"{target.slug}: region {named[-1].label} runs past residue {target.aa}")

    out: list[dict] = []

    # Cut into pieces, each with its own numbering, or one chain.
    pieces = target.cleaved and sum(1 for region in named if region.kept) >= 2
    # One chain is counted from the first residue after the leading run of
    # removed regions that starts at residue 1: a signal peptide, an initiator
    # methionine. Nothing removed at the front leaves it counted from 1.
    mature = 1
    for region in named:
        if region.start != mature or region.kept:
            break
        mature = region.end + 1

    def origin_of(start: int, kept: bool, own: bool) -> int:
        if not kept:
            return 1
        if pieces:
            return start if own else 1
        return mature

    def filler(start: int, end: int, before: Region | None, after: Region | None) -> None:
        if start > end:
            return
        if target.cleaved and before is not None and after is not None:
            label = f"{before.short} / {after.short} cleavage site"
        elif target.cleaved and (before is None) != (after is None):
            # Past the last named piece of a cut precursor, or before the first,
            # is what gets trimmed off it, named as the gene page names it
            # (R3.6). Named for the chain instead, glucagon's trailing RK read
            # "removed with Proglucagon".
            label = "N-terminal extension" if before is None else "C-terminal extension"
        else:
            label = target.chain_label
        # A gap in a precursor is the basic residues a protease cuts at, and
        # they go with the cut. A gap in a folded chain is just more of that
        # chain, and stays.
        out.append({
            "label": label, "short": "", "start": start, "end": end,
            "origin": origin_of(start, not target.cleaved, own=False),
            "kept": not target.cleaved,
        })

    cursor = 1
    previous: Region | None = None
    for region in named:
        filler(cursor, region.start - 1, previous, region)
        out.append({
            "label": region.label, "short": region.short,
            "start": region.start, "end": region.end,
            "origin": origin_of(region.start, region.kept, own=True),
            "kept": region.kept,
        })
        cursor = region.end + 1
        previous = region
    filler(cursor, target.aa, previous, None)

    covered = sum(r["end"] - r["start"] + 1 for r in out)
    if covered != target.aa or out[0]["start"] != 1 or out[-1]["end"] != target.aa:
        raise ValueError(f"{target.slug}: partition covers {covered} of {target.aa} residues")
    return out
