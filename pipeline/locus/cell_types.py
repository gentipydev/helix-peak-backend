"""The Human Protein Atlas's cell types, by its own names, version 25.1.

The zoom's path goes from an organ into one kind of cell, and the cell must be
one that lives in that organ. The Atlas's single cell type reading names a
kind of cell, never where it lives, so this module says, once for every gene:

- `CLASSES`: the fifteen classes the Atlas sorts its cell types into, which
  the app draws a cell's shape from;
- `CELL_TYPES`: each of its 154 single cell types, with its class as the
  Atlas gives it and the consensus tissues (`tissues.TISSUES`) it lives in,
  the first its home. A cell type found along nerves or in connective tissue
  everywhere, or in a tissue the consensus reading does not sample (the
  lacrimal gland, the conjunctiva), has no home;
- `PAIR_CELL_CLASSES`: the class of a cell type the tissue cell type reading
  names that the single cell vocabulary does not (`Beta cells`,
  `Skeletal myocytes`);
- `ANUCLEATE`: the two cell types with no nucleus, and the cell that still
  has one, where it is found;
- `SUBCELLULAR` and `SECRETOME`: the Atlas's words for where in a cell a
  protein is and where it is secreted to.

The names and classes are the Atlas's search vocabularies and its single cell
type page as served for 25.1. The homes are general biology, not any gene's.
The Atlas writes some names in lower case in some answers (`monocytes`,
`cDC`), so every lookup ignores case.
"""

from __future__ import annotations

import re

CLASSES = (
    "Neuronal cells", "Glial cells", "Endocrine cells",
    "Squamous epithelial cells", "Pigment cells", "Ciliated cells",
    "Specialized epithelial cells", "Glandular epithelial cells", "Germ cells",
    "Trophoblast cells", "Muscle cells", "Endothelial and mural cells",
    "Mesenchymal cells", "Blood and immune cells",
    "Stem and proliferating cells",
)

_N, _G, _E = "Neuronal cells", "Glial cells", "Endocrine cells"
_SQ, _PI, _CI = "Squamous epithelial cells", "Pigment cells", "Ciliated cells"
_SP, _GL, _GE = ("Specialized epithelial cells", "Glandular epithelial cells",
                 "Germ cells")
_TR, _MU, _EN = "Trophoblast cells", "Muscle cells", "Endothelial and mural cells"
_ME, _BL, _ST = ("Mesenchymal cells", "Blood and immune cells",
                 "Stem and proliferating cells")

# Each single cell type: its class, and the consensus tissues it lives in,
# its home first.
CELL_TYPES: dict[str, tuple[str, tuple[str, ...]]] = {
    "Brain excitatory neurons": (_N, ("brain",)),
    "Brain inhibitory neurons": (_N, ("brain",)),
    "Retinal amacrine cells": (_N, ("retina",)),
    "Retinal horizontal cells": (_N, ("retina",)),
    "Retinal ganglion cells": (_N, ("retina",)),
    "Retinal bipolar cells": (_N, ("retina",)),
    "Rod photoreceptor cells": (_N, ("retina",)),
    "Cone photoreceptor cells": (_N, ("retina",)),
    "Other brain neurons": (_N, ("brain",)),
    "Astrocytes": (_G, ("brain",)),
    "Bergmann glia": (_G, ("brain",)),
    "Müller glia": (_G, ("retina",)),
    "Pituicytes/FSCs": (_G, ("pituitary gland",)),
    "Oligodendrocyte progenitor cells": (_G, ("brain",)),
    "Oligodendrocytes": (_G, ("brain",)),
    "Schwann cells": (_G, ()),
    "Microglia": (_G, ("brain",)),
    "Adrenal cortex cells": (_E, ("adrenal gland",)),
    "Adrenal medulla cells": (_E, ("adrenal gland",)),
    "Corticotrophs": (_E, ("pituitary gland",)),
    "Gonadotrophs": (_E, ("pituitary gland",)),
    "Lactotrophs": (_E, ("pituitary gland",)),
    "Somatotrophs": (_E, ("pituitary gland",)),
    "Thyrotrophs": (_E, ("pituitary gland",)),
    "Neuroendocrine cells": (_E, ("intestine", "stomach", "lung")),
    "Pancreatic islet cells": (_E, ("pancreas",)),
    "Leydig cells": (_E, ("testis",)),
    "Esophageal apical cells": (_SQ, ("esophagus",)),
    "Esophageal suprabasal cells": (_SQ, ("esophagus",)),
    "Suprabasal keratinocytes": (_SQ, ("skin",)),
    "Melanocytes": (_PI, ("skin",)),
    "Retinal pigment epithelial cells": (_PI, ("retina",)),
    "Ependymal cells": (_CI, ("brain",)),
    "Choroid plexus epithelial cells": (_CI, ("choroid plexus",)),
    "Respiratory ciliated cells": (_CI, ("lung",)),
    "Respiratory deuterosomal cells": (_CI, ("lung",)),
    "Epididymal efferent duct ciliated cells": (_CI, ("epididymis",)),
    "Fallopian tube ciliated cells": (_CI, ("fallopian tube",)),
    "Endometrial ciliated cells": (_CI, ("endometrium",)),
    "Alveolar cells type 1": (_SP, ("lung",)),
    "Alveolar cells type 2": (_SP, ("lung",)),
    "Transitional alveolar cells": (_SP, ("lung",)),
    "Ocular epithelial cells": (_SP, ("retina",)),
    "Mesothelial cells": (_SP, ("lung", "adipose tissue")),
    "Urothelial cells": (_SP, ("urinary bladder",)),
    "Epicardial cells": (_SP, ("heart muscle",)),
    "Prostatic hillock cells": (_SP, ("prostate",)),
    "Respiratory ionocytes": (_SP, ("lung",)),
    "Salivary ionocytes": (_SP, ("salivary gland",)),
    "Salivary myoepithelial cells": (_SP, ("salivary gland",)),
    "Paneth cells": (_SP, ("intestine",)),
    "Tuft cells": (_SP, ("intestine", "stomach", "lung")),
    "Breast myoepithelial cells": (_SP, ("breast",)),
    "Medullary thymic epithelial cells": (_SP, ("lymphoid tissue",)),
    "Salivary duct cells": (_SP, ("salivary gland",)),
    "Cholangiocytes": (_SP, ("liver", "gallbladder")),
    "Pancreatic duct cells": (_SP, ("pancreas",)),
    "Renal collecting duct intercalated cells": (_SP, ("kidney",)),
    "Renal collecting duct principal cells": (_SP, ("kidney",)),
    "Renal connecting tubule cells": (_SP, ("kidney",)),
    "Epididymal efferent duct absorptive cells": (_SP, ("epididymis",)),
    "Hepatocytes": (_SP, ("liver",)),
    "Podocytes": (_SP, ("kidney",)),
    "Proximal tubule cells": (_SP, ("kidney",)),
    "Loop of henle epithelial cells": (_SP, ("kidney",)),
    "Papillary tip epithelial cells": (_SP, ("kidney",)),
    "Distal convoluted tubule cells": (_SP, ("kidney",)),
    "Sertoli cells": (_SP, ("testis",)),
    "Epididymal clear cells": (_SP, ("epididymis",)),
    "Epididymal principal cells": (_SP, ("epididymis",)),
    "Granulosa cells": (_SP, ("ovary",)),
    "Respiratory basal cells": (_SP, ("lung",)),
    "Esophageal basal cells": (_SP, ("esophagus",)),
    "Salivary basal cells": (_SP, ("salivary gland",)),
    "Basal keratinocytes": (_SP, ("skin",)),
    "Epididymal basal cells": (_SP, ("epididymis",)),
    "Basal prostatic cells": (_SP, ("prostate",)),
    "Lacrimal acinar cells": (_GL, ()),
    "Salivary acinar cells": (_GL, ("salivary gland",)),
    "Pancreatic acinar cells": (_GL, ("pancreas",)),
    "Foveolar cells": (_GL, ("stomach",)),
    "Mucous neck cells": (_GL, ("stomach",)),
    "Parietal cells": (_GL, ("stomach",)),
    "Gastric chief cells": (_GL, ("stomach",)),
    "Enterocytes": (_GL, ("intestine",)),
    "Colonocytes": (_GL, ("intestine",)),
    "Submucosal glandular cells": (_GL, ("lung", "esophagus")),
    "Prostatic glandular cells": (_GL, ("prostate",)),
    "Prostatic club cells": (_GL, ("prostate",)),
    "Endometrial glandular cells": (_GL, ("endometrium",)),
    "Conjunctival goblet cells": (_GL, ()),
    "Goblet cells": (_GL, ("intestine",)),
    "Respiratory secretory cells": (_GL, ("lung",)),
    "Breast secretory cells": (_GL, ("breast",)),
    "Breast lactating cells": (_GL, ("breast",)),
    "Breast hormone-responsive cells": (_GL, ("breast",)),
    "Fallopian secretory cells": (_GL, ("fallopian tube",)),
    "Endometrial secretory cells": (_GL, ("endometrium",)),
    "Endometrial luminal cells": (_GL, ("endometrium",)),
    "Undifferentiated spermatogonia": (_GE, ("testis",)),
    "Differentiating spermatogonia": (_GE, ("testis",)),
    "Early primary spermatocytes": (_GE, ("testis",)),
    "Late primary spermatocytes": (_GE, ("testis",)),
    "Early spermatids": (_GE, ("testis",)),
    "Late spermatids": (_GE, ("testis",)),
    "Oocytes": (_GE, ("ovary",)),
    "Cytotrophoblasts": (_TR, ("placenta",)),
    "Migrating cytotrophoblasts": (_TR, ("placenta",)),
    "Syncytiotrophoblasts": (_TR, ("placenta",)),
    "Extravillous trophoblasts": (_TR, ("placenta",)),
    "Cardiomyocytes": (_MU, ("heart muscle",)),
    "Myonuclei": (_MU, ("skeletal muscle", "tongue")),
    "Smooth muscle cells": (_MU, ("smooth muscle",)),
    "Vascular endothelial cells": (_EN, ("blood vessel",)),
    "Lymphatic endothelial cells": (_EN, ("lymphoid tissue",)),
    "Pericytes": (_EN, ("blood vessel",)),
    "Vascular smooth muscle cells": (_EN, ("blood vessel",)),
    "Adipocytes": (_ME, ("adipose tissue", "breast")),
    "Fibro-adipogenic progenitors": (_ME, ("skeletal muscle",)),
    "Fibroblasts": (_ME, ()),
    "Hepatic stellate cells": (_ME, ("liver",)),
    "Peritubular myoid cells": (_ME, ("testis",)),
    "Ovarian stromal cells": (_ME, ("ovary",)),
    "Endometrial stromal cells": (_ME, ("endometrium",)),
    "Decidual stromal cells": (_ME, ("placenta",)),
    "Thymic myoid cells": (_ME, ("lymphoid tissue",)),
    "Mast cells": (_BL, ()),
    "Neutrophil progenitors": (_BL, ("bone marrow",)),
    "Neutrophils": (_BL, ("bone marrow",)),
    "Monocyte progenitors": (_BL, ("bone marrow",)),
    "Monocytes": (_BL, ("bone marrow",)),
    "Macrophages": (_BL, ()),
    "Kupffer cells": (_BL, ("liver",)),
    "Hofbauer cells": (_BL, ("placenta",)),
    "T-cells": (_BL, ("lymphoid tissue",)),
    "Thymocytes": (_BL, ("lymphoid tissue",)),
    "Innate lymphoid cells": (_BL, ("lymphoid tissue",)),
    "B-cells": (_BL, ("lymphoid tissue",)),
    "Plasma cells": (_BL, ("lymphoid tissue", "bone marrow")),
    "CDC": (_BL, ("lymphoid tissue",)),
    "PDCs": (_BL, ("lymphoid tissue",)),
    "NK-cells": (_BL, ("lymphoid tissue",)),
    "Megakaryocyte progenitors": (_BL, ("bone marrow",)),
    "Megakaryocyte-Erythroid progenitors": (_BL, ("bone marrow",)),
    "Megakaryocytes": (_BL, ("bone marrow",)),
    "Platelets": (_BL, ("bone marrow",)),
    "Erythrocyte progenitors": (_BL, ("bone marrow",)),
    "Erythrocytes": (_BL, ("bone marrow",)),
    "Gastric progenitor cells": (_ST, ("stomach",)),
    "Enteric transient amplifying cells": (_ST, ("intestine",)),
    "Pituitary stem cells": (_ST, ("pituitary gland",)),
    "Enteric stem cells": (_ST, ("intestine",)),
    "Myosatellite cells": (_ST, ("skeletal muscle",)),
    "Hematopoietic stem cells": (_ST, ("bone marrow",)),
}

# The tissue cell type reading's own names that the single cell vocabulary
# lacks, with their class. A name it shares with that vocabulary
# (`Hepatocytes`, `Adipocytes`) takes the class there.
PAIR_CELL_CLASSES = {
    "Alpha cells": _E,
    "Beta cells": _E,
    "Delta cells": _E,
    "Somatotropes": _E,
    "Corticotropes": _E,
    "Gonadotropes": _E,
    "Lactotropes": _E,
    "Thyrotropes": _E,
    "Thyroid glandular cells": _GL,
    "Skeletal myocytes": _MU,
    "Erythroid cells": _BL,
    "Exocrine glandular cells": _GL,
    "Eccrine sweat gland cells": _GL,
    "Minor salivary glandular cells": _GL,
    "Breast glandular cells": _GL,
    "Colon enterocytes": _GL,
    "Spermatogonia": _GE,
    "Adipose progenitor cells": _ME,
    "Endothelial cells": _EN,
    "Mitotic cells": _ST,
}

# The cell types with no nucleus, and the cell that still has one: red cells
# lose theirs as they mature in the marrow; platelets are fragments shed by
# megakaryocytes and never had one.
ANUCLEATE = {
    "Erythrocytes": {"cell": "erythroblasts", "place": "bone marrow"},
    "Platelets": {"cell": "megakaryocytes", "place": "bone marrow"},
}

# Where in a cell the Atlas finds a protein, by its own names.
SUBCELLULAR = (
    "Acrosome", "Actin filaments", "Aggresome", "Annulus", "Basal body",
    "Calyx", "Cell Junctions", "Centriolar satellite", "Centrosome",
    "Cleavage furrow", "Connecting piece", "Cytokinetic bridge",
    "Cytoplasmic bodies", "Cytosol", "End piece", "Endoplasmic reticulum",
    "Endosomes", "Equatorial segment", "Flagellar centriole",
    "Focal adhesion sites", "Golgi apparatus", "Intermediate filaments",
    "Kinetochore", "Lipid droplets", "Lysosomes", "Microtubule ends",
    "Microtubules", "Mid piece", "Midbody", "Midbody ring", "Mitochondria",
    "Mitotic chromosome", "Mitotic spindle", "Nuclear bodies",
    "Nuclear membrane", "Nuclear speckles", "Nucleoli",
    "Nucleoli fibrillar center", "Nucleoli rim", "Nucleoplasm",
    "Perinuclear theca", "Peroxisomes", "Plasma membrane", "Primary cilium",
    "Primary cilium tip", "Primary cilium transition zone", "Principal piece",
    "Rods & Rings", "Vesicles",
)

# Where the Atlas finds a protein secreted to, by its own names.
SECRETOME = (
    "Immunoglobulin genes", "Intracellular and membrane",
    "Secreted - unknown location", "Secreted in brain",
    "Secreted in female reproductive system",
    "Secreted in male reproductive system", "Secreted in other tissues",
    "Secreted to blood", "Secreted to digestive system",
    "Secreted to extracellular matrix",
)

_BY_FOLD = {name.casefold(): name for name in CELL_TYPES}
_PAIR_BY_FOLD = {name.casefold(): name for name in PAIR_CELL_CLASSES}
# A pair's cell type can name its tissue after it: `Mitotic cells (Stomach)`.
_QUALIFIER = re.compile(r"\s*\([^)]*\)$")


def single(name: str) -> str | None:
    """The single cell type [name] is, as the vocabulary writes it."""
    return _BY_FOLD.get(name.strip().casefold())


def homes(name: str) -> tuple[str, ...]:
    """The consensus tissues a single cell type lives in, its home first."""
    found = single(name)
    return CELL_TYPES[found][1] if found else ()


def cell_class(name: str) -> str | None:
    """The class of a cell type either reading names, or None for a name
    neither vocabulary has."""
    bare = _QUALIFIER.sub("", name.strip())
    found = single(name) or single(bare)
    if found:
        return CELL_TYPES[found][0]
    pair = _PAIR_BY_FOLD.get(bare.casefold())
    return PAIR_CELL_CLASSES[pair] if pair else None


def anucleate(name: str) -> dict | None:
    """The cell with a nucleus a zoom lands in instead, for a cell type with
    none."""
    found = single(name)
    return dict(ANUCLEATE[found]) if found in ANUCLEATE else None
