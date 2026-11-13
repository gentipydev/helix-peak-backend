"""The Human Protein Atlas's tissues, by its own names, version 25.1.

The Atlas names a tissue in two vocabularies, and the zoom's path needs both
to agree on one place:

- `TISSUES`: the tissues its consensus RNA reading names, the keys of a gene's
  `RNA tissue specific nTPM` (`bone marrow`, `choroid plexus`). A tissue it
  sampled more than once can carry a sample number (`stomach 1`), which names
  the sample, not another tissue.
- `PAIR_TISSUES`: the tissues its tissue cell type reading names, the part
  before ` - ` in a gene's `RNA tissue cell type enrichment`
  (`Adipose subcutaneous - Adipocytes`), each with the consensus tissue it is.

Both lists are the Atlas's search vocabularies as its site served them for
25.1 (`tissue_category_rnatissue` on the single cell type page and
`ce_tissue` on the tissue cell type page). General anatomy, not any gene's.
"""

from __future__ import annotations

import re

# The consensus tissues, lower case, as a gene's reading names them.
TISSUES = (
    "adipose tissue", "adrenal gland", "blood vessel", "bone marrow", "brain",
    "breast", "cervix", "choroid plexus", "endometrium", "epididymis",
    "esophagus", "fallopian tube", "gallbladder", "heart muscle", "intestine",
    "kidney", "liver", "lung", "lymphoid tissue", "ovary", "pancreas",
    "parathyroid gland", "pituitary gland", "placenta", "prostate", "retina",
    "salivary gland", "seminal vesicle", "skeletal muscle", "skin",
    "smooth muscle", "stomach", "testis", "thyroid gland", "tongue",
    "urinary bladder", "vagina",
)

# The tissue cell type reading's tissues, and the consensus tissue each is.
# Its two adipose depots are one consensus tissue; its colon is the
# consensus intestine; its spleen is lymphoid tissue, which the consensus
# reading groups the spleen into.
PAIR_TISSUES = {
    "Adipose subcutaneous": "adipose tissue",
    "Adipose visceral": "adipose tissue",
    "Adrenal gland": "adrenal gland",
    "Breast": "breast",
    "Colon": "intestine",
    "Heart muscle": "heart muscle",
    "Kidney": "kidney",
    "Liver": "liver",
    "Lung": "lung",
    "Minor Salivary Gland": "salivary gland",
    "Pancreas": "pancreas",
    "Pituitary gland": "pituitary gland",
    "Prostate": "prostate",
    "Skeletal muscle": "skeletal muscle",
    "Skin": "skin",
    "Spleen": "lymphoid tissue",
    "Stomach": "stomach",
    "Testis": "testis",
    "Thyroid gland": "thyroid gland",
}

# The tissues the consensus level table names (`rna_tissue_consensus.tsv`) that
# its gene readings group under one of `TISSUES`: the regions of the brain
# and of the spinal cord as the brain, the parts of the gut as the intestine,
# and the lymphoid organs as lymphoid tissue, as the Atlas groups them. The
# table's other 34 names are `TISSUES` as they are. A group's level is the
# highest of its members'.
CONSENSUS_GROUPS = {
    "amygdala": "brain", "basal ganglia": "brain", "cerebellum": "brain",
    "cerebral cortex": "brain", "hippocampal formation": "brain",
    "hypothalamus": "brain", "midbrain": "brain", "spinal cord": "brain",
    "colon": "intestine", "duodenum": "intestine", "rectum": "intestine",
    "small intestine": "intestine",
    "appendix": "lymphoid tissue", "lymph node": "lymphoid tissue",
    "spleen": "lymphoid tissue", "thymus": "lymphoid tissue", "tonsil": "lymphoid tissue",
}

_SAMPLE = re.compile(r"\s+\d+$")


def tissue_name(raw: str) -> str | None:
    """The consensus tissue a reading's name is, without its sample number,
    or None for a name the Atlas's list does not have."""
    name = _SAMPLE.sub("", raw.strip()).lower()
    return name if name in TISSUES else None


def consensus_group(raw: str) -> str | None:
    """The tissue of `TISSUES` a consensus level's tissue counts toward: itself,
    or the group the Atlas puts it in; None for a name neither list has."""
    return tissue_name(raw) or CONSENSUS_GROUPS.get(raw.strip().lower())


def pair_tissue(raw: str) -> str | None:
    """The consensus tissue a tissue cell type pair's tissue is, or None."""
    return PAIR_TISSUES.get(raw.strip())
