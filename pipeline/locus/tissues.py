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

_SAMPLE = re.compile(r"\s+\d+$")


def tissue_name(raw: str) -> str | None:
    """The consensus tissue a reading's name is, without its sample number,
    or None for a name the Atlas's list does not have."""
    name = _SAMPLE.sub("", raw.strip()).lower()
    return name if name in TISSUES else None


def pair_tissue(raw: str) -> str | None:
    """The consensus tissue a tissue cell type pair's tissue is, or None."""
    return PAIR_TISSUES.get(raw.strip())
