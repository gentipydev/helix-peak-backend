"""The assay track behind each AVI feature: one Atlas query per substitution.

Ported from the AVI skill's `_find_adata_max` and `_resolve_track_attributions`
(`.claude/skills/alphagenome-variant-impact-score/scripts/alphagenome_atlas_avi.py`),
which this venv cannot import: it needs polars at import time. A winner is the
skill's reconstruction, the track with the largest absolute effect for that
scorer; the AVI model reports only the feature's attribution, not its input.

A reply (genes x tracks for twelve scorers, ~130 KB) is reduced on arrival and
never stored. What is kept per feature is the scorer, the track's assay,
biosample, ontology term, gene and target (factor, mark or splice site), and
the assay's own signed value.
"""
from __future__ import annotations

import time

# Feature -> its one scorer. MERGED_SPLICING is three, compared at these scales.
SCORERS = {
    "MAX_ABS_RNA_SEQ": "RNA_SEQ",
    "MAX_ABS_ATAC": "ATAC",
    "MAX_ABS_DNASE": "DNASE",
    "MAX_ABS_CHIP_TF": "CHIP_TF",
    "MAX_ABS_CHIP_HISTONE": "CHIP_HISTONE",
    "MAX_ABS_CAGE": "CAGE",
    "MAX_ABS_PROCAP": "PROCAP",
    "MAX_ABS_POLYADENYLATION": "POLYADENYLATION",
    "MAX_ABS_CONTACT_MAPS": "CONTACT_MAPS",
}
SPLICING = (("SPLICE_SITES", 1.0), ("SPLICE_SITE_USAGE", 1.0), ("SPLICE_JUNCTIONS", 0.2))
REQUESTED = (*SCORERS.values(), *(name for name, _ in SPLICING))
RETRIED = ("RESOURCE_EXHAUSTED", "UNAVAILABLE", "DEADLINE_EXCEEDED",
           "failed to connect", "Broken pipe")


def _text(row, column: str) -> str:
    value = row.get(column, "") if row is not None else ""
    return "" if value is None or value != value else str(value)


def strongest(adata, *, absolute: bool) -> tuple[list[str], float] | None:
    """The (track, value) with the largest |value| (or value), as the skill finds it."""
    import numpy as np

    if adata is None or adata.X is None or adata.X.size == 0:
        return None
    values = np.nan_to_num(np.asarray(adata.X, dtype=np.float64))
    if values.ndim == 1:
        values = values.reshape(1, -1)
    flat = int(np.argmax(np.abs(values) if absolute else values))
    g, t = (int(i) for i in np.unravel_index(flat, values.shape))
    track = adata.var.iloc[t] if len(adata.var) > t else None
    gene = adata.obs.iloc[g] if len(adata.obs) > g else None
    target = (_text(track, "transcription_factor") or _text(track, "histone_mark")
              or ("" if _text(track, "Assay title") else _text(track, "name")))
    return [
        _text(track, "Assay title"),
        _text(track, "biosample_name"),
        _text(track, "ontology_curie"),
        _text(gene, "gene_name"),
        target,
    ], float(values[g, t])


def reduce(replies) -> dict[str, list]:
    """Feature -> [scorer, assay, biosample, ontology, gene, target, value]."""
    out = {}
    for feature, scorer in SCORERS.items():
        hit = strongest(replies.get(scorer), absolute=True)
        if hit is not None:
            out[feature] = [scorer, *hit[0], hit[1]]
    splice = []
    for scorer, scale in SPLICING:
        hit = strongest(replies.get(scorer), absolute=False)
        if hit is not None:
            splice.append((hit[1] * scale, [scorer, *hit[0], hit[1]]))
    if splice:
        out["MERGED_SPLICING"] = max(splice, key=lambda s: s[0])[1]
    return out


def winners(client, variant) -> dict[str, list]:
    """One substitution's winners, retried while the Atlas is busy."""
    for attempt in range(8):
        try:
            return reduce(client.query_variant(variant, requested_scorers=list(REQUESTED)))
        except Exception as error:
            if attempt == 7 or not any(word in str(error) for word in RETRIED):
                raise
            time.sleep(min(5 * 2**attempt, 60))
    raise AssertionError("unreachable")
