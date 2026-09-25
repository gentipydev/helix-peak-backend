"""The pipeline package, as the web service's test suite sees it.

`pipeline/` is the bake tools, moved here from the client in Phase 3 of
HANDOFF-ONDEMAND.md. Nothing under `app/` imports it yet, so these hold the
package to the two things the move promised: it imports from the repository
root on the service's own Python, and the curated rows it seeds from are the
twenty `targets.py` bakes.
"""

import copy
import importlib

import pytest

from pipeline import check_assets, seed_catalog
from pipeline.targets import TARGETS


@pytest.mark.parametrize("module", [
    "pipeline.paths",
    "pipeline.targets",
    "pipeline.check_assets",
    "pipeline.seed_catalog",
    "pipeline.upload_tracks",
    "pipeline.fetch_tracks",
    "pipeline.impact.check_explanations",
])
def test_the_tools_import_from_the_repository_root(module):
    importlib.import_module(module)


def test_the_curated_rows_are_the_twenty_targets():
    curated = seed_catalog.curated_rows()
    assert sorted(curated) == sorted(t.slug for t in TARGETS)
    assert sorted(row["order"] for row in curated.values()) == list(range(len(TARGETS)))
    for target in TARGETS:
        row = check_assets.curated_catalog()[target.slug]
        assert (row["gene"], row["uniprot"], row["accession"]) == (
            target.gene, target.uniprot, target.source.accession)


def _served(target):
    """The row `/protein/{slug}` serves for a target seeded from these tables."""
    curated = seed_catalog.curated_rows()[target.slug]
    row = seed_catalog.protein_row(target, curated)
    return {
        "gene": row["gene"],
        "uniprot": row["uniprot"],
        "accession": row["accession"],
        "transcript_id": row["transcript_id"],
        "protein_id": row["protein_id"],
        "mature_peptides": row["mature_peptides"],
        "display": row["display"],
        "summary": row["summary"],
        "chain": row["chain_name"],
        "facts": {key: row[key] for key in ("residues", "exons", "chains", "bridges")},
        "regions": row["regions"],
        "disulfides": row["disulfides"],
        "chains": (row["structure"] or {}).get("chains", []),
        "structure": (row["structure"] or {}).get("chrome"),
        "catalog_order": row["catalog_order"],
        "tracks": {kind: "ready" for kind in
                   ("constraint", "impact", "clinvar", "impact_explanations")},
    }


@pytest.fixture
def against(monkeypatch):
    """`check_against` over served rows this test controls, and its findings."""
    def run(served):
        monkeypatch.setattr(check_assets, "service_rows", lambda base_url: served)
        monkeypatch.setattr(check_assets, "problems", [])
        check_assets.check_against("http://test", check_assets.curated_catalog())
        return check_assets.problems
    return run


def test_a_faithful_service_passes(against):
    assert against({t.slug: _served(t) for t in TARGETS}) == []


def test_the_fold_prose_is_compared_whole(against):
    served = {t.slug: _served(t) for t in TARGETS}
    edited = copy.deepcopy(served["insulin"])
    edited["structure"]["sentence"] += " Edited."
    served["insulin"] = edited
    assert against(served) == ["insulin: served structure.sentence is not the curated row's"]


def test_the_reading_order_is_compared(against):
    served = {t.slug: _served(t) for t in TARGETS}
    served["prion"]["catalog_order"] = 0
    found = against(served)
    assert len(found) == 1 and found[0].startswith("prion: served catalog_order 0")
