"""The locus bake's pinned inputs: the snapshot, made and read, offline.

A snapshot is built here from the locus fixtures, made small
(`testing.py`). Nothing fetches.
One test reads the real pinned snapshot, and skips where it has not been built
(`snapshot.py`), as the tests that need stored tracks skip until
`fetch_tracks.py` has filled `pipeline/data/`.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.locus import bake_locus, inputs, snapshot  # noqa: E402
from pipeline.locus.testing import RELEASES, atlas_genes, downloads, tables, text  # noqa: E402
from pipeline.paths import DATA  # noqa: E402


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """A snapshot written and read back: its path, its digest, and it read."""
    folder = downloads(tmp_path_factory.mktemp("downloads"))
    path, digest = snapshot.write(text(folder), tmp_path_factory.mktemp("snapshot"))
    return path, digest, inputs.read(path, digest)


def test_the_same_sources_make_the_same_snapshot(tmp_path):
    folder = downloads(tmp_path / "downloads")
    first, digest = snapshot.write(text(folder), tmp_path / "a")
    second, again = snapshot.write(text(folder), tmp_path / "b")
    assert digest == again
    assert first.read_bytes() == second.read_bytes()
    assert first.name == f"locus-inputs.{digest[:12]}.jsonl.gz"
    # The digest is the text's, so another compressor still matches it.
    assert inputs.read(first, digest).sha256 == digest


def test_another_date_is_another_snapshot(tmp_path):
    folder = downloads(tmp_path / "downloads")
    assert text(folder, "2026-11-13") != text(folder, "2026-11-14")


def test_any_other_snapshot_is_refused(built, tmp_path):
    path, digest, _ = built
    with pytest.raises(inputs.InputsUnavailable, match="not the pinned"):
        inputs.read(path, "f" * 64)
    with pytest.raises(inputs.InputsUnavailable, match="snapshot.py"):
        inputs.read(tmp_path / "nowhere.jsonl.gz", digest)
    # A file that is not one at all.
    junk = tmp_path / "junk.jsonl.gz"
    junk.write_bytes(gzip.compress(b'{"kind": "something else", "format": 1}\n'))
    with pytest.raises(inputs.InputsUnavailable):
        inputs.read(junk, None)


def test_a_gene_reads_as_the_atlas_answered_it(built):
    _, _, read = built
    for gene in atlas_genes():
        found = read.atlas_gene(gene["Ensembl"])
        assert found == {name: gene.get(name) for name in snapshot.ATLAS_FIELDS}
        assert bake_locus.expression_of(found) == bake_locus.expression_of(gene)
    assert read.atlas_gene("ENSG00000000000") is None


def test_its_levels_are_the_tables_own(built):
    _, _, read = built
    tissues, cells = read.levels("ENSG00000254647")
    assert set(tissues) == {"colon", "pancreas", "small intestine"}
    assert tissues["small intestine"] == 0
    assert cells == {"enterocytes": 2.5, "pancreatic islet cells": 105}
    assert read.levels("ENSG00000000000") is None


def test_it_finds_a_gene_by_uniprot_and_by_hgnc(built):
    _, _, read = built
    assert read.by_uniprot("P01308") == ("ENSG00000254647",)
    assert read.by_uniprot("Q00000") == ()
    assert read.hgnc_ensembl("HGNC:9187") == "ENSG00000181222"
    assert read.hgnc_ensembl("HGNC:999999") is None
    assert read.hgnc_ensembl(None) is None


def test_it_places_a_scaffold_by_its_refseq_accession(built):
    _, _, read = built
    (kir,) = read.placements("NT_187669.1")
    assert (kir["parent_acc"], kir["ori"], kir["parent_start"]) == ("NC_000019.10", "-", "54698377")
    assert read.placements("NW_000000000.1") == []


def test_a_table_from_it_is_the_table_as_read(built):
    _, _, read = built
    for track, table in tables().items():
        kept = bake_locus.table_from(read, track)
        assert kept.sha256 == table.sha256
        assert kept.provenance() == table.provenance()


def test_it_names_its_sources(built):
    _, _, read = built
    assert (read.atlas["version"], read.atlas["site"]) == ("25.1", snapshot.HPA_SITE)
    assert read.mane_release == "v1.5"
    assert read.sources["hgnc"]["date"] == snapshot.HGNC_DATE
    assert read.sources["placement"]["assembly"] == "GRCh38.p14"
    atlas = bake_locus.atlas_from(read)
    assert atlas.provenance("ENSG00000254647")["url"] == \
        "https://v25.proteinatlas.org/ENSG00000254647.json"


def test_it_refuses_another_atlas_release(tmp_path):
    folder = downloads(tmp_path / "downloads")
    page = RELEASES.replace("version 25.1", "version 26.0", 1)
    with pytest.raises(ValueError, match="26.0"):
        snapshot.snapshot_text(folder=folder, release_page=page, tables=tables(),
                               created="2026-11-13")


def test_it_is_looked_for_where_it_is_named(tmp_path, monkeypatch):
    monkeypatch.setenv(inputs.ENV, str(tmp_path))
    assert inputs.default_path("ab" * 32) == tmp_path / "locus-inputs.abababababab.jsonl.gz"
    named = tmp_path / "elsewhere.jsonl.gz"
    monkeypatch.setenv(inputs.ENV, str(named))
    assert inputs.default_path() == named
    monkeypatch.delenv(inputs.ENV)
    assert inputs.default_path().parent == inputs.FOLDER


_PINNED = inputs.default_path()


@pytest.mark.skipif(not _PINNED.is_file(), reason="the pinned snapshot is not built here")
def test_the_pinned_snapshot_gives_every_stored_curated_reading():
    """Every curated protein's stored `expression`, as the snapshot gives it:
    the Atlas 25.1 the twenty and those after them were baked from."""
    pinned = inputs.load()
    assert pinned.sha256 == inputs.INPUTS_SHA256
    assert pinned.created == inputs.INPUTS_CREATED
    assert pinned.atlas["version"] == "25.1"
    stored = sorted((DATA / "assets" / "locus").glob("*_locus.json"))
    if not stored:
        pytest.skip("no stored locus tracks: run pipeline/fetch_tracks.py")
    for path in stored:
        track = json.loads(path.read_bytes())
        expression = dict(track["expression"])
        ensembl = expression.pop("ensembl_gene")
        assert bake_locus.expression_of(pinned.atlas_gene(ensembl)) == expression, path.name


def test_the_snapshot_folder_is_ignored_by_git():
    ignore = (BACKEND / ".gitignore").read_text(encoding="utf-8").splitlines()
    assert "pipeline/data/" in ignore
    assert inputs.FOLDER.is_relative_to(BACKEND / "pipeline" / "data")


def test_the_bake_fetches_nothing_but_uniprot():
    """The Atlas, UCSC and MANE come from the snapshot alone: the bake names
    no fetch of them."""
    source = Path(bake_locus.__file__).read_text(encoding="utf-8")
    main = source[source.index("def main("):]
    assert "fetch_table(" not in main and "current_summary" not in source
    assert "proteinatlas.org/about" not in source
