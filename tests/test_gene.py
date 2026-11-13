"""Tests for GET /gene/{id}/{gene}.

Every test patches Entrez.efetch against the saved NG_007114 record, so the
default run makes no network call.
"""

from urllib.error import URLError

from tests.ncbi_errors import (
    MISSING_ACCESSION_BODY,
    UNPARSEABLE_ID_BODY,
    http_error,
)


def test_returns_the_insulin_gene(client, genbank_text, efetch_returning):
    efetch_returning(genbank_text)

    response = client.get("/gene/NG_007114/INS")

    assert response.status_code == 200
    body = response.json()
    assert body["gene"] == "INS"
    assert body["location"] == {"start": 4986, "end": 6416, "strand": 1}
    assert len(body["sequence"]) == 1431
    assert [exon["number"] for exon in body["exons"]] == [1, 2, 3]
    assert body["protein"]["product"] == "insulin preproprotein"
    assert [peptide["product"] for peptide in body["peptides"]] == [
        "insulin B chain",
        "C-peptide",
        "insulin A chain",
    ]


def test_response_is_a_json_object(client, genbank_text, efetch_returning):
    """The Flutter client rejects a bare top-level array."""
    efetch_returning(genbank_text)

    assert isinstance(client.get("/gene/NG_007114/INS").json(), dict)


def test_does_not_leak_the_overlapping_readthrough(
    client, genbank_text, efetch_returning
):
    """INS-IGF2 shares INS's start coordinate, so a sloppy filter would show it."""
    efetch_returning(genbank_text)

    text = str(client.get("/gene/NG_007114/INS").json())
    assert "INS-IGF2" not in text
    assert "insulin, isoform 2 precursor" not in text
    assert "tyrosine" not in text


def test_id_and_gene_are_both_free_parameters(client, genbank_text, efetch_returning):
    calls = efetch_returning(genbank_text)

    body = client.get("/gene/NG_007114/TH").json()

    assert calls["id"] == "NG_007114"
    assert body["gene"] == "TH"
    assert body["location"] == {"start": 1, "end": 2266, "strand": 1}


def test_absent_gene_is_404(client, genbank_text, efetch_returning):
    efetch_returning(genbank_text)

    response = client.get("/gene/NG_007114/BRCA1")

    assert response.status_code == 404
    detail = response.json()["detail"]
    assert "BRCA1" in detail
    assert "NG_007114" in detail


def test_unparseable_id_is_404_not_an_error_body(client, efetch_returning):
    """NCBI answers HTTP 200 with plain text here; SeqIO rejects it."""
    efetch_returning(UNPARSEABLE_ID_BODY)

    response = client.get("/gene/NOT_A_REAL_ID/INS")

    assert response.status_code == 404
    assert "NOT_A_REAL_ID" in response.json()["detail"]


def test_missing_accession_is_404(client, efetch_raising):
    efetch_raising(http_error(400, MISSING_ACCESSION_BODY))

    response = client.get("/gene/NG_999999/INS")

    assert response.status_code == 404
    assert "NG_999999" in response.json()["detail"]


def test_upstream_server_error_is_502(client, efetch_raising):
    efetch_raising(http_error(500, b"Internal Server Error"))

    response = client.get("/gene/NG_007114/INS")

    assert response.status_code == 502
    detail = response.json()["detail"]
    assert "500" in detail
    assert "Internal Server Error" in detail


def test_unreachable_ncbi_is_502(client, efetch_raising):
    efetch_raising(URLError("Name or service not known"))

    response = client.get("/gene/NG_007114/INS")

    assert response.status_code == 502
    assert "Could not reach NCBI" in response.json()["detail"]


def test_timeout_is_502(client, efetch_raising):
    """TimeoutError is not a URLError -- it lands on the bare OSError clause."""
    efetch_raising(TimeoutError("timed out"))

    response = client.get("/gene/NG_007114/INS")

    assert response.status_code == 502
    assert "Could not reach NCBI" in response.json()["detail"]


# -------------------------------------------- what the route may cost (2026-11-13)
#
# The route reads NCBI for whoever asks, and the app no longer reads it, so it
# is bounded: one record of a gene's size, a size limit on what is read, and
# so many NCBI reads an hour.

import pytest  # noqa: E402

from app import record_cache  # noqa: E402
from app.config import settings  # noqa: E402


@pytest.mark.parametrize("id", [
    "NG_007114,NG_012232",       # efetch would read several records
    "NC_000001.11",              # a chromosome whole
    "NT_187361.1",               # a contig
    "NW_025791820.1",            # a scaffold
    "AC_000133.1",               # an alternate assembly's chromosome
    "NG_007114.1:1-500",         # a range
    "ng_007114",                 # not as NCBI spells one
    "NG_007114.1234",
    "%2e%2e",
])
def test_an_id_that_is_not_one_genes_record_never_reaches_ncbi(client, efetch_raising, id):
    efetch_raising(AssertionError("NCBI was asked for " + id))
    response = client.get("/gene/{}/INS".format(id))
    assert response.status_code == 404
    assert response.json()["detail"].startswith("No GenBank record found for id")


@pytest.mark.parametrize("id", ["NG_007114", "NG_007114.1", "NM_000207.3", "M10051", "AF123456.2"])
def test_one_record_of_a_genes_size_is_read(client, genbank_text, efetch_returning, id):
    calls = efetch_returning(genbank_text)
    client.get("/gene/{}/INS".format(id))
    assert calls["id"] == id


def test_a_record_longer_than_the_service_reads_is_refused_unread(
        client, genbank_text, efetch_returning, monkeypatch):
    monkeypatch.setattr(settings, "max_record_bytes", 1000)
    efetch_returning(genbank_text)

    response = client.get("/gene/NG_007114/INS")

    assert response.status_code == 422
    assert response.json()["detail"] == \
        "Record 'NG_007114' is longer than this service reads (1,000 characters)."


def test_ncbi_is_read_so_many_times_an_hour(client, genbank_text, efetch_returning, monkeypatch):
    monkeypatch.setattr(settings, "gene_fetches_per_hour", 2)
    record_cache.reads.clear()
    efetch_returning(genbank_text)
    try:
        said = [client.get("/gene/NG_007114/INS").status_code for _ in range(3)]
        refused = client.get("/gene/NG_007114/INS")
    finally:
        record_cache.reads.clear()

    assert said == [200, 200, 429]
    assert refused.json()["detail"] == (
        "This service reads NCBI 2 times an hour at most, and this hour's are taken. "
        "Ask again later.")


def test_a_genes_name_is_bounded(client, efetch_raising):
    efetch_raising(AssertionError("NCBI was asked"))
    assert client.get("/gene/NG_007114/" + "A" * 65).status_code == 422
