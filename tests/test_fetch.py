"""Tests for GET /fetch/{id}.

Every test patches Entrez.efetch, so the default run makes no network call.
"""

import email.message
import io
from urllib.error import HTTPError, URLError

import pytest

from app import entrez_client

# NCBI's two real failure shapes, captured from the live API. Neither is an
# empty body, so neither is caught by a naive "if not content" check.
UNPARSEABLE_ID_BODY = (
    "Error: F a i l e d  t o  u n d e r s t a n d  i d :  N O T _ A _ R E A L _ I D \n\n\n"
)
MISSING_ACCESSION_BODY = (
    b"Error: CEFetchPApplication::proxy_stream(): Error: "
    b"F a i l e d  t o  r e t r i e v e  s e q u e n c e :  N G _ 9 9 9 9 9 9 \n\n"
)


def _http_error(code, body):
    return HTTPError(
        url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        code=code,
        msg="Bad Request",
        hdrs=email.message.Message(),
        fp=io.BytesIO(body),
    )


@pytest.fixture
def efetch_returning(monkeypatch):
    """Patch Entrez.efetch to return the given text, recording its arguments."""

    def _install(text):
        calls = {}

        def fake_efetch(**kwargs):
            calls.update(kwargs)
            return io.StringIO(text)

        monkeypatch.setattr(entrez_client.Entrez, "efetch", fake_efetch)
        return calls

    return _install


@pytest.fixture
def efetch_raising(monkeypatch):
    """Patch Entrez.efetch to raise the given exception."""

    def _install(exception):
        def fake_efetch(**kwargs):
            raise exception

        monkeypatch.setattr(entrez_client.Entrez, "efetch", fake_efetch)

    return _install


def test_returns_raw_genbank_text(client, genbank_text, efetch_returning):
    efetch_returning(genbank_text)

    response = client.get("/fetch/NG_007114")

    assert response.status_code == 200
    body = response.json()
    assert body["id"] == "NG_007114"
    content = body["content"]

    assert content.startswith(
        "LOCUS       NG_007114               8416 bp    DNA     linear   PRI"
    )
    # The INS CDS block. NG_007114 also carries the adjacent TH gene, so assert
    # on INS-specific strings rather than on any CDS.
    assert '     CDS             join(5224..5410,6198..6343)' in content
    assert '/gene="INS"' in content
    assert '/product="insulin preproprotein"' in content
    assert '/protein_id="NP_000198.1"' in content
    # Returned as-is, terminator and all -- no parsing, no reformatting.
    assert content == genbank_text


def test_query_is_hardcoded_except_for_id(client, genbank_text, efetch_returning):
    calls = efetch_returning(genbank_text)

    client.get("/fetch/NG_007114")

    assert calls == {
        "db": "nucleotide",
        "id": "NG_007114",
        "rettype": "gb",
        "retmode": "text",
    }


def test_unparseable_id_is_404_not_an_error_body(client, efetch_returning):
    """NCBI answers HTTP 200 with an error body here; we must not pass it on."""
    efetch_returning(UNPARSEABLE_ID_BODY)

    response = client.get("/fetch/NOT_A_REAL_ID")

    assert response.status_code == 404
    assert "NOT_A_REAL_ID" in response.json()["detail"]


def test_missing_accession_is_404(client, efetch_raising):
    """A well-formed but unknown accession comes back as an upstream 400."""
    efetch_raising(_http_error(400, MISSING_ACCESSION_BODY))

    response = client.get("/fetch/NG_999999")

    assert response.status_code == 404
    assert "NG_999999" in response.json()["detail"]


def test_upstream_server_error_is_502_with_ncbi_message(client, efetch_raising):
    efetch_raising(_http_error(500, b"Internal Server Error"))

    response = client.get("/fetch/NG_007114")

    assert response.status_code == 502
    detail = response.json()["detail"]
    assert "500" in detail
    assert "Internal Server Error" in detail


def test_unreachable_ncbi_is_502(client, efetch_raising):
    efetch_raising(URLError("Name or service not known"))

    response = client.get("/fetch/NG_007114")

    assert response.status_code == 502
    assert "Could not reach NCBI" in response.json()["detail"]


def test_timeout_is_502(client, efetch_raising):
    efetch_raising(TimeoutError("timed out"))

    response = client.get("/fetch/NG_007114")

    assert response.status_code == 502
    assert "Could not reach NCBI" in response.json()["detail"]


def test_health(client):
    assert client.get("/health").json() == {"status": "ok"}
