"""Tests for the catalog endpoints.

No database is involved. A fake pool stands in for psycopg so the policy --
what a page returns, how a search ranks, what an unreadable catalog reports,
and how a track becomes a URL -- is testable without a socket, the same trade
`test_record_cache.py` makes. The SQL itself is exercised against the real
database by running the service.
"""

import pytest

from app import catalog, db, tracks
from app.config import settings


def _protein_row(
    slug="insulin", display="Insulin", gene="INS", uniprot="P01308",
    accession="NG_007114", summary="The hormone that clears glucose.",
    residues=110, exons=3, chains=3, bridges=3, catalog_order=0,
):
    return (slug, display, gene, uniprot, accession, summary,
            residues, exons, chains, bridges, catalog_order)


_DETAIL_TAIL = (
    None,                       # chain_name
    True,                       # mature_peptides
    "NM_000207.3",              # transcript_id
    "NP_000198.1",              # protein_id
    [{"label": "B chain", "short": "B", "start": 25, "end": 54, "origin": 1,
      "kept": True}],
    [[31, 96]],
    {"chrome": {"pdb": "3I40", "modelled": None, "label": "the hormone",
                "count": 51, "unit": "residues", "sentence": "Two chains.",
                "semantics": "Drag to turn it."},
     "chains": [{"node": "chainA", "tint": "mature3"},
                {"node": "chainB", "tint": "mature1"}]},
    {"prose": "hand"},          # provenance
    1,                          # resolver_version
)


class FakeCursor:
    def __init__(self, rows):
        self._rows = rows

    def fetchone(self):
        return self._rows[0] if self._rows else None

    def fetchall(self):
        return self._rows


class FakePool:
    """Answers each statement by the shape of its opening keyword."""

    def __init__(self, proteins=None, track_rows=None, failing=False):
        self.proteins = proteins or []
        self.track_rows = track_rows or []
        self.failing = failing
        self.statements = []

    def connection(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def execute(self, sql, params=()):
        if self.failing:
            raise OSError("connection closed")
        self.statements.append((sql, params))
        flat = " ".join(sql.split())
        if flat.startswith("select 1 from protein"):
            return FakeCursor([(1,)] if self.proteins else [])
        if "from protein_track" in flat and flat.startswith("select slug"):
            return FakeCursor(self.track_rows)
        if "from protein_track" in flat:
            return FakeCursor(self.track_rows)
        if "with matched" in flat:
            # The search adds a trailing `tier` column the summary ignores.
            return FakeCursor([row + (0,) for row in self.proteins])
        if "from protein p" in flat:
            return FakeCursor(self.proteins)
        return FakeCursor([])


@pytest.fixture
def fake_pool(monkeypatch):
    def _install(**kwargs):
        pool = FakePool(**kwargs)
        monkeypatch.setattr(db, "pool", pool)
        return pool

    return _install


@pytest.fixture
def no_pool(monkeypatch):
    monkeypatch.setattr(db, "pool", None)


# --- the catalog is unavailable, never empty --------------------------------


def test_no_database_reports_unavailable_not_an_empty_catalog(client, no_pool):
    """An unreadable catalog must not answer that nothing exists.

    There is nothing to read through to here, so an empty answer would tell a
    client that a protein does not exist when the truth is that nobody could
    look.
    """
    for path in ("/catalog", "/catalog/search?q=ins", "/protein/insulin",
                 "/protein/insulin/tracks"):
        assert client.get(path).status_code == 503, path


def test_a_broken_database_reports_unavailable(client, fake_pool):
    fake_pool(failing=True)
    assert client.get("/catalog").status_code == 503


def test_the_gene_route_is_unaffected_by_the_catalog(client, no_pool, genbank_text,
                                                     efetch_returning):
    """`/gene` keeps answering with no database, exactly as it did before."""
    efetch_returning(genbank_text)
    assert client.get("/gene/NG_007114/INS").status_code == 200


# --- pages ------------------------------------------------------------------


def test_a_page_carries_the_reading_order(client, fake_pool):
    fake_pool(proteins=[_protein_row()])
    body = client.get("/catalog").json()
    assert body["proteins"][0]["slug"] == "insulin"
    assert body["proteins"][0]["catalog_order"] == 0
    assert body["proteins"][0]["facts"] == {
        "residues": 110, "exons": 3, "chains": 3, "bridges": 3
    }
    assert body["next"] is None


def test_a_page_is_keyed_by_slug_so_a_cursor_stays_stable(client, fake_pool):
    """Ordered by the column the cursor compares, or rows shift under a client."""
    pool = fake_pool(proteins=[_protein_row()])
    client.get("/catalog?cursor=insulin")
    page_sql = " ".join(pool.statements[0][0].split())
    assert "order by p.slug" in page_sql
    assert "p.slug > %s" in page_sql


def test_a_full_page_offers_a_cursor(client, fake_pool):
    fake_pool(proteins=[_protein_row(slug="a"), _protein_row(slug="b")])
    body = client.get("/catalog?limit=1").json()
    assert len(body["proteins"]) == 1
    assert body["next"] == "a"


def test_the_catalog_is_the_list_and_built_proteins_are_not_on_it(client, fake_pool):
    """A protein built on demand has a row with no `catalog_order`.

    It is found through `/proteins/suggest`; the list stays the one someone
    chose, and so does searching it.
    """
    pool = fake_pool(proteins=[_protein_row()])
    client.get("/catalog")
    client.get("/catalog/search?q=ins")
    flat = [" ".join(sql.split()) for sql, _ in pool.statements]
    (page_sql,) = [sql for sql in flat if "from protein p" in sql and "order by p.slug" in sql]
    (search_sql,) = [sql for sql in flat if sql.startswith("with matched")]
    assert "where p.catalog_order is not null" in page_sql
    assert "where p.catalog_order is not null" in search_sql


# --- search -----------------------------------------------------------------


def test_search_matches_the_dart_tiers(client, fake_pool):
    """The five whole names `_closeness` compares, and its four tiers."""
    pool = fake_pool(proteins=[_protein_row()])
    client.get("/catalog/search?q=INS")
    sql, params = pool.statements[0]
    assert params["needle"] == "ins"
    assert params["prefix"] == "ins%"
    assert params["infix"] == "%ins%"
    assert params["whole"] == ["display", "gene", "slug", "uniprot", "accession"]
    assert "order by tier, p.catalog_order nulls last" in " ".join(sql.split())


def test_a_query_is_a_query_not_a_pattern(client, fake_pool):
    """Someone typing "p53%" is looking for a protein, not writing a LIKE."""
    pool = fake_pool(proteins=[])
    client.get("/catalog/search?q=p53%25")
    params = pool.statements[0][1]
    assert params["prefix"] == "p53\\%%"
    assert params["infix"] == "%p53\\%%"


def test_an_empty_query_returns_the_catalog(client, fake_pool):
    fake_pool(proteins=[_protein_row()])
    body = client.get("/catalog/search?q=  ").json()
    assert [p["slug"] for p in body["proteins"]] == ["insulin"]
    assert body["candidates"] == []


# --- one protein ------------------------------------------------------------


def test_detail_carries_what_targets_py_held(client, fake_pool):
    fake_pool(proteins=[_protein_row() + _DETAIL_TAIL])
    body = client.get("/protein/insulin").json()
    assert body["transcript_id"] == "NM_000207.3"
    assert body["disulfides"] == [[31, 96]]
    assert body["regions"][0]["short"] == "B"
    assert body["structure"]["pdb"] == "3I40"
    # The tints run in the order the record lists the mature peptides, which
    # is why insulin's B chain is first and its A chain third.
    assert body["chains"] == [{"node": "chainA", "tint": "mature3"},
                              {"node": "chainB", "tint": "mature1"}]
    assert body["provenance"] == {"prose": "hand"}


def test_an_unknown_slug_is_a_404_not_a_503(client, fake_pool):
    fake_pool(proteins=[])
    assert client.get("/protein/nosuchthing").status_code == 404


# --- tracks -----------------------------------------------------------------


def _track_row(kind="constraint", state="ready", reason=None,
               bucket="tracks", path="constraint/insulin.json.gz"):
    return (kind, state, reason, bucket, path, 161352, "abc123", "gzip",
            "json", {"model": "facebook/esm2_t33_650M_UR50D"})


def test_a_ready_track_becomes_a_public_url(client, fake_pool, monkeypatch):
    monkeypatch.setattr(settings, "supabase_url", "https://ref.supabase.co")
    fake_pool(proteins=[_protein_row()], track_rows=[_track_row()])
    body = client.get("/protein/insulin/tracks").json()
    constraint = body["constraint"]
    assert constraint["state"] == "ready"
    assert constraint["url"] == (
        "https://ref.supabase.co/storage/v1/object/public/tracks/"
        "constraint/insulin.json.gz"
    )
    assert constraint["sha256"] == "abc123"
    assert constraint["content_encoding"] == "gzip"
    assert constraint["provenance"]["model"] == "facebook/esm2_t33_650M_UR50D"


def test_every_kind_is_present_and_a_missing_row_is_absent(client, fake_pool,
                                                           monkeypatch):
    """A missing row and an explicit "nothing is coming" say the same thing."""
    monkeypatch.setattr(settings, "supabase_url", "https://ref.supabase.co")
    fake_pool(proteins=[_protein_row()], track_rows=[_track_row()])
    body = client.get("/protein/insulin/tracks").json()
    assert set(body) == set(tracks.KINDS)
    assert body["clinvar"]["state"] == "absent"
    assert body["clinvar"]["url"] is None


def test_pending_and_refused_are_distinct_states(client, fake_pool, monkeypatch):
    """R9.3: an unbaked gene and one whose bake is arriving are not one state."""
    monkeypatch.setattr(settings, "supabase_url", "https://ref.supabase.co")
    fake_pool(
        proteins=[_protein_row()],
        track_rows=[
            _track_row(kind="clinvar", state="pending", bucket=None, path=None),
            _track_row(kind="structure", state="refused", bucket=None, path=None,
                       reason="No entry covers the mature chain."),
        ],
    )
    body = client.get("/protein/insulin/tracks").json()
    assert body["clinvar"]["state"] == "pending"
    assert body["clinvar"]["reason"] is None
    assert body["structure"]["state"] == "refused"
    assert body["structure"]["reason"] == "No entry covers the mature chain."


def test_a_ready_track_this_service_cannot_address_is_a_fault(client, fake_pool,
                                                              monkeypatch):
    """Reporting ready without a URL would hand back a track nobody can fetch."""
    monkeypatch.setattr(settings, "supabase_url", None)
    fake_pool(proteins=[_protein_row()], track_rows=[_track_row()])
    body = client.get("/protein/insulin/tracks").json()
    assert body["constraint"]["state"] == "refused"
    assert body["constraint"]["url"] is None


def test_tracks_does_not_read_the_protein_row(client, fake_pool, monkeypatch):
    """The most-polled endpoint takes the cheapest existence check."""
    monkeypatch.setattr(settings, "supabase_url", "https://ref.supabase.co")
    pool = fake_pool(proteins=[_protein_row()], track_rows=[_track_row()])
    client.get("/protein/insulin/tracks")
    assert " ".join(pool.statements[0][0].split()).startswith("select 1 from protein")


def test_models_and_tracks_live_in_different_buckets():
    assert tracks.bucket_for("structure") == settings.models_bucket
    assert tracks.bucket_for("clinvar") == settings.tracks_bucket
