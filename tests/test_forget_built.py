"""`scripts/forget_built.py`: which objects are a protein's, and, against a real Postgres, that
deleting a protein built on demand takes its rows, its requests and its objects, and nothing
of any other protein's.

The Postgres half is opt-in, as `pipeline/resolver/test_worker_pg.py` is:

    RESOLVER_TEST_DATABASE_URL=postgresql://postgres@127.0.0.1:54329/postgres \\
        .venv/bin/pytest tests/test_forget_built.py
"""

from __future__ import annotations

import json
import os
import sys
import uuid
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[1]
for path in (str(BACKEND), str(BACKEND / "scripts")):
    if path not in sys.path:
        sys.path.insert(0, path)

import forget_built  # noqa: E402
from forget_built import Refused, objects_of, owned  # noqa: E402

ADMIN_URL = os.environ.get("RESOLVER_TEST_DATABASE_URL")
needs_postgres = pytest.mark.skipif(
    not ADMIN_URL, reason="set RESOLVER_TEST_DATABASE_URL to delete against Postgres")

MIGRATIONS = sorted((BACKEND / "migrations").glob("0*.sql"))


def test_an_object_is_a_proteins_by_its_whole_name_one_folder_deep():
    assert owned("b2m", "record/b2m.2b96c1622fb4.json")
    assert owned("b2m", "structure/b2m.4ad25d77e05f.glb")
    assert owned("hla-a", "fold_audit/hla-a.99cbee03e8e2.json")
    # Another protein whose name starts or ends with it.
    assert not owned("b2m", "record/b2mx.2b96c1622fb4.json")
    assert not owned("b2m", "record/ab2m.2b96c1622fb4.json")
    # Not the bakes' own shape.
    assert not owned("b2m", "b2m.2b96c1622fb4.json")
    assert not owned("b2m", "x/record/b2m.2b96c1622fb4.json")
    assert not owned("b2m", "record/b2m.2b96c1622fb4.json.bak")
    assert not owned("b2m", "record/b2m.json")


def test_its_rows_name_each_tracks_object_the_glb_and_every_generations():
    tracks = [
        {"bucket": "tracks", "object_path": "record/b2m.2b96c1622fb4.json", "provenance": {}},
        {"bucket": "models", "object_path": "structure/b2m.77272b304683.fsceneb",
         "provenance": {"glb": {"path": "structure/b2m.4ad25d77e05f.glb"}}},
        # A track with nothing stored names nothing.
        {"bucket": None, "object_path": None, "provenance": {}},
    ]
    generations = [
        {"structure_track": {"bucket": "models", "object_path": "structure/b2m.edc63541580a.fsceneb",
                             "provenance": {"glb": {"path": "structure/b2m.faba9e1f8751.glb"}}},
         "folding_track": None, "audit": None},
        {"structure_track": {"bucket": "models", "object_path": "structure/b2m.77272b304683.fsceneb",
                             "provenance": {"glb": {"path": "structure/b2m.4ad25d77e05f.glb"}}},
         "folding_track": {"bucket": "tracks", "object_path": "folding/b2m.26b8108b5cbb.json"},
         "audit": {"bucket": "tracks", "path": "fold_audit/b2m.99cbee03e8e2.json"}},
    ]
    assert objects_of("b2m", tracks, generations) == {
        ("tracks", "record/b2m.2b96c1622fb4.json"),
        ("models", "structure/b2m.77272b304683.fsceneb"),
        ("models", "structure/b2m.4ad25d77e05f.glb"),
        ("models", "structure/b2m.edc63541580a.fsceneb"),
        ("models", "structure/b2m.faba9e1f8751.glb"),
        ("tracks", "folding/b2m.26b8108b5cbb.json"),
        ("tracks", "fold_audit/b2m.99cbee03e8e2.json"),
    }


def test_a_row_naming_another_proteins_object_refuses_the_lot():
    tracks = [{"bucket": "tracks", "object_path": "record/insulin.2b96c1622fb4.json",
               "provenance": {}}]
    with pytest.raises(Refused, match="insulin"):
        objects_of("b2m", tracks, [])


def test_genes_or_all_never_both():
    with pytest.raises(SystemExit):
        forget_built.main(["B2M", "--all"])
    with pytest.raises(SystemExit):
        forget_built.main([])


# --------------------------------------------------------------- against Postgres


class FakeStorage:
    """Storage as `forget_built.Storage` reads it: folders, names under them, and deletes."""

    def __init__(self, held):
        self.held = set(held)  # {(bucket, path)}
        self.deleted = []

    def folders(self, bucket):
        return sorted({p.split("/")[0] for b, p in self.held if b == bucket})

    def objects(self, bucket, folder, start):
        return sorted(p for b, p in self.held
                      if b == bucket and p.startswith(folder + "/" + start))

    def delete(self, bucket, paths):
        self.deleted.extend((bucket, p) for p in paths)
        self.held -= {(bucket, p) for p in paths}


@pytest.fixture(scope="module")
def database_url():
    import psycopg
    from psycopg.conninfo import conninfo_to_dict, make_conninfo

    name = "helixpeak_forget_{}".format(uuid.uuid4().hex[:8])
    with psycopg.connect(ADMIN_URL, autocommit=True) as admin:
        admin.execute("create database {}".format(name))
    url = make_conninfo(**{**conninfo_to_dict(ADMIN_URL), "dbname": name})
    with psycopg.connect(url, autocommit=True) as conn:
        for migration in MIGRATIONS:
            conn.execute(migration.read_text())
    yield url
    with psycopg.connect(ADMIN_URL, autocommit=True) as admin:
        admin.execute("drop database {} with (force)".format(name))


def _protein(conn, slug, gene, order=None):
    conn.execute(
        """insert into protein (slug, gene, uniprot, accession, display, summary, residues,
               exons, chains, bridges, regions, disulfides, provenance, resolver_version,
               catalog_order)
           values (%s, %s, 'P00000', 'NG_000000.1', %s, '', 100, 1, 1, 0, '[]', '[]', '{}',
                   %s, %s)""",
        (slug, gene, gene, 0 if order is not None else 1, order))
    conn.execute("insert into protein_alias (slug, alias, kind) values (%s, %s, 'gene')",
                 (slug, gene))


def _track(conn, slug, kind, bucket, path, provenance=None):
    conn.execute(
        """insert into protein_track (slug, kind, state, bucket, object_path, format, provenance)
           values (%s, %s, 'ready', %s, %s, 'json', %s)""",
        (slug, kind, bucket, path, _jsonb(provenance or {})))


def _jsonb(value):
    from psycopg.types.json import Jsonb

    return Jsonb(value)


@pytest.fixture
def conn(database_url):
    from pipeline.resolver import store

    with store.connect(database_url) as connection:
        connection.execute(
            "truncate resolve_request, bake_job, protein_track, protein_alias, protein, "
            "fold_selection, fold_generation restart identity cascade")
        # A curated protein, with what a curated protein has.
        _protein(connection, "insulin", "INS", order=0)
        _track(connection, "insulin", "record", "tracks", "record/insulin.aaaaaaaaaaaa.json")
        # Two built on demand: B2M with every table filled, TTR beside it.
        _protein(connection, "b2m", "B2M")
        _track(connection, "b2m", "record", "tracks", "record/b2m.2b96c1622fb4.json")
        _track(connection, "b2m", "folding", "tracks", "folding/b2m.26b8108b5cbb.json")
        _track(connection, "b2m", "structure", "models", "structure/b2m.77272b304683.fsceneb",
               {"glb": {"path": "structure/b2m.4ad25d77e05f.glb"}})
        job = connection.execute(
            "insert into bake_job (slug, kind, state) values ('b2m', 'structure', 'done') "
            "returning id").fetchone()[0]
        selection = connection.execute(
            """insert into fold_selection (slug, job_id, policy, inputs, record_sha256, outcome,
                   candidate, entry)
               values ('b2m', %s, 1, 'x', 'y', 'chosen', '1K5N', '1K5N') returning id""",
            (job,)).fetchone()[0]
        connection.execute(
            """insert into fold_generation (slug, source, state, structure_track, active)
               values ('b2m', 'alphafold', 'ready', %s, false)""",
            (_jsonb({"bucket": "models", "object_path": "structure/b2m.edc63541580a.fsceneb",
                     "provenance": {"glb": {"path": "structure/b2m.faba9e1f8751.glb"}}}),))
        connection.execute(
            """insert into fold_generation (slug, source, state, job_id, selection_id,
                   structure_track, folding_track, audit, active)
               values ('b2m', 'pdb', 'ready', %s, %s, %s, %s, %s, true)""",
            (job, selection,
             _jsonb({"bucket": "models", "object_path": "structure/b2m.77272b304683.fsceneb",
                     "provenance": {"glb": {"path": "structure/b2m.4ad25d77e05f.glb"}}}),
             _jsonb({"bucket": "tracks", "object_path": "folding/b2m.26b8108b5cbb.json"}),
             _jsonb({"bucket": "tracks", "path": "fold_audit/b2m.99cbee03e8e2.json"})))
        connection.execute(
            "insert into resolve_request (gene, uniprot, slug, state) "
            "values ('B2M', 'P61769', 'b2m', 'done')")
        _protein(connection, "ttr", "TTR")
        _track(connection, "ttr", "record", "tracks", "record/ttr.111111111111.json")
        connection.execute(
            "insert into resolve_request (gene, uniprot, slug, state) "
            "values ('TTR', 'P02766', 'ttr', 'done')")
        yield connection


def _storage():
    return FakeStorage({
        ("tracks", "record/insulin.aaaaaaaaaaaa.json"),
        ("tracks", "record/b2m.2b96c1622fb4.json"),
        ("tracks", "folding/b2m.26b8108b5cbb.json"),
        ("tracks", "fold_audit/b2m.99cbee03e8e2.json"),
        # Left by a re-bake: no row names it.
        ("tracks", "constraint/b2m.0123456789ab.json"),
        ("models", "structure/b2m.77272b304683.fsceneb"),
        ("models", "structure/b2m.4ad25d77e05f.glb"),
        ("models", "structure/b2m.edc63541580a.fsceneb"),
        ("models", "structure/b2m.faba9e1f8751.glb"),
        # Another protein whose name starts with it.
        ("tracks", "record/b2mx.222222222222.json"),
        ("tracks", "record/ttr.111111111111.json"),
    })


def _count(conn, table, where, *params):
    return conn.execute("select count(*) from {} where {}".format(table, where), params).fetchone()[0]


@needs_postgres
def test_a_built_protein_goes_whole_and_nothing_else_does(conn):
    storage = _storage()
    said = []
    assert forget_built.forget(conn, storage, ["B2M"], False, False, said.append) == 1

    for table in ("protein", "protein_alias", "protein_track", "bake_job", "fold_selection",
                  "fold_generation"):
        assert _count(conn, table, "slug = 'b2m'") == 0, table
    assert _count(conn, "resolve_request", "gene = 'B2M'") == 0
    assert not {o for o in storage.held if "/b2m." in o[1]}
    assert ("tracks", "constraint/b2m.0123456789ab.json") in storage.deleted
    # The curated protein and the other built one are as they were.
    assert _count(conn, "protein", "slug in ('insulin', 'ttr')") == 2
    assert _count(conn, "protein_track", "slug in ('insulin', 'ttr')") == 2
    assert _count(conn, "resolve_request", "gene = 'TTR'") == 1
    assert storage.held == {
        ("tracks", "record/insulin.aaaaaaaaaaaa.json"),
        ("tracks", "record/b2mx.222222222222.json"),
        ("tracks", "record/ttr.111111111111.json"),
    }
    assert said[-1] == "Deleted 8 objects."


@needs_postgres
def test_all_is_every_built_protein_and_never_a_curated_one(conn):
    storage = _storage()
    assert forget_built.forget(conn, storage, [], True, False, lambda line: None) == 2
    assert [r[0] for r in conn.execute("select slug from protein")] == ["insulin"]
    assert _count(conn, "resolve_request", "true") == 0
    assert storage.held == {
        ("tracks", "record/insulin.aaaaaaaaaaaa.json"),
        ("tracks", "record/b2mx.222222222222.json"),
    }


@needs_postgres
def test_a_dry_run_deletes_nothing(conn):
    storage = _storage()
    said = []
    assert forget_built.forget(conn, storage, [], True, True, said.append) == 0
    assert _count(conn, "protein", "true") == 3
    assert storage.deleted == []
    assert said[-1] == "Dry run: nothing deleted (2 proteins, 9 objects)."


@needs_postgres
def test_a_curated_or_unknown_gene_is_refused_before_anything_is_written(conn):
    storage = _storage()
    with pytest.raises(Refused, match="curated, never deleted here: INS"):
        forget_built.forget(conn, storage, ["B2M", "INS"], False, False, lambda line: None)
    with pytest.raises(Refused, match="no protein for BRCA1"):
        forget_built.forget(conn, storage, ["BRCA1"], False, False, lambda line: None)
    assert _count(conn, "protein", "true") == 3
    assert storage.deleted == []


@needs_postgres
def test_a_build_under_way_is_refused_and_nothing_is_written(conn):
    conn.execute("insert into bake_job (slug, kind, state) values ('ttr', 'constraint', 'running')")
    storage = _storage()
    with pytest.raises(Refused, match="constraint bake of ttr running"):
        forget_built.forget(conn, storage, [], True, False, lambda line: None)
    assert _count(conn, "protein", "true") == 3
    assert _count(conn, "resolve_request", "true") == 2
    assert storage.deleted == []
