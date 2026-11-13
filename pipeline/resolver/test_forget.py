"""`forget.purge_discarded`: the objects of a protein the service discarded go,
and nothing else does -- not another protein's, not one a row names now, not
one made after the discard (the gene built again since), and not one whose
time storage does not say. Offline: a fake connection and a fake storage.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

from pipeline.resolver import forget

AT = datetime(2026, 11, 13, 16, 41, 30, tzinfo=timezone.utc)


def _stamp(when: datetime) -> str:
    return when.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"


class Conn:
    """Answers the purge's statements: whether 0014 is there, the discards not
    purged, what a protein's rows name now, and the purge's own update."""

    def __init__(self, discards=(), ready=True, tracks=(), generations=()):
        self.discards = list(discards)
        self.ready = ready
        self.tracks = list(tracks)
        self.generations = list(generations)
        self.purged = []

    def execute(self, sql, params=None):
        flat = " ".join(sql.split())
        if flat.startswith("select to_regclass('discard_request')"):
            return Rows([(self.ready,)])
        if flat.startswith("select id, slug, requested_at from discard_request"):
            return Rows(self.discards)
        if flat.startswith("select bucket, object_path, provenance from protein_track"):
            return Rows(self.tracks)
        if flat.startswith("select structure_track, folding_track, audit from fold_generation"):
            return Rows(self.generations)
        if flat.startswith("update discard_request set purged_at = now()"):
            self.purged.append(params)
            return Rows([])
        raise AssertionError(f"unexpected statement: {flat}")


class Rows(list):
    def fetchone(self):
        return self[0] if self else None

    def fetchall(self):
        return list(self)


class Storage:
    """Folders, and in each the objects with when they were made; what is deleted."""

    def __init__(self, objects, failing=False):
        # {(bucket, path): created_at stamp or None}
        self.made = dict(objects)
        self.failing = failing
        self.deleted = []

    def folders(self, bucket):
        if self.failing:
            raise OSError("storage unreachable")
        return sorted({path.split("/")[0] for b, path in self.made if b == bucket})

    def stamped(self, bucket, folder, start):
        return [(path, forget._made(stamp)) for (b, path), stamp in sorted(self.made.items())
                if b == bucket and path.startswith(folder + "/" + start)]

    def delete(self, bucket, paths):
        self.deleted.append((bucket, list(paths)))


BEFORE = _stamp(AT - timedelta(minutes=20))
AFTER = _stamp(AT + timedelta(minutes=5))


def test_a_discarded_proteins_objects_go_and_the_discard_is_marked_purged():
    conn = Conn(discards=[(3, "rbp4", AT)])
    storage = Storage({
        ("tracks", "record/rbp4.6e89b1dd7c5f.json"): BEFORE,
        ("tracks", "impact/rbp4.c6aa14908e05.json"): BEFORE,
        ("models", "structure/rbp4.b6da3fc2d242.fsceneb"): BEFORE,
        ("models", "structure/rbp4.f09c3651b1be.glb"): BEFORE,
        # Another protein's, whose name starts the same, and one not named as a bake's.
        ("tracks", "record/rbp4x.0123456789ab.json"): BEFORE,
        ("tracks", "record/rbp4.notes.json"): BEFORE,
    })
    assert forget.purge_discarded(conn, storage) == [{"slug": "rbp4", "objects": 4}]
    assert storage.deleted == [
        ("tracks", ["impact/rbp4.c6aa14908e05.json", "record/rbp4.6e89b1dd7c5f.json"]),
        ("models", ["structure/rbp4.b6da3fc2d242.fsceneb", "structure/rbp4.f09c3651b1be.glb"]),
    ]
    assert conn.purged == [(4, 3)]


def test_a_build_made_after_the_discard_keeps_every_object():
    # The gene built again: its new objects are newer than the discard, and its
    # rows name them; an object whose time storage does not say is left too.
    conn = Conn(discards=[(3, "rbp4", AT)],
                tracks=[("tracks", "record/rbp4.aaaaaaaaaaaa.json", {})])
    storage = Storage({
        ("tracks", "record/rbp4.6e89b1dd7c5f.json"): BEFORE,
        ("tracks", "record/rbp4.aaaaaaaaaaaa.json"): BEFORE,
        ("tracks", "impact/rbp4.bbbbbbbbbbbb.json"): AFTER,
        ("tracks", "clinvar/rbp4.cccccccccccc.json"): None,
    })
    assert forget.purge_discarded(conn, storage) == [{"slug": "rbp4", "objects": 1}]
    assert storage.deleted == [("tracks", ["record/rbp4.6e89b1dd7c5f.json"])]


def test_nothing_is_done_before_0014_or_with_nothing_discarded():
    for conn in (Conn(ready=False, discards=[(3, "rbp4", AT)]), Conn()):
        storage = Storage({("tracks", "record/rbp4.6e89b1dd7c5f.json"): BEFORE})
        assert forget.purge_discarded(conn, storage) == []
        assert storage.deleted == [] and conn.purged == []


def test_storage_out_of_reach_leaves_the_discard_for_the_next_sweep():
    conn = Conn(discards=[(3, "rbp4", AT)])
    storage = Storage({("tracks", "record/rbp4.6e89b1dd7c5f.json"): BEFORE}, failing=True)
    assert forget.purge_discarded(conn, storage) == []
    assert conn.purged == []


def test_storage_times_are_read_whatever_their_fraction_and_zone():
    assert forget._made("2026-11-13T16:15:29.054Z") == datetime(
        2026, 11, 13, 16, 15, 29, 54000, tzinfo=timezone.utc)
    assert forget._made("2026-11-13T16:15:29Z") == datetime(
        2026, 11, 13, 16, 15, 29, tzinfo=timezone.utc)
    assert forget._made("2026-11-13T17:15:29.1234567+01:00") == datetime(
        2026, 11, 13, 16, 15, 29, 123456, tzinfo=timezone.utc)
    assert forget._made(None) is None and forget._made("yesterday") is None


def test_a_proteins_locus_is_forgotten_with_it():
    # The locus a built protein is given (`worker.make_locus`) is stored as the
    # uploader names every track, so its row and its name both make it the
    # protein's, and `forget_built.py` deletes it with the rest.
    locus = {"slug": "b2m", "kind": "locus", "bucket": "tracks",
             "object_path": "locus/b2m.0123456789ab.json"}
    assert forget.owned("b2m", "locus/b2m.0123456789ab.json")
    assert not forget.owned("b2m", "locus/b2mx.0123456789ab.json")
    assert forget.objects_of("b2m", [locus], []) == {("tracks", "locus/b2m.0123456789ab.json")}
    # A refused locus has no object to forget.
    refused = {"slug": "b2m", "kind": "locus", "bucket": None, "object_path": None}
    assert forget.objects_of("b2m", [refused], []) == set()
