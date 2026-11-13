"""Which objects in storage are a protein's, and deleting a discarded protein's.

Shared by `scripts/forget_built.py`, which deletes proteins built on demand by
hand, rows and objects, and by the resolver's sweep (`worker.sweep`), which
deletes the objects of each protein the service discarded (`resolves.discard`,
`discard_request`, 0014). The service deletes a discarded protein's rows at
once but holds no storage key, so its objects wait here for the next sweep.

An object is a protein's only by its name, `<folder>/<slug>.<sha12>.<ext>`
(`owned`), in the `tracks` and `models` buckets, as every bake writes it.

A discarded protein's objects go only where no row names them now and they
were made before the discard (`purge_discarded`): the same gene built again
after its discard takes the same slug, and keeps every object of its new build.
"""

from __future__ import annotations

import json
import logging
import re
import urllib.request
from datetime import datetime, timezone
from typing import Dict, Iterable, List, Optional, Set, Tuple

BUCKETS = ("tracks", "models")

Object = Tuple[str, str]  # (bucket, path)

log = logging.getLogger("resolver-worker")


class Refused(Exception):
    """Nothing was written: the sentence says why."""


def owned(slug: str, path: str) -> bool:
    """Whether `path` is an object of `slug`'s: `<folder>/<slug>.<sha12>.<ext>`, exactly.

    The name is the protein's own, never a longer one that starts with it
    (`b2m` is not `b2mx`), and sits one folder deep, as every bake writes it."""
    return re.fullmatch(r"[a-z0-9_-]+/" + re.escape(slug) + r"\.[0-9a-f]{12}\.[a-z]+", path) is not None


def _named(holder: Optional[dict], key: str = "object_path") -> List[Object]:
    """The objects one row of a track, as stored or as a generation keeps it, names: its
    own and, for a structure, the `.glb` it was compiled from, in the same bucket."""
    if not holder or not holder.get("bucket"):
        return []
    found = []
    if holder.get(key):
        found.append((holder["bucket"], holder[key]))
    glb = ((holder.get("provenance") or {}).get("glb") or {}).get("path")
    if glb:
        found.append((holder["bucket"], glb))
    return found


def objects_of(slug: str, tracks: Iterable[dict], generations: Iterable[dict]) -> Set[Object]:
    """Every object `slug`'s rows name. Raises Refused on one that is not its own, which
    would be another protein's object."""
    found = set()  # type: Set[Object]
    for track in tracks:
        found.update(_named(track))
    for generation in generations:
        found.update(_named(generation.get("structure_track")))
        found.update(_named(generation.get("folding_track")))
        found.update(_named(generation.get("audit"), key="path"))
    strangers = sorted("{}/{}".format(*o) for o in found if not owned(slug, o[1]))
    if strangers:
        raise Refused("{}'s rows name objects that are not its own: {}".format(
            slug, ", ".join(strangers)))
    return found


def _made(stamp: Optional[str]) -> Optional[datetime]:
    """Storage's `created_at`, `2026-11-13T16:15:29.054Z`, as a time; None where it
    has none or one this cannot read, which no purge then takes."""
    if not stamp:
        return None
    found = re.fullmatch(r"(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d)(?:\.(\d+))?(Z|[+-]\d\d:\d\d)?", stamp)
    if found is None:
        return None
    whole, fraction, zone = found.groups()
    micro = (fraction or "0")[:6].ljust(6, "0")
    offset = "+00:00" if zone in (None, "Z") else zone
    return datetime.fromisoformat("{}.{}{}".format(whole, micro, offset)).astimezone(timezone.utc)


class Storage:
    """Lists and deletes objects in Supabase storage, with the service key."""

    def __init__(self, url: str, key: str, uploads=None):
        if uploads is None:
            from pipeline import upload_tracks

            uploads = upload_tracks.Storage(url, key)
        self._uploads = uploads
        self._folders = {}  # type: Dict[str, List[str]]

    @classmethod
    def over(cls, uploads) -> "Storage":
        """Listing over a client the caller holds already (`store.TrackStorage`)."""
        return cls("", "", uploads=uploads)

    def _list(self, bucket: str, prefix: str, search: str = "") -> List[dict]:
        found, offset = [], 0  # type: List[dict], int
        while True:
            body = json.dumps({"prefix": prefix, "search": search, "limit": 100, "offset": offset,
                               "sortBy": {"column": "name", "order": "asc"}}).encode()
            request = urllib.request.Request(
                "{}/object/list/{}".format(self._uploads.base, bucket), data=body, method="POST",
                headers={**self._uploads.headers, "Content-Type": "application/json"})
            with urllib.request.urlopen(request, timeout=120) as response:
                page = json.loads(response.read())
            found.extend(page)
            if len(page) < 100:
                return found
            offset += 100

    def folders(self, bucket: str) -> List[str]:
        if bucket not in self._folders:
            self._folders[bucket] = [o["name"] for o in self._list(bucket, "") if o.get("id") is None]
        return self._folders[bucket]

    def objects(self, bucket: str, folder: str, start: str) -> List[str]:
        """The names in `folder` of `bucket` that start with `start`, as paths."""
        return [path for path, _made_at in self.stamped(bucket, folder, start)]

    def stamped(self, bucket: str, folder: str,
                start: str) -> List[Tuple[str, Optional[datetime]]]:
        """As `objects`, each with when storage says it was made."""
        return [("{}/{}".format(folder, o["name"]), _made(o.get("created_at")))
                for o in self._list(bucket, folder + "/", start) if o.get("id") is not None]

    def delete(self, bucket: str, paths: List[str]) -> None:
        self._uploads.delete(bucket, paths)


def lister(storage):
    """`storage` where it lists and deletes as `Storage` does already, else a
    `Storage` over it: the worker holds the uploader's client, which only puts."""
    return storage if hasattr(storage, "stamped") else Storage.over(storage)


def stored(storage, slug: str) -> Set[Object]:
    """Every object in storage named as `slug`'s, whatever row names it or none."""
    found = set()  # type: Set[Object]
    for bucket in BUCKETS:
        for folder in storage.folders(bucket):
            found.update((bucket, path) for path in storage.objects(bucket, folder, slug + ".")
                         if owned(slug, path))
    return found


def named(conn, slug: str) -> Set[Object]:
    """Every object `slug`'s rows name now."""
    tracks = [{"bucket": bucket, "object_path": path, "provenance": provenance}
              for bucket, path, provenance in conn.execute(
                  "select bucket, object_path, provenance from protein_track where slug = %s",
                  (slug,))]
    generations = [{"structure_track": s, "folding_track": f, "audit": a}
                   for s, f, a in conn.execute(
                       "select structure_track, folding_track, audit from fold_generation "
                       "where slug = %s", (slug,))]
    return objects_of(slug, tracks, generations)


def discards_ready(conn) -> bool:
    """Whether the database holds the discards to purge (0014)."""
    return conn.execute("select to_regclass('discard_request') is not null").fetchone()[0]


_UNPURGED = """
select id, slug, requested_at from discard_request
where purged_at is null and slug is not null
order by requested_at, id
limit %s
"""
_PURGED = """
update discard_request set purged_at = now(), objects = %s where id = %s
"""


def purge_discarded(conn, storage, limit: int = 10) -> List[dict]:
    """Delete the objects of proteins the service discarded, oldest discard first.

    Of each one's objects, only those that no row names now and that storage
    says were made before the discard: the gene built again since keeps its
    new build whole. Each discard is marked purged with how many went. One that
    could not be purged (storage unreachable) is said and left for the next
    sweep; nothing raised here keeps a sweep from resolving.

    `storage` lists and deletes as `Storage` does. Returns each discard purged,
    `{"slug", "objects"}`.
    """
    if not discards_ready(conn):
        return []
    done = []  # type: List[dict]
    for job, slug, asked_at in conn.execute(_UNPURGED, (limit,)).fetchall():
        try:
            kept = named(conn, slug)
            gone = sorted(
                (bucket, path)
                for bucket in BUCKETS
                for folder in storage.folders(bucket)
                for path, made in storage.stamped(bucket, folder, slug + ".")
                if owned(slug, path) and (bucket, path) not in kept
                and made is not None and made < asked_at)
            for bucket in BUCKETS:
                paths = [path for b, path in gone if b == bucket]
                for start in range(0, len(paths), 100):
                    storage.delete(bucket, paths[start:start + 100])
            conn.execute(_PURGED, (len(gone), job))
        except Exception as exc:  # noqa: BLE001 - said, and tried again next sweep
            log.warning("%s: its discarded objects could not be deleted now (%s: %s); "
                        "the next sweep tries again", slug, type(exc).__name__, exc)
            continue
        log.info("%s: discarded, %d object(s) deleted from storage", slug, len(gone))
        done.append({"slug": slug, "objects": len(gone)})
    return done
