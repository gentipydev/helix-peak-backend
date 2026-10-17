"""`fetch_tracks.py --check`, offline: the service's answers are stubbed and the
data directory is a temporary one, so nothing is fetched or written.
"""

from __future__ import annotations

import hashlib
import sys
import urllib.error
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import fetch_tracks  # noqa: E402
from pipeline.targets import BY_SLUG  # noqa: E402

INSULIN, LYSOZYME = BY_SLUG["insulin"], BY_SLUG["lysozyme"]


def _digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _stored(monkeypatch, held: dict):
    def stored_objects(_service, target):
        answer = held[target.slug]
        if isinstance(answer, Exception):
            raise answer
        return answer

    monkeypatch.setattr(fetch_tracks, "stored_objects", stored_objects)


def _write(data: Path, kind: str, target, payload: bytes) -> None:
    path = data / fetch_tracks.asset_path(kind, target)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def test_a_file_the_same_as_stored_passes_and_nothing_is_written(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(fetch_tracks, "DATA", tmp_path)
    _write(tmp_path, "record", INSULIN, b"record")
    _write(tmp_path, "structure", INSULIN, b"model")
    _stored(monkeypatch, {"insulin": {"record": ("u", _digest(b"record")),
                                      "structure": ("u", _digest(b"model"))}})
    before = sorted(p.relative_to(tmp_path) for p in tmp_path.rglob("*"))
    assert fetch_tracks.check("s", [INSULIN], fetch_tracks.KINDS) == 0
    assert sorted(p.relative_to(tmp_path) for p in tmp_path.rglob("*")) == before
    assert "2 the same as stored, 0 differ, 0 missing, 0 unread" in capsys.readouterr().err


def test_a_file_that_differs_or_is_missing_fails(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(fetch_tracks, "DATA", tmp_path)
    _write(tmp_path, "record", INSULIN, b"re-baked")
    _stored(monkeypatch, {"insulin": {"record": ("u", _digest(b"record")),
                                      "folding": ("u", _digest(b"fold"))}})
    assert fetch_tracks.check("s", [INSULIN], fetch_tracks.KINDS) == 1
    said = capsys.readouterr().err
    assert "insulin record: differs from its row's sha256" in said
    assert "insulin folding: no local file" in said
    # Only the kinds asked for.
    assert fetch_tracks.check("s", [INSULIN], ("folding",)) == 1
    assert fetch_tracks.check("s", [INSULIN], ("structure",)) == 0


def test_a_row_not_in_the_catalog_yet_is_said_and_passed_by(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(fetch_tracks, "DATA", tmp_path)
    _write(tmp_path, "record", LYSOZYME, b"record")
    unseeded = urllib.error.HTTPError("u", 404, "Not Found", None, None)
    _stored(monkeypatch, {"insulin": unseeded,
                          "lysozyme": {"record": ("u", _digest(b"record"))}})
    assert fetch_tracks.check("s", [INSULIN, LYSOZYME], fetch_tracks.KINDS) == 0
    assert "insulin: not in the catalog yet" in capsys.readouterr().err
    # Any other failure to read is a failure.
    _stored(monkeypatch, {"insulin": urllib.error.HTTPError("u", 502, "Bad", None, None),
                          "lysozyme": {}})
    assert fetch_tracks.check("s", [INSULIN, LYSOZYME], fetch_tracks.KINDS) == 1


def test_a_fetch_passes_a_row_not_in_the_catalog_yet(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(fetch_tracks, "DATA", tmp_path)
    _stored(monkeypatch, {"insulin": urllib.error.HTTPError("u", 404, "Not Found", None, None)})
    assert fetch_tracks.fetch("s", [INSULIN], fetch_tracks.KINDS) == 0
    assert "insulin: not in the catalog yet" in capsys.readouterr().err
