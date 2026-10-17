"""`sync_fixtures.py`, offline: it adds a curated row to the app's fixtures and
never changes what the app has committed. The app is a throwaway git checkout
and the tracks are stand-ins, so nothing is fetched; the last test holds the
local build of the twenty's rows to the app's own fixtures, where both are here.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

BACKEND = Path(__file__).resolve().parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import sync_fixtures  # noqa: E402
from pipeline.paths import CLIENT, DATA  # noqa: E402
from pipeline.targets import BY_SLUG, TARGETS  # noqa: E402

INSULIN, LYSOZYME = BY_SLUG["insulin"], BY_SLUG["lysozyme"]
KINDS = ["record", "constraint", "impact", "clinvar", "structure", "impact_explanations"]


def _row(slug: str, order: int) -> dict:
    return {"slug": slug, "display": slug.title(), "catalog_order": order,
            "tracks": {kind: "ready" for kind in KINDS}}


def _git(root: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True)


@pytest.fixture
def app(tmp_path):
    """An app checkout whose fixtures hold insulin, committed."""
    root = tmp_path / "app"
    (root / "test" / "fixtures" / "mock").mkdir(parents=True)
    (root / sync_fixtures.CATALOG).write_bytes(
        sync_fixtures.encode({"proteins": [_row("insulin", 0)], "next": None}))
    (root / "test" / "fixtures" / "mock" / "gene_ins.json").write_text("insulin's record")
    _git(root, "init", "-q")
    _git(root, "-c", "user.email=t@t", "-c", "user.name=t", "add", ".")
    _git(root, "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "fixtures")
    return root


def _stand_ins(monkeypatch, rows: dict, held: dict):
    monkeypatch.setattr(sync_fixtures, "live_row", lambda _s, target, kinds: {
        **rows[target.slug], "tracks": {k: "ready" for k in kinds}})
    monkeypatch.setattr(sync_fixtures, "payloads", lambda target, _src, _s: held[target.slug])


def test_a_new_row_and_its_files_are_added_after_the_others(app, monkeypatch):
    _stand_ins(monkeypatch, {"lysozyme": _row("lysozyme", 3)},
               {"lysozyme": {"record": b"lysozyme's record", "locus": b"its place",
                             "folding": b"its fold", "structure": b"glTF"}})
    assert sync_fixtures.main(["--target", "lysozyme", "--app", str(app)]) == 0
    catalog = json.loads((app / sync_fixtures.CATALOG).read_text())
    assert [row["slug"] for row in catalog["proteins"]] == ["insulin", "lysozyme"]
    assert catalog["next"] is None
    assert (app / "test/fixtures/mock/gene_lyz.json").read_bytes() == b"lysozyme's record"
    assert (app / "test/features/zoom/fixtures/lysozyme_locus.json").read_bytes() == b"its place"
    assert (app / "test/shared/folding/fixtures/lysozyme_folding.json").read_bytes() == b"its fold"
    assert (app / "test/fixtures/models/lysozyme.glb").read_bytes() == b"glTF"
    # Its own formatting, so the diff is the row and nothing else.
    assert (app / sync_fixtures.CATALOG).read_bytes() == sync_fixtures.encode(catalog)
    # Written again, uncommitted, it may change; and checking writes nothing.
    assert sync_fixtures.main(["--target", "lysozyme", "--app", str(app), "--check"]) == 0


def test_a_committed_row_or_file_is_never_changed(app, monkeypatch, capsys):
    _stand_ins(monkeypatch, {"insulin": {**_row("insulin", 0), "display": "Insulin!"}},
               {"insulin": {"record": b"insulin's record"}})
    before = (app / sync_fixtures.CATALOG).read_bytes()
    assert sync_fixtures.main(["--target", "insulin", "--app", str(app)]) == 1
    assert "the app has committed another row for it" in capsys.readouterr().err
    assert (app / sync_fixtures.CATALOG).read_bytes() == before

    _stand_ins(monkeypatch, {"insulin": _row("insulin", 0)},
               {"insulin": {"record": b"another record"}})
    assert sync_fixtures.main(["--target", "insulin", "--app", str(app)]) == 1
    assert "gene_ins.json: the app has committed other bytes" in capsys.readouterr().err
    assert (app / "test/fixtures/mock/gene_ins.json").read_text() == "insulin's record"


def test_a_check_finds_what_differs_and_writes_nothing(app, monkeypatch, capsys):
    _stand_ins(monkeypatch, {"lysozyme": _row("lysozyme", 3)},
               {"lysozyme": {"record": b"lysozyme's record"}})
    assert sync_fixtures.main(["--target", "lysozyme", "--app", str(app), "--check"]) == 1
    said = capsys.readouterr().err
    assert "differs: test/fixtures/mock/gene_lyz.json" in said
    assert not (app / "test/fixtures/mock/gene_lyz.json").exists()


@pytest.mark.skipif(not (CLIENT / sync_fixtures.CATALOG).exists()
                    or not all((DATA / t.mock_asset).exists() for t in TARGETS[:20]),
                    reason="needs the app's checkout and the stored tracks (fetch_tracks.py)")
def test_the_twentys_rows_built_here_are_the_apps_own():
    app = sync_fixtures.App(CLIENT)
    catalog = json.loads(app.committed(sync_fixtures.CATALOG))
    kinds = sync_fixtures.fixture_kinds(catalog)
    for fixture in catalog["proteins"]:
        built = sync_fixtures.local_row(BY_SLUG[fixture["slug"]], kinds)
        assert built == fixture, fixture["slug"]
        assert json.dumps(built) == json.dumps(fixture), f"{fixture['slug']}: key order"
