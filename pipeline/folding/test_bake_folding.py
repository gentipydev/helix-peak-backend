"""The folding bake and its check, offline.

The entries are in the repository. Tests that build a payload need the stored
record and model, and skip until `fetch_tracks.py --kind record --kind
structure` has put them in `pipeline/data/`. They place the payload in the
frame the stored model's bridges give exactly, so no PyMOL is needed here; the
bake itself finds the same frame with the structure bake's own export.
"""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import numpy as np
import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import fetch_tracks, upload_tracks  # noqa: E402
from pipeline.folding import bake_folding, check_folding  # noqa: E402
from pipeline.folding.bake_folding import Frame  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.structure.frame import fit_size, model_frame  # noqa: E402
from pipeline.structure.glb import read_glb  # noqa: E402
from pipeline.structure.pdb import ca_atoms, chain_residues  # noqa: E402
from pipeline.targets import BY_SLUG  # noqa: E402


def _stored(slug: str) -> bool:
    target = BY_SLUG[slug]
    return (DATA / target.mock_asset).exists() and (DATA / target.structure_asset).exists()


needs = pytest.mark.skipif(
    not all(_stored(s) for s in ("prion", "app", "insulin", "glucagon")),
    reason="run fetch_tracks.py --kind record --kind structure first")


def _track(slug: str, frame: Frame | None = None) -> dict:
    """The payload, in the exact frame the stored model's bridges give."""
    target = BY_SLUG[slug]
    glb = (DATA / target.structure_asset).read_bytes()
    if frame is None:
        fit, method = model_frame(target.structure, read_glb(glb))
        assert method.startswith("disulfide")
        frame = Frame(centre=fit.centre, length=fit.length, method=method, matched=None)
    record = json.loads((DATA / target.mock_asset).read_bytes())
    return bake_folding.payload(target, record, glb, frame)


# -- reading an entry against its record ---------------------------------------


def test_the_entrys_numbering_is_found_by_its_letters():
    translation = "M" + "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAV"
    residues = chain_residues(BY_SLUG["hemoglobin"].structure, "B")[:40]
    # 2DN1 numbers beta from Val1, after the initiator methionine: precursor 2.
    assert bake_folding.numbering_offset(residues, translation) == 1


def test_a_repeat_is_read_as_the_first_of_its_kind():
    one = chain_residues(BY_SLUG["ubiquitin"].structure, "A")
    letters = "".join(r.letter for r in one)
    assert bake_folding.numbering_offset(one, letters * 3 + "C") == 0


def test_the_mature_chain_is_the_chain_not_a_domain():
    # Myoglobin's region table names the globin fold (2-148), a domain; the
    # chain is everything after the initiator methionine.
    assert bake_folding.mature_span(BY_SLUG["myoglobin"], [3, 150]) == (2, 154)
    insulin = BY_SLUG["insulin"]
    assert bake_folding.mature_span(insulin, [90, 110]) == (90, 110)
    assert bake_folding.mature_span(insulin, [25, 54]) == (25, 54)
    assert bake_folding.mature_span(BY_SLUG["tnf"], [81, 233]) == (77, 233)
    with pytest.raises(ValueError):
        bake_folding.mature_span(insulin, [40, 100])      # no one chain holds both


# -- the payload -----------------------------------------------------------------


@needs
def test_the_prions_first_hundred_residues_have_no_place():
    (chain,) = _track("prion")["chains"]
    states = {r["n"]: r["state"] for r in chain["residues"]}
    assert (chain["first"], chain["last"]) == (23, 230)
    assert len(chain["residues"]) == 208
    assert states[23] == "absent"          # before what was crystallised
    assert all(states[n] == "disordered" for n in range(24, 117))
    assert all(states[n] == "ordered" for n in range(117, 226))
    assert all(states[n] == "disordered" for n in range(226, 231))
    assert sum(1 for s in states.values() if s != "ordered") == 99


@needs
def test_what_was_never_crystallised_is_absent_not_disordered():
    (chain,) = _track("app")["chains"]
    states = [r["state"] for r in chain["residues"]]
    assert (chain["first"], chain["last"], len(states)) == (18, 770, 753)
    assert states.count("absent") == 580 and states[-580:] == ["absent"] * 580
    assert states.count("disordered") == 11


@needs
def test_an_ordered_residue_is_the_entrys_ca_in_the_models_frame():
    track = _track("insulin")
    target = BY_SLUG["insulin"]
    placed = np.array([r["ca"] for c in track["chains"] for r in c["residues"]])
    frame = track["frame"]
    back = placed * frame["length_angstrom"] + np.array(frame["centre_angstrom"])
    assert np.abs(back - ca_atoms(target.structure)).max() < 1e-3
    labels = {r["ss"] for c in track["chains"] for r in c["residues"]}
    assert labels == {"helix", "coil"}


@needs
def test_the_letters_are_the_records_and_the_entrys_are_noted():
    (b_chain,) = [c for c in _track("insulin")["chains"] if c["node"] == "chainB"]
    # 3I40 carries Ala at B30 where the gene makes Thr: the letter is the gene's.
    assert b_chain["residues"][-1]["aa"] == "T"
    assert b_chain["entry_differs"] == [{"n": 54, "entry": "A", "record": "T"}]


@needs
def test_one_residue_a_line_and_it_reads_back():
    track = _track("prion")
    text = bake_folding.encode(track).decode()
    assert json.loads(text) == track
    lines = [line.strip() for line in text.splitlines()]
    assert '{"n": 23, "aa": "K", "state": "absent"},' in lines
    assert sum(1 for line in lines if line.startswith('{"n": ')) == 208


# -- the check -------------------------------------------------------------------


@needs
def test_a_right_track_passes():
    for slug in ("prion", "app", "insulin"):
        assert check_folding.problems_of(BY_SLUG[slug], _track(slug)) == [], slug


@needs
def test_a_dropped_residue_is_caught():
    track = _track("prion")
    del track["chains"][0]["residues"][150]
    assert any("mature chain" in p for p in check_folding.problems_of(BY_SLUG["prion"], track))


@needs
def test_a_hole_in_the_chain_is_caught():
    track = _track("prion")
    residue = track["chains"][0]["residues"][150]
    for key in ("ca", "ss"):
        residue.pop(key)
    residue["state"] = "absent"
    found = check_folding.problems_of(BY_SLUG["prion"], track)
    assert any("absent residue between" in p for p in found)


@needs
def test_a_wrong_letter_is_caught():
    track = _track("prion")
    track["chains"][0]["residues"][10]["aa"] = "W"
    assert any("record's" in p for p in check_folding.problems_of(BY_SLUG["prion"], track))


@needs
def test_a_moved_frame_is_caught():
    track = copy.deepcopy(_track("insulin"))
    track["frame"]["centre_angstrom"][0] += 0.5
    found = check_folding.problems_of(BY_SLUG["insulin"], track)
    assert any("not the entry's in this frame" in p for p in found)


@needs
def test_a_fitted_frame_is_not_close_enough():
    # Glucagon has no bridges, and a CA fit to its one straight helix comes out
    # 3.9% short: its ends sit up to 1.5 A from the model's. The check holds
    # every helix CA to its ribbon, so the fitted track fails where the bake's
    # own frame passes.
    target = BY_SLUG["glucagon"]
    glb = (DATA / target.structure_asset).read_bytes()
    meshes = read_glb(glb)
    vertices = np.vstack([m.positions for m in meshes.values()])
    fit = fit_size(vertices, ca_atoms(target.structure))
    fitted = _track("glucagon", Frame(fit.centre, fit.length, "fitted", None))
    assert any("off its ribbon" in p for p in check_folding.problems_of(target, fitted))


# -- where it lives ---------------------------------------------------------------


def test_the_fetcher_and_the_uploader_agree_on_where_it_lives():
    prion = BY_SLUG["prion"]
    assert fetch_tracks.asset_path("folding", prion) == bake_folding.folding_asset(prion)
    assert upload_tracks.asset_of("folding", prion) == DATA / bake_folding.folding_asset(prion)
    assert "folding" in fetch_tracks.KINDS


_BAKED = [t for t in bake_folding.FOLDING_TARGETS if (DATA / bake_folding.folding_asset(t)).exists()]


@pytest.mark.skipif(not _BAKED or not all(_stored(t.slug) for t in _BAKED),
                    reason="run bake_folding.py and fetch_tracks.py first")
@pytest.mark.parametrize("target", _BAKED, ids=lambda t: t.slug)
def test_every_baked_track_checks_out(target):
    assert check_folding.problems_of(target) == []
