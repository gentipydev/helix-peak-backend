"""The audio bake, its file and its check, offline.

The mapping is tested on its own tables. Tests that score a protein need its
stored record, folding, constraint and ClinVar tracks, and skip until
`fetch_tracks.py` and `bake_folding.py` have put them in `pipeline/data/`.
Encoding and listening need PyAV (`requirements.txt`) and skip without it.
"""

from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

import numpy as np
import pytest

BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline import fetch_tracks, upload_tracks  # noqa: E402
from pipeline.audio import bake_audio as bake  # noqa: E402
from pipeline.audio import check_audio, m4a  # noqa: E402
from pipeline.audio.bake_audio import Score  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import BY_SLUG  # noqa: E402

try:
    import av  # noqa: F401
    HAS_AV = True
except ImportError:
    HAS_AV = False

needs_av = pytest.mark.skipif(not HAS_AV, reason="pip install -r pipeline/audio/requirements.txt")


def _stored(slug: str) -> bool:
    target = BY_SLUG[slug]
    return all((DATA / p).exists() for p in (
        target.mock_asset, bake.folding_asset(target), target.constraint_asset,
        bake.clinvar_asset(target)))


needs_inputs = pytest.mark.skipif(
    not all(_stored(s) for s in ("insulin", "glucagon")),
    reason="run fetch_tracks.py --kind record --kind constraint --kind clinvar "
           "and folding/bake_folding.py first")


def _score(slug: str) -> Score:
    target = BY_SLUG[slug]
    record, folding, constraint, clinvar, _ = bake.inputs(target)
    return bake.score(target, record, folding, constraint, clinvar)


# -- pitch: hydropathy on a pentatonic scale -----------------------------------


def test_more_hydrophobic_is_never_lower():
    letters = sorted(bake.HYDROPATHY, key=bake.HYDROPATHY.get)
    pitches = [bake.pitch_of(letter) for letter in letters]
    assert pitches == sorted(pitches)


def test_every_note_is_a_step_of_c_major_pentatonic():
    assert {p % 12 for p in bake.STEPS} == {0, 2, 4, 7, 9}
    assert (bake.STEPS[0], bake.STEPS[-1]) == (60, 84)
    assert all(bake.pitch_of(letter) in bake.STEPS for letter in bake.HYDROPATHY)


def test_each_residue_has_its_note():
    # Arginine, the most hydrophilic, is middle C; isoleucine and valine the
    # top. No residue's hydropathy falls nearest E4 or D5.
    assert {letter: bake.pitch_of(letter) for letter in "RKDENQHPSTWYGAMCFLVI"} == {
        "R": 60, "K": 62, "D": 62, "E": 62, "N": 62, "Q": 62, "H": 62, "P": 67,
        "S": 69, "T": 69, "W": 69, "Y": 69, "G": 72, "A": 76, "M": 76, "C": 79,
        "F": 79, "L": 81, "V": 84, "I": 84,
    }


def test_the_scale_is_equal_tempered_from_a_440():
    assert bake.hertz(69) == 440.0
    assert bake.hertz(60) == pytest.approx(261.626, abs=1e-3)


# -- the tempo -----------------------------------------------------------------


def test_a_residue_lasts_an_eighth_of_a_second_up_to_960():
    assert bake.note_samples(110) == bake.NOTE_SAMPLES == 2000
    assert bake.note_samples(960) == 2000
    assert bake.note_samples(961) < 2000


def test_no_piece_runs_past_two_minutes():
    for residues in (1, 960, 961, 1480, 3685, 4800):
        assert bake.note_samples(residues) * residues <= 120 * bake.SAMPLE_RATE


def test_a_protein_whose_notes_would_be_clicks_is_refused():
    assert bake.note_samples(4800) == bake.SHORTEST_NOTE_SAMPLES
    with pytest.raises(ValueError, match="click"):
        bake.note_samples(34350)            # titin


def test_onsets_are_one_note_apart_from_zero_to_the_nearest_millisecond():
    assert [bake.onset_ms(i, 2000) for i in range(4)] == [0, 125, 250, 375]
    # Dystrophin's 521-sample notes: 8 x 521 / 16 is 260.5, which rounds up.
    assert bake.onset_ms(8, 521) == 261
    assert bake.onset_ms(3685, 521) == 119993


# -- timbre, loudness and the accent from their tracks -------------------------


def _folding(*residues: dict) -> dict:
    return {"chains": [{"residues": list(residues)}]}


def test_timbre_is_the_folding_tracks_structure_and_nothing_else():
    folding = _folding(
        {"n": 2, "state": "ordered", "ss": "helix"},
        {"n": 3, "state": "ordered", "ss": "strand"},
        {"n": 4, "state": "ordered", "ss": "coil"},
        {"n": 5, "state": "disordered"},
        {"n": 6, "state": "absent"},
    )
    assert bake.timbres_from(folding, 7) == "-HEC---"


def test_a_protein_with_no_folding_track_is_one_timbre():
    assert bake.timbres_from(None, 5) == "-----"


def test_a_residue_two_chains_place_is_refused_not_guessed():
    twice = {"chains": [{"residues": [{"n": 1, "state": "ordered", "ss": "helix"}]},
                        {"residues": [{"n": 1, "state": "ordered", "ss": "strand"}]}]}
    with pytest.raises(ValueError, match="twice"):
        bake.timbres_from(twice, 3)


def test_loudness_is_conservation_to_two_places():
    constraint = {"sequence": "MA", "positions": [
        {"index": 0, "conservation": 0.93464}, {"index": 1, "conservation": 0.176825}]}
    assert bake.loudness_from(constraint, "MA") == (0.93, 0.18)
    with pytest.raises(ValueError):
        bake.loudness_from(constraint, "MK")


def test_an_unscored_protein_plays_every_note_at_one_level():
    assert bake.loudness_from(None, "MA") is None
    assert bake.gain_db(None) == bake.FLAT_DB == -6.0
    flat = bake.synthesise(Score("AAAA", 2000, (76,) * 4, "----", None, ())).reshape(4, 2000)
    assert all(np.array_equal(flat[0], note) for note in flat[1:])
    # Flat is the middle of the scored range: a conservation of one half.
    halfway = bake.synthesise(Score("AAAA", 2000, (76,) * 4, "----", (0.5,) * 4, ()))
    assert np.array_equal(halfway.reshape(4, 2000), flat)


def test_the_most_conserved_is_twelve_decibels_over_the_least():
    loud = bake.synthesise(Score("A", 2000, (76,), "H", (1.0,), ()))
    quiet = bake.synthesise(Score("A", 2000, (76,), "H", (0.0,), ()))
    ratio = 20 * np.log10(np.sqrt(np.mean(loud.astype(float) ** 2) / np.mean(quiet.astype(float) ** 2)))
    assert ratio == pytest.approx(12.0, abs=1e-3)


def test_ticks_are_the_residues_clinvar_has_records_at():
    clinvar = {"protein_sequence": "MAL", "variants": [
        {"residue": 2}, {"residue": None}, {"residue": 2}, {"residue": 3}]}
    assert bake.accents_from(clinvar, "MAL") == (2, 3)
    assert bake.accents_from(None, "MAL") == ()
    with pytest.raises(ValueError):
        bake.accents_from({"protein_sequence": "MAL", "variants": [{"residue": 4}]}, "MAL")


# -- the sound -----------------------------------------------------------------


def test_nothing_clips():
    for n in (bake.NOTE_SAMPLES, bake.SHORTEST_NOTE_SAMPLES):
        for timbre in bake.TIMBRES:
            for pitch in bake.STEPS:
                note = bake.synthesise(Score("A", n, (pitch,), timbre.code, (1.0,), (1,)))
                assert np.abs(note).max() < 0.89


def test_every_note_falls_silent_before_the_next():
    for timbre in bake.TIMBRES:
        note = bake.synthesise(Score("A", 2000, (72,), timbre.code, (1.0,), (1,)))
        assert note[0] == 0.0
        assert note[-1] == 0.0


def test_every_partial_is_under_the_band_the_encoder_keeps():
    played = Score("I" * 8, 2000, (84,) * 8, "E" * 8, (1.0,) * 8, ())
    spectrum = np.abs(np.fft.rfft(bake.synthesise(played).astype(float)))
    freqs = np.fft.rfftfreq(8 * 2000, 1 / bake.SAMPLE_RATE)
    above = spectrum[freqs > bake.BANDWIDTH_HZ + 100].max()
    assert above < 1e-3 * spectrum.max()


def test_the_tick_is_clear_of_every_harmonic_a_note_plays():
    for timbre in bake.TIMBRES:
        for pitch in bake.STEPS:
            f = bake.hertz(pitch)
            for k, a in enumerate(timbre.partials, 1):
                if a and k * f < bake.BANDWIDTH_HZ:
                    assert abs(k * f - bake.TICK_HZ) > 170


def test_every_timbre_is_as_loud_as_the_others():
    levels = []
    for timbre in bake.TIMBRES:
        wave = bake.wave(timbre, 72, 16000)
        levels.append(np.sqrt(np.mean(wave ** 2)))
    assert max(levels) - min(levels) < 1e-3


# -- the file ------------------------------------------------------------------


@needs_av
def test_the_file_is_one_mono_aac_stream_that_starts_on_the_first_note():
    played = Score("MALWR" * 20, 2000, tuple(bake.pitch_of(a) for a in "MALWR" * 20),
                   "HHEEC" * 20, tuple([0.3, 0.9, 0.5, 1.0, 0.0] * 20), (1, 7, 50))
    audio = bake.encode(bake.synthesise(played))
    facts = m4a.facts(audio)
    assert (facts.brand, facts.codec, facts.object_type) == ("M4A ", "mp4a", m4a.AAC_LC)
    assert (facts.sample_rate, facts.channels) == (16000, 1)
    assert facts.priming == 1024
    assert facts.samples == played.samples == 200000
    heard = check_audio.listen(m4a.with_map(audio, {}), {
        "sequence": played.sequence, "tempo": {"note_samples": 2000}, "pitch": list(played.pitch),
        "timbre": played.timbre, "loudness": list(played.loudness), "accent": list(played.accent)})
    assert heard.lag == 0
    assert heard.pitches == 1.0
    assert heard.timbres >= 0.95


@needs_av
def test_the_same_samples_make_the_same_bytes():
    played = Score("MALW" * 10, 2000, tuple(bake.pitch_of(a) for a in "MALW" * 10),
                   "-" * 40, None, ())
    pcm = bake.synthesise(played)
    assert bake.encode(pcm) == bake.encode(pcm)


@needs_av
def test_the_map_rides_last_and_reads_back():
    audio = bake.encode(np.zeros(4096, dtype=np.float32))
    carried = m4a.with_map(audio, {"slug": "insulin", "onset_ms": [0, 125]})
    assert m4a.read_map(carried) == {"slug": "insulin", "onset_ms": [0, 125]}
    assert m4a.carries_map_last(carried)
    assert m4a.without_map(carried) == audio
    assert m4a.facts(carried) == m4a.facts(audio)
    with pytest.raises(ValueError, match="already"):
        m4a.with_map(carried, {})
    with pytest.raises(ValueError, match="0 maps"):
        m4a.read_map(audio)


def test_a_box_that_overruns_its_file_is_refused():
    with pytest.raises(ValueError, match="claims"):
        m4a.boxes(b"\x00\x00\x00\x40ftypM4A \x00\x00\x00\x00")


# -- the proteins --------------------------------------------------------------


@needs_inputs
def test_insulins_signal_peptide_core_is_its_highest_run():
    played = _score("insulin")
    runs = [np.mean(played.pitch[i:i + 10]) for i in range(len(played.pitch) - 9)]
    # LLPLLALLAL, the hydrophobic core of MALWMRLLPLLALLALWGPDPAAA.
    assert int(np.argmax(runs)) + 1 == 7
    assert played.sequence[6:16] == "LLPLLALLAL"


@needs_inputs
def test_insulin_is_timbred_where_its_fold_is_solved():
    played = _score("insulin")
    # The fold is the B chain (25-54) and the A chain (90-110); the signal
    # peptide and the C peptide between them have no structure.
    assert set(played.timbre[:24]) == {"-"}
    assert set(played.timbre[56:88]) == {"-"}
    assert "H" in played.timbre[24:54] and "H" in played.timbre[89:110]


@needs_inputs
def test_glucagon_ticks_at_its_three_clinvar_residues():
    assert len(_score("glucagon").accent) == 3


@needs_inputs
def test_a_missing_input_is_an_error_not_a_plainer_piece(tmp_path, monkeypatch):
    monkeypatch.setattr(bake, "DATA", tmp_path)
    with pytest.raises(FileNotFoundError, match="fetch_tracks.py --kind record"):
        bake.inputs(BY_SLUG["insulin"])


@needs_inputs
def test_an_unfolded_unscored_protein_is_one_flat_timbre():
    target = dataclasses.replace(BY_SLUG["insulin"], structure=None, scored=False,
                                 clinvar_available=False)
    record, _, _, _, sources = bake.inputs(target)
    assert sources["folding"] is sources["constraint"] is sources["clinvar"] is None
    played = bake.score(target, record, None, None, None)
    assert set(played.timbre) == {"-"}
    assert played.loudness is None and played.accent == ()


# -- the check -----------------------------------------------------------------


@pytest.fixture(scope="module")
def insulin_track():
    if not (HAS_AV and _stored("insulin")):
        pytest.skip("needs PyAV and insulin's stored tracks")
    return bake.track_of(BY_SLUG["insulin"])


def _remapped(blob: bytes, change) -> bytes:
    carried = json.loads(json.dumps(m4a.read_map(blob)))
    change(carried)
    return m4a.with_map(m4a.without_map(blob), carried)


def test_a_true_track_passes(insulin_track):
    blob, carried = insulin_track
    assert check_audio.problems_of(BY_SLUG["insulin"], blob) == []
    assert carried["residues"] == len(carried["onset_ms"]) == 110


@pytest.mark.parametrize("change, finding", [
    (lambda c: c["onset_ms"].pop(), "the timing map has 109 onsets, the protein 110"),
    (lambda c: c["onset_ms"].__setitem__(5, 626), "not one note apart"),
    (lambda c: c["pitch"].__setitem__(0, 84), "hydropathy"),
    (lambda c: c.__setitem__("timbre", "H" + c["timbre"][1:]), "folding track"),
    (lambda c: c["loudness"].__setitem__(0, 0.5), "conservation"),
    (lambda c: c["accent"].pop(), "ClinVar"),
    (lambda c: c.__setitem__("slug", "glucagon"), "slug"),
    (lambda c: c["sources"]["constraint"].__setitem__("sha256", "0" * 64),
     "another constraint track"),
    (lambda c: c["mapping"]["pitch"].__setitem__("higher", "less hydrophobic"), "mapping"),
    (lambda c: c["audio"].__setitem__("samples", 1), "the map says samples"),
])
def test_the_check_catches_a_damaged_map(insulin_track, change, finding):
    blob, _ = insulin_track
    found = check_audio.problems_of(BY_SLUG["insulin"], _remapped(blob, change), hear=False)
    assert any(finding in problem for problem in found), found


def test_the_check_hears_audio_that_does_not_play_its_map(insulin_track):
    blob, carried = insulin_track
    played = bake.Score(carried["sequence"], 2000, tuple(reversed(carried["pitch"])),
                        carried["timbre"], tuple(carried["loudness"]), tuple(carried["accent"]))
    wrong = m4a.with_map(bake.encode(bake.synthesise(played)), carried)
    found = check_audio.problems_of(BY_SLUG["insulin"], wrong)
    assert any("keep their pitch" in problem for problem in found), found


def test_the_check_hears_audio_that_starts_late(insulin_track):
    blob, carried = insulin_track
    played = bake.Score(carried["sequence"], 2000, tuple(carried["pitch"]), carried["timbre"],
                        tuple(carried["loudness"]), tuple(carried["accent"]))
    pcm = bake.synthesise(played)
    late = np.concatenate([np.zeros(1024, dtype=np.float32), pcm[:-1024]])
    found = check_audio.problems_of(BY_SLUG["insulin"],
                                    m4a.with_map(bake.encode(late), carried))
    assert any("starts 1024 samples off" in problem for problem in found), found


def test_a_file_with_no_map_is_not_a_track(insulin_track):
    blob, _ = insulin_track
    found = check_audio.problems_of(BY_SLUG["insulin"], m4a.without_map(blob))
    assert found and "not a track" in found[0]


# -- storage -------------------------------------------------------------------


def test_the_fetcher_and_the_uploader_agree_on_where_it_lives():
    insulin = BY_SLUG["insulin"]
    assert "audio" in fetch_tracks.KINDS
    assert fetch_tracks.asset_path("audio", insulin) == bake.audio_asset(insulin)
    assert upload_tracks.asset_of("audio", insulin) == DATA / bake.audio_asset(insulin)


def test_the_uploader_names_the_mapping_and_refuses_damage(insulin_track):
    blob, carried = insulin_track
    provenance = upload_tracks.validate("audio", BY_SLUG["insulin"], blob)
    assert set(provenance) == {"audio", "tempo", "mapping", "sources", "schema_version",
                               "built_by"}
    assert provenance["mapping"]["pitch"]["property"] == "hydropathy"
    with pytest.raises(ValueError):
        upload_tracks.validate("audio", BY_SLUG["insulin"],
                               _remapped(blob, lambda c: c["onset_ms"].pop()))


@pytest.mark.parametrize("target", [t.slug for t in bake.AUDIO_TARGETS])
@pytest.mark.skipif(not HAS_AV or not (DATA / "assets/audio").is_dir(),
                    reason="run bake_audio.py --all first")
def test_every_baked_track_checks_out(target):
    if not (DATA / bake.audio_asset(BY_SLUG[target])).exists():
        pytest.skip(f"{target} not baked")
    assert check_audio.problems_of(BY_SLUG[target]) == []
