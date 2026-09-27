"""Bake `audio`: each protein as sound, one note a residue, and the map between them.

    .venv/Scripts/python -m pip install -r pipeline/audio/requirements.txt
    .venv/Scripts/python pipeline/fetch_tracks.py --kind record --kind constraint --kind clinvar
    .venv/Scripts/python pipeline/audio/bake_audio.py --all        # after bake_folding.py
    .venv/Scripts/python pipeline/audio/bake_audio.py --target insulin

The lab's Listen plays a protein residue by residue with a playhead over its
grid. This track is what it plays: per protein, an audio file of one note a
residue, and a timing map from each residue to the millisecond its note starts.

**The mapping.** Four properties, each on its own channel of the note:

- Pitch is **hydropathy**, Kyte and Doolittle's scale, quantised to the nearest
  step of a C major pentatonic scale over two octaves from middle C: more
  hydrophobic is higher, so a signal peptide or a transmembrane helix is a high
  run. A pentatonic scale has no step that clashes with another.
- Timbre is **secondary structure**, from the `folding` track: a helix is a held
  reed, a strand a plucked string, a coil a soft bell. A residue the fold does
  not place (outside the solved chains, or loose in them) is a hollow swell,
  and a protein with no folding track is that one timbre throughout.
- Loudness is **conservation**, from the ESM-2 `constraint` track: the most
  conserved residue at full level, the least twelve decibels under it. A
  protein with no constraint track plays every note at the middle of the range.
- A short tick opens the note of any residue **a ClinVar record** is reported
  at, from the `clinvar` snapshot, whatever the record's classification.

The mapping is arbitrary, chosen for teaching; the properties behind it are
measured. The map carries the mapping (the scale, the notes, the timbres and
what each sounds like) so the app says what it plays from the file itself.

**The tempo.** A residue lasts an eighth of a second, eight a second, up to 960
residues. A longer protein plays faster, so that no piece runs past two minutes:
CFTR at about twelve notes a second and dystrophin at about thirty-one. The
tempo is the size: a note a residue at one tempo would make dystrophin's
3,685 residues over seven minutes, and over a megabyte, where every other
protein is a few hundred kilobytes at most. A protein that could only fit by
notes under 25 ms, clicks rather than pitches, is refused: over 4,800 residues.

**The file** is AAC-LC in an MP4 container (`.m4a`), mono, 16 kHz, 20 kbit/s,
with the encoder's low-pass at 4.5 kHz. Every partial the mapping plays sits
under 4.4 kHz, and at this rate each channel survives the codec: decoded, 98
to 100% of each protein's notes keep their pitch and their timbre, and their
level to within a decibel (`check_audio.py` holds every file to it). The
container's edit list skips the encoder's priming, so position
zero is the first note on every player (see `m4a.py`), and the timing map is
exact. The map rides inside the file, in a `uuid` box after the audio, so the
track is one object and one digest, and the file still plays anywhere.

A bake is deterministic: the same inputs give the same bytes (ffmpeg's
bit-exact flags keep its version out of the file).

Which proteins get the track: all twenty (`AUDIO_TARGETS`), recorded here,
never in targets.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.audio import m4a  # noqa: E402
from pipeline.folding.bake_folding import FOLDING_TARGETS, folding_asset  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import BY_SLUG, TARGETS, Target  # noqa: E402

AUDIO_TARGETS = TARGETS

SCHEMA_VERSION = 1

# -- the file ------------------------------------------------------------------

SAMPLE_RATE = 16_000
BITRATE = 20_000
# The encoder's low-pass, and the band every partial is kept under.
CUTOFF_HZ = 4_500
BANDWIDTH_HZ = 4_400
# AAC-LC codes 1,024 samples a frame; the encoder is fed that many at a time.
AAC_FRAME = 1_024

# -- the tempo -----------------------------------------------------------------

# A residue an eighth of a second: 2,000 samples at 16 kHz.
NOTE_SAMPLES = 2_000
# No piece longer than two minutes.
LONGEST_SAMPLES = 120 * SAMPLE_RATE
# No note shorter than 25 ms: under that a note is a click, not a pitch. A
# protein that could only fit two minutes that way -- over 4,800 residues,
# like titin's 34,350 -- is refused rather than played as noise.
SHORTEST_NOTE_SAMPLES = 400

# -- pitch ---------------------------------------------------------------------

# Kyte and Doolittle's hydropathy index (J Mol Biol 1982;157:105-132, table 2).
HYDROPATHY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}
# C major pentatonic, two octaves up from middle C, as MIDI note numbers.
STEPS = (60, 62, 64, 67, 69, 72, 74, 76, 79, 81, 84)
_LOWEST, _HIGHEST = min(HYDROPATHY.values()), max(HYDROPATHY.values())

# -- timbre --------------------------------------------------------------------


@dataclass(frozen=True)
class Timbre:
    code: str                      # its letter in the map's `timbre` string
    structure: str                 # what the folding track calls it
    sound: str                     # what it sounds like, for the app to say
    partials: tuple[float, ...]    # each harmonic's amplitude, the fundamental's first
    attack: float                  # the part of a note spent rising
    decay: float | None            # e-folding time as a part of a note; None holds it


# DSSP's letters for a helix, a strand and a coil, and `-` where nothing is known.
TIMBRES = (
    Timbre("H", "helix", "a held reed", (1.0, 0.5, 0.33, 0.25), 0.06, None),
    Timbre("E", "strand", "a plucked string", (1.0, 0.5, 0.33, 0.25, 0.2, 0.16), 0.02, 0.30),
    Timbre("C", "coil", "a soft bell", (1.0, 0.15), 0.08, 0.60),
    Timbre("-", "none", "a hollow swell", (1.0, 0.0, 0.11, 0.0, 0.04), 0.25, None),
)
BY_CODE = {t.code: t for t in TIMBRES}
BY_STRUCTURE = {t.structure: t for t in TIMBRES}
NONE = BY_STRUCTURE["none"]

# Every note falls to silence over its last part, so no note clicks into the next.
RELEASE = 0.15

# -- loudness ------------------------------------------------------------------

# The least conserved residue plays this far under the most conserved; a
# protein with no constraint track plays every note at the middle.
QUIETEST_DB = -12.0
FLAT_DB = -6.0

# -- the accent ----------------------------------------------------------------

# A tick between the harmonics any note carries: none lies within 170 Hz of it.
TICK_HZ = 3_700.0
TICK_LEVEL = 0.25
TICK_DECAY_S = 0.004
TICK_S = 0.020

# The loudest waveform any timbre makes is a strand's, 1.874 times its RMS at
# its peak; at this level a note at full loudness stays under 0.64 of full
# scale and a tick adds at most 0.25, so nothing clips (`test_nothing_clips`).
LEVEL = 0.34


def audio_asset(target: Target) -> str:
    return f"assets/audio/{target.slug}.m4a"


# -- the map -------------------------------------------------------------------


def note_samples(residues: int) -> int:
    """How long each residue's note is, in samples: an eighth of a second, or
    shorter where that would run the piece past two minutes."""
    if residues < 1:
        raise ValueError("a protein of no residues")
    samples = min(NOTE_SAMPLES, LONGEST_SAMPLES // residues)
    if samples < SHORTEST_NOTE_SAMPLES:
        raise ValueError(f"{residues:,} residues in two minutes are notes of "
                         f"{samples * 1000 / SAMPLE_RATE:.1f} ms; under "
                         f"{SHORTEST_NOTE_SAMPLES * 1000 // SAMPLE_RATE} ms a note is a click")
    return samples


def step_of(letter: str) -> int:
    """The step of the scale nearest a residue's hydropathy, lowest first."""
    fraction = (HYDROPATHY[letter] - _LOWEST) / (_HIGHEST - _LOWEST)
    return math.floor(fraction * (len(STEPS) - 1) + 0.5)


def pitch_of(letter: str) -> int:
    """A residue's note, as a MIDI note number."""
    return STEPS[step_of(letter)]


def hertz(pitch: int) -> float:
    return 440.0 * 2.0 ** ((pitch - 69) / 12)


def onset_ms(residue: int, samples_per_note: int) -> int:
    """When the note of residue `residue` (0-based) starts, to the nearest ms."""
    return (2 * residue * samples_per_note * 1000 + SAMPLE_RATE) // (2 * SAMPLE_RATE)


def timbres_from(folding: dict | None, residues: int) -> str:
    """One code a residue: the structure the folding track places it in, else `-`.

    With no folding track at all, every residue is `-`, the one timbre. A
    residue two chains both place is refused rather than guessed between.
    """
    codes = [NONE.code] * residues
    if folding is None:
        return "".join(codes)
    for chain in folding["chains"]:
        for residue in chain["residues"]:
            if residue["state"] != "ordered":
                continue
            n = residue["n"]
            if not 1 <= n <= residues:
                raise ValueError(f"the folding track places residue {n} of {residues}")
            if codes[n - 1] != NONE.code:
                raise ValueError(f"the folding track places residue {n} twice")
            codes[n - 1] = BY_STRUCTURE[residue["ss"]].code
    return "".join(codes)


def loudness_from(constraint: dict | None, sequence: str) -> tuple[float, ...] | None:
    """Each residue's conservation, to two places; None with no constraint track."""
    if constraint is None:
        return None
    if constraint.get("sequence") != sequence:
        raise ValueError("the constraint track scored another sequence")
    positions = constraint["positions"]
    if [p["index"] for p in positions] != list(range(len(sequence))):
        raise ValueError("the constraint track does not score every residue once, in order")
    return tuple(round(float(p["conservation"]), 2) for p in positions)


def accents_from(clinvar: dict | None, sequence: str) -> tuple[int, ...]:
    """The residues at least one ClinVar record is reported at, in order."""
    if clinvar is None:
        return ()
    if clinvar.get("protein_sequence") != sequence:
        raise ValueError("the ClinVar snapshot places records on another sequence")
    marked = {v["residue"] for v in clinvar["variants"] if v.get("residue") is not None}
    outside = sorted(n for n in marked if not 1 <= n <= len(sequence))
    if outside:
        raise ValueError(f"ClinVar records at residues {outside[:3]}, outside the protein")
    return tuple(sorted(marked))


@dataclass(frozen=True)
class Score:
    """What each residue plays, before any sound: the map's own arrays."""

    sequence: str
    note_samples: int
    pitch: tuple[int, ...]
    timbre: str
    loudness: tuple[float, ...] | None
    accent: tuple[int, ...]

    @property
    def residues(self) -> int:
        return len(self.sequence)

    @property
    def samples(self) -> int:
        return self.residues * self.note_samples

    def onsets(self) -> list[int]:
        return [onset_ms(i, self.note_samples) for i in range(self.residues)]


def score(target: Target, record: dict, folding: dict | None,
          constraint: dict | None, clinvar: dict | None) -> Score:
    if record.get("gene") != target.gene:
        raise ValueError(f"{target.slug}: the record is {record.get('gene')}'s")
    sequence = record["protein"]["translation"]
    if len(sequence) != target.aa:
        raise ValueError(f"{target.slug}: the record translates {len(sequence)} "
                         f"residues, the table says {target.aa}")
    unknown = sorted(set(sequence) - set(HYDROPATHY))
    if unknown:
        raise ValueError(f"{target.slug}: no hydropathy for {unknown}")
    return Score(
        sequence=sequence,
        note_samples=note_samples(len(sequence)),
        pitch=tuple(pitch_of(letter) for letter in sequence),
        timbre=timbres_from(folding, len(sequence)),
        loudness=loudness_from(constraint, sequence),
        accent=accents_from(clinvar, sequence),
    )


# -- the sound -----------------------------------------------------------------


def gain_db(level: float | None) -> float:
    """A note's level in decibels under full: from its conservation, or flat."""
    return FLAT_DB if level is None else QUIETEST_DB * (1.0 - level)


def wave(timbre: Timbre, pitch: int, n: int) -> np.ndarray:
    """`n` samples of a timbre at a pitch, as loud (RMS) as a sine of amplitude 1.

    Only the harmonics under the band the encoder keeps are played, so what
    the file holds is what was synthesised.
    """
    f = hertz(pitch)
    t = np.arange(n) / SAMPLE_RATE
    kept = [(k, a) for k, a in enumerate(timbre.partials, 1) if a and k * f < BANDWIDTH_HZ]
    norm = math.sqrt(sum(a * a for _, a in kept) / 2)
    return sum(a * np.sin(2 * np.pi * k * f * t) for k, a in kept) / norm


def envelope(timbre: Timbre, n: int) -> np.ndarray:
    """A note's shape: its rise, its decay if it has one, and its fall to silence."""
    shape = np.ones(n)
    rise = max(1, round(timbre.attack * n))
    shape[:rise] = np.arange(rise) / rise
    if timbre.decay is not None:
        shape *= np.exp(-np.arange(n) / (timbre.decay * n))
    fall = max(2, round(RELEASE * n))
    shape[-fall:] *= np.linspace(1.0, 0.0, fall)
    return shape


def tick(n: int) -> np.ndarray:
    """The accent: a short, fast-decaying tone at the start of a note."""
    k = min(n, round(TICK_S * SAMPLE_RATE))
    t = np.arange(k) / SAMPLE_RATE
    out = np.zeros(n)
    out[:k] = TICK_LEVEL * np.sin(2 * np.pi * TICK_HZ * t) * np.exp(-t / TICK_DECAY_S)
    return out


def synthesise(played: Score) -> np.ndarray:
    """The whole piece, one note a residue, as float32 samples at SAMPLE_RATE."""
    n = played.note_samples
    shapes: dict[tuple[int, str], np.ndarray] = {}
    accent = tick(n)
    accented = set(played.accent)
    out = np.zeros(played.samples)
    for i, (pitch, code) in enumerate(zip(played.pitch, played.timbre)):
        key = (pitch, code)
        if key not in shapes:
            timbre = BY_CODE[code]
            shapes[key] = LEVEL * wave(timbre, pitch, n) * envelope(timbre, n)
        level = None if played.loudness is None else played.loudness[i]
        note = shapes[key] * 10.0 ** (gain_db(level) / 20.0)
        if i + 1 in accented:
            note = note + accent
        out[i * n:(i + 1) * n] = note
    return out.astype(np.float32)


def encode(pcm: np.ndarray) -> bytes:
    """AAC-LC in MP4, mono, at SAMPLE_RATE and BITRATE, through PyAV.

    `faststart` puts the index before the audio, and needs a file it can read
    back. `bitexact` keeps ffmpeg's version out of the file, so the same
    samples always make the same bytes. A movie timescale of the sample rate
    lets the edit list count samples, not milliseconds.
    """
    import av  # the bake's one requirement beyond numpy: requirements.txt

    with tempfile.TemporaryDirectory() as work:
        path = Path(work) / "audio.m4a"
        with av.open(str(path), "w", format="ipod", container_options={
            "movflags": "+faststart",
            "fflags": "+bitexact",
            "movie_timescale": str(SAMPLE_RATE),
        }) as container:
            stream = container.add_stream("aac", rate=SAMPLE_RATE, layout="mono", options={
                "flags": "+bitexact",
                "cutoff": str(CUTOFF_HZ),
            })
            stream.bit_rate = BITRATE
            for start in range(0, len(pcm), AAC_FRAME):
                frame = av.AudioFrame.from_ndarray(
                    np.ascontiguousarray(pcm[start:start + AAC_FRAME]).reshape(1, -1),
                    format="fltp", layout="mono")
                frame.sample_rate = SAMPLE_RATE
                frame.pts = start
                for packet in stream.encode(frame):
                    container.mux(packet)
            for packet in stream.encode(None):
                container.mux(packet)
        return path.read_bytes()


# -- the track -----------------------------------------------------------------


def mapping(folding: bool, constraint: bool, clinvar: bool) -> dict:
    """Which property drives which channel, and how: what the app says it plays."""
    return {
        "pitch": {
            "property": "hydropathy",
            "scale": "Kyte & Doolittle 1982",
            "higher": "more hydrophobic",
            "notes": "C major pentatonic, C4 to C6",
            "steps": list(STEPS),
            "hydropathy": dict(HYDROPATHY),
        },
        "timbre": {
            "property": "secondary structure",
            "source": "folding" if folding else None,
            "timbres": [{"code": t.code, "structure": t.structure, "sound": t.sound}
                        for t in TIMBRES],
        },
        "loudness": {
            "property": "conservation",
            "source": "constraint" if constraint else None,
            "louder": "more conserved",
            "range_db": [QUIETEST_DB, 0.0],
            "flat_db": FLAT_DB,
        },
        "accent": {
            "property": "a ClinVar record reported at the residue",
            "source": "clinvar" if clinvar else None,
            "sound": "a short tick as the note starts",
        },
    }


def _source(relative: str, blob: bytes, **fields) -> dict:
    return {"path": relative, "sha256": hashlib.sha256(blob).hexdigest(), **fields}


def clinvar_asset(target: Target) -> str:
    return f"assets/clinvar/{target.slug}_clinvar.json"


def inputs(target: Target) -> tuple[dict, dict | None, dict | None, dict | None, dict]:
    """The record and whichever of the folding, constraint and ClinVar tracks
    this protein has, with what they are: (record, folding, constraint,
    clinvar, sources). A track the protein should have and the disk does not
    hold is an error, never a quieter or plainer piece."""
    wanted = {
        "record": target.mock_asset,
        "folding": folding_asset(target) if target in FOLDING_TARGETS else None,
        "constraint": target.constraint_asset if target.scored else None,
        "clinvar": clinvar_asset(target) if target.clinvar_available else None,
    }
    loaded: dict[str, dict | None] = {}
    sources: dict[str, dict | None] = {}
    for name, relative in wanted.items():
        if relative is None:
            loaded[name] = sources[name] = None
            continue
        path = DATA / relative
        if not path.exists():
            how = ("run folding/bake_folding.py" if name == "folding"
                   else f"run fetch_tracks.py --kind {name}")
            raise FileNotFoundError(f"{path}: {how} first")
        blob = path.read_bytes()
        data = json.loads(blob)
        loaded[name] = data
        extra = {
            "record": lambda d: {"accession": target.source.accession},
            "folding": lambda d: {"pdb": d["pdb"]},
            "constraint": lambda d: {"model": d["model"]},
            "clinvar": lambda d: {"retrieved_at": d["retrieved_at"]},
        }[name](data)
        sources[name] = _source(relative, blob, **extra)
    return loaded["record"], loaded["folding"], loaded["constraint"], loaded["clinvar"], sources


def the_map(target: Target, played: Score, sources: dict, facts: m4a.Facts) -> dict:
    """What the file plays, residue by residue, and what it was made from."""
    n = played.note_samples
    return {
        "slug": target.slug,
        "gene": target.gene,
        "uniprot": target.uniprot,
        "schema_version": SCHEMA_VERSION,
        "residues": played.residues,
        "sequence": played.sequence,
        "audio": {
            "format": "m4a",
            "codec": "AAC-LC",
            "sample_rate": facts.sample_rate,
            "channels": facts.channels,
            "bitrate": BITRATE,
            "cutoff_hz": CUTOFF_HZ,
            "priming": facts.priming,
            "samples": facts.samples,
            "duration_ms": onset_ms(played.residues, n),
        },
        "tempo": {
            "note_samples": n,
            "note_ms": round(n * 1000 / SAMPLE_RATE, 4),
            "longest_s": LONGEST_SAMPLES // SAMPLE_RATE,
            "rule": ("a residue lasts 125 ms, eight a second, unless that would run "
                     "the piece past two minutes; then all of it plays faster"),
        },
        "onset_ms": played.onsets(),
        "pitch": list(played.pitch),
        "timbre": played.timbre,
        "loudness": None if played.loudness is None else list(played.loudness),
        "accent": list(played.accent),
        "mapping": mapping(sources["folding"] is not None, sources["constraint"] is not None,
                           sources["clinvar"] is not None),
        "sources": sources,
        "built_by": "pipeline/audio/bake_audio.py",
    }


def track_of(target: Target) -> tuple[bytes, dict]:
    """The file for one protein, the map inside it and all: (bytes, map)."""
    record, folding, constraint, clinvar, sources = inputs(target)
    played = score(target, record, folding, constraint, clinvar)
    audio = encode(synthesise(played))
    facts = m4a.facts(audio)
    if facts.samples != played.samples:
        raise ValueError(f"{target.slug}: the edit list plays {facts.samples} samples, "
                         f"the piece is {played.samples}")
    carried = the_map(target, played, sources, facts)
    return m4a.with_map(audio, carried), carried


def bake(target: Target) -> dict:
    blob, carried = track_of(target)
    destination = DATA / audio_asset(target)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(blob)
    return carried


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if not args.all and not args.target:
        parser.error("name --target or --all")
    unknown = sorted(set(args.target or []) - set(BY_SLUG))
    if unknown:
        parser.error(f"no such target: {unknown}")
    chosen = AUDIO_TARGETS if args.all else tuple(BY_SLUG[s] for s in args.target)
    print(f"{'slug':<16}{'residues':>9}{'note ms':>9}{'per s':>7}{'seconds':>9}"
          f"{'audio B':>10}{'map B':>9}{'file B':>10}")
    total = 0
    largest = (0, "")
    for target in chosen:
        carried = bake(target)
        blob = (DATA / audio_asset(target)).read_bytes()
        audio = len(m4a.without_map(blob))
        size = len(blob)
        total += size
        largest = max(largest, (size, target.slug))
        tempo = carried["tempo"]
        print(f"{target.slug:<16}{carried['residues']:>9,}{tempo['note_ms']:>9.1f}"
              f"{SAMPLE_RATE / tempo['note_samples']:>7.1f}"
              f"{carried['audio']['duration_ms'] / 1000:>9.1f}"
              f"{audio:>10,}{size - audio:>9,}{size:>10,}")
    print(f"{len(chosen)} tracks, {total:,} B; the largest is {largest[1]}'s, {largest[0]:,} B",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
