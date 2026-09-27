"""Check every baked `audio` track against its protein, its inputs and its own sound.

    .venv/Scripts/python pipeline/audio/check_audio.py                   # exits 1 on any failure
    .venv/Scripts/python pipeline/audio/check_audio.py --target insulin

For each target, the file under `pipeline/data/assets/audio/`:

- is an .m4a of one mono AAC-LC stream at 16 kHz whose edit list skips the
  encoder's priming and plays exactly the piece: one note a residue, each the
  tempo's length;
- carries one map, as its last box, naming its protein;
- has a timing map as long as the protein is: one onset per residue, the first
  at zero, each one note after the last, and the last note ending where the
  file does;
- plays what its inputs say: each pitch the hydropathy of the record's residue,
  each timbre the folding track's structure there, each loudness the constraint
  track's conservation, each tick a residue the ClinVar snapshot has a record
  at. The inputs are read from disk and held to the digests the map names;
- names its mapping and every source it was made from;
- sounds like its map: decoded, it starts on the first note, and note by note
  its pitch, level, timbre and ticks are those of the same map synthesised
  afresh (`listen`).

Exits 1 with the problems listed. Listening decodes, so it needs PyAV
(`requirements.txt`); without it the check says so and fails.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.audio import bake_audio as bake  # noqa: E402
from pipeline.audio import m4a  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import BY_SLUG, Target  # noqa: E402

# What the decoded file must keep of its map, note by note (measured on all
# twenty: see the README). A note's pitch is the frequency its spectrum peaks
# at, its level its RMS, its timbre the nearest of the four as the fresh
# synthesis sounds them, and a tick the energy around TICK_HZ as it starts.
PITCHES_KEPT = 0.97          # of notes whose pitch is the fresh synthesis's, within 3%
LEVEL_KEPT_DB = 1.5          # the mean difference in a note's level
TIMBRES_KEPT = 0.95          # of notes heard as their own timbre
TICKS_HEARD_DB = 10.0        # how much louder a note's opening is at TICK_HZ with a tick


def record_of(target: Target) -> tuple[dict, bytes] | None:
    path = DATA / target.mock_asset
    if not path.exists():
        return None
    blob = path.read_bytes()
    return json.loads(blob), blob


# -- listening -----------------------------------------------------------------


def decode(blob: bytes) -> np.ndarray:
    """The file as a player hears it: its edit list applied, float samples."""
    import av  # requirements.txt

    chunks = []
    with av.open(io.BytesIO(m4a.without_map(blob)), format="mp4") as container:
        for frame in container.decode(audio=0):
            chunks.append(frame.to_ndarray().astype(np.float64).reshape(-1))
    return np.concatenate(chunks) if chunks else np.zeros(0)


def _features(note: np.ndarray) -> tuple[float, float, np.ndarray, float]:
    """(pitch Hz, level dB, shape, tick dB) of one note's samples."""
    n = len(note)
    spectrum = np.abs(np.fft.rfft(note * np.hanning(n), n=8192))
    freqs = np.fft.rfftfreq(8192, 1 / bake.SAMPLE_RATE)
    band = (freqs > 200) & (freqs < 1200)
    pitch = float(freqs[band][np.argmax(spectrum[band])])
    level = 20 * np.log10(np.sqrt(np.mean(note ** 2)) + 1e-9)
    # The timbre: how the note's energy falls over its four quarters, and how
    # strong its second to fourth harmonics are against its fundamental.
    total = np.sum(note ** 2) + 1e-12
    q = n // 4
    quarters = [10 * np.log10((np.sum(note[j * q:(j + 1) * q] ** 2) + 1e-12) / total)
                for j in range(4)]

    def at(f: float) -> float:
        k = int(np.argmin(np.abs(freqs - f)))
        return float(spectrum[max(0, k - 3):k + 4].max()) + 1e-9

    harmonics = [20 * np.log10(at(pitch * k) / at(pitch)) if pitch * k < bake.SAMPLE_RATE / 2
                 else -120.0 for k in (2, 3, 4)]
    head = note[:min(n, round(bake.TICK_S * bake.SAMPLE_RATE))]
    opening = np.abs(np.fft.rfft(head * np.hanning(len(head)), n=4096))
    near = np.fft.rfftfreq(4096, 1 / bake.SAMPLE_RATE)
    window = (near > bake.TICK_HZ - 250) & (near < bake.TICK_HZ + 250)
    tick = 20 * np.log10(float(opening[window].max()) + 1e-9)
    return pitch, level, np.array(quarters + harmonics), tick


@dataclass(frozen=True)
class Heard:
    lag: int                # samples the decoded audio sits after the fresh synthesis
    pitches: float          # of notes at the fresh synthesis's pitch
    level_db: float         # mean difference in a note's level
    timbres: float          # of notes heard as their own timbre
    ticks_db: float | None  # ticked notes' openings over the others', or None without both


def listen(blob: bytes, carried: dict) -> Heard:
    """The decoded file against its own map synthesised afresh, note by note."""
    played = bake.Score(
        sequence=carried["sequence"],
        note_samples=carried["tempo"]["note_samples"],
        pitch=tuple(carried["pitch"]),
        timbre=carried["timbre"],
        loudness=None if carried["loudness"] is None else tuple(carried["loudness"]),
        accent=tuple(carried["accent"]),
    )
    fresh = bake.synthesise(played).astype(np.float64)
    heard = decode(blob)

    # Where the decoded audio lines up with the fresh synthesis: at zero if the
    # edit list skipped the priming, at 1,024 if it did not.
    reference = fresh[:bake.SAMPLE_RATE // 2]
    best, lag = -np.inf, 0
    for shift in range(-bake.AAC_FRAME - 64, bake.AAC_FRAME + 65):
        a = heard[max(0, shift):max(0, shift) + len(reference) - 1200]
        b = reference[max(0, -shift):max(0, -shift) + len(reference) - 1200]
        m = min(len(a), len(b))
        score = float(np.dot(a[:m], b[:m]))
        if score > best:
            best, lag = score, shift
    heard = heard[lag:] if lag >= 0 else np.concatenate([np.zeros(-lag), heard])
    heard = np.concatenate([heard, np.zeros(max(0, len(fresh) - len(heard)))])

    n = played.note_samples
    ours = [_features(fresh[i * n:(i + 1) * n]) for i in range(played.residues)]
    theirs = [_features(heard[i * n:(i + 1) * n]) for i in range(played.residues)]
    pitches = float(np.mean([abs(t[0] - o[0]) <= 0.03 * o[0] for o, t in zip(ours, theirs)]))
    level = float(np.mean([abs(t[1] - o[1]) for o, t in zip(ours, theirs)]))

    centres = {code: np.mean([o[2] for o, c in zip(ours, played.timbre) if c == code], axis=0)
               for code in set(played.timbre)}
    own = [min(centres, key=lambda code: float(np.sum((t[2] - centres[code]) ** 2))) == c
           for t, c in zip(theirs, played.timbre)]

    accented = set(played.accent)
    on = [t[3] for i, t in enumerate(theirs, 1) if i in accented]
    off = [t[3] for i, t in enumerate(theirs, 1) if i not in accented]
    ticks = float(np.median(on) - np.median(off)) if on and off else None
    return Heard(lag=lag, pitches=pitches, level_db=level, timbres=float(np.mean(own)),
                 ticks_db=ticks)


def heard_problems(where: str, blob: bytes, carried: dict) -> list[str]:
    try:
        heard = listen(blob, carried)
    except ImportError:
        return [f"{where}: cannot listen without PyAV; pip install -r pipeline/audio/requirements.txt"]
    out = []
    if heard.lag != 0:
        out.append(f"{where}: the audio starts {heard.lag} samples off its first note "
                   f"(an edit list not skipping the priming is 1,024)")
    if heard.pitches < PITCHES_KEPT:
        out.append(f"{where}: {heard.pitches:.1%} of notes keep their pitch, under {PITCHES_KEPT:.0%}")
    if heard.level_db > LEVEL_KEPT_DB:
        out.append(f"{where}: notes' levels are {heard.level_db:.2f} dB off, over {LEVEL_KEPT_DB}")
    if heard.timbres < TIMBRES_KEPT:
        out.append(f"{where}: {heard.timbres:.1%} of notes sound as their timbre, "
                   f"under {TIMBRES_KEPT:.0%}")
    if heard.ticks_db is not None and heard.ticks_db < TICKS_HEARD_DB:
        out.append(f"{where}: a tick adds {heard.ticks_db:.1f} dB to a note's opening, "
                   f"under {TICKS_HEARD_DB}")
    return out


# -- the file and its map ------------------------------------------------------


def _held_to(where: str, name: str, source: dict | None, relative: str | None) -> tuple[list[str], dict | None]:
    """The input the map names, read from disk and held to the digest it names."""
    if relative is None:
        if source is not None:
            return [f"{where}: names a {name} source this protein has no track for"], None
        return [], None
    if source is None:
        return [f"{where}: names no {name} source, and the protein has a {name} track"], None
    if source.get("path") != relative:
        return [f"{where}: its {name} source is {source.get('path')!r}, not {relative!r}"], None
    path = DATA / relative
    if not path.exists():
        return [f"{where}: no {name} track on disk to hold it to ({path})"], None
    blob = path.read_bytes()
    if hashlib.sha256(blob).hexdigest() != source.get("sha256"):
        return [f"{where}: was made from another {name} track than the one on disk"], None
    return [], json.loads(blob)


def problems_of(target: Target, blob: bytes | None = None, *, hear: bool = True) -> list[str]:
    """Everything wrong with one protein's track: the file's, or `blob`."""
    where = target.slug
    if blob is None:
        path = DATA / bake.audio_asset(target)
        if not path.exists():
            return [f"{where}: not baked ({path})"]
        blob = path.read_bytes()

    try:
        facts = m4a.facts(blob)
        carried = m4a.read_map(blob)
    except (ValueError, KeyError, IndexError, json.JSONDecodeError) as problem:
        return [f"{where}: not a track: {problem}"]

    out: list[str] = []
    if not m4a.carries_map_last(blob):
        out.append(f"{where}: the map is not the file's last box")
    for key, want in (("slug", target.slug), ("gene", target.gene), ("uniprot", target.uniprot),
                      ("schema_version", bake.SCHEMA_VERSION), ("residues", target.aa),
                      ("built_by", "pipeline/audio/bake_audio.py")):
        if carried.get(key) != want:
            out.append(f"{where}: {key} is {carried.get(key)!r}, expected {want!r}")

    # The file: what a decoder is told.
    for label, got, want in (("brand", facts.brand, "M4A "), ("codec", facts.codec, "mp4a"),
                             ("audio object type", facts.object_type, m4a.AAC_LC),
                             ("sample rate", facts.sample_rate, bake.SAMPLE_RATE),
                             ("channels", facts.channels, 1)):
        if got != want:
            out.append(f"{where}: {label} {got!r}, expected {want!r}")
    if facts.priming <= 0:
        out.append(f"{where}: the edit list skips no priming, so the first note is late")
    audio = carried.get("audio") or {}
    for key, got in (("sample_rate", facts.sample_rate), ("channels", facts.channels),
                     ("priming", facts.priming), ("samples", facts.samples)):
        if audio.get(key) != got:
            out.append(f"{where}: the map says {key} {audio.get(key)!r}, the file {got!r}")

    # The timing map: as long as the protein, and the file.
    onsets = carried.get("onset_ms")
    tempo = carried.get("tempo") or {}
    n = tempo.get("note_samples")
    if not isinstance(onsets, list) or len(onsets) != target.aa:
        count = len(onsets) if isinstance(onsets, list) else None
        return out + [f"{where}: the timing map has {count} onsets, the protein {target.aa} residues"]
    if n != bake.note_samples(target.aa):
        out.append(f"{where}: notes of {n} samples, the tempo gives {bake.note_samples(target.aa)}")
        return out
    if onsets != [bake.onset_ms(i, n) for i in range(target.aa)]:
        out.append(f"{where}: the onsets are not one note apart from zero")
    if facts.samples != target.aa * n:
        out.append(f"{where}: the file plays {facts.samples} samples, "
                   f"{target.aa} notes of {n} are {target.aa * n}")
    if audio.get("duration_ms") != bake.onset_ms(target.aa, n):
        out.append(f"{where}: the last note does not end where the file does")

    # The channels, each from its input.
    sequence = carried.get("sequence")
    for key in ("pitch", "timbre", "loudness", "accent"):
        if key not in carried:
            out.append(f"{where}: no {key}")
    if out:
        return out
    if not isinstance(sequence, str) or len(sequence) != target.aa:
        return out + [f"{where}: the sequence is not the protein's length"]
    if len(carried["pitch"]) != target.aa or len(carried["timbre"]) != target.aa:
        return out + [f"{where}: a channel is not one entry per residue"]
    if carried["loudness"] is not None and len(carried["loudness"]) != target.aa:
        return out + [f"{where}: loudness is not one entry per residue"]

    sources = carried.get("sources") or {}
    read = record_of(target)
    if read is None:
        return out + [f"{where}: no stored record to hold it to; run fetch_tracks.py --kind record"]
    record, record_blob = read
    if sequence != record["protein"]["translation"]:
        out.append(f"{where}: the sequence is not the record's translation")
    if (sources.get("record") or {}).get("sha256") != hashlib.sha256(record_blob).hexdigest():
        out.append(f"{where}: was made from another record than the one on disk")
    if carried["pitch"] != [bake.pitch_of(letter) for letter in sequence]:
        out.append(f"{where}: a pitch is not its residue's hydropathy on the scale")

    wanted = {
        "folding": bake.folding_asset(target) if target in bake.FOLDING_TARGETS else None,
        "constraint": target.constraint_asset if target.scored else None,
        "clinvar": bake.clinvar_asset(target) if target.clinvar_available else None,
    }
    inputs = {}
    for name, relative in wanted.items():
        found, data = _held_to(where, name, sources.get(name), relative)
        out += found
        inputs[name] = data
    try:
        if wanted["folding"] is None or inputs["folding"] is not None:
            if carried["timbre"] != bake.timbres_from(inputs["folding"], target.aa):
                out.append(f"{where}: a timbre is not the folding track's structure")
        if wanted["constraint"] is None or inputs["constraint"] is not None:
            want = bake.loudness_from(inputs["constraint"], sequence)
            if (carried["loudness"] if carried["loudness"] is None else tuple(carried["loudness"])) != want:
                out.append(f"{where}: a loudness is not the constraint track's conservation")
        if wanted["clinvar"] is None or inputs["clinvar"] is not None:
            if tuple(carried["accent"]) != bake.accents_from(inputs["clinvar"], sequence):
                out.append(f"{where}: the ticks are not the residues ClinVar has records at")
    except (ValueError, KeyError) as problem:
        out.append(f"{where}: {problem}")

    want_mapping = bake.mapping(wanted["folding"] is not None, wanted["constraint"] is not None,
                                wanted["clinvar"] is not None)
    if carried.get("mapping") != want_mapping:
        out.append(f"{where}: the mapping it names is not the one it was baked with")

    if hear and not out:
        out += heard_problems(where, blob, carried)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--target", action="append", default=None)
    args = parser.parse_args()
    unknown = sorted(set(args.target or []) - set(BY_SLUG))
    if unknown:
        parser.error(f"no such target: {unknown}")
    chosen = (tuple(BY_SLUG[s] for s in args.target) if args.target
              else bake.AUDIO_TARGETS)
    problems: list[str] = []
    for target in chosen:
        found = problems_of(target)
        problems += found
        if not found:
            blob = (DATA / bake.audio_asset(target)).read_bytes()
            heard = listen(blob, m4a.read_map(blob))
            ticks = "   none" if heard.ticks_db is None else f"{heard.ticks_db:5.1f} dB"
            print(f"{target.slug:<16} ok  lag {heard.lag:+d}  pitch {heard.pitches:6.1%}  "
                  f"level {heard.level_db:4.2f} dB  timbre {heard.timbres:6.1%}  ticks {ticks}")
    for problem in problems:
        print(f"  {problem}", file=sys.stderr)
    print(f"{len(chosen)} tracks checked, {len(problems)} problems", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
