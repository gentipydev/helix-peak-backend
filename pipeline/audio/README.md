# Baking `audio`: each protein as sound, for Listen

The lab's Listen plays a protein one note a residue, with a playhead moving
over its grid. It needs, per protein, an audio file and a timing map from each
residue to the millisecond its note starts, so the playhead can follow the
player's own position rather than a clock of its own. Nothing else carries
either.

This is a separate track kind, `audio`. Nothing is added to `targets.py`,
`curated/catalog.json` or the `protein` table, and no other track is touched:
the bake only reads the record, folding, constraint and ClinVar tracks.

```sh
.venv/Scripts/python -m pip install -r pipeline/audio/requirements.txt   # PyAV
.venv/Scripts/python pipeline/fetch_tracks.py --kind record --kind constraint --kind clinvar
.venv/Scripts/python pipeline/audio/bake_audio.py --all      # writes pipeline/data/assets/audio/
.venv/Scripts/python pipeline/audio/check_audio.py           # every check below; exits 1 on any failure
.venv/Scripts/python pipeline/upload_tracks.py --kind audio --dry-run
```

The folding tracks are the folding bake's (`folding/bake_folding.py`), until
they are uploaded and `fetch_tracks.py --kind folding` can bring them back.

Upload only when asked: `DATABASE_URL`, `SUPABASE_URL` and `SUPABASE_SERVICE_KEY`
set, `migrations/0008_audio.sql` applied (after `0007`), then the same command
without `--dry-run`. Objects go to the `tracks` bucket as
`audio/<slug>.<sha12>.m4a`, `audio/mp4`: the first non-JSON object in that
bucket, so if the bucket was given a list of allowed types, `audio/mp4` has to
join it first. The playbook named the migration `0007_audio`; `0007` is the
locus track's.

Which proteins get the track: all twenty (`AUDIO_TARGETS` in `bake_audio.py`).

## The mapping

One note a residue, and four properties, each on a channel of its own:

| channel | property | from | how |
|---|---|---|---|
| pitch | hydropathy | Kyte & Doolittle 1982 | quantised to the nearest step of C major pentatonic from C4 to C6; more hydrophobic is higher |
| timbre | secondary structure | the `folding` track | helix a held reed, strand a plucked string, coil a soft bell; anything the fold does not place a hollow swell |
| loudness | conservation | the ESM-2 `constraint` track | the most conserved residue at full level, the least 12 dB under it |
| accent | a ClinVar record at the residue | the `clinvar` snapshot | a short tick as the note starts, whatever the record's classification |

The mapping is arbitrary, chosen for teaching; the properties behind it are
measured. So a signal peptide's hydrophobic core sounds as a high run
(insulin's LLPLLALLAL, residues 7-16, is its highest ten notes), a helix holds
its notes where a strand plucks them, and a residue the language model finds
tolerant of any substitution plays quietly.

What a protein does not have, the piece says plainly: a residue outside the
solved fold (a signal peptide, a C peptide, dystrophin's 3,447 unsolved
residues) is the hollow swell; a protein with no folding track at all is that
one timbre throughout; one with no constraint track plays every note at the
middle of the range, -6 dB; one with no ClinVar snapshot has no ticks. A track
a protein should have and the disk does not hold is an error, never a plainer
piece.

Each timbre is as loud as the others (the same RMS), every partial it plays
lies under 4.4 kHz, and the tick, at 3.7 kHz, lies at least 170 Hz from any
harmonic a note carries, so a tick is never mistaken for a note.

## The tempo, and the size

A residue lasts 125 ms, eight a second, up to 960 residues. A longer protein
plays faster so that no piece runs past two minutes: CFTR at about twelve
notes a second, dystrophin at about thirty-one. A protein that could only fit
by notes under 25 ms, clicks rather than pitches, is refused: over 4,800
residues (titin's 34,350 would be).

This is the one kind that could be large, and the lab's cache is 100 MB
(`LabScope.budget`). At one tempo for all, dystrophin's 3,685 residues would
play for seven minutes and 41 seconds, over a megabyte, and CFTR's for three;
with the two-minute ceiling the largest file is dystrophin's at 377 KB, and
all twenty come to 2.25 MB, 2.2% of the lab's budget. Everything the bake
decides about size is two constants: `NOTE_SAMPLES` and `LONGEST_SAMPLES`.

## What the file is

An `.m4a`: AAC-LC in an MP4 container, mono, 16 kHz, 20 kbit/s, the encoder's
low-pass at 4.5 kHz. The rate was measured, not guessed: on insulin,
hemoglobin, CFTR and dystrophin, encoded, decoded and compared note by note
with the synthesis, AAC at 16 kbit/s lost the timbre of up to one note in nine;
20 kbit/s with the low-pass kept pitch and timbre on 98 to 100% of notes and
levels within a decibel; MP3 at 16 kbit/s did as well in fewer bytes. MP4 won
anyway, on timing: an AAC encoder
emits 1,024 samples of priming, and the container's edit list says to skip
them, which ExoPlayer, AVFoundation and ffmpeg all honour, so position zero on
every player is the first note. MP3 carries its delay in an encoder tag not
every player reads, and a playhead would lead the sound by up to 70 ms there.

The timing map rides inside the file, in a `uuid` box after the audio: the
standard's own place for data it does not define (ISO/IEC 14496-12). Players
skip it, so the file plays anywhere as what it is; the app finds it by its
identifier (`MAP_UUID` in `m4a.py`). The track is one object and one digest.
`m4a.py` reads the boxes with the standard library, so the check and the
uploader need no decoder to hold the file to its map.

A bake is deterministic on one machine: re-baked here, all twenty came out
byte for byte the same, and ffmpeg's bit-exact flags keep its version out of
the file. Another machine's floating point can move a sample, and with it the
bytes, as the structure bakes' can, so across machines the proof is
`check_audio.py` (the file sounds like its map), not the digest.

## What the map holds

One JSON object per protein, inside the file:

| field | what |
|---|---|
| `slug`, `gene`, `uniprot`, `schema_version`, `residues`, `sequence` | the protein, and the record's translation it plays |
| `audio` | format, codec, sample rate, channels, bit rate, cut-off, the priming the edit list skips, the samples it plays, and `duration_ms` |
| `tempo` | `note_samples`, `note_ms`, the two-minute ceiling and the rule |
| `onset_ms` | **the timing map**: one integer per residue, the millisecond its note starts, to the nearest; residue `i` (from 1) is `onset_ms[i - 1]` |
| `pitch` | one MIDI note number per residue |
| `timbre` | one letter per residue: `H` helix, `E` strand, `C` coil, `-` none |
| `loudness` | one conservation per residue, to two places, or null for a protein with no constraint track |
| `accent` | the residues a ClinVar record is reported at, in order |
| `mapping` | which property drives which channel, the scale and its steps, the hydropathy of each residue, each timbre and what it sounds like, the loudness range, and which track each channel came from (null where the protein has none) |
| `sources` | the record, folding, constraint and ClinVar files it was made from: path, sha256, and the accession, PDB entry, model or retrieval date |

The row's provenance is the map less its arrays.

## What the check holds

`check_audio.py`, and the uploader before it writes, hold each track to:

- its file: one mono AAC-LC stream at 16 kHz, an `M4A ` file whose edit list
  skips the priming and plays exactly the piece, one note a residue at the
  tempo's length;
- its map: one, the file's last box, naming its protein;
- **the timing map's length: one onset per residue, exactly the protein's
  residue count**, the first at zero, each one note after the last, the last
  note ending where the file does;
- its inputs: each pitch the hydropathy of the record's residue, each timbre
  the folding track's structure, each loudness the constraint track's
  conservation, each tick a residue ClinVar has a record at, the inputs read
  from disk and held to the digests the map names;
- its mapping: the one it was baked with;
- its sound: decoded, the audio starts on its first note (lag 0), and note by
  note it keeps its map, synthesised afresh: at least 97% of pitches, levels
  within 1.5 dB on average, at least 95% of notes heard as their own timbre,
  and a tick adding at least 10 dB to a note's opening at 3.7 kHz.

## What it found

| slug | residues | note | per second | seconds | helix | strand | coil | none | ticks | audio B | map B | file B |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| insulin | 110 | 125.0 ms | 8.0 | 13.8 | 30 | 0 | 21 | 59 | 65 | 36,565 | 4,025 | 40,590 |
| hemoglobin | 147 | 125.0 ms | 8.0 | 18.4 | 127 | 0 | 18 | 2 | 146 | 48,429 | 4,896 | 53,325 |
| myoglobin | 154 | 125.0 ms | 8.0 | 19.2 | 132 | 0 | 17 | 5 | 23 | 50,672 | 4,615 | 55,287 |
| p53 | 393 | 125.0 ms | 8.0 | 49.1 | 24 | 64 | 106 | 199 | 391 | 128,000 | 9,767 | 137,767 |
| lysozyme | 148 | 125.0 ms | 8.0 | 18.5 | 55 | 8 | 67 | 18 | 45 | 48,753 | 4,579 | 53,332 |
| relaxin | 185 | 125.0 ms | 8.0 | 23.1 | 42 | 4 | 5 | 134 | 40 | 60,563 | 5,161 | 65,724 |
| oxytocin | 125 | 125.0 ms | 8.0 | 15.6 | 0 | 0 | 9 | 116 | 9 | 41,286 | 4,095 | 45,381 |
| somatotropin | 217 | 125.0 ms | 8.0 | 27.1 | 104 | 0 | 82 | 31 | 94 | 71,054 | 5,888 | 76,942 |
| ubiquitin | 229 | 125.0 ms | 8.0 | 28.6 | 16 | 33 | 27 | 153 | 14 | 74,741 | 5,679 | 80,420 |
| dystrophin | 3,685 | 32.6 ms | 30.7 | 120.0 | 167 | 3 | 68 | 3,447 | 3,038 | 302,275 | 75,165 | 377,440 |
| vasopressin | 164 | 125.0 ms | 8.0 | 20.5 | 0 | 0 | 9 | 155 | 86 | 53,780 | 4,977 | 58,757 |
| glucagon | 180 | 125.0 ms | 8.0 | 22.5 | 28 | 0 | 1 | 151 | 3 | 59,035 | 4,968 | 64,003 |
| app | 770 | 125.0 ms | 8.0 | 96.2 | 27 | 58 | 77 | 608 | 293 | 249,507 | 15,413 | 264,920 |
| cftr | 1,480 | 81.1 ms | 12.3 | 120.0 | 807 | 91 | 241 | 341 | 1,398 | 310,672 | 31,777 | 342,449 |
| erythropoietin | 193 | 125.0 ms | 8.0 | 24.1 | 105 | 6 | 55 | 27 | 40 | 63,358 | 5,319 | 68,677 |
| leptin | 167 | 125.0 ms | 8.0 | 20.9 | 95 | 0 | 35 | 37 | 53 | 54,900 | 4,917 | 59,817 |
| tnf | 233 | 125.0 ms | 8.0 | 29.1 | 5 | 78 | 70 | 80 | 16 | 76,165 | 5,818 | 81,983 |
| sod1 | 154 | 125.0 ms | 8.0 | 19.2 | 17 | 58 | 78 | 1 | 110 | 50,646 | 4,884 | 55,530 |
| amylase | 511 | 125.0 ms | 8.0 | 63.9 | 134 | 85 | 277 | 15 | 12 | 166,138 | 10,215 | 176,353 |
| prion | 253 | 125.0 ms | 8.0 | 31.6 | 71 | 13 | 25 | 144 | 115 | 82,577 | 6,516 | 89,093 |

2,247,790 bytes for all twenty. Decoded, every one starts on its first note;
the worst any keeps of its map is dystrophin's 98.9% of pitches, at 33 ms a
note, and SOD1's 98.1% of timbres.

- **Dystrophin** plays for two minutes at 31 notes a second, and all but 238 of
  its residues in the hollow swell: its only solved stretch is the actin-binding
  domain (1DXX). 3,038 of its residues tick.
- **Hemoglobin** ticks at 146 of 147 residues and **p53** at 391 of 393: almost
  every position of either has a ClinVar record. Glucagon and amylase tick at
  3 and 12.
- **Oxytocin** and **vasopressin** are nine placed residues each, the hormone,
  in a precursor of 125 and 164: the rest is the swell.
- **CFTR**'s 341 unplaced residues are mostly its regulatory region, which the
  fold page already says is too mobile to resolve.
