"""The parts of an .m4a the audio track relies on, read and written in plain Python.

The bake encodes with PyAV: ffmpeg's AAC encoder and its MP4 muxer. Everything
the track is then held to -- that the file is one mono AAC-LC stream at the
rate its map was cut for, that its edit list starts playback where the first
note starts, how many samples it plays -- is written in the file's boxes
(ISO/IEC 14496-12 and 14496-14), and so is the map itself. None of it needs a
decoder: the check and the uploader read the file with the standard library,
and the app reads the map the same way.

**The map rides in a `uuid` box**, the standard's own place for data it does
not define. A reader skips a box it does not know, so every player plays the
file as the audio it is, and the app finds the box by its identifier. The box
goes after everything the muxer wrote, so no offset in the file moves.

**The edit list is what makes the map exact.** An AAC encoder emits 1024
samples of priming before the first real one; the muxer's edit list says to
skip them, and ExoPlayer, AVFoundation and ffmpeg all honour it, so position
zero on any of them is the first note's first sample. That is why the track is
AAC in MP4 and not MP3, whose delay rides in an encoder tag not every player
reads.
"""

from __future__ import annotations

import json
import struct
from dataclasses import dataclass

# The identifier of the box that carries the map: uuid5 of the URL namespace
# and `https://helix-peak-backend.onrender.com/tracks/audio/map`. It names what
# the box holds, not a version of it; the map's `schema_version` does that.
MAP_UUID = bytes.fromhex("6a05e1b031335606aa35c3378b607c32")

# AAC's sampling frequency index (ISO/IEC 14496-3, 1.6.3.4).
_FREQUENCIES = (96000, 88200, 64000, 48000, 44100, 32000, 24000, 22050,
                16000, 12000, 11025, 8000, 7350)

# MPEG-4 audio in an elementary stream descriptor, and AAC's low-complexity
# profile in the AudioSpecificConfig: what every phone decodes.
MPEG4_AUDIO = 0x40
AAC_LC = 2


@dataclass(frozen=True)
class Box:
    kind: bytes
    start: int      # where its header begins
    size: int       # header included
    header: int     # 8, 16 with a 64-bit size, and 16 more for a `uuid`'s identifier

    @property
    def end(self) -> int:
        return self.start + self.size

    def body(self, blob: bytes) -> bytes:
        return blob[self.start + self.header:self.end]


def boxes(blob: bytes, start: int = 0, end: int | None = None) -> list[Box]:
    """The boxes laid end to end in `blob[start:end]`, refusing any that overrun it."""
    end = len(blob) if end is None else end
    found: list[Box] = []
    at = start
    while at < end:
        if end - at < 8:
            raise ValueError(f"{end - at} stray bytes at {at}")
        size, kind = struct.unpack_from(">I4s", blob, at)
        header = 8
        if size == 1:
            if end - at < 16:
                raise ValueError(f"the {kind!r} box at {at} is cut short")
            size = struct.unpack_from(">Q", blob, at + 8)[0]
            header = 16
        elif size == 0:
            size = end - at
        if kind == b"uuid":
            header += 16
        if size < header or at + size > end:
            raise ValueError(f"the {kind!r} box at {at} claims {size} bytes")
        found.append(Box(kind, at, size, header))
        at += size
    return found


def find(blob: bytes, *path: bytes) -> Box:
    """The one box at `path` (`moov`, `trak`, `mdia` ...), refusing none or several."""
    start, end = 0, len(blob)
    box = None
    for depth, kind in enumerate(path):
        matches = [b for b in boxes(blob, start, end) if b.kind == kind]
        if len(matches) != 1:
            where = "/".join(p.decode("latin1") for p in path[:depth + 1])
            raise ValueError(f"{len(matches)} {where} boxes, not one")
        box = matches[0]
        start, end = box.start + box.header, box.end
    if box is None:
        raise ValueError("no path")
    return box


@dataclass(frozen=True)
class Facts:
    """What a decoder is told about the audio, read from the boxes alone."""

    brand: str          # the file type's major brand: `M4A ` for an audio-only MP4
    codec: str          # the sample entry: `mp4a`
    object_type: int    # the AudioSpecificConfig's audio object type: 2, AAC-LC
    sample_rate: int
    channels: int
    bitrate: int        # the elementary stream's average bit rate, as the muxer wrote it
    priming: int        # samples the edit list skips before the first one it plays
    samples: int        # samples the edit list plays


def _timescale(body: bytes) -> tuple[int, int]:
    """(timescale, duration) from an `mvhd` or `mdhd` body, either version."""
    if body[0] == 1:
        return struct.unpack_from(">IQ", body, 20)
    return struct.unpack_from(">II", body, 12)


def _descriptor(data: bytes, at: int) -> tuple[int, int, int]:
    """(tag, where its payload starts, its length): an MPEG-4 descriptor header."""
    tag = data[at]
    at += 1
    length = 0
    for _ in range(4):
        byte = data[at]
        at += 1
        length = (length << 7) | (byte & 0x7F)
        if not byte & 0x80:
            break
    return tag, at, length


def _audio_specific_config(esds: bytes) -> tuple[int, int, int, int, int]:
    """What an `esds` body tells a decoder.

    (object type indication, average bit rate) from the DecoderConfigDescriptor,
    then (audio object type, sampling frequency, channel configuration) from the
    AudioSpecificConfig inside it.
    """
    tag, at, length = _descriptor(esds, 4)          # past the full box's version and flags
    if tag != 0x03:
        raise ValueError("the esds holds no ES descriptor")
    flags = esds[at + 2]
    at += 3
    if flags & 0x80:
        at += 2
    if flags & 0x40:
        at += 1 + esds[at]
    if flags & 0x20:
        at += 2
    tag, at, length = _descriptor(esds, at)
    if tag != 0x04:
        raise ValueError("the ES descriptor holds no decoder config")
    indication = esds[at]
    bitrate = struct.unpack_from(">I", esds, at + 9)[0]
    tag, at, length = _descriptor(esds, at + 13)
    if tag != 0x05 or length < 2:
        raise ValueError("the decoder config holds no AudioSpecificConfig")
    bits = int.from_bytes(esds[at:at + length], "big")
    width = 8 * length
    object_type = bits >> (width - 5)
    index = (bits >> (width - 9)) & 0xF
    if object_type == 31 or index == 15:
        raise ValueError("an escaped AudioSpecificConfig; this track never writes one")
    channels = (bits >> (width - 13)) & 0xF
    if index >= len(_FREQUENCIES):
        raise ValueError(f"sampling frequency index {index}")
    return indication, bitrate, object_type, _FREQUENCIES[index], channels


def facts(blob: bytes) -> Facts:
    """The audio's type, codec, rate, channels and played length, from its boxes."""
    top = boxes(blob)
    if not top or top[0].kind != b"ftyp":
        raise ValueError("not an MP4: it does not open with a file type box")
    brand = top[0].body(blob)[:4].decode("latin1")

    movie_scale, _ = _timescale(find(blob, b"moov", b"mvhd").body(blob))
    media_scale, _ = _timescale(find(blob, b"moov", b"trak", b"mdia", b"mdhd").body(blob))

    elst = find(blob, b"moov", b"trak", b"edts", b"elst").body(blob)
    version = elst[0]
    count = struct.unpack_from(">I", elst, 4)[0]
    if count != 1:
        raise ValueError(f"{count} edits, not one")
    if version == 1:
        segment, media_time, rate = struct.unpack_from(">QqI", elst, 8)
    else:
        segment, media_time, rate = struct.unpack_from(">IiI", elst, 8)
    if rate != 0x00010000:
        raise ValueError("the edit plays at another rate than 1")
    if media_time < 0:
        raise ValueError("the edit is empty")

    stsd = find(blob, b"moov", b"trak", b"mdia", b"minf", b"stbl", b"stsd").body(blob)
    if struct.unpack_from(">I", stsd, 4)[0] != 1:
        raise ValueError("more than one sample description")
    entry = boxes(stsd, 8)[0]
    body = entry.body(stsd)
    # An AudioSampleEntry: 6 reserved bytes and a data reference index, then 8
    # reserved, channel count, sample size, 4 more, and the rate as 16.16.
    channels = struct.unpack_from(">H", body, 16)[0]
    esds = [b for b in boxes(stsd, entry.start + entry.header + 28, entry.end)
            if b.kind == b"esds"]
    if len(esds) != 1:
        raise ValueError(f"{len(esds)} elementary stream descriptors, not one")
    indication, bitrate, object_type, rate_hz, configured = _audio_specific_config(
        esds[0].body(stsd))
    if indication != MPEG4_AUDIO:
        raise ValueError(f"object type indication {indication:#x}, not MPEG-4 audio")
    if configured != channels:
        raise ValueError(f"the sample entry says {channels} channels, the stream {configured}")
    if media_scale != rate_hz:
        raise ValueError(f"the track counts {media_scale} a second, the stream is {rate_hz} Hz")
    return Facts(
        brand=brand,
        codec=entry.kind.decode("latin1"),
        object_type=object_type,
        sample_rate=rate_hz,
        channels=channels,
        bitrate=bitrate,
        priming=media_time,
        # The segment is counted in the movie's timescale; a sample in the media's.
        samples=segment * media_scale // movie_scale,
    )


def _map_boxes(blob: bytes) -> list[Box]:
    return [b for b in boxes(blob)
            if b.kind == b"uuid" and blob[b.start + b.header - 16:b.start + b.header] == MAP_UUID]


def carries_map_last(blob: bytes) -> bool:
    """Whether the map is the file's last box, where no muxer offset can move it."""
    top = boxes(blob)
    maps = _map_boxes(blob)
    return len(maps) == 1 and top[-1] == maps[0]


def read_map(blob: bytes) -> dict:
    """The map the file carries, refusing a file with none or with two."""
    maps = _map_boxes(blob)
    if len(maps) != 1:
        raise ValueError(f"{len(maps)} maps in the file, not one")
    found = json.loads(maps[0].body(blob))
    if not isinstance(found, dict):
        raise ValueError("the map is not an object")
    return found


def with_map(m4a: bytes, the_map: dict) -> bytes:
    """The file with `the_map` appended as its last box."""
    if _map_boxes(m4a):
        raise ValueError("the file already carries a map")
    body = json.dumps(the_map, ensure_ascii=False, separators=(",", ":")).encode()
    return m4a + struct.pack(">I4s", 8 + 16 + len(body), b"uuid") + MAP_UUID + body


def without_map(blob: bytes) -> bytes:
    """The file as the muxer wrote it: every map box taken out."""
    maps = _map_boxes(blob)
    return b"".join(blob[b.start:b.end] for b in boxes(blob) if b not in maps)
