"""Bake `trafficking`: where in the cell each protein goes, as UniProt has it.

    .venv/Scripts/python pipeline/trafficking/bake_trafficking.py --all
    .venv/Scripts/python pipeline/trafficking/bake_trafficking.py --target cftr

The app works out a protein's route through the cell from its sequence
features, and four of them are in the records already: the signal peptide, the
disulfides, the cuts, and the prion's GPI-anchor signal, all regions of the
precursor. The fifth is not. A region table says what a precursor is cut
into, not where any of it sits, so no record says whether a stretch crosses a
membrane, and nineteen of the twenty routes stop there. This track carries
what the records cannot, per protein, read from its UniProt entry:

- `transmembrane`: each `Transmembrane` feature, as a span in precursor
  numbering with UniProt's description ("Helical; Signal-anchor for type II
  membrane protein") and its evidence codes. An empty list means UniProt
  annotates none, which is not the same as not having looked.
- `gpi_anchor`: the residue a `GPI-anchor` lipidation is attached to, and the
  signal after it that is cut off in its place, or null.
- `location`: every place the entry's subcellular location comments name,
  with the topology and orientation UniProt gives each, and the molecule a
  comment is about where it is about one (an isoform, a chain cut from it).

Each payload names its protein, its length and the entry's own versions, and
the release that served it and the day it was fetched. So two bakes within
one release differ in `retrieved` alone.

Which proteins get the track: all twenty (`TRAFFICKING_TARGETS`), recorded
here, never in targets.py.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.paths import DATA  # noqa: E402
from pipeline.targets import BY_SLUG, TARGETS, Target  # noqa: E402
from pipeline.uniprot import Entry, fetch_entry  # noqa: E402

# Every protein gets the track: each has a route to draw.
TRAFFICKING_TARGETS = TARGETS

SCHEMA_VERSION = 1

# The lipidation UniProt writes for a GPI anchor begins with this, and goes on
# to name the residue: "GPI-anchor amidated serine".
_GPI = "GPI-anchor"


def trafficking_asset(target: Target) -> str:
    return f"assets/trafficking/{target.slug}_trafficking.json"


def _evidence(node: dict) -> list[str]:
    return sorted({e["evidenceCode"] for e in node.get("evidences", []) if "evidenceCode" in e})


def _span(feature: dict) -> tuple[int, int]:
    """A feature's first and last residue, refusing one with an unknown end."""
    start = feature["location"]["start"].get("value")
    end = feature["location"]["end"].get("value")
    if start is None or end is None:
        raise ValueError(f"{feature['type']} with an unknown end: {feature['location']}")
    return int(start), int(end)


def transmembrane(body: dict) -> list[dict]:
    spans = []
    for feature in body.get("features", []):
        if feature["type"] != "Transmembrane":
            continue
        start, end = _span(feature)
        spans.append({
            "start": start,
            "end": end,
            "description": feature.get("description", ""),
            "evidence": _evidence(feature),
        })
    return sorted(spans, key=lambda s: s["start"])


def gpi_anchor(body: dict) -> dict | None:
    """The GPI-anchored residue, and everything after it, which is the signal.

    The transamidase cuts the chain after that residue and attaches the anchor
    in the signal's place, so the signal runs from the next residue to the end.
    """
    residues = body["sequence"]["length"]
    anchors = [
        f for f in body.get("features", [])
        if f["type"] == "Lipidation" and f.get("description", "").startswith(_GPI)
    ]
    if not anchors:
        return None
    if len(anchors) > 1:
        raise ValueError(f"{len(anchors)} GPI anchors")
    site, end = _span(anchors[0])
    if site != end:
        raise ValueError(f"a GPI anchor on {site}-{end}, not on one residue")
    return {
        "site": site,
        "residue": anchors[0]["description"],
        "signal": {"start": site + 1, "end": residues},
        "evidence": _evidence(anchors[0]),
    }


def location(body: dict) -> list[dict]:
    found = []
    for comment in body.get("comments", []):
        if comment.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for place in comment.get("subcellularLocations", []):
            found.append({
                "molecule": comment.get("molecule"),
                "location": place["location"]["value"],
                "topology": (place.get("topology") or {}).get("value"),
                "orientation": (place.get("orientation") or {}).get("value"),
                "evidence": _evidence(place["location"]),
            })
    return found


def payload(target: Target, entry: Entry) -> dict:
    """The track for one protein, from the entry UniProt served for it."""
    body = entry.body
    if body.get("primaryAccession") != target.uniprot:
        raise ValueError(f"{target.slug}: asked for {target.uniprot}, "
                         f"UniProt answered {body.get('primaryAccession')}")
    if not entry.release or not entry.release_date:
        raise ValueError(f"{target.slug}: UniProt did not say which release served it")
    audit = body.get("entryAudit", {})
    return {
        "slug": target.slug,
        "gene": target.gene,
        "uniprot": target.uniprot,
        "residues": body["sequence"]["length"],
        "source": "UniProtKB",
        "release": entry.release,
        "release_date": entry.release_date,
        "retrieved": entry.retrieved,
        "entry_version": audit.get("entryVersion"),
        "sequence_version": audit.get("sequenceVersion"),
        "schema_version": SCHEMA_VERSION,
        "transmembrane": transmembrane(body),
        "gpi_anchor": gpi_anchor(body),
        "location": location(body),
        "built_by": "pipeline/trafficking/bake_trafficking.py",
    }


def encode(track: dict) -> bytes:
    return (json.dumps(track, indent=2, ensure_ascii=False) + "\n").encode()


def bake(target: Target) -> dict:
    track = payload(target, fetch_entry(target.uniprot))
    destination = DATA / trafficking_asset(target)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(encode(track))
    return track


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
    chosen = TRAFFICKING_TARGETS if args.all else tuple(BY_SLUG[s] for s in args.target)
    for target in chosen:
        track = bake(target)
        spans = track["transmembrane"]
        gpi = track["gpi_anchor"]
        places = sorted({p["location"] for p in track["location"] if p["molecule"] is None})
        print(f"{target.slug:<16} {track['residues']:>5} aa  "
              f"{len(spans):>2} span(s)  {'GPI ' + str(gpi['site']) if gpi else 'no GPI':<8}  "
              f"{'; '.join(places) or '-'}  ({track['release']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
