"""Bake every exact-allele AVI attribution for insulin, the contributions showcase.

Schema 2 carries, for each substitution the AVI asset scores, the Atlas's raw AVI
score, all 18 signed feature attributions (not a renormalized percentage) and,
for the assay features, the track with the largest absolute effect
(`track_winners.py`). Raw window responses checkpoint outside mobile assets;
the per-substitution track replies are reduced on arrival and only their
winners checkpoint.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime as dt
import hashlib
import importlib.metadata
import json
import math
import os
import pathlib
import statistics
import subprocess
import sys
import time
import urllib.parse

BACKEND = pathlib.Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.impact.bake_impact import BASES, COMPLEMENT, BakeError, phred_of, windows  # noqa: E402
from pipeline.impact import track_winners  # noqa: E402
from pipeline.paths import DATA, SKILLS  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
SCORER = "AVI_SCORE_FEATURE_IMPORTANCE"
PILOT = ("insulin",)
# Substitutions carry no insertion or deletion; a non-zero flag means a misread row.
INDELS = ("IS_INSERTION", "IS_DELETION")
# raw - sum(attributions) must hold this still to be stated as one baseline.
BASELINE_SPREAD = 1e-3
WORKERS = 3


def finite(values: list[float]) -> list[float]:
    if not values or not all(math.isfinite(v) for v in values):
        raise BakeError("Non-finite or empty attribution row")
    return values


def query(client, chromosome: str, low: int, high: int) -> dict:
    from alphagenome.data import genome
    import numpy as np

    for attempt in range(8):
        try:
            result = client.query_interval(
                genome.Interval(chromosome, low - 1, high),
                requested_scorers=["AVI_SCORE", SCORER],
                max_workers=4, progress_bar=False,
            )
            break
        except Exception as error:
            if attempt == 7 or not any(word in str(error) for word in track_winners.RETRIED):
                raise
            delay = min(5 * 2**attempt, 60)
            print(f"  retrying Atlas in {delay}s", flush=True)
            time.sleep(delay)
    fi = result[SCORER]
    avi = result["AVI_SCORE"]
    names = list(fi.var["name"] if "name" in fi.var else fi.var_names)
    scores = {
        str(v): (round(phred_of(float(q)), 1), float(x))
        for v, q, x in zip(avi.obs["variant"], np.ravel(avi.layers["quantiles"]),
                           np.ravel(avi.X), strict=True)
    }
    rows = []
    for v, row in zip(fi.obs["variant"], fi.X, strict=True):
        values = [float(n) for n in np.ravel(row)]
        if len(values) != len(names) or str(v) not in scores:
            raise BakeError("Attribution and AVI rows do not align")
        phred, raw = scores[str(v)]
        finite([raw, *values])  # gate before checkpointing
        rows.append([int(v.position), v.reference_bases, v.alternate_bases, phred, raw, values])
    return {"features": names, "rows": rows}


def atlas_template() -> str:
    skill = SKILLS / "alphagenome-atlas-website-links"
    variant = "chr1:1:A>C"
    url = subprocess.check_output(
        ["uv", "run", "scripts/alphagenome_atlas_links.py", "variant", variant],
        cwd=skill, text=True,
    ).strip()
    parts = urllib.parse.urlsplit(url)
    query = urllib.parse.parse_qsl(parts.query)
    if dict(query).get("q") != variant:
        raise BakeError("Unexpected Atlas link format")
    encoded = urllib.parse.urlencode([(k, "{variant}" if k == "q" else v) for k, v in query])
    return urllib.parse.urlunsplit(parts._replace(query=encoded)).replace("%7Bvariant%7D", "{variant}")


def gather_winners(client, chromosome: str, keys: list[tuple[int, str, str]],
                   checkpoint: pathlib.Path) -> dict[tuple[int, str, str], dict]:
    """Every substitution's track winners, resumed from the JSONL checkpoint."""
    from alphagenome.data import genome

    done = {}
    if checkpoint.exists():
        for line in checkpoint.read_text().splitlines():
            entry = json.loads(line)
            done[tuple(entry["key"])] = entry["winners"]
    todo = [k for k in keys if k not in done]
    print(f"  track winners: {len(done):,} checkpointed, {len(todo):,} to ask", flush=True)
    with checkpoint.open("a") as out, ThreadPoolExecutor(WORKERS) as pool:
        futures = {
            pool.submit(track_winners.winners, client,
                        genome.Variant(chromosome, position, ref, alt)): (position, ref, alt)
            for position, ref, alt in todo
        }
        for n, future in enumerate(as_completed(futures), 1):
            key = futures[future]
            done[key] = future.result()
            out.write(json.dumps({"key": list(key), "winners": done[key],
                                  "at": dt.datetime.now(dt.timezone.utc).isoformat()},
                                 separators=(",", ":")) + "\n")
            out.flush()
            if n % 250 == 0 or n == len(todo):
                print(f"  track winners: {n:,}/{len(todo):,}", flush=True)
    return done


def bake(slug: str, client, template: str) -> None:
    source = DATA / f"assets/impact/{slug}_avi.json"
    raw_bytes = source.read_bytes()
    impact = json.loads(raw_bytes)
    digest = hashlib.sha256(raw_bytes).hexdigest()
    mapping = {
        r["local"] + i: r["genomic"] + r["step"] * i
        for r in impact["runs"] for i in range(r["length"])
    }
    wanted = {mapping[int(p)] for p in impact["positions"]}
    plan = windows(wanted, chunk=1024, merge_gap=32)
    cache_dir = HERE / "attribution_raw" / digest
    cache_dir.mkdir(parents=True, exist_ok=True)
    gathered = {}
    features = None
    for i, (low, high) in enumerate(plan, 1):
        checkpoint = cache_dir / f"{low}-{high}.json"
        payload = json.loads(checkpoint.read_text()) if checkpoint.exists() else None
        if payload is None or any(len(row) != 6 for row in payload["rows"]):
            # Absent, or written before the raw score was kept.
            payload = query(client, impact["chromosome"], low, high)
            payload["retrieved_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
            checkpoint.write_text(json.dumps(payload, separators=(",", ":")))
            time.sleep(1.5)
        if features is not None and features != payload["features"]:
            raise BakeError("Feature metadata changed during bake")
        features = payload["features"]
        for position, ref, alt, phred, raw, values in payload["rows"]:
            gathered[position, ref, alt] = phred, raw, values
        print(f"{slug}: {i}/{len(plan)} windows", flush=True)

    indels = [features.index(name) for name in INDELS]
    keys = {}
    for local, existing in impact["positions"].items():
        ref = impact["sequence"][int(local) - impact["start"]]
        for alt, expected in zip((b for b in BASES if b != ref), existing, strict=True):
            gr, ga = (ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)) if impact["complemented"] else (ref, alt)
            key = mapping[int(local)], gr, ga
            if key not in gathered:
                raise BakeError(f"No exact attribution for {key}")
            phred, raw, values = gathered[key]
            if phred != expected:
                raise BakeError(f"AVI changed at {key}: {expected} -> {phred}; rebake AVI first")
            if any(values[i] != 0 for i in indels):
                raise BakeError(f"Insertion/deletion attribution on a substitution at {key}")
            keys[local, alt] = key

    found = gather_winners(client, impact["chromosome"], sorted(set(keys.values())),
                           cache_dir / "track_winners.jsonl")

    tracks: list[list[str]] = []
    index: dict[tuple[str, ...], int] = {}
    offsets = []
    positions = {}
    for local, existing in impact["positions"].items():
        ref = impact["sequence"][int(local) - impact["start"]]
        rows = []
        for alt in (b for b in BASES if b != ref):
            phred, raw, values = gathered[keys[local, alt]]
            offsets.append(raw - math.fsum(values))
            won = []
            for feature, (*track, value) in sorted(found[keys[local, alt]].items(),
                                                   key=lambda item: features.index(item[0])):
                if tuple(track) not in index:
                    index[tuple(track)] = len(tracks)
                    tracks.append(track)
                won.append([features.index(feature), index[tuple(track)], round(value, 6)])
            rows.append([phred, round(raw, 6), [round(v, 6) for v in values], won])
        positions[local] = rows

    spread = max(offsets) - min(offsets)
    offset = statistics.median(offsets)
    print(f"  raw - sum: median {offset:+.6f}, spread {spread:.6f}", flush=True)
    asset = {
        "schema_version": 2,
        **{k: impact[k] for k in ("gene", "uniprot", "accession", "assembly", "annotation", "chromosome", "transcript", "start", "sequence", "complemented", "orientation", "runs")},
        "scorer": SCORER,
        "scope": "across_atlas_genes_and_biosamples",
        "units": "raw_score_attribution",
        "selection": "all",
        "alt_order": "ACGT minus wildtype",
        "impact_sha256": digest,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "client_version": importlib.metadata.version("alphagenome"),
        "atlas_release": None,
        "atlas_url_template": template,
        # Stated only where raw - sum holds still; the offsets are recorded either way.
        "baseline": round(offset, 6) if spread <= BASELINE_SPREAD else None,
        "raw_minus_sum": {"median": round(offset, 6), "min": round(min(offsets), 6),
                          "max": round(max(offsets), 6)},
        "track_scorers": list(track_winners.REQUESTED),
        "features": features,
        "tracks": tracks,
        "positions": positions,
    }
    out = DATA / f"assets/impact_explanations/{slug}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    temp = out.with_suffix(".tmp")
    temp.write_text(json.dumps(asset, separators=(",", ":")) + "\n")
    temp.replace(out)
    print(f"Wrote {out.name}: {len(positions) * 3:,} exact alternatives, {len(tracks):,} tracks, "
          f"{out.stat().st_size / 1024:.0f} KiB", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", choices=PILOT, action="append")
    args = parser.parse_args()
    # Same credential source as the existing offline bake; never serialized.
    import dotenv
    from alphagenome.atlas import atlas
    dotenv.load_dotenv(pathlib.Path.home() / ".env")
    template = atlas_template()
    client = atlas.create(os.environ["ALPHAGENOME_API_KEY"])
    for slug in args.target or PILOT:
        bake(slug, client, template)


if __name__ == "__main__":
    main()
