"""Bake exact-allele AVI attributions for the first three supported genes.

Raw responses checkpoint outside mobile assets. Mobile carries the three largest
absolute contributions, with their signs, not a renormalized percentage.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.metadata
import json
import math
import os
import pathlib
import subprocess
import sys
import time
import urllib.parse

BACKEND = pathlib.Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.impact.bake_impact import BASES, COMPLEMENT, BakeError, phred_of, windows  # noqa: E402
from pipeline.paths import DATA, SKILLS  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
SCORER = "AVI_SCORE_FEATURE_IMPORTANCE"
PILOT = ("insulin", "hemoglobin", "cftr")


def compact(values: list[float]) -> list[list[float | int]]:
    if not values or not all(math.isfinite(v) for v in values):
        raise BakeError("Non-finite or empty attribution row")
    order = sorted(range(len(values)), key=lambda i: (-abs(values[i]), i))
    return [[i, round(values[i], 7)] for i in order[:3] if abs(values[i]) >= 0.00000005]


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
            if attempt == 7 or not any(word in str(error) for word in (
                "RESOURCE_EXHAUSTED", "UNAVAILABLE", "DEADLINE_EXCEEDED",
                "failed to connect", "Broken pipe",
            )):
                raise
            delay = min(5 * 2**attempt, 60)
            print(f"  retrying Atlas in {delay}s", flush=True)
            time.sleep(delay)
    fi = result[SCORER]
    avi = result["AVI_SCORE"]
    names = list(fi.var["name"] if "name" in fi.var else fi.var_names)
    scores = {
        str(v): round(phred_of(float(q)), 1)
        for v, q in zip(avi.obs["variant"], np.ravel(avi.layers["quantiles"]), strict=True)
    }
    rows = []
    for v, row in zip(fi.obs["variant"], fi.X, strict=True):
        values = [float(n) for n in np.ravel(row)]
        if len(values) != len(names) or str(v) not in scores:
            raise BakeError("Attribution and AVI rows do not align")
        compact(values)  # finite-value gate before checkpointing
        rows.append([int(v.position), v.reference_bases, v.alternate_bases,
                     scores[str(v)], values])
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


def bake(slug: str, client, template: str) -> None:
    source = DATA / f"assets/impact/{slug}_avi.json"
    raw = source.read_bytes()
    impact = json.loads(raw)
    digest = hashlib.sha256(raw).hexdigest()
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
        if checkpoint.exists():
            payload = json.loads(checkpoint.read_text())
        else:
            payload = query(client, impact["chromosome"], low, high)
            payload["retrieved_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
            checkpoint.write_text(json.dumps(payload, separators=(",", ":")))
            time.sleep(1.5)
        if features is not None and features != payload["features"]:
            raise BakeError("Feature metadata changed during bake")
        features = payload["features"]
        for position, ref, alt, phred, values in payload["rows"]:
            gathered[position, ref, alt] = phred, values
        print(f"{slug}: {i}/{len(plan)} windows", flush=True)

    positions = {}
    for local, existing in impact["positions"].items():
        ref = impact["sequence"][int(local) - impact["start"]]
        rows = []
        for alt, expected in zip((b for b in BASES if b != ref), existing, strict=True):
            gr, ga = (ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)) if impact["complemented"] else (ref, alt)
            key = mapping[int(local)], gr, ga
            if key not in gathered:
                raise BakeError(f"No exact attribution for {key}")
            phred, values = gathered[key]
            if phred != expected:
                raise BakeError(f"AVI changed at {key}: {expected} -> {phred}; rebake AVI first")
            rows.append([phred, compact(values)])
        positions[local] = rows

    asset = {
        "schema_version": 1,
        **{k: impact[k] for k in ("gene", "uniprot", "accession", "assembly", "annotation", "chromosome", "transcript", "start", "sequence", "complemented", "orientation", "runs")},
        "scorer": SCORER,
        "scope": "across_atlas_genes_and_biosamples",
        "units": "raw_score_attribution",
        "selection": "top_3_absolute_signed",
        "alt_order": "ACGT minus wildtype",
        "impact_sha256": digest,
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "client_version": importlib.metadata.version("alphagenome"),
        "atlas_release": None,
        "atlas_url_template": template,
        "features": features,
        "positions": positions,
    }
    out = DATA / f"assets/impact_explanations/{slug}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    temp = out.with_suffix(".tmp")
    temp.write_text(json.dumps(asset, separators=(",", ":")) + "\n")
    temp.replace(out)
    print(f"Wrote {out.name}: {len(positions) * 3:,} exact alternatives, {out.stat().st_size / 1024:.0f} KiB", flush=True)


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
