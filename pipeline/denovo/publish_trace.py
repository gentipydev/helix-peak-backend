"""Publish the recorded training run the app plays, and record where it landed.

    python3 pipeline/denovo/publish_trace.py --dry-run
    python3 pipeline/denovo/publish_trace.py --dry-run --against https://helix-peak-backend.onrender.com
    DATABASE_URL=... SUPABASE_URL=... SUPABASE_SERVICE_KEY=... \
        python3 pipeline/denovo/publish_trace.py

Run from the backend's `.venv`, as `upload_tracks.py` is: the recorder needs
torch, this needs psycopg, and neither environment has both. The payload is
the file `trace.py` wrote, held to `check_trace.validate` first, so a run that
would fail on the way out never becomes a row.

The bytes go to the public tracks bucket under a path that carries their
digest, `peptide-lm/<name>.<sha12>.json`, stored plain as every track is (the
CDN compresses in transit), and the `lm_run` row (migration 0013) is pointed at
them. The object it replaces is deleted only once the row no longer names it.
`GET /peptide-lm/runs/<name>` then serves the row.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, Optional, Sequence

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.denovo.check_trace import DEFAULT, Refused, validate  # noqa: E402
from pipeline.upload_tracks import Storage  # noqa: E402

BUCKET = "tracks"
FOLDER = "peptide-lm"

_SUMMARY = ("schema", "peptides", "targets", "parameters", "steps", "checkpoints", "loss", "floor")

_UPSERT = """
insert into lm_run (name, bucket, object_path, bytes, sha256, format, provenance, updated_at)
values (%(name)s, %(bucket)s, %(object_path)s, %(bytes)s, %(sha256)s, %(format)s,
        %(provenance)s, now())
on conflict (name) do update set
    bucket = excluded.bucket, object_path = excluded.object_path, bytes = excluded.bytes,
    sha256 = excluded.sha256, format = excluded.format, provenance = excluded.provenance,
    updated_at = now()
"""


def row_for(payload: bytes, summary: Dict, commit: Optional[str] = None,
            made: Optional[str] = None) -> Dict:
    """The `lm_run` row naming `payload`. What made it and when go here, not in
    the payload, so the same run recorded again is the same bytes."""
    digest = hashlib.sha256(payload).hexdigest()
    provenance = {key: summary[key] for key in _SUMMARY}
    provenance.update(summary["provenance"])
    if commit:
        provenance["commit"] = commit
    if made:
        provenance["made"] = made
    return {
        "name": summary["name"],
        "bucket": BUCKET,
        "object_path": "%s/%s.%s.json" % (FOLDER, summary["name"], digest[:12]),
        "bytes": len(payload),
        "sha256": digest,
        "format": "json",
        "provenance": provenance,
    }


def served(service: str, name: str) -> Optional[Dict]:
    """The row the service names for `name`, or None where it names none."""
    url = "%s/peptide-lm/runs/%s" % (service.rstrip("/"), name)
    request = urllib.request.Request(url, headers={"User-Agent": "helix-peek-denovo"})
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            return json.load(response)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        raise


def _commit() -> Optional[str]:
    try:
        return subprocess.run(
            ["git", "-C", str(BACKEND), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True,
        ).stdout.strip() or None
    except (OSError, subprocess.CalledProcessError):
        return None


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Publish the recorded training run.")
    parser.add_argument("--trace", type=Path, default=DEFAULT)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--against", default=None, metavar="SERVICE",
                        help="with --dry-run: whether the service already serves these bytes")
    arguments = parser.parse_args(argv if argv is not None else sys.argv[1:])
    if arguments.against and not arguments.dry_run:
        parser.error("--against compares a dry run with what is served; add --dry-run")

    payload = arguments.trace.read_bytes()
    try:
        summary = validate(payload)
    except Refused as exc:
        print("%s: REFUSED -- %s" % (arguments.trace, exc), file=sys.stderr)
        return 1
    row = row_for(payload, summary, commit=_commit(),
                  made=datetime.date.today().isoformat())
    print("  %-10s %12s B  %s/%s" % (row["name"], format(row["bytes"], ","), BUCKET,
                                     row["object_path"]))
    print("  %d steps, checkpoints %s, loss %.4f -> %.4f"
          % (summary["steps"], ", ".join(str(s) for s in summary["checkpoints"]),
             summary["loss"][0], summary["loss"][1]))
    if arguments.against:
        found = served(arguments.against, row["name"])
        if found is None:
            print("  %s serves no run %r" % (arguments.against, row["name"]))
        elif found.get("sha256") == row["sha256"]:
            print("  %s already serves these bytes" % arguments.against)
        else:
            print("  %s serves other bytes (%s)"
                  % (arguments.against, str(found.get("sha256"))[:12]))
    if arguments.dry_run:
        return 0

    url, key = os.environ.get("SUPABASE_URL"), os.environ.get("SUPABASE_SERVICE_KEY")
    database = os.environ.get("DATABASE_URL")
    if not (url and key and database):
        raise SystemExit("SUPABASE_URL, SUPABASE_SERVICE_KEY and DATABASE_URL must be set.")

    import psycopg
    from psycopg.types.json import Jsonb

    storage = Storage(url, key)
    with psycopg.connect(database) as conn:
        with conn.cursor() as cur:
            old = cur.execute("select object_path from lm_run where name = %s",
                              (row["name"],)).fetchone()
            storage.put(BUCKET, row["object_path"], payload, "application/json")
            cur.execute(_UPSERT, {**row, "provenance": Jsonb(row["provenance"])})
        conn.commit()

    # Only after the row points elsewhere: nothing is deleted while it is still
    # what a client would be told to fetch.
    if old and old[0] != row["object_path"]:
        storage.delete(BUCKET, [old[0]])
        print("  removed %s" % old[0], file=sys.stderr)
    print("published %s: %s" % (row["name"], row["sha256"]), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
