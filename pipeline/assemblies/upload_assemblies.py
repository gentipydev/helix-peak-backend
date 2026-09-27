"""Upload an assembly's morph, and the models beside it, and record where they landed.

    .venv/Scripts/python pipeline/assemblies/upload_assemblies.py --dry-run
    DATABASE_URL=... SUPABASE_URL=... SUPABASE_SERVICE_KEY=... \
        .venv/Scripts/python pipeline/assemblies/upload_assemblies.py

What `upload_tracks.py` does for a protein, for an assembly: validated first
(the whole of `check_assembly.py`), stored plain at a path carrying the digest,
`morph/<slug>.<sha12>.json` in `tracks` and each state's `.glb` in `models`
beside it, named in the row's provenance. The row is `assembly_track`
(migrations/0006_assemblies.sql, which must be applied, with the `assembly`
row seeded first). Upload only when asked, with `--dry-run` first.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.assemblies.assemblies import ASSEMBLIES  # noqa: E402
from pipeline.assemblies.bake_assembly import model_asset, morph_asset  # noqa: E402
from pipeline.assemblies.check_assembly import problems_of  # noqa: E402
from pipeline.paths import DATA  # noqa: E402
from pipeline.upload_tracks import Storage  # noqa: E402

_UPSERT = """
insert into assembly_track
    (slug, kind, state, reason, bucket, object_path, bytes, sha256,
     content_encoding, format, provenance, updated_at)
values (%(slug)s, 'morph', 'ready', null, 'tracks', %(object_path)s, %(bytes)s,
        %(sha256)s, null, 'json', %(provenance)s, now())
on conflict (slug, kind) do update set
    state = 'ready', reason = null, bucket = excluded.bucket,
    object_path = excluded.object_path, bytes = excluded.bytes,
    sha256 = excluded.sha256, content_encoding = null,
    format = excluded.format, provenance = excluded.provenance, updated_at = now()
"""


def planned() -> list[dict]:
    rows = []
    for assembly in ASSEMBLIES:
        path = DATA / morph_asset(assembly)
        if not path.exists():
            continue
        found = problems_of(assembly)
        if found:
            print(f"  {assembly.slug}: REFUSED -- {'; '.join(found)}", file=sys.stderr)
            continue
        payload = path.read_bytes()
        track = json.loads(payload)
        digest = hashlib.sha256(payload).hexdigest()
        models = []
        for state in assembly.states:
            glb = (DATA / model_asset(assembly, state.name)).read_bytes()
            sha = hashlib.sha256(glb).hexdigest()
            models.append((state.name, f"morph/{assembly.slug}.{state.name}.{sha[:12]}.glb", glb, sha))
        rows.append({
            "slug": assembly.slug,
            "object_path": f"morph/{assembly.slug}.{digest[:12]}.json",
            "bytes": len(payload),
            "sha256": digest,
            "payload": payload,
            "models": models,
            "provenance": {
                "states": [{"name": s["name"], "pdb": s["pdb"], "ligand": s["ligand"]}
                           for s in track["states"]],
                "superposition": {k: v for k, v in track["frame"]["superposition"].items()
                                  if k not in ("rotation", "shift_angstrom")},
                "models": {name: {"path": p, "sha256": sha, "bytes": len(glb)}
                           for name, p, glb, sha in models},
                "schema_version": track["schema_version"],
                "built_by": track["built_by"],
            },
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    rows = planned()
    for row in rows:
        print(f"  {row['slug']:<16} {row['bytes']:>10,} B  tracks/{row['object_path']}")
        for _, path, glb, _ in row["models"]:
            print(f"  {row['slug']:<16} {len(glb):>10,} B  models/{path}")
    print(f"\n{len(rows)} assembly morph(s)", file=sys.stderr)
    if args.dry_run:
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
            for row in rows:
                storage.put("tracks", row["object_path"], row["payload"], "application/json")
                for _, path, glb, _ in row["models"]:
                    storage.put("models", path, glb, "model/gltf-binary")
                cur.execute(_UPSERT, {**{k: row[k] for k in ("slug", "object_path", "bytes", "sha256")},
                                      "provenance": Jsonb(row["provenance"])})
                print(f"  uploaded {row['slug']}", file=sys.stderr)
        conn.commit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
