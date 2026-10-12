"""Dry run of the experimental-structure choice. Writes nothing but its report.

    .venv/bin/python scripts/dry_fold.py --out /tmp/folds TTR INSL3     # resolved here, as the worker would
    .venv/bin/python scripts/dry_fold.py --out /tmp/folds --curated     # the twenty, as biological cases
    .venv/bin/python scripts/dry_fold.py --out /tmp/folds --built       # every protein built on demand
    .venv/bin/python scripts/dry_fold.py --out /tmp/folds --built --bake   # and bake what is chosen

For each protein it builds what a structure is chosen for (`fold_job.inputs_of`)
the way the worker will: a gene named here is resolved exactly as the resolver
does (`resolve.resolve`, against live NCBI and UniProt), one of the twenty from
its `targets.py` row and its stored record, and a built protein from its live
row and stored record. Then it runs the choice (`fold_job.choose`) against
live RCSB and PDBe, and writes `report.json` and `report.md` into `--out`.

`--bake` goes on to what the worker will do with a choice: each candidate
whose file passes is baked (`structure/experimental.py`, in the structure
bake's environment, a process a protein) and held to the app's own check, and
the report gains what the bake measured. What is baked stays under `--out`.

Nothing is queued, stored or uploaded. The database is read in read-only
transactions, through `DATABASE_URL` in `--env` (this checkout's `.env`, else
the main checkout's). The only other writes are the record builder's GenBank
cache and the bakers' data directory, in a temporary directory removed when the
run ends. Stored records are read from `--records` (`pipeline/data/`, filled by
`fetch_tracks.py`).
"""

from __future__ import annotations

import argparse
import atexit
import json
import os
import shutil
import sys
import tempfile
import time
import traceback
import urllib.request
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BACKEND))


def _env_file(given: str | None) -> Path:
    for candidate in ([Path(given)] if given else
                      [BACKEND / ".env", BACKEND.parent / "helix-peek-backend" / ".env"]):
        if candidate.exists():
            return candidate
    raise SystemExit("No .env with DATABASE_URL and NCBI_EMAIL; pass --env.")


def main(arguments: list) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("genes", nargs="*")
    parser.add_argument("--curated", action="store_true", help="the twenty")
    parser.add_argument("--built", action="store_true", help="every protein built on demand")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--env", default=None)
    parser.add_argument("--records", type=Path,
                        default=BACKEND.parent / "helix-peek-backend" / "pipeline" / "data",
                        help="where the twenty's stored records are (fetch_tracks.py)")
    parser.add_argument("--limit", type=int, default=None, help="entry files examined at most")
    parser.add_argument("--no-recovery", action="store_true")
    parser.add_argument("--only", action="append", default=None, metavar="SLUG",
                        help="of --curated or --built, only these; repeatable")
    parser.add_argument("--jobs", type=Path, default=None,
                        help="proteins as an earlier run resolved them (its jobs.pickle), "
                             "in place of resolving them again")
    parser.add_argument("--bake", action="store_true",
                        help="bake each choice and hold it to the app's check")
    parser.add_argument("--structure-python", type=Path, default=None,
                        help="the structure bake's Python (pipeline/structure/venv)")
    parser.add_argument("--app", type=Path, default=BACKEND.parent / "helix-peek")
    parser.add_argument("--draws", type=int, default=None, help="candidates baked at most")
    args = parser.parse_args(arguments)

    from dotenv import dotenv_values

    env = dotenv_values(_env_file(args.env))
    # Named before any baker is imported: each reads its directory once, then.
    scratch = Path(tempfile.mkdtemp(prefix="helixpeek-dry-fold-"))
    atexit.register(shutil.rmtree, scratch, True)
    os.environ["NCBI_EMAIL"] = env["NCBI_EMAIL"]
    os.environ["HELIXPEEK_GB_CACHE"] = str(scratch / "genbank")
    os.environ["HELIXPEEK_DATA"] = str(scratch / "data")

    from pipeline import uniprot
    from pipeline.structure import fold_job, fold_report

    import pickle

    # (name, how, target, translation, record, row) or (name, how, error): the
    # record's bytes and the protein's row are what a bake needs beside the rest.
    jobs = []
    if args.jobs:
        jobs = pickle.loads(args.jobs.read_bytes())
    elif args.curated:
        from pipeline import seed_catalog
        from pipeline.targets import TARGETS

        curated = seed_catalog.curated_rows()
        for target in TARGETS:
            record = args.records / target.mock_asset
            if not record.exists():
                jobs.append((target.slug, "curated", f"no stored record at {record}"))
                continue
            stored = record.read_bytes()
            translation = json.loads(stored)["protein"]["translation"]
            jobs.append((target.slug, "curated", target, translation, stored,
                         seed_catalog.protein_row(target, curated[target.slug])))
    if not args.jobs and (args.built or args.genes):
        jobs += _from_database(env, args.genes, args.built)
    args.out.mkdir(parents=True, exist_ok=True)
    if not args.jobs:
        (args.out / "jobs.pickle").write_bytes(pickle.dumps(jobs))

    if args.only:
        jobs = [job for job in jobs if job[1] == "resolved here" or job[0] in args.only]

    report = []
    for job in jobs:
        name, how = job[0], job[1]
        print(f"{name} ({how})", flush=True)
        if len(job) == 3:
            report.append({"name": name, "how": how, "error": job[2]})
            print(f"  {job[2]}", flush=True)
            continue
        target, translation = job[2], job[3]
        started = time.monotonic()
        if args.bake:
            if len(job) < 6:
                raise SystemExit("These jobs were resolved without their records; run without "
                                 "--jobs to resolve them again.")
            said = _baked(args, name, how, target, job[4], job[5])
            said.setdefault("seconds", round(time.monotonic() - started, 1))
            report.append(said)
            baked = said.get("baked") or {}
            print(f"  {said.get('outcome', said.get('error'))} "
                  f"{(said.get('chosen') or {}).get('entry', '')} "
                  + (f"scene {baked['scene_bytes']:,} B, bind {baked.get('bind_ms')} ms "
                     if baked else "") + f"in {said['seconds']} s", flush=True)
            continue
        try:
            entry = uniprot.fetch_entry(target.uniprot)
            from pipeline.structure import pdb_search

            inputs = fold_job.inputs_of(target, translation, entry.body,
                                        pdb_search.pfam_domains(target.uniprot))
            kwargs = {"recovery": not args.no_recovery}
            if args.limit:
                kwargs["limit"] = args.limit
            choice = fold_job.choose(inputs, **kwargs)
            said = fold_report.describe(name, how, target, inputs, choice, entry.release)
        except Exception as exc:  # noqa: BLE001 -- the report says what broke, and goes on
            said = {"name": name, "how": how, "error": f"{type(exc).__name__}: {exc}",
                    "trace": traceback.format_exc()}
        said["seconds"] = round(time.monotonic() - started, 1)
        report.append(said)
        print(f"  {said.get('outcome', said.get('error'))} {said.get('chosen', {}).get('entry', '')} "
              f"in {said['seconds']} s", flush=True)

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "report.json").write_text(json.dumps(report, indent=1, ensure_ascii=False,
                                                     default=str) + "\n")
    (args.out / "report.md").write_text(fold_report.markdown(report))
    print(f"{len(report)} protein(s): {args.out / 'report.md'}")
    return 0


def _baked(args, name: str, how: str, target, record: bytes, protein: dict) -> dict:
    """One protein's choice, baked and checked, by `structure/experimental.py`
    in the structure bake's environment: its part of the report."""
    import subprocess

    from pipeline.structure import experimental

    python = args.structure_python or next(
        (candidate for candidate in (
            BACKEND / "pipeline" / "structure" / "venv" / "bin" / "python",
            BACKEND.parent / "helix-peek-backend" / "pipeline" / "structure" / "venv" / "bin"
            / "python") if candidate.exists()), None)
    if python is None:
        raise SystemExit("No structure bake environment (pipeline/structure/venv); "
                         "pass --structure-python.")
    folder = args.out / "baked" / name
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    (folder / "record.json").write_bytes(record)
    job = {"target": experimental.target_json(target), "record": "record.json",
           "row": experimental.served_row(protein)}
    if how == "curated" and target.structure is not None:
        job["hand_pick"] = target.structure.pdb
    (folder / "job.json").write_text(json.dumps(job, ensure_ascii=False), encoding="utf-8")
    ran = subprocess.run(
        [str(python), "-u", str(BACKEND / "pipeline" / "structure" / "experimental.py"),
         "--job", str(folder / "job.json"), "--out", str(folder), "--app", str(args.app),
         "--describe", name, how]
        + (["--limit", str(args.limit)] if args.limit else [])
        + (["--draws", str(args.draws)] if args.draws else []),
        cwd=str(BACKEND), stdin=subprocess.DEVNULL, capture_output=True, text=True)
    (folder / "stdout.txt").write_text(ran.stdout)
    (folder / "stderr.txt").write_text(ran.stderr)
    described = folder / "report.json"
    if ran.returncode != 0 or not described.exists():
        return {"name": name, "how": how,
                "error": f"the bake exited with status {ran.returncode}: "
                         f"{' '.join(ran.stderr.split())[-400:]}",
                "trace": ran.stderr[-4000:]}
    return json.loads(described.read_text(encoding="utf-8"))


def _from_database(env: dict, genes: list, built: bool) -> list:
    """Jobs for genes resolved here, and for every protein built on demand, from
    the live rows, read in read-only transactions."""
    import psycopg

    from pipeline.resolver import store
    from pipeline.resolver.resolve import Refused, resolve, target_of

    jobs = []
    with psycopg.connect(env["DATABASE_URL"], prepare_threshold=None) as conn:
        with conn.transaction():
            conn.execute("set transaction read only")
            mane = store.mane_release(conn)
            index = []
            for gene in genes:
                hit = conn.execute("select uniprot from protein_index where gene = %s and buildable",
                                   (gene,)).fetchone()
                index.append((gene, store.index_row(conn, hit[0], gene) if hit else None))
            rows = []
            if built:
                for (slug,) in conn.execute(
                        "select slug from protein where catalog_order is null and "
                        "resolver_version > 0 order by resolved_at").fetchall():
                    rows.append((slug, store.protein(conn, slug),
                                 store.ready_track(conn, slug, "record")))
    from pipeline import uniprot

    for gene, row in index:
        if row is None:
            jobs.append((gene, "resolved here", "no buildable index entry"))
            continue
        try:
            resolution = resolve(row, uniprot.fetch_entry(row.uniprot), mane_release=mane)
        except Refused as exc:
            jobs.append((gene, "resolved here", f"the resolver refuses it: {exc}"))
            continue
        translation = json.loads(resolution.record)["protein"]["translation"]
        jobs.append((gene, "resolved here", resolution.target, translation, resolution.record,
                     resolution.protein))
    public = env["SUPABASE_URL"].rstrip("/") + "/storage/v1/object/public"
    for slug, protein, record in rows:
        if record is None:
            jobs.append((slug, "built", "no ready record"))
            continue
        with urllib.request.urlopen(f"{public}/{record['bucket']}/{record['object_path']}",
                                    timeout=120) as response:
            stored = response.read()
        translation = json.loads(stored)["protein"]["translation"]
        jobs.append((slug, "built", target_of(protein), translation, stored, protein))
    return jobs


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
