"""Dry run of the experimental-structure choice. Writes nothing but its report.

    .venv/bin/python scripts/dry_fold.py --out /tmp/folds TTR INSL3     # resolved here, as the worker would
    .venv/bin/python scripts/dry_fold.py --out /tmp/folds --curated     # the twenty, as biological cases
    .venv/bin/python scripts/dry_fold.py --out /tmp/folds --built       # every protein built on demand

For each protein it builds what a structure is chosen for (`fold_job.inputs_of`)
the way the worker will: a gene named here is resolved exactly as the resolver
does (`resolve.resolve`, against live NCBI and UniProt), one of the twenty from
its `targets.py` row and its stored record, and a built protein from its live
row and stored record. Then it runs the choice (`fold_job.choose`) against
live RCSB and PDBe, and writes `report.json` and `report.md` into `--out`.

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
    from pipeline.structure import fold_job

    import pickle

    jobs = []           # (name, how, target, translation) or (name, how, error)
    if args.jobs:
        jobs = pickle.loads(args.jobs.read_bytes())
    elif args.curated:
        from pipeline.targets import TARGETS

        for target in TARGETS:
            record = args.records / target.mock_asset
            if not record.exists():
                jobs.append((target.slug, "curated", f"no stored record at {record}"))
                continue
            translation = json.loads(record.read_bytes())["protein"]["translation"]
            jobs.append((target.slug, "curated", target, translation))
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
        try:
            entry = uniprot.fetch_entry(target.uniprot)
            from pipeline.structure import pdb_search

            inputs = fold_job.inputs_of(target, translation, entry.body,
                                        pdb_search.pfam_domains(target.uniprot))
            kwargs = {"recovery": not args.no_recovery}
            if args.limit:
                kwargs["limit"] = args.limit
            choice = fold_job.choose(inputs, **kwargs)
            said = describe(name, how, target, inputs, choice, entry.release)
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
    (args.out / "report.md").write_text(markdown(report))
    print(f"{len(report)} protein(s): {args.out / 'report.md'}")
    return 0


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
        jobs.append((gene, "resolved here", resolution.target, translation))
    public = env["SUPABASE_URL"].rstrip("/") + "/storage/v1/object/public"
    for slug, protein, record in rows:
        if record is None:
            jobs.append((slug, "built", "no ready record"))
            continue
        with urllib.request.urlopen(f"{public}/{record['bucket']}/{record['object_path']}",
                                    timeout=120) as response:
            translation = json.loads(response.read())["protein"]["translation"]
        jobs.append((slug, "built", target_of(protein), translation))
    return jobs


def describe(name: str, how: str, target, inputs, choice, release) -> dict:
    """One protein's part of the report: what was searched, ranked, examined and
    chosen, and the words the fold page would carry."""
    from pipeline.structure import normalize, pdb_policy, pdb_words
    from pipeline.structure.fold_job import candidate_id

    ranking = choice.ranking
    said = {
        "name": name, "how": how, "gene": inputs.gene, "accession": inputs.accession,
        "display": inputs.display, "residues": len(inputs.translation), "uniprot_release": release,
        "products": [{"label": p.label, "start": p.start, "end": p.end} for p in inputs.products],
        "domains": [{"kind": d.kind, "name": d.name, "start": d.start, "end": d.end,
                     "source": d.source} for d in inputs.domains],
        "outcome": choice.outcome,
        "search": choice.discovery.snapshot, "recovered": list(choice.recovered),
        "recovery_refusals": choice.recovery_refusals,
        "candidates": len(ranking.candidates), "rejected": ranking.reasons(),
        "rejection_examples": _examples(ranking.rejections),
        "top": [_row(i + 1, c) for i, c in enumerate(_distinct(ranking.candidates)[:8])],
        "examined": [{"candidate": candidate_id(e.candidate),
                      "passed": e.examined is not None,
                      "reason": None if e.refused is None else e.refused.reason,
                      "detail": None if e.refused is None else e.refused.detail}
                     for e in choice.examinations],
        "timings": choice.seconds,
    }
    said["otherwise"] = alternatives(inputs, choice)
    if how == "curated" and target.structure is not None:
        said["historical"] = _historical(target.structure.pdb, ranking)
    winner = choice.winner
    if winner is not None:
        candidate, examined = winner.candidate, winner.examined
        others = [c for c in ranking.candidates if c.entry != candidate.entry]
        said["chosen"] = {**_row(ranking.candidates.index(candidate) + 1, candidate),
                          "why": pdb_policy.why(candidate, others[0] if others else None),
                          "examined": normalize.summary(examined)}
        lone = candidate.picks[0].entity if len(candidate.picks) == 1 else None
        said["words"] = pdb_words.words(inputs, candidate, examined,
                                        choice.entries[candidate.entry],
                                        choice.entities.get(lone) if lone else None)
    return said


def alternatives(inputs, choice) -> dict:
    """What would rank first under each setting the policy could have had
    instead, for the report's reader to weigh. None of these is examined
    against its file: each is the head of a ranking, not a choice.

    - `quality_first`: the experiment's quality before the residues that are
      not the gene's (a sharper mutant before a coarser wild type).
    - `measured_first`: a region measured residue by residue before any that
      is not, whatever its method (an NMR model before a crystal nobody
      deposited the data of).
    - `map_fit_30`: an EM chain refused where 30% of the residues drawn fit
      their map poorly, not 50%.
    """
    from pipeline.structure import pdb_policy

    ranked = choice.ranking.candidates

    def head(key):
        return _row(1, min(ranked, key=key)) if ranked else None

    # The key's parts: scope, chains together, coverage, precursor, fidelity,
    # no resolution, unmeasured, resolution, fit, entry, chains.
    out = {
        "quality_first": head(lambda c: c.key()[:4] + c.key()[5:9] + (c.key()[4],) + c.key()[9:]),
        "measured_first": head(lambda c: c.key()[:5] + c.key()[6:]),
    }
    held = pdb_policy.Q_SHARE
    pdb_policy.Q_SHARE = 0.30
    try:
        strict = pdb_policy.rank(inputs, list(choice.entities.values())).candidates
    finally:
        pdb_policy.Q_SHARE = held
    out["map_fit_30"] = _row(1, strict[0]) if strict else None
    return out


def _distinct(candidates):
    """The best candidate of each entry, in rank order."""
    seen, out = set(), []
    for c in candidates:
        if c.entry not in seen:
            seen.add(c.entry)
            out.append(c)
    return out


def _row(rank: int, c) -> dict:
    return {"rank": rank, "entry": c.entry, "scope": c.scope,
            "chains": [f"{p.asym}:{p.first}-{p.last} ({p.placed})" for p in c.picks],
            "domains": list(c.domains), "method": c.method, "resolution": c.resolution,
            "useful": c.useful, "coverage": round(c.coverage, 3), "share": round(c.share, 3),
            "shortfall": c.shortfall,
            "precursor": c.precursor, "mismatches": c.mismatches, "allowed": c.allowed,
            "tier": c.tier, "support": None if c.support is None else round(c.support, 3),
            "measures": c.measures}


def _examples(rejections) -> dict:
    out: dict = {}
    for r in rejections:
        held = out.setdefault(r.reason, [])
        if len(held) < 4:
            held.append(f"{r.entity}{'/' + r.chain if r.chain else ''}: {r.detail}")
    return out


def _historical(pdb: str, ranking) -> dict:
    """Where the entry chosen by hand for one of the twenty stands under the rule."""
    for i, c in enumerate(ranking.candidates):
        if c.entry.upper() == pdb.upper():
            return {"entry": pdb, "rank": i + 1, "of": len(ranking.candidates)}
    reasons = sorted({r.reason + ": " + r.detail for r in ranking.rejections
                      if r.entity.upper().startswith(pdb.upper() + "_")})
    return {"entry": pdb, "rank": None, "rejected": reasons[:4]}


def markdown(report: list) -> str:
    lines = ["# Experimental folds: dry run", "",
             "| protein | outcome | entry | method | scope | placed / spanned | sentence |",
             "|---|---|---|---|---|---|---|"]
    for p in report:
        if "error" in p:
            lines.append(f"| {p['name']} | error | | | | | {p['error'][:80]} |")
            continue
        chosen, words = p.get("chosen") or {}, (p.get("words") or {}).get("chrome") or {}
        examined = chosen.get("examined") or {}
        placed = sum(c["placed"] for c in examined.get("chains", []))
        listed = sum(c["listed"] for c in examined.get("chains", []))
        lines.append(f"| {p['name']} | {p['outcome']} | {chosen.get('entry', '')} | "
                     f"{chosen.get('method', '')} {chosen.get('resolution') or ''} | "
                     f"{chosen.get('scope', '')} | {f'{placed}/{listed}' if chosen else ''} | "
                     f"{words.get('sentence', '')} |")
    for p in report:
        lines += ["", f"## {p['name']} ({p['how']})", ""]
        if "error" in p:
            lines += ["```", p.get("trace") or p["error"], "```"]
            continue
        lines.append(f"{p['display']}, UniProt {p['accession']} ({p['uniprot_release']}), "
                     f"{p['residues']:,} residues. Products: "
                     + "; ".join(f"{x['label']} {x['start']}–{x['end']}" for x in p["products"])
                     + f". Domains annotated: {len(p['domains'])}.")
        s = p["search"]
        lines.append(f"Search: {s['entities']} entities (release {s['release']}), recovered "
                     f"{len(p['recovered'])}; {p['candidates']} eligible; rejected "
                     + (", ".join(f"{k} {v}" for k, v in p["rejected"].items()) or "none")
                     + f". Timings {p['timings']}, {p['seconds']} s.")
        if p.get("historical"):
            h = p["historical"]
            lines.append(f"Chosen by hand: {h['entry']}: "
                         + (f"rank {h['rank']} of {h['of']}" if h.get("rank")
                            else f"not eligible ({'; '.join(h.get('rejected') or ['not found'])})"))
        lines += ["", "| # | entry | scope | chains | method | Å | placed | short | precursor | "
                      "mism. | tier | support |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in p["top"]:
            lines.append(f"| {r['rank']} | {r['entry']} | {r['scope']} | {' '.join(r['chains'])} | "
                         f"{r['method']} | {r['resolution'] or ''} | {r['useful']} | "
                         f"{r['shortfall']} | {'yes' if r['precursor'] else ''} | "
                         f"{r['mismatches']}/{r['allowed']} | {r['tier']} | {r['support']} |")
        refused = [e for e in p["examined"] if not e["passed"]]
        if refused:
            lines += ["", "Refused on examination: " + "; ".join(
                f"{e['candidate']} ({e['reason']}: {e['detail']})" for e in refused)]
        if p.get("chosen"):
            c, w = p["chosen"], p["words"]
            e = c["examined"]
            lines += ["", f"**Chosen: {c['entry']}** {', '.join(c['chains'])}: {c['why']}. "
                          f"Bridges {[(b['a'], b['b']) for b in e['bridges']]}"
                          f"{', dropped ' + str(e['bridges_dropped']) if e['bridges_dropped'] else ''}; "
                          f"helices {e['helix']}, strands {e['strand']} "
                          f"({e['secondary_method'] or 'PyMOL to assign'}); modified "
                          f"{len(e['modified'])}; mismatches {e['mismatches']}.",
                      "", f"- label: {w['chrome']['label']}; count {w['chrome']['count']}; "
                          f"modelled {w['chrome']['modelled']}; tints "
                          + ", ".join(f"{x['node']}={x['tint']}" for x in w["chains"]),
                      f"- sentence: {w['chrome']['sentence']}",
                      f"- semantics: {w['chrome']['semantics']}"]
        if p["rejection_examples"]:
            lines += ["", "<details><summary>rejections, examples</summary>", ""]
            for reason, examples in p["rejection_examples"].items():
                lines.append(f"- {reason}: " + "; ".join(examples))
            lines.append("</details>")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
