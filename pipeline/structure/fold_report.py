"""What a search for a protein's experimental structure found, as a report.

`describe` is one protein's part of the dry run's report (`scripts/dry_fold.py`),
`markdown` the whole report as a page, and `selection` what a fold's audit
keeps of the search that chose it. All three read a `fold_job.Choice` and write
nothing.
"""

from __future__ import annotations

import hashlib
import json

from pipeline.structure import normalize, pdb_policy, pdb_words
from pipeline.structure.fold_job import candidate_id


def describe(name: str, how: str, target, inputs, choice, release) -> dict:
    """One protein's part of the report: what was searched, ranked, examined and
    chosen, the words the fold page would carry, and what its bake measured
    where it was baked."""
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
                      "passed": e.passed,
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
        built = winner.built
        if built is not None:
            said["baked"] = {**built.measured,
                             "representation": built.provenance["representation"],
                             "sampling": built.provenance["sampling"],
                             "secondary_structure": built.provenance["secondary_structure"],
                             "residues": built.provenance["residues"],
                             "frame": built.provenance["frame"]}
    return said


def selection(inputs, choice, release) -> dict:
    """What an audit keeps of a search: what was searched and when, the rule,
    every rejection by its reason, each candidate examined and how it fared,
    and why the winner beat the next. Enough to say why a protein has the
    fold it has, or none."""
    ranking = choice.ranking
    said = {
        "policy": pdb_policy.POLICY_VERSION,
        "inputs": inputs_digest(inputs),
        "uniprot_release": release,
        "outcome": choice.outcome,
        "search": choice.discovery.snapshot,
        "recovered": list(choice.recovered),
        "recovery_refusals": choice.recovery_refusals,
        "candidates": len(ranking.candidates),
        "rejected": ranking.reasons(),
        "rejection_examples": _examples(ranking.rejections),
        "examined": [{"candidate": candidate_id(e.candidate), "passed": e.passed,
                      "reason": None if e.refused is None else e.refused.reason,
                      "detail": None if e.refused is None else e.refused.detail}
                     for e in choice.examinations],
        "top": [_row(i + 1, c) for i, c in enumerate(_distinct(ranking.candidates)[:8])],
    }
    # What a resumed search skips: every candidate turned down, in this attempt
    # and those it went on from.
    said["refused"] = choice.refused
    said["resumed"] = len(choice.resumed)
    winner = choice.winner
    if winner is not None:
        others = [c for c in ranking.candidates if c.entry != winner.candidate.entry]
        said["candidate"] = candidate_id(winner.candidate)
        said["chosen"] = {**_row(ranking.candidates.index(winner.candidate) + 1, winner.candidate),
                          "why": pdb_policy.why(winner.candidate, others[0] if others else None)}
    return said


def inputs_digest(inputs) -> str:
    """What a choice was made for, as one digest: the protein's letters, its
    products, what is removed, its domains and its disulfides. A search that
    found nothing suitable is remembered for these inputs under this policy,
    and for no others."""
    said = {"accession": inputs.accession, "translation": inputs.translation,
            "products": [[p.label, p.start, p.end] for p in inputs.products],
            "removed": [list(r) for r in inputs.removed],
            "domains": [[d.kind, d.name, d.start, d.end, d.source] for d in inputs.domains],
            "disulfides": [list(pair) for pair in inputs.disulfides]}
    return hashlib.sha256(json.dumps(said, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


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
             "| protein | outcome | entry | method | scope | placed / spanned | sentence | "
             "scene | bake | bind |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for p in report:
        if "error" in p:
            lines.append(f"| {p['name']} | error | | | | | {p['error'][:80]} | | | |")
            continue
        chosen, words = p.get("chosen") or {}, (p.get("words") or {}).get("chrome") or {}
        examined = chosen.get("examined") or {}
        placed = sum(c["placed"] for c in examined.get("chains", []))
        listed = sum(c["listed"] for c in examined.get("chains", []))
        baked = p.get("baked") or {}
        lines.append(f"| {p['name']} | {p['outcome']} | {chosen.get('entry', '')} | "
                     f"{chosen.get('method', '')} {chosen.get('resolution') or ''} | "
                     f"{chosen.get('scope', '')} | {f'{placed}/{listed}' if chosen else ''} | "
                     f"{words.get('sentence', '')} | "
                     + (f"{baked['scene_bytes']:,} B | {baked['bake_seconds']} s | "
                        f"{baked.get('bind_ms', '')} ms |" if baked else "| | |"))
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
        if p.get("baked"):
            b = p["baked"]
            lines += ["", f"Baked as a {b['representation']} at sampling {b['sampling']} "
                          f"({b['secondary_structure']}): scene {b['scene_bytes']:,} B, "
                          f"{b.get('vertices', 0):,} vertices, {b.get('rings', 0):,} rings; "
                          f"residues {b['residues']}; bake {b['bake_seconds']} s, the app's "
                          f"check {b['check_seconds']} s (bind {b.get('bind_ms')} ms, 480 frames "
                          f"{b.get('frames_ms')} ms, slowest {b.get('slowest_frame_us')} us, "
                          f"largest step {b.get('largest_step', 0):.2f} A)."]
        if p["rejection_examples"]:
            lines += ["", "<details><summary>rejections, examples</summary>", ""]
            for reason, examples in p["rejection_examples"].items():
                lines.append(f"- {reason}: " + "; ".join(examples))
            lines.append("</details>")
    return "\n".join(lines) + "\n"
