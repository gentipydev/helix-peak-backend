from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

HERE = Path(__file__).resolve().parent
if str(HERE.parents[1]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[1]))

from pipeline.denovo.sequences import clean, is_modelable

DATA = HERE / "data"
CORPUS = DATA / "peptides.txt"
PROVENANCE = DATA / "provenance.json"

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
USER_AGENT = "helix-peek-denovo (offline learning project)"

NOT_TOXIN = "NOT (keyword:KW-0800)"
WHOLE_QUERY = (
    "(reviewed:true) AND (length:[10 TO 50]) AND (fragment:false) " + NOT_TOXIN
)
PRODUCT_QUERY = "(reviewed:true) AND (ft_peptide:*) " + NOT_TOXIN
PRODUCT_FIELDS = ("accession", "sequence", "ft_peptide", "ft_propep")

MATURE = "Peptide"
PROPEPTIDE = "Propeptide"


def fetch(query: str, fmt: str, fields: Sequence[str] = (), timeout: int = 600) -> Tuple[str, str]:
    params = {"query": query, "format": fmt}
    if fields:
        params["fields"] = ",".join(fields)
    url = STREAM + "?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        release = response.headers.get("X-UniProt-Release", "")
        return response.read().decode("utf-8"), release


def parse_fasta(text: str) -> List[str]:
    sequences, current = [], []
    for line in text.splitlines():
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
        elif line.strip():
            current.append(line.strip())
    if current:
        sequences.append("".join(current))
    return sequences


def _exact(place: Dict) -> bool:
    modifier = place.get("modifier")
    return modifier in (None, "EXACT") and isinstance(place.get("value"), int)


def products_from_entry(entry: Dict, kinds: Sequence[str] = (MATURE, PROPEPTIDE)) -> List[Tuple[str, str]]:
    sequence = (entry.get("sequence") or {}).get("value") or ""
    found = []
    for feature in entry.get("features", []):
        kind = feature.get("type")
        if kind not in kinds:
            continue
        location = feature.get("location") or {}
        start, end = location.get("start") or {}, location.get("end") or {}
        if not (_exact(start) and _exact(end)):
            continue
        first, last = start["value"], end["value"]
        if first < 1 or last > len(sequence) or last < first:
            continue
        found.append((kind, sequence[first - 1 : last]))
    return found


def products_from_json(text: str, kinds: Sequence[str] = (MATURE, PROPEPTIDE)) -> List[Tuple[str, str]]:
    entries = json.loads(text).get("results", [])
    found = []
    for entry in entries:
        found.extend(products_from_entry(entry, kinds=kinds))
    return found


def write_corpus(sequences: Sequence[str], path: Path = CORPUS) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "\n".join(sequences) + "\n"
    path.write_text(body)
    return hashlib.sha256(body.encode()).hexdigest()


def main(argv: Sequence[str] = ()) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=CORPUS)
    parser.add_argument(
        "--no-propeptides",
        action="store_true",
        help="keep only PEPTIDE products, not the propeptides cut off beside them",
    )
    args = parser.parse_args(argv or sys.argv[1:])

    kinds = (MATURE,) if args.no_propeptides else (MATURE, PROPEPTIDE)

    print("whole entries  : %s" % WHOLE_QUERY)
    fasta, release = fetch(WHOLE_QUERY, "fasta")
    whole = parse_fasta(fasta)
    print("  fetched %d entries (UniProt release %s)" % (len(whole), release or "?"))

    print("mature products: %s" % PRODUCT_QUERY)
    raw, _ = fetch(PRODUCT_QUERY, "json", PRODUCT_FIELDS)
    products = products_from_json(raw, kinds=kinds)
    by_kind = {kind: sum(1 for k, _ in products if k == kind) for kind in kinds}
    print("  fetched %d products (%s)" % (
        len(products),
        ", ".join("%s %d" % (k, n) for k, n in by_kind.items()),
    ))

    product_sequences = [s for _, s in products]
    in_range = [s for s in product_sequences if is_modelable(s)]
    print("  %d of them are 10-50 canonical residues" % len(in_range))

    corpus = clean(list(whole) + product_sequences)
    print(
        "\ncorpus: %d whole + %d products -> %d after range, alphabet and exact duplicates"
        % (len(whole), len(product_sequences), len(corpus))
    )

    digest = write_corpus(corpus, args.out)
    PROVENANCE.write_text(json.dumps({
        "fetched": date.today().isoformat(),
        "uniprot_release": release,
        "source": "https://rest.uniprot.org/uniprotkb",
        "queries": {"whole_entries": WHOLE_QUERY, "mature_products": PRODUCT_QUERY},
        "product_kinds": list(kinds),
        "counts": {
            "whole_entries": len(whole),
            "products": len(product_sequences),
            "products_by_kind": by_kind,
            "products_in_range": len(in_range),
            "corpus": len(corpus),
        },
        "corpus_sha256": digest,
    }, indent=2) + "\n")
    print("wrote %s (sha256 %s)" % (args.out, digest[:16]))
    print("wrote %s" % PROVENANCE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
