"""Ask a running service for suggestions, and check the ones that must hold.

    .venv/bin/python scripts/check_suggest.py https://helix-peak-backend.onrender.com
    .venv/bin/python scripts/check_suggest.py http://localhost:8000

The unit tests fake the database, so this is where the ranking SQL meets the
real index: each query below states what a biologist typing it would expect
to see, and the run fails if any of them does not hold. It then times a warm
round of queries end to end.
"""

import json
import statistics
import sys
import time
import urllib.parse
import urllib.request
from typing import Callable, List, Tuple


def _ask(base: str, query: str, limit: int = 12) -> Tuple[dict, float]:
    url = "{}/proteins/suggest?{}".format(
        base.rstrip("/"), urllib.parse.urlencode({"q": query, "limit": limit}))
    started = time.perf_counter()
    with urllib.request.urlopen(url, timeout=90) as response:
        body = json.loads(response.read())
    return body, time.perf_counter() - started


def _genes(body: dict) -> List[str]:
    return [s["gene"] or s["uniprot"] for s in body["suggestions"]]


def first(gene: str, status: str = None) -> Callable[[dict], bool]:
    def check(body):
        top = body["suggestions"][:1]
        return bool(top) and top[0]["gene"] == gene and (status is None or top[0]["status"] == status)
    check.__doc__ = "{} first{}".format(gene, ", " + status if status else "")
    return check


def within(genes: List[str], n: int) -> Callable[[dict], bool]:
    def check(body):
        return set(genes) <= set(_genes(body)[:n])
    check.__doc__ = "{} in the top {}".format(", ".join(genes), n)
    return check


def status_of(gene: str, status: str) -> Callable[[dict], bool]:
    def check(body):
        return any(s["gene"] == gene and s["status"] == status for s in body["suggestions"])
    check.__doc__ = "{} is {}".format(gene, status)
    return check


GOLDEN = [
    ("insulin", first("INS", "listed")),
    ("INS", first("INS", "listed")),
    ("P01308", first("INS")),
    ("insul", first("INS")),
    ("insuln", within(["INS"], 3)),
    ("p53", first("TP53", "listed")),
    ("P04637", first("TP53")),
    ("hemoglobin", first("HBB", "listed")),
    ("hemoglobin", within(["HBA1", "HBA2"], 10)),
    ("hba", within(["HBA1", "HBA2"], 5)),
    ("brca", within(["BRCA1", "BRCA2"], 5)),
    # MANE's titin is isoform Q8WZ42-12, so it is unavailable before the exon
    # budget (R2.4) is ever asked.
    ("ttn", first("TTN", "unavailable")),
    ("glucagon-like peptide 1", within(["GCG"], 3)),
    ("gpr89b", first("GPR89B", "buildable")),
    ("TP53BP2", status_of("TP53BP2", "unavailable")),
    ("α-actinin", within(["ACTN1"], 5)),
]

TIMED = ["ins", "insu", "insul", "p5", "p53", "hem", "hemo", "kin", "kinase", "brc",
         "a", "zz", "receptor", "insuln", "tp53bp", "cftr", "dystro", "amyl", "glp",
         "P04637"]


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    base = sys.argv[1]
    _ask(base, "warm up")

    failures = 0
    for query, check in GOLDEN:
        body, _ = _ask(base, query)
        held = check(body)
        failures += not held
        print("{}  {:<28} {:<34} {}".format(
            "ok  " if held else "FAIL", repr(query), check.__doc__,
            " ".join(_genes(body)[:6])))

    timings = []
    for _ in range(3):
        for query in TIMED:
            timings.append(_ask(base, query)[1])
    timings.sort()
    p95 = timings[int(len(timings) * 0.95) - 1]
    print("\n{} warm queries: median {:.0f} ms, p95 {:.0f} ms, slowest {:.0f} ms".format(
        len(timings), statistics.median(timings) * 1000, p95 * 1000, timings[-1] * 1000))
    print("release: {}".format(_ask(base, "ins")[0]["release"]))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
