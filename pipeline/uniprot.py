"""UniProt, as the bakers read it: one entry at a time, from rest.uniprot.org.

Two bakers read the same entries. `mock/build_gene_record.py` checks each
record's protein against the entry's canonical sequence, and
`trafficking/bake_trafficking.py` reads its membrane topology, GPI anchor and
subcellular location. Both fetch through `fetch_entry`, so there is one place
that knows the URL, and one that knows which release answered.
"""

from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass
from datetime import date, datetime, timezone

ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{}.json"

_MONTHS = ("January", "February", "March", "April", "May", "June", "July",
           "August", "September", "October", "November", "December")


@dataclass(frozen=True)
class Entry:
    """One UniProtKB entry as the service answered it, and when."""

    accession: str
    body: dict
    # The release that served it, as its `X-UniProt-Release` header names it
    # ("2026_03"), and that release's date from `X-UniProt-Release-Date`, in
    # ISO form. None where the service left a header off.
    release: str | None
    release_date: str | None
    # The day the entry was fetched, in UTC, ISO form.
    retrieved: str


def release_day(header: str | None) -> str | None:
    """`X-UniProt-Release-Date` ("02-September-2026") as an ISO date.

    Read against a fixed list of English month names, not `strptime`'s `%B`,
    which follows the locale the bake happens to run in.
    """
    if not header:
        return None
    day, month, year = header.strip().split("-")
    return date(int(year), _MONTHS.index(month) + 1, int(day)).isoformat()


def fetch_entry(uniprot: str) -> Entry:
    with urllib.request.urlopen(ENTRY_URL.format(uniprot), timeout=60) as response:
        body = json.load(response)
        headers = response.headers
    return Entry(
        accession=uniprot,
        body=body,
        release=headers.get("X-UniProt-Release"),
        release_date=release_day(headers.get("X-UniProt-Release-Date")),
        retrieved=datetime.now(timezone.utc).date().isoformat(),
    )


def canonical_sequence(uniprot: str) -> str:
    return fetch_entry(uniprot).body["sequence"]["value"]
