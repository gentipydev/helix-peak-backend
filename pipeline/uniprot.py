"""UniProt, as the bakers read it: one entry at a time, from rest.uniprot.org.

Moved here from `mock/build_gene_record.py`, which checks each record's protein
against the entry's canonical sequence, so that a second baker reading the same
entries shares the fetch rather than copying it.
"""

from __future__ import annotations

import json
import urllib.request


def canonical_sequence(uniprot: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot}.json"
    with urllib.request.urlopen(url, timeout=60) as response:
        return json.load(response)["sequence"]["value"]
