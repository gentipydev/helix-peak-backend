"""MANE, as the tools read it: the current release's summary, from NCBI.

`scripts/load_protein_index.py` reads it through `current_summary`, to join
every reviewed human UniProt entry to its MANE Select transcripts, so there is
one place that knows where the release is listed and what its summary file is
called. `locus/bake_locus.py` makes the same join with MANE v1.5 as its pinned
snapshot holds it (`locus/inputs.py`), never the current release.
"""

import gzip
import io
import re
from typing import Callable, Dict, Tuple

from app.protein_index import parse_mane

MANE_DIR = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"

# A fetch: the body at a URL, and its response headers.
Get = Callable[[str], Tuple[bytes, Dict[str, str]]]


def current_summary(get: Get) -> Tuple[Dict[str, dict], str]:
    """MANE Select rows keyed by versionless ENST, and the release ("v1.5")."""
    listing, _ = get(MANE_DIR)
    names = re.findall(r"MANE\.GRCh38\.v([0-9.]+)\.summary\.txt\.gz", listing.decode())
    if not names:
        raise SystemExit("No MANE summary in {}".format(MANE_DIR))
    mane_release = "v" + names[0]
    body, _ = get(MANE_DIR + "MANE.GRCh38.{}.summary.txt.gz".format(mane_release))
    mane = parse_mane(io.StringIO(gzip.decompress(body).decode()))
    return mane, mane_release
