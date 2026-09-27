"""MANE, as the tools read it: the current release's summary, from NCBI.

Two tools read it: `scripts/load_protein_index.py`, which joins every reviewed
human UniProt entry to its MANE Select transcripts, and `locus/bake_locus.py`,
which places each of the twenty on its chromosome by the same join. Both fetch
through `current_summary`, so there is one place that knows where the release
is listed and what its summary file is called.
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
