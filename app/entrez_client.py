"""Thin wrapper around the Entrez call this service makes.

Text and parsing are separate so the cache has something to store: the record
goes to the database as the text NCBI sent, and is parsed on the way out of
both paths.
"""

import io

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord


def fetch_genbank_text(id: str) -> str:
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
    try:
        return handle.read()
    finally:
        handle.close()


def parse_genbank_text(text: str) -> SeqRecord:
    """Raises ValueError when ``text`` is not a GenBank record.

    NCBI answers an unusable id with plain text and an HTTP 200, so this is
    the check that tells a real record from a rejection.
    """
    return SeqIO.read(io.StringIO(text), "genbank")


def fetch_genbank_record(id: str) -> SeqRecord:
    """Fetch and parse in one step, with no cache involved."""
    return parse_genbank_text(fetch_genbank_text(id))
