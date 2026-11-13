"""Thin wrapper around the Entrez call this service makes.

Text and parsing are separate so the cache has something to store: the record
goes to the database as the text NCBI sent, and is parsed on the way out of
both paths.
"""

import io
from typing import Optional

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord


class RecordTooLarge(Exception):
    """NCBI's record is longer than the service reads. Not a ValueError: that
    is how a record that is not one is told (`parse_genbank_text`)."""

    def __init__(self, id: str, limit: int):
        super().__init__("Record {!r} is longer than {:,} characters.".format(id, limit))
        self.id = id
        self.limit = limit


def fetch_genbank_text(id: str, limit: Optional[int] = None) -> str:
    """The record as NCBI sends it. With `limit`, never more than that is read:
    one character past it raises RecordTooLarge, and the rest stays unread, so
    a whole chromosome (hundreds of MB) never reaches this process's memory."""
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
    try:
        if limit is None:
            return handle.read()
        text = handle.read(limit + 1)
        if len(text) > limit:
            raise RecordTooLarge(id, limit)
        return text
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
