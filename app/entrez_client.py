"""Thin wrapper around the Entrez call this service makes."""

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord


def fetch_genbank_record(id: str) -> SeqRecord:
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
    try:
        return SeqIO.read(handle, "genbank")
    finally:
        handle.close()
