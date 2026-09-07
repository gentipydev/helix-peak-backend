"""Thin wrapper around the single Entrez call this service makes."""

from Bio import Entrez


def fetch_genbank(id: str) -> str:
    """Return the raw GenBank flat file for ``id`` from NCBI's nucleotide db.

    Returns the upstream text verbatim, including NCBI's own error bodies --
    classifying those is the router's job.
    """
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
    content = handle.read()
    handle.close()
    return content
