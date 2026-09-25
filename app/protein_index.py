"""The protein index: how its rows and search terms are made from UniProt and
MANE, and how a query is normalised to meet them.

Pure functions only. `scripts/load_protein_index.py` downloads the sources and
writes the tables; `app/suggest.py` reads them. Both normalise through
``normalize`` here, which is the whole contract between a stored term and a
typed query: if the two ever normalised differently, a prefix that should match
would fall between them without an error anywhere.
"""

import re
import unicodedata
from typing import Dict, Iterable, List, Optional, Set, Tuple

# Spelled out before the ASCII fold, which would otherwise drop them: an
# "α-actinin" has to be findable as "alpha actinin".
_GREEK = {
    "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon",
    "ζ": "zeta", "η": "eta", "θ": "theta", "ι": "iota", "κ": "kappa",
    "λ": "lambda", "μ": "mu", "ν": "nu", "ξ": "xi", "ο": "omicron", "π": "pi",
    "ρ": "rho", "σ": "sigma", "ς": "sigma", "τ": "tau", "υ": "upsilon",
    "φ": "phi", "χ": "chi", "ψ": "psi", "ω": "omega",
}

_NOT_ALNUM = re.compile(r"[^a-z0-9]+")

# Kinds, in the order a search ranks them after exact matches. A word is its
# own kind because only a whole term may match exactly.
ACCESSION, SYMBOL, SYNONYM, NAME, OTHER_NAME, WORD = 0, 1, 2, 3, 4, 5

# Words that would make a name's word terms match half the index. A name is
# still a term whole, so "protein kinase C" stays findable by its prefix.
_STOPWORDS = frozenset({
    "and", "for", "from", "like", "the", "with", "protein", "probable",
    "putative", "uncharacterized", "containing", "domain", "family",
    "isoform", "subunit", "type", "member", "chain", "homolog",
})

# A word shorter than this is only matched as part of a whole name.
_MIN_WORD = 3


def normalize(text: str) -> str:
    """Lower-case ASCII letters and digits, separated by single spaces."""
    lowered = text.lower()
    for letter, spelled in _GREEK.items():
        if letter in lowered:
            lowered = lowered.replace(letter, " {} ".format(spelled))
    folded = unicodedata.normalize("NFKD", lowered).encode("ascii", "ignore").decode("ascii")
    return _NOT_ALNUM.sub(" ", folded).strip()


def bounds(needle: str) -> Tuple[str, str]:
    """The range of terms that start with ``needle``, under "C" collation.

    A normalised needle ends in a letter or digit, and the character after it
    is still a single byte, so the upper bound is exact: every term that starts
    with the needle sorts below it, and nothing else in between does.
    """
    if not needle:
        raise ValueError("an empty needle has no range")
    return needle, needle[:-1] + chr(ord(needle[-1]) + 1)


def versionless(accession: str) -> str:
    """ENST00000381330.5 -> ENST00000381330."""
    return accession.split(".", 1)[0]


# ----------------------------------------------------------------- sources


def parse_mane(lines: Iterable[str]) -> Dict[str, dict]:
    """MANE's summary file, as MANE Select rows keyed by versionless ENST."""
    rows = iter(lines)
    header = next(rows).lstrip("#").rstrip("\n").split("\t")
    found = {}
    for line in rows:
        values = line.rstrip("\n").split("\t")
        if len(values) != len(header):
            continue
        row = dict(zip(header, values))
        if row["MANE_status"] != "MANE Select":
            continue
        found[versionless(row["Ensembl_nuc"])] = row
    return found


def parse_refseqgene(lines: Iterable[str]) -> Dict[str, str]:
    """LRG_RefSeqGene, as the RefSeqGene record each RefSeq protein sits on.

    Any category counts. DMD's MANE protein, NP_003997.2, is on NG_012232.1
    only as "aligned: Selected"; the "reference standard" row there names the
    older NP_003997.1. Where a protein is on more than one record, the
    reference-standard row wins.
    """
    found = {}  # protein -> (rank, record)
    for line in lines:
        if line.startswith("#"):
            continue
        values = line.rstrip("\n").split("\t")
        if len(values) < 10:
            continue
        record, protein, category = values[3], values[7], values[9]
        if not record or not protein:
            continue
        rank = 0 if category == "reference standard" else 1
        held = found.get(protein)
        if held is None or rank < held[0]:
            found[protein] = (rank, record)
    return {protein: record for protein, (_, record) in found.items()}


# ----------------------------------------------------------------- an entry


def _value(node: Optional[dict]) -> Optional[str]:
    return (node or {}).get("value")


def _names(description: dict) -> Tuple[str, List[str]]:
    """The recommended name, and every other name UniProt gives the protein.

    The others are the recommended name's short forms, the alternative names
    and theirs, CD antigen and INN names, and the names of what the precursor
    is cut into or made of -- glucagon-like peptide 1 is found under GCG.
    """
    recommended = description.get("recommendedName") or {}
    name = _value(recommended.get("fullName"))
    if name is None:
        submitted = (description.get("submissionNames") or [{}])[0]
        name = _value(submitted.get("fullName")) or ""

    others = []

    def collect(block):
        others.append(_value(block.get("fullName")))
        others.extend(_value(short) for short in block.get("shortNames") or [])

    collect(recommended)
    for block in description.get("alternativeNames") or []:
        collect(block)
    for key in ("cdAntigenNames", "innNames"):
        others.extend(_value(node) for node in description.get(key) or [])
    for key in ("contains", "includes"):
        for part in description.get(key) or []:
            collect(part.get("recommendedName") or {})
            for block in part.get("alternativeNames") or []:
                collect(block)
    return name, [other for other in others if other and other != name]


def genes_of(entry: dict) -> Tuple[List[str], List[str]]:
    """The entry's gene names, in UniProt's order, and every synonym."""
    names, synonyms = [], []
    for gene in entry.get("genes") or []:
        named = _value(gene.get("geneName"))
        if named:
            names.append(named)
        for key in ("synonyms", "orfNames"):
            synonyms.extend(_value(node) for node in gene.get(key) or [])
    return names, [synonym for synonym in synonyms if synonym]


def _existence(entry: dict) -> int:
    """"1: Evidence at protein level" -> 1."""
    stated = entry.get("proteinExistence") or ""
    head = stated.split(":", 1)[0].strip()
    return int(head) if head.isdigit() else 5


def _xrefs(entry: dict, database: str) -> List[dict]:
    return [x for x in entry.get("uniProtKBCrossReferences") or []
            if x.get("database") == database]


def _property(xref: dict, key: str) -> Optional[str]:
    for prop in xref.get("properties") or []:
        if prop.get("key") == key:
            return prop.get("value")
    return None


def rows_for_entry(
    entry: dict,
    mane: Dict[str, dict],
    refseqgene: Dict[str, str],
    mane_release: str,
) -> List[dict]:
    """One UniProt entry's rows: one per MANE gene, or one under its own name.

    UniProt leaves the isoform off a MANE cross-reference when MANE's protein
    is the canonical sequence, and names it when it is not -- 1,006 of 18,601
    in release 2026_03, 28 of them tagged "-1" because their canonical is
    another isoform. A tagged row is not buildable yet: UniProt numbers its
    features on the canonical sequence, and the resolver reads regions,
    disulfides and names from those features.
    """
    accession = entry["primaryAccession"]
    name, _ = _names(entry.get("proteinDescription") or {})
    gene_names, _ = genes_of(entry)
    base = {
        "uniprot": accession,
        "name": name,
        "length": int((entry.get("sequence") or {}).get("length") or 0),
        "annotation_score": int(round(float(entry.get("annotationScore") or 0))),
        "existence": _existence(entry),
        "gene_id": None, "hgnc_id": None,
        "refseq_nuc": None, "refseq_prot": None,
        "ensembl_nuc": None, "ensembl_prot": None, "mane_isoform": None,
        "chrom_acc": None, "chrom_start": None, "chrom_end": None,
        "chrom_strand": None, "refseqgene": None,
    }
    own_gene = gene_names[0] if gene_names else ""
    gene_ids = [x.get("id") for x in _xrefs(entry, "GeneID")]
    hgnc_ids = [x.get("id") for x in _xrefs(entry, "HGNC")]

    rows = {}
    for xref in _xrefs(entry, "MANE-Select"):
        found = mane.get(versionless(xref.get("id") or ""))
        if found is None:
            row = dict(base, gene=own_gene, buildable=False,
                       unavailable_reason="MANE {} has no transcript {}.".format(
                           mane_release, xref.get("id")))
            rows.setdefault(row["gene"], row)
            continue
        isoform = xref.get("isoformId")
        row = dict(
            base,
            gene=found["symbol"],
            gene_id=int(found["NCBI_GeneID"].split(":")[-1]),
            hgnc_id=found["HGNC_ID"] or None,
            refseq_nuc=found["RefSeq_nuc"],
            refseq_prot=found["RefSeq_prot"],
            ensembl_nuc=found["Ensembl_nuc"],
            ensembl_prot=found["Ensembl_prot"],
            mane_isoform=isoform,
            chrom_acc=found["GRCh38_chr"],
            chrom_start=int(found["chr_start"]),
            chrom_end=int(found["chr_end"]),
            chrom_strand=1 if found["chr_strand"] == "+" else -1,
            refseqgene=refseqgene.get(found["RefSeq_prot"]),
            buildable=isoform is None,
            unavailable_reason=None if isoform is None else (
                "MANE Select encodes isoform {}, and UniProt numbers its "
                "features on the canonical sequence.".format(isoform)),
        )
        rows[row["gene"]] = row

    if not rows:
        reason = ("No MANE Select transcript." if own_gene
                  else "UniProt names no gene, and there is no MANE Select transcript.")
        rows[own_gene] = dict(base, gene=own_gene, buildable=False,
                              unavailable_reason=reason)
        if len(gene_ids) == 1:
            rows[own_gene]["gene_id"] = int(gene_ids[0])
        if len(hgnc_ids) == 1:
            rows[own_gene]["hgnc_id"] = hgnc_ids[0]
    return list(rows.values())


def terms_for_row(entry: dict, row: dict, mane_name: Optional[str] = None) -> Set[Tuple[str, int]]:
    """Every (term, kind) a search can find this row by, normalised.

    The row's own gene is its symbol. The entry's other gene names are its
    synonyms, so "HBA2" still finds HBA1's row -- ranked below HBA2's own.
    """
    name, others = _names(entry.get("proteinDescription") or {})
    gene_names, synonyms = genes_of(entry)
    terms = set()

    def add(text, kind):
        term = normalize(text or "")
        if term:
            terms.add((term, kind))

    add(entry["primaryAccession"], ACCESSION)
    add(row["gene"], SYMBOL)
    for other in gene_names + synonyms:
        if normalize(other) != normalize(row["gene"]):
            add(other, SYNONYM)

    for text, kind in [(name, NAME), (mane_name, NAME)] + [(other, OTHER_NAME) for other in others]:
        add(text, kind)
        terms.update(_words(text))
    return terms


def _words(text: Optional[str]) -> Set[Tuple[str, int]]:
    """The words a name can be found by, apart from the name itself."""
    return {
        (word, WORD)
        for word in normalize(text or "").split(" ")
        if len(word) >= _MIN_WORD and word not in _STOPWORDS
    }


def curated_terms(display: str, slug: str, accession: str) -> Set[Tuple[str, int]]:
    """What a listed protein is also called: the app's own name for it.

    "hemoglobin" is the slug HBB has always had, so it has to find HBB exactly,
    ahead of every other hemoglobin.
    """
    terms = set()
    for text, kind in ((display, NAME), (slug, NAME), (accession, ACCESSION)):
        term = normalize(text)
        if term:
            terms.add((term, kind))
    return terms | _words(display)
