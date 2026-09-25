"""Build the protein index from UniProt and MANE, and write it to Supabase.

    .venv/bin/python scripts/load_protein_index.py --dry-run
    .venv/bin/python scripts/load_protein_index.py [--report PATH]

`DATABASE_URL` comes from `.env`, as it does for the service. Run it again
whenever UniProt (every eight weeks) or MANE publishes a release.

It downloads three things, none of them large:

- every reviewed human UniProt entry, as JSON, so names arrive structured and
  nothing here parses a name out of a string of parentheses;
- MANE's summary for the current release: the transcript, protein, gene and
  GRCh38 span for each MANE Select gene;
- NCBI's LRG_RefSeqGene table: which RefSeqGene record each RefSeq protein is
  annotated on.

The listed proteins -- the `protein` rows with a `catalog_order` -- are checked
against the rows built for their genes, and lend them the app's own names:
"hemoglobin" is HBB's slug and has to find HBB exactly.

Every row of `protein_index` and `protein_index_term` is replaced in one
transaction, with plain deletes rather than a truncate, so a search running
meanwhile reads the old index until the new one is whole. Nothing is written
when a check fails, and nothing at all with `--dry-run`.
"""

import argparse
import collections
import gzip
import io
import json
import os
import re
import sys
import time
import urllib.request
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import psycopg  # noqa: E402

from app.config import settings  # noqa: E402
from app.protein_index import (  # noqa: E402
    curated_terms,
    genes_of,
    parse_mane,
    parse_refseqgene,
    rows_for_entry,
    terms_for_row,
    versionless,
)

UNIPROT = (
    "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=json"
    "&query=%28organism_id%3A9606%29%20AND%20%28reviewed%3Atrue%29"
    "&fields=accession,gene_names,protein_name,length,annotation_score,"
    "protein_existence,xref_mane-select,xref_hgnc,xref_geneid"
)
MANE_DIR = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
REFSEQGENE = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene"

_AGENT = "helixpeek-index-loader"

_ROW_COLUMNS = (
    "uniprot", "gene", "name", "length", "annotation_score", "existence",
    "gene_id", "hgnc_id", "refseq_nuc", "refseq_prot", "ensembl_nuc",
    "ensembl_prot", "mane_isoform", "chrom_acc", "chrom_start", "chrom_end",
    "chrom_strand", "refseqgene", "buildable", "unavailable_reason",
)

_CATALOG = """
select slug, gene, uniprot, display, accession, protein_id, transcript_id
from protein
where catalog_order is not null
order by catalog_order
"""


def _get(url: str) -> Tuple[bytes, Dict[str, str]]:
    request = urllib.request.Request(url, headers={"User-Agent": _AGENT})
    with urllib.request.urlopen(request, timeout=300) as response:
        return response.read(), {k.lower(): v for k, v in response.headers.items()}


def download() -> Tuple[List[dict], str, Dict[str, dict], str, Dict[str, str]]:
    """UniProt entries and release, MANE rows and release, RefSeqGene map."""
    started = time.time()
    body, headers = _get(UNIPROT)
    entries = json.loads(gzip.decompress(body))["results"]
    uniprot_release = headers.get("x-uniprot-release", "unknown")
    print("UniProt {}: {:,} entries in {:.0f} s".format(
        uniprot_release, len(entries), time.time() - started))

    listing, _ = _get(MANE_DIR)
    names = re.findall(r"MANE\.GRCh38\.v([0-9.]+)\.summary\.txt\.gz", listing.decode())
    if not names:
        raise SystemExit("No MANE summary in {}".format(MANE_DIR))
    mane_release = "v" + names[0]
    body, _ = _get(MANE_DIR + "MANE.GRCh38.{}.summary.txt.gz".format(mane_release))
    mane = parse_mane(io.StringIO(gzip.decompress(body).decode()))
    print("MANE {}: {:,} MANE Select transcripts".format(mane_release, len(mane)))

    body, _ = _get(REFSEQGENE)
    refseqgene = parse_refseqgene(io.StringIO(body.decode()))
    print("LRG_RefSeqGene: {:,} proteins on a RefSeqGene record".format(len(refseqgene)))
    return entries, uniprot_release, mane, mane_release, refseqgene


def build(entries, mane, mane_release, refseqgene):
    """Every row and every term, before the catalog is laid over them."""
    rows, terms = [], {}
    own_genes = {}
    for entry in entries:
        # Synonyms count: B7ZAQ6 is "GPHRA" to UniProt and GPR89A to MANE,
        # and says so only among its synonyms.
        names, synonyms = genes_of(entry)
        own_genes[entry["primaryAccession"]] = {name.lower() for name in names + synonyms}
        for row in rows_for_entry(entry, mane, refseqgene, mane_release):
            found = mane.get(versionless(row["ensembl_nuc"] or ""))
            rows.append(row)
            terms[(row["uniprot"], row["gene"])] = terms_for_row(
                entry, row, found["name"] if found else None)

    # Two entries can list the same MANE transcript: B7ZAQ6 and P0CG08, Golgi
    # pH regulators A and B, are one sequence, and each lists both genes'
    # transcripts. A build is keyed by its gene, so one row keeps it: the entry
    # that names the gene, then the better annotated one.
    duplicates = []
    by_gene = collections.defaultdict(list)
    for row in rows:
        if row["buildable"]:
            by_gene[row["gene"].lower()].append(row)
    for gene, group in by_gene.items():
        if len(group) < 2:
            continue
        group.sort(key=lambda row: (
            gene not in own_genes[row["uniprot"]],
            -row["annotation_score"],
            row["uniprot"],
        ))
        keeper = group[0]
        for other in group[1:]:
            other["buildable"] = False
            other["unavailable_reason"] = (
                "UniProt's entry for {}'s MANE Select protein is {}.".format(
                    keeper["gene"], keeper["uniprot"]))
            duplicates.append((keeper["gene"], keeper["uniprot"], other["uniprot"]))

    # A losing row is only noise where its entry kept a gene of its own: a
    # search for GPR89B would offer P0CG08, and then B7ZAQ6 saying it is not
    # GPR89B's. It stays, unbuildable, only where it is all its entry has.
    rows_of = collections.Counter(row["uniprot"] for row in rows)
    kept = []
    for row in rows:
        lost = (not row["buildable"]
                and row["unavailable_reason"].startswith("UniProt's entry for "))
        if lost and rows_of[row["uniprot"]] > 1:
            rows_of[row["uniprot"]] -= 1
            del terms[(row["uniprot"], row["gene"])]
            continue
        kept.append(row)
    return kept, terms, duplicates, own_genes


def overlay_catalog(catalog, rows, terms) -> Tuple[List[str], List[str]]:
    """Check the listed proteins against their rows, and lend them their names.

    Returns (failures, notes). A failure stops the write: a listed protein
    with no buildable row would vanish from its own search, and a slug a build
    would take from another protein would open the wrong walk.
    """
    failures, notes = [], []
    buildable = {row["gene"]: row for row in rows if row["buildable"]}
    slugs = {slug: gene for (slug, gene, *_rest) in catalog}

    for slug, gene, uniprot, display, accession, protein_id, transcript_id in catalog:
        row = buildable.get(gene)
        if row is None:
            failures.append("{} ({}): no buildable row for the gene".format(slug, gene))
            continue
        if row["uniprot"] != uniprot:
            notes.append("{}: catalog says {}, the index row is {}".format(
                slug, uniprot, row["uniprot"]))
        if protein_id and row["refseq_prot"] != protein_id:
            notes.append("{}: catalog protein {}, MANE {}".format(
                slug, protein_id, row["refseq_prot"]))
        if transcript_id and row["refseq_nuc"] != transcript_id:
            notes.append("{}: catalog transcript {}, MANE {}".format(
                slug, transcript_id, row["refseq_nuc"]))
        if accession.startswith("NG_"):
            if versionless(row["refseqgene"] or "") != versionless(accession):
                notes.append("{}: catalog record {}, RefSeqGene for MANE's protein {}".format(
                    slug, accession, row["refseqgene"]))
        elif row["refseqgene"]:
            notes.append("{}: read from a slice of {}, though {} now has MANE's protein".format(
                slug, accession, row["refseqgene"]))
        terms[(row["uniprot"], row["gene"])] |= curated_terms(display, slug, accession)

    for gene, row in buildable.items():
        owner = slugs.get(gene.lower())
        if owner is not None and owner != gene:
            failures.append("{} would take slug {!r}, which is {}'s".format(
                gene, gene.lower(), owner))
    return failures, notes


def report(rows, terms, duplicates, failures, notes, releases) -> str:
    # Reasons grouped by their wording, with the accession, transcript or gene
    # each one names taken out.
    reasons = collections.Counter(
        re.sub(r"entry for \S+'s", "entry for …'s",
               re.sub(r"[A-Z][0-9][A-Z0-9]{3}[0-9]([A-Z][A-Z0-9]{2}[0-9])?(-[0-9]+)?"
                      r"|ENST[0-9.]+", "…",
                      row["unavailable_reason"]))
        for row in rows if not row["buildable"])
    lines = [
        "Protein index · UniProt {} · MANE {}".format(*releases),
        "",
        "{:,} entries, {:,} rows, {:,} buildable, {:,} search terms".format(
            len({row["uniprot"] for row in rows}), len(rows),
            sum(1 for row in rows if row["buildable"]),
            sum(len(found) for found in terms.values())),
        "{:,} buildable rows read from a RefSeqGene record, the rest from a slice".format(
            sum(1 for row in rows if row["buildable"] and row["refseqgene"])),
        "",
        "Not buildable, by reason:",
    ]
    lines += ["  {:>6,}  {}".format(count, reason) for reason, count in reasons.most_common()]
    lines += ["", "Genes two entries claim ({}):".format(len(duplicates))]
    lines += ["  {} kept {}, not {}".format(*found) for found in duplicates]
    lines += ["", "Listed proteins: {} failures, {} notes".format(len(failures), len(notes))]
    lines += ["  FAIL " + failure for failure in failures]
    lines += ["  note " + note for note in notes]
    return "\n".join(lines) + "\n"


def write(conn, rows, terms, releases) -> None:
    with conn.transaction():
        conn.execute("delete from protein_index_term")
        conn.execute("delete from protein_index")
        with conn.cursor() as cur:
            with cur.copy("copy protein_index ({}) from stdin".format(
                    ", ".join(_ROW_COLUMNS))) as copy:
                for row in rows:
                    copy.write_row(tuple(row[column] for column in _ROW_COLUMNS))
            with cur.copy("copy protein_index_term (uniprot, gene, term, kind) from stdin") as copy:
                for (uniprot, gene), found in terms.items():
                    for term, kind in sorted(found):
                        copy.write_row((uniprot, gene, term, kind))
        conn.execute(
            """
            insert into protein_index_release (only_row, uniprot, mane, entries, buildable, loaded_at)
            values (true, %s, %s, %s, %s, now())
            on conflict (only_row) do update set
                uniprot = excluded.uniprot, mane = excluded.mane,
                entries = excluded.entries, buildable = excluded.buildable,
                loaded_at = excluded.loaded_at
            """,
            (releases[0], releases[1], len({row["uniprot"] for row in rows}),
             sum(1 for row in rows if row["buildable"])),
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--dry-run", action="store_true",
                        help="build and check everything, write nothing")
    parser.add_argument("--report", help="also write the report to this file")
    arguments = parser.parse_args()

    if not settings.database_url:
        raise SystemExit("DATABASE_URL is not set; the catalog check needs it.")

    entries, uniprot_release, mane, mane_release, refseqgene = download()
    rows, terms, duplicates, _ = build(entries, mane, mane_release, refseqgene)

    with psycopg.connect(settings.database_url, prepare_threshold=None) as conn:
        catalog = conn.execute(_CATALOG).fetchall()
        failures, notes = overlay_catalog(catalog, rows, terms)
        text = report(rows, terms, duplicates, failures, notes,
                      (uniprot_release, mane_release))
        print("\n" + text)
        if arguments.report:
            with open(arguments.report, "w") as handle:
                handle.write(text)
        if failures:
            print("Not written: {} listed-protein check(s) failed.".format(len(failures)))
            return 1
        if arguments.dry_run:
            print("Dry run: nothing written.")
            return 0
        started = time.time()
        write(conn, rows, terms, (uniprot_release, mane_release))
    print("Written in {:.1f} s.".format(time.time() - started))
    return 0


if __name__ == "__main__":
    sys.exit(main())
