"""Dry run of the resolver against live NCBI and UniProt. Writes nothing.

    .venv/bin/python scripts/dry_resolve.py RBP4 SEC61G
    .venv/bin/python scripts/dry_resolve.py --sample      # a fixed sample of 20

For each gene it reads the buildable `protein_index` row (a read-only
transaction), runs `resolve.resolve` exactly as the worker does, and puts the
record through the uploader's gate. No request, row or object is created: the
only thing written is the builder's GenBank cache, in a temporary directory
removed when the run ends. ESM-2 is not run, so it cannot say whether a
protein's scores will clear the alignment gate.
"""

from __future__ import annotations

import atexit
import json
import os
import shutil
import sys
import tempfile
import time
import traceback
from pathlib import Path

BACKEND = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BACKEND))

from dotenv import dotenv_values  # noqa: E402

# Named before any baker is imported: each reads its directory once, then.
_scratch = Path(tempfile.mkdtemp(prefix="helixpeek-dry-resolve-"))
atexit.register(shutil.rmtree, _scratch, True)
_env = dotenv_values(BACKEND / ".env")
os.environ["NCBI_EMAIL"] = _env["NCBI_EMAIL"]
os.environ["HELIXPEEK_GB_CACHE"] = str(_scratch / "genbank")
os.environ["HELIXPEEK_DATA"] = str(_scratch / "data")

import psycopg  # noqa: E402

from pipeline import uniprot, upload_tracks  # noqa: E402
from pipeline.resolver import store  # noqa: E402
from pipeline.resolver.resolve import Refused, resolve  # noqa: E402

SAMPLE_SLICE = 14
SAMPLE_REFSEQGENE = 6


def rows_for(conn, genes: list[str]):
    mane = store.mane_release(conn)
    found = []
    for gene in genes:
        hit = conn.execute(
            "select uniprot from protein_index where gene = %s and buildable", (gene,)).fetchone()
        if hit is None:
            print(f"{gene}: no buildable index entry")
            continue
        found.append(store.index_row(conn, hit[0], gene))
    return mane, found


def sample(conn) -> list[str]:
    """A fixed sample: the order is a hash of the gene, so it is the same every run."""
    def pick(where: str, count: int) -> list[str]:
        return [r[0] for r in conn.execute(
            f"select gene from protein_index where buildable and gene <> '' and {where} "
            "and chrom_end - chrom_start < 400000 order by md5(gene) limit %s", (count,)).fetchall()]
    return pick("refseqgene is null", SAMPLE_SLICE) + pick("refseqgene is not null", SAMPLE_REFSEQGENE)


def main(genes: list[str]) -> int:
    with psycopg.connect(_env["DATABASE_URL"], prepare_threshold=None) as conn:
        with conn.transaction():
            conn.execute("set transaction read only")
            if "--sample" in genes:
                genes = [g for g in genes if g != "--sample"] + sample(conn)
            mane, rows = rows_for(conn, genes)

    tally = {"ok": 0, "refused": 0, "error": 0}
    for row in rows:
        path = "RefSeqGene" if row.refseqgene else "slice"
        started = time.perf_counter()
        try:
            resolution = resolve(row, uniprot.fetch_entry(row.uniprot), mane_release=mane)
            provenance = upload_tracks.validate("record", resolution.target, resolution.record)
        except Refused as exc:
            tally["refused"] += 1
            print(f"REFUSED {row.gene:10s} {path:10s} {row.length:5d} aa  {' '.join(str(exc).split())[:230]}")
            continue
        except ValueError as exc:
            tally["refused"] += 1
            print(f"GATE    {row.gene:10s} {path:10s} {row.length:5d} aa  {' '.join(str(exc).split())[:230]}")
            continue
        except Exception as exc:  # noqa: BLE001 -- these are the ones to read
            tally["error"] += 1
            where = traceback.extract_tb(exc.__traceback__)[-1]
            print(f"ERROR   {row.gene:10s} {path:10s} {row.length:5d} aa  {type(exc).__name__}: "
                  f"{' '.join(str(exc).split())[:200]}  [{Path(where.filename).name}:{where.lineno} {where.name}]")
            continue
        tally["ok"] += 1
        target, protein = resolution.target, resolution.protein
        record = json.loads(resolution.record)
        source = target.source
        cut = f"{source.accession}" + (f":{source.seq_start}-{source.seq_stop}" if source.seq_start else "")
        regions = " ".join(f"{r['short'] or '·'}[{r['start']}-{r['end']}{'' if r['kept'] else ' removed'}]"
                           for r in protein["regions"])
        print(f"ok      {row.gene:10s} {path:10s} {row.length:5d} aa  {cut}  "
              f"exons {protein['exons']}, record {len(resolution.record):,} B"
              f"{', introns shortened' if record.get('intron_compression') or record.get('compressed') else ''}, "
              f"variants {len(protein['provenance']['uniprot_variants'])}, bridges {protein['bridges']}, "
              f"cleaved {target.cleaved}; {regions}; gate ok ({provenance['accession']}) "
              f"[{time.perf_counter() - started:.1f}s]")
    print(f"-- {tally} of {len(rows)}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    sys.exit(main(sys.argv[1:]))
