# helix-peek-backend

A minimal FastAPI service that fetches a GenBank record from NCBI and returns one
gene lifted out of it. It exists so the Helix Peek Flutter app has a real backend
to talk to.

The endpoint wraps a single Biopython call:

```python
Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
```

`db`, `rettype` and `retmode` are hardcoded; `id` is the only input to it.

## The endpoint

```
GET /gene/{id}/{gene}    one gene from a GenBank record, structured
```

### `GET /gene/{id}/{gene}`

Parses the FEATURES table and returns just the requested gene: its span, sequence,
exons, transcript, protein, and peptides. Coordinates are **1-based inclusive**,
matching what NCBI shows for the record rather than Biopython's 0-based half-open
internals. `join(...)` locations keep their parts in `segments`, so introns are not
silently spliced out.

The payload is only what the app draws. A length a client can compute from what is
already here — a span's `end - start + 1`, a translation's length — is not sent, so
there is never a second copy of a number to disagree with the first.

```console
$ curl localhost:8000/gene/NG_007114/INS
{
  "gene": "INS",
  "location": {"start": 4986, "end": 6416, "strand": 1},
  "sequence": "AGCCCTCCAGGACAGGCTGCATCAGAAGAGG...",
  "transcript": {"segments": [{"start": 4986, "end": 5027}, ...]},
  "protein": {
    "product": "insulin preproprotein",
    "translation": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT...",
    "segments": [{"start": 5224, "end": 5410}, {"start": 6198, "end": 6343}]
  },
  "exons": [{"number": 1, "start": 4986, "end": 5027}, ...],
  "signal_peptide": {"product": null, "segments": [...], "translation": "MALWMRLLPLLALLALWGPDPAAA"},
  "proprotein": {"product": "proinsulin", "segments": [...], "translation": "FVNQHLCGSHLVEA..."},
  "peptides": [
    {"product": "insulin B chain", "segments": [...], "translation": "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"},
    ...
  ]
}
```

Both `id` and `gene` are free parameters: which record to read and which gene to
keep are the client's choice. The Flutter app defaults to `NG_007114` / `INS`.

#### Why the gene filter is an exact match

`NG_007114` is a RefSeqGene *region*, not a single gene. It carries three: `TH`
(1..2266), `INS` (4986..6416) and the `INS-IGF2` readthrough (4986..>8416). `INS`
and `INS-IGF2` **start at the same base and share a signal peptide**, so no
coordinate window separates them — `efetch` with `seq_start=4986&seq_stop=6416`
still returns five `INS-IGF2` features. Selection is therefore an exact comparison
against the `/gene` qualifier; `"INS-IGF2".startswith("INS")` is why a prefix or
substring test would be wrong.

There is also `GET /health` → `{"status": "ok"}`, which checks the service is up
without spending an NCBI request.

### AVI contribution details

`GET /gene/{id}/{gene}/impact-explanations` serves the versioned, generated
attribution payload for INS, HBB and CFTR. It reads the same JSON that the mobile
app bundles in mock mode, using exact accession/gene identity, with no NCBI or
Atlas calls. Missing coverage returns 404; malformed generated data returns 503.

`IMPACT_EXPLANATIONS_DIR` defaults to the companion checkout's
`helix-peek/assets/impact_explanations`. Package or mount that generated directory
and configure the path when deploying. The endpoint does not require an
AlphaGenome API key. [Baking and the contract](../helix-peek/docs/avi-contributions.md).

## Setup

`NCBI_EMAIL` is required. NCBI rejects and eventually blocks traffic with no
contact address, so the app refuses to start without it rather than failing at
request time.

```bash
cp .env.example .env
# edit .env and set NCBI_EMAIL to a real address
```

### The cache database

`DATABASE_URL` is optional. Unset, the service reads through to NCBI for every
request; the pool is simply never created. It is opened in the app lifespan
rather than at import, so a wrong or unreachable database is reported on
`GET /health/db` instead of killing the worker before it can say why.

On Supabase, take the **connection pooler** URI (Project Settings → Database),
not the direct `db.<ref>.supabase.co` one: direct hosts resolve to IPv6 only,
which Render cannot reach. The URI needs `?sslmode=require`.

```
DATABASE_URL=postgresql://postgres.<ref>:<password>@aws-<n>-<region>.pooler.supabase.com:5432/postgres?sslmode=require
```

Port 5432 is the session pooler, which is what a long-lived server with its own
client-side pool wants. Port 6543 is the transaction pooler, for serverless
callers; it disallows prepared statements, so it needs psycopg configured for
that before it will work here.

`GET /health/db` returns `unconfigured`, `ok` with the server version, or
`error` with just the exception class. The endpoint is public and psycopg names
the host, port and user in its connection errors, so the full message goes to
the log instead of the response body.

#### What is cached

Whole GenBank records, keyed by accession -- not by `(accession, gene)`. The
expensive step is the NCBI round trip for a record; `extract_gene` is pure and
costs ~0.1 ms, so one cached row serves every gene in it.

What is stored is the text NCBI returned, not a parsed payload. Reparsing takes
~2 ms against a 0.5-3 s fetch, and raw text keeps the cache independent of both
the parser and the response schema: a fix to `genbank_parser` changes what
clients see on the next request without invalidating a thing.

The table is created at startup (`create table if not exists`) -- one table, no
history to migrate, and a deploy that cannot half-apply.

`CACHE_TTL_DAYS` defaults to 30. A record for a given accession is effectively
immutable, so the TTL is about picking up NCBI's own corrections rather than
about staleness.

Every cache operation is best-effort. A database that is missing, asleep or
broken is logged and read through to NCBI, so it costs latency and nothing
else. Nothing unparseable is ever written: NCBI answers a bad id with plain
text and HTTP 200, and the parse happens before the write.

### Run locally

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements-dev.txt
.venv/bin/uvicorn app.main:app --reload
```

Verified on macOS system Python 3.9.6; the code avoids 3.10+ syntax so it runs on
both that and the 3.12 in the Docker image.

### Run with Docker

```bash
docker compose up --build
```

`.env` is read via `env_file`, so create it first.

### Tests

```bash
.venv/bin/pytest
```

Every test patches `Entrez.efetch` against a saved fixture
(`tests/fixtures/ng_007114.gb`, a real NCBI response), so the default run makes no
network call. The parser tests skip HTTP entirely — `extract_gene` takes an
already-parsed record and does no I/O.

## Error handling

NCBI never signals "no such record" with an empty body, and it does not use HTTP
status codes consistently. Both real failure shapes are mapped to `404`:

| `id` | What NCBI actually returns | We return |
|---|---|---|
| `NG_007114` | `200`, the GenBank record | `200` |
| `NOT_A_REAL_ID_XYZ` (unparseable) | **`200`**, body `Error: F a i l e d  t o  u n d e r s t a n d  i d : ...` | `404` |
| `NG_999999` (well-formed, unknown) | **`400`**, `Error: ... F a i l e d  t o  r e t r i e v e  s e q u e n c e : ...` | `404` |
| — | connection refused, DNS failure, timeout | `502`, upstream message in `detail` |
| — | any other upstream status (5xx) | `502`, upstream message in `detail` |

(The spaced-out letters are literal. NCBI really sends `F a i l e d`.)

There is one more 404: the record parsed, but no feature in it carries that gene
name.

A valid response is recognised *positively* — `SeqIO.read` raises `ValueError` on
a body that is not a GenBank record. Checking for an empty body instead would pass
NCBI's error text straight through to the client as a `200`. Any upstream `400` is
treated as "not found" because `id` is the only part of the request that can be
invalid.

Requests time out after `NCBI_TIMEOUT_SECONDS` (default 20). `Entrez.efetch` takes
no timeout argument and urllib has no default, so this is enforced with
`socket.setdefaulttimeout` at startup.

## Suggestions: every reviewed human protein

`GET /proteins/suggest?q=insul&limit=12` answers as someone types. The twenty
listed proteins are among the answers, tagged `listed`, but so is every other
reviewed human protein, each saying what the app can do with it: `ready`
(built earlier), `buildable`, or `unavailable` with a `reason`.

It reads `protein_index`, built from UniProt's reviewed human entries joined to
MANE Select, one row per (entry, gene) because some entries are made by two
genes (P69905 by HBA1 and HBA2). Search terms are normalised to lower-case
ASCII and stored under "C" collation, so every prefix is a plain btree range
that works under psycopg's prepared, generic plans too. A trigram index answers
typos, but only when no prefix matches at all. Server-side, a query takes 6-30 ms.

To build it, or to refresh it after a UniProt or MANE release:

```sh
.venv/bin/python scripts/apply_migration.py migrations/0002_protein_index.sql   # once
.venv/bin/python scripts/load_protein_index.py --dry-run                         # report only
.venv/bin/python scripts/load_protein_index.py                                   # replace the index
.venv/bin/python scripts/check_suggest.py https://helix-peak-backend.onrender.com
```

The loader refuses to write if a listed protein has no buildable row, or if a
build would take a listed protein's slug. `apply_migration.py` is `psql -f` for a
machine without psql.

## Rate limiting

Not implemented here, and not needed yet: Biopython already sleeps ~0.37s between
calls to stay under NCBI's 3 requests/second limit for keyless clients. Add a real
limiter once this serves more than one developer.

## CORS

`allow_origins=["*"]`, so `flutter run -d chrome` works. That is fine for a local
single-user dev service and should be tightened before this is exposed anywhere.

## Deliberately not built

Auth is not here yet, and neither is `/summary`. Parsing covers one gene per
request — there is no endpoint that returns every gene in a record, because no
client needs one yet.

## Layout

```
app/
  main.py             FastAPI app: Entrez.email, socket timeout, CORS, router,
                      pool lifespan, /health and /health/db
  config.py           pydantic-settings; NCBI_EMAIL required, DATABASE_URL optional
  db.py               the psycopg connection pool; sync, to match the read path
  entrez_client.py    the Entrez.efetch wrapper; fetch and parse kept separate
  record_cache.py     read-through cache of whole records, keyed by accession
  genbank_parser.py   extract_gene(record, gene) -- pure, no I/O
  protein_index.py    how index rows and search terms are made, and normalize()
  suggest.py          /proteins/suggest: ranked prefix tiers, near misses last
  schemas.py          pydantic response models
  router.py           the endpoints and the 404/502/503 mapping
migrations/           applied by hand, in order
pipeline/             the bake tools, moved from helix-peek/tool/ in Phase 3; see
                      pipeline/README.md. Nothing under app/ imports it yet
  targets.py          the twenty curated proteins' bake table
  curated/catalog.json  their hand-written prose, facts and tints
  fetch_tracks.py     storage's tracks back into pipeline/data/, sha256-checked
  check_assets.py     the tracks against each other, the tables and the service
  upload_tracks.py    validated tracks into storage, and their protein_track rows
  seed_catalog.py     the protein and protein_alias rows; --check diffs the live ones
  mock/ constraint/ impact/ clinvar/ structure/   one baker each
scripts/
  apply_migration.py      psql -f, for a machine without psql
  load_protein_index.py   UniProt + MANE + LRG_RefSeqGene -> protein_index
  check_suggest.py        golden queries and warm timings against a running service
tests/
  test_gene.py            /gene against a mocked Entrez.efetch
  test_health.py          the liveness check
  test_genbank_parser.py  the extractor, straight off the fixture
  test_config.py          startup fails without NCBI_EMAIL
  test_health_db.py       the readiness contract, without opening a socket
  test_record_cache.py    hit, miss, eviction and degradation, against a fake pool
  test_protein_index.py   normalising, MANE and RefSeqGene parsing, rows and terms
  test_suggest.py         statuses, the short and near-miss rules, against a fake pool
  test_pipeline.py        the pipeline imports, the curated rows, check_against's comparisons
  conftest.py             shared fixtures, incl. the efetch patches
  ncbi_errors.py          real NCBI failure bodies
  fixtures/ng_007114.gb
```
