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

## Rate limiting

Not implemented here, and not needed yet: Biopython already sleeps ~0.37s between
calls to stay under NCBI's 3 requests/second limit for keyless clients. Add a real
limiter once this serves more than one developer.

## CORS

`allow_origins=["*"]`, so `flutter run -d chrome` works. That is fine for a local
single-user dev service and should be tightened before this is exposed anywhere.

## Deliberately not built

`/search` and `/summary` are not here yet, and neither is caching or auth. Parsing
covers one gene per request — there is no endpoint that returns every gene in a
record, because no client needs one yet.

## Layout

```
app/
  main.py             FastAPI app: Entrez.email, socket timeout, CORS, router,
                      pool lifespan, /health and /health/db
  config.py           pydantic-settings; NCBI_EMAIL required, DATABASE_URL optional
  db.py               the psycopg connection pool; sync, to match the read path
  entrez_client.py    the Entrez.efetch wrapper, returning a parsed record
  genbank_parser.py   extract_gene(record, gene) -- pure, no I/O
  schemas.py          pydantic response models for /gene
  router.py           the endpoint and the 404/502 mapping
tests/
  test_gene.py            /gene against a mocked Entrez.efetch
  test_health.py          the liveness check
  test_genbank_parser.py  the extractor, straight off the fixture
  test_config.py          startup fails without NCBI_EMAIL
  test_health_db.py       the readiness contract, without opening a socket
  conftest.py             shared fixtures, incl. the efetch patches
  ncbi_errors.py          real NCBI failure bodies
  fixtures/ng_007114.gb
```
