# helix-peak-backend

A minimal FastAPI service with exactly one endpoint: fetch a GenBank record from
NCBI and return its raw text. It exists so the HelixPeak Flutter app has a real
backend to talk to — the smallest slice that works end to end.

It wraps a single Biopython call:

```python
Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
```

## The endpoint

```
GET /fetch/{id}
```

`db`, `rettype` and `retmode` are hardcoded; `id` is the only input.

```console
$ curl localhost:8000/fetch/NG_007114
{
  "id": "NG_007114",
  "content": "LOCUS       NG_007114               8416 bp    DNA     linear   PRI 06-SEP-2026\nDEFINITION  Homo sapiens insulin (INS), RefSeqGene on chromosome 11.\n..."
}
```

The GenBank text is returned verbatim. Parsing the FEATURES table into structured
JSON is a later phase.

There is also `GET /health` → `{"status": "ok"}`, which checks the service is up
without spending an NCBI request.

## Setup

`NCBI_EMAIL` is required. NCBI rejects and eventually blocks traffic with no
contact address, so the app refuses to start without it rather than failing at
request time.

```bash
cp .env.example .env
# edit .env and set NCBI_EMAIL to a real address
```

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
network call.

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

So a valid response is recognised *positively*, by its `LOCUS` header — checking
for an empty body would pass NCBI's error text straight through to the client as
a `200`. Any upstream `400` is treated as "not found" because `id` is the only
part of the request that can be invalid.

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

`/search` and `/summary` are not here yet, and neither is GenBank parsing, caching,
or auth. This is intentionally the smallest thing that actually works: Flutter hits
it, gets real NCBI data back.

## Layout

```
app/
  main.py           FastAPI app: Entrez.email, socket timeout, CORS, router
  config.py         pydantic-settings; NCBI_EMAIL required
  entrez_client.py  the Entrez.efetch wrapper, nothing else
  router.py         GET /fetch/{id} and the 404/502 mapping
tests/
  test_fetch.py     endpoint tests against a mocked Entrez.efetch
  test_config.py    startup fails without NCBI_EMAIL
  fixtures/ng_007114.gb
```
