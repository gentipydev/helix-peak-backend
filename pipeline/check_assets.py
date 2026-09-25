"""Check the baked tracks against each other, and against the curated rows.

Several bakes write the tracks from one table, `targets.py`, and
`curated/catalog.json` holds the hand-written half of each row -- what
`protein_catalog.dart` held until Phase 3 of HANDOFF-ONDEMAND.md. Nothing but
this notices when they come apart: a constraint track baked before a gene
record was rebuilt still parses, still loads, and quietly stops colouring the
page it is for.

    python3 pipeline/check_assets.py
    python3 pipeline/check_assets.py --against https://helix-peak-backend.onrender.com

The offline checks read the files under `pipeline/data/` (`paths.DATA`): a
bake's own output before `upload_tracks.py` sends it, or storage's copy after
`fetch_tracks.py` brings it back. Against storage's copy they prove the stored
tracks agree with each other and with the tables; only a bake's output, or
`fetch_tracks.py --verify-client`, compares storage with something else.

Exits non-zero on the first disagreement, with both sides named.

`--against` checks the copy the service serves: its catalog rows, and the
record and scene containers they name. It reads `targets.py` and the curated
file directly rather than going through `seed_catalog.py`, so a seeder that
writes the wrong thing is caught rather than confirmed -- the same reason
`constraint/verify_cpu.py` re-derives its numbers instead of calling the scorer.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import struct
import sys

BACKEND = Path(__file__).resolve().parents[1]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))

from pipeline.impact.check_explanations import validate as validate_explanations  # noqa: E402
from pipeline.paths import CURATED, DATA  # noqa: E402
from pipeline.targets import TARGETS, Target, partition  # noqa: E402

problems: list[str] = []


def fail(message: str) -> None:
    problems.append(message)


def glb_nodes(path: Path) -> set[str]:
    """The node names in a binary glTF, read straight out of its JSON chunk."""
    blob = path.read_bytes()
    magic, _, _ = struct.unpack_from("<III", blob, 0)
    if magic != 0x46546C67:
        raise ValueError(f"{path} is not a .glb")
    length, kind = struct.unpack_from("<II", blob, 12)
    if kind != 0x4E4F534A:
        raise ValueError(f"{path} does not start with a JSON chunk")
    document = json.loads(blob[20 : 20 + length])
    return {n["name"] for n in document.get("nodes", []) if n.get("name")}


def curated_catalog() -> dict[str, dict]:
    """The curated rows, keyed by slug, in the shape the checks below compare.

    Until Phase 3 these were read off `protein_catalog.dart` with regexes, and
    the row kept the names that reading gave it -- `tints`, `modelled`, the four
    booleans the seeded track maps stood for -- so every comparison below is
    still the one it was. The families a protein has are `tracks` in the file:
    a family listed there is one the protein is baked with.
    """
    rows: dict[str, dict] = {}
    for protein in json.loads(CURATED.read_text())["proteins"]:
        chrome = protein.get("structure") or {}
        kinds = set(protein.get("tracks") or [])
        modelled = chrome.get("modelled")
        rows[protein["slug"]] = {
            "slug": protein["slug"],
            "catalog_order": protein["catalog_order"],
            "display": protein["display"],
            "gene": protein["gene"],
            "uniprot": protein["uniprot"],
            "accession": protein["accession"],
            "summary": protein["summary"],
            "chain": protein.get("chain"),
            "facts": protein["facts"],
            "chains": [chain["node"] for chain in protein["chains"]],
            "tints": [(chain["node"], chain["tint"]) for chain in protein["chains"]],
            "chrome": chrome,
            "pdb": chrome.get("pdb"),
            "modelled": tuple(modelled) if modelled else None,
            "count": chrome.get("count"),
            "scored": "constraint" in kinds,
            "impact_scored": "impact" in kinds,
            "clinvar_available": "clinvar" in kinds,
            "impact_explanations": "impact_explanations" in kinds,
        }
    return rows


def service_rows(base_url: str) -> dict[str, dict]:
    """`GET /protein/{slug}` for every target, keyed by slug."""
    import urllib.error
    import urllib.request

    found: dict[str, dict] = {}
    for target in TARGETS:
        url = f"{base_url.rstrip('/')}/protein/{target.slug}"
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                found[target.slug] = json.loads(response.read())
        except urllib.error.HTTPError as exc:
            fail(f"{target.slug}: {url} returned {exc.code}")
        except (urllib.error.URLError, OSError, ValueError) as exc:
            raise SystemExit(f"Could not read {url}: {exc}")
    return found


def check_scenes(base_url: str) -> None:
    """The container the service names, against the node names the app paints.

    This is `structure_view.dart`'s `_paint` contract, and it is the one thing
    the offline suite can no longer see: the `.glb` on disk proves what the bake
    made, not what was compiled and uploaded from it. A container whose nodes
    were renamed or dropped draws an unpainted molecule on a phone -- since
    Phase 5, one that throws rather than draws -- and nothing else would notice.

    Names are looked for as raw bytes. `.fsceneb` is a binary container and this
    is deliberately not a parser: a substring miss is a real miss, and a hit is
    the name being in there somewhere, which is as much as a grep can promise
    and enough to catch a bake that renamed a chain.
    """
    import urllib.error
    import urllib.request

    for target in TARGETS:
        url = f"{base_url.rstrip('/')}/protein/{target.slug}/tracks"
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                tracks = json.loads(response.read())
        except (urllib.error.URLError, OSError, ValueError) as exc:
            raise SystemExit(f"Could not read {url}: {exc}")
        row = tracks.get("structure") or {}
        if row.get("state") != "ready" or not row.get("url"):
            fail(f"{target.slug}: the service has no structure to fetch ({row.get('state')})")
            continue
        if row.get("format") != "fsceneb":
            fail(f"{target.slug}: structure format is {row.get('format')!r}, not 'fsceneb'")
        try:
            with urllib.request.urlopen(row["url"], timeout=60) as response:
                blob = response.read()
        except (urllib.error.URLError, OSError) as exc:
            raise SystemExit(f"Could not read {row['url']}: {exc}")
        digest = hashlib.sha256(blob).hexdigest()
        if digest != row.get("sha256"):
            fail(f"{target.slug}: the container does not match the sha256 the row carries")
        wanted = {c.node for c in target.structure.chains} | (
            {"bonds"} if target.structure.bonds else set()
        )
        missing = {n for n in wanted if n.encode() not in blob}
        if missing:
            fail(f"{target.slug}: the served container has no node named {sorted(missing)}")


def check_records(base_url: str) -> None:
    """The gene record the service names, against the one the bake wrote.

    The walk reads its record from storage since Phase 2 of HANDOFF-ONDEMAND.md,
    so what is served has to be the file every offline check above was run on,
    byte for byte: the bake's gates held for those bytes and no others.
    """
    import urllib.error
    import urllib.request

    for target in TARGETS:
        url = f"{base_url.rstrip('/')}/protein/{target.slug}/tracks"
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                tracks = json.loads(response.read())
        except (urllib.error.URLError, OSError, ValueError) as exc:
            raise SystemExit(f"Could not read {url}: {exc}")
        row = tracks.get("record") or {}
        if row.get("state") != "ready" or not row.get("url"):
            fail(f"{target.slug}: the service has no record to fetch ({row.get('state')})")
            continue
        try:
            with urllib.request.urlopen(row["url"], timeout=60) as response:
                blob = response.read()
        except (urllib.error.URLError, OSError) as exc:
            raise SystemExit(f"Could not read {row['url']}: {exc}")
        if hashlib.sha256(blob).hexdigest() != row.get("sha256"):
            fail(f"{target.slug}: the record does not match the sha256 the row carries")
        if blob != (DATA / target.mock_asset).read_bytes():
            fail(f"{target.slug}: the served record is not {target.mock_asset}")


def check_against(base_url: str, catalog: dict[str, dict]) -> None:
    """The service's catalog row against the two tables it was seeded from.

    Every field the app reads has to survive the move. A field that quietly
    arrives null is a protein page that draws one thing less than it used to,
    and nothing else in the suite would notice.
    """
    served = service_rows(base_url)
    for target in TARGETS:
        row = served.get(target.slug)
        if row is None:
            continue
        curated = catalog[target.slug]
        where = f"{target.slug}: served"
        source = target.source

        expected = {
            "gene": target.gene,
            "uniprot": target.uniprot,
            "accession": source.accession,
            "transcript_id": source.transcript_id,
            "protein_id": source.protein_id,
            "mature_peptides": target.mature_peptides,
            "display": curated["display"],
            "summary": curated["summary"],
            "chain": curated.get("chain"),
        }
        for key, want in expected.items():
            if row.get(key) != want:
                fail(f"{where} {key} is {row.get(key)!r}, table says {want!r}")

        if row.get("facts") != curated["facts"]:
            fail(f"{where} facts {row.get('facts')} != catalog {curated['facts']}")

        want_regions = partition(target)
        if row.get("regions") != want_regions:
            fail(f"{where} {len(row.get('regions') or [])} regions, "
                 f"partition() gives {len(want_regions)}")

        want_bonds = [list(pair) for pair in target.disulfides]
        if row.get("disulfides") != want_bonds:
            fail(f"{where} disulfides {row.get('disulfides')} != {want_bonds}")

        # The node names are the contract with structure_view.dart, and the
        # tint order is the order the record lists the mature peptides in.
        want_chains = [{"node": node, "tint": tint} for node, tint in curated["tints"]]
        if row.get("chains") != want_chains:
            fail(f"{where} chains {row.get('chains')} != catalog {want_chains}")

        chrome = row.get("structure") or {}
        if chrome.get("pdb") != target.structure.pdb:
            fail(f"{where} pdb {chrome.get('pdb')!r} != {target.structure.pdb!r}")
        want_modelled = list(curated["modelled"]) if curated["modelled"] else None
        if chrome.get("modelled") != want_modelled:
            fail(f"{where} modelled {chrome.get('modelled')} != {want_modelled}")
        if chrome.get("count") != curated["count"]:
            fail(f"{where} count {chrome.get('count')} != {curated['count']}")
        # The prose is compared whole, not only held non-empty: this file is where
        # it is edited now, and a seed that dropped an edit would otherwise pass.
        for key in ("label", "unit", "sentence", "semantics"):
            if not chrome.get(key):
                fail(f"{where} structure.{key} is empty")
            elif chrome.get(key) != curated["chrome"].get(key):
                fail(f"{where} structure.{key} is not the curated row's")
        if row.get("catalog_order") != curated["catalog_order"]:
            fail(f"{where} catalog_order {row.get('catalog_order')} != "
                 f"curated {curated['catalog_order']}")

        # The four booleans, as the four states that replaced them. Every one
        # of the twenty is baked, so a track that is neither ready nor absent
        # means the seed and the table disagree about what exists.
        states = row.get("tracks") or {}
        for kind, baked in (("constraint", target.scored),
                            ("impact", target.impact_scored),
                            ("clinvar", target.clinvar_available),
                            ("impact_explanations", curated["impact_explanations"])):
            if kind not in states:
                fail(f"{where} has no {kind} track state")
            elif states[kind] == "refused" and baked:
                fail(f"{where} {kind} is refused but the table says it is baked")

    if not problems:
        print(f"{len(served)} catalog rows on {base_url} agree with targets.py "
              f"and curated/catalog.json.")


# The standard code, which every protein here uses: no selenocysteine, no
# alternative start. `TCA` is serine, `TAA` a stop.
_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = {
    a + b + c: _AMINO[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}


def translate(bases: str) -> str:
    return "".join(CODONS.get(bases[i : i + 3], "X") for i in range(0, len(bases) - 2, 3))


def positions(mock: dict, segments: list[dict]) -> list[int]:
    """Every position of `segments`, in transcript order."""
    reverse = mock["location"].get("strand") == -1
    out: list[int] = []
    for segment in sorted(segments, key=lambda s: s["start"], reverse=reverse):
        span = range(segment["start"], segment["end"] + 1)
        out.extend(reversed(span) if reverse else span)
    return out


def bases(mock: dict, at: list[int]) -> str:
    """The record's own bases at `at`. A minus-strand `sequence` is already
    reverse-complemented (R2.1), so it is indexed from the far end there."""
    start, end = mock["location"]["start"], mock["location"]["end"]
    reverse = mock["location"].get("strand") == -1
    sequence = mock["sequence"]
    return "".join(sequence[end - p] if reverse else sequence[p - start] for p in at)


def check_record(mock: dict, where: str) -> None:
    """What the bake proves once, proved again from the baked file alone.

    The record's pieces are read from positions, and the page draws them from
    positions, so every piece is checked the way the page will read it:

    - the CDS translates to the protein (R2.2);
    - the signal peptide, the proprotein and every peptide are whole codons of
      one stretch of that CDS, and translate to the protein at their own
      offset. A neighbouring gene's feature — LTA's peptides beside TNF, the
      INS-IGF2 readthrough beside insulin — fails here, which is what re-proves
      the backend's exact `/gene` match for every record rather than for the
      one it has a test fixture of;
    - every intron starts and ends as an intron can: GT-AG, GC-AG or AT-AC
      (R2.3). Shortened introns keep their own ends, so this holds for
      dystrophin's too.
    """
    protein = mock.get("protein")
    if not protein:
        fail(f"{where}: the record has no protein")
        return
    translation = protein["translation"]
    cds = positions(mock, protein["segments"])
    if translate(bases(mock, cds)).rstrip("*") != translation:
        fail(f"{where}: the CDS does not translate to the record's protein")
        return
    index = {p: i for i, p in enumerate(cds)}

    pieces = [("the signal peptide", mock.get("signal_peptide")), ("the proprotein", mock.get("proprotein"))]
    pieces += [(p.get("product") or "a peptide", p) for p in mock["peptides"]]
    for name, piece in pieces:
        if not piece:
            continue
        at = positions(mock, piece["segments"])
        offsets = [index.get(p) for p in at]
        if None in offsets or offsets != list(range(offsets[0], offsets[0] + len(at))) or offsets[0] % 3 or len(at) % 3:
            fail(f"{where}: {name} is not whole codons of one stretch of the CDS")
            continue
        first = offsets[0] // 3
        expected = translation[first : first + len(at) // 3]
        if translate(bases(mock, at)) != expected or piece.get("translation") != expected:
            fail(f"{where}: {name} does not translate to residues {first + 1}-{first + len(at) // 3}")

    transcript = mock.get("transcript") or {}
    segments = transcript.get("segments") or mock["exons"]
    reverse = mock["location"].get("strand") == -1
    ordered = sorted(segments, key=lambda s: s["start"], reverse=reverse)
    for number, (before, after) in enumerate(zip(ordered, ordered[1:]), 1):
        if reverse:
            intron = list(range(before["start"] - 1, after["end"], -1))
        else:
            intron = list(range(before["end"] + 1, after["start"]))
        if len(intron) < 4:
            continue
        ends = (bases(mock, intron[:2]), bases(mock, intron[-2:]))
        if ends not in {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}:
            fail(f"{where}: intron {number} runs {ends[0]}…{ends[1]}, which no spliceosome cuts")


def check(target: Target, catalog: dict[str, dict]) -> None:
    where = target.slug

    mock_path = DATA / target.mock_asset
    constraint_path = DATA / target.constraint_asset
    impact_path = DATA / target.impact_asset
    model_path = DATA / target.structure_asset
    required = [mock_path, model_path]
    if target.scored:
        required.append(constraint_path)
    if target.impact_scored:
        required.append(impact_path)
    for path in required:
        if not path.exists():
            fail(f"{where}: missing {path.relative_to(DATA)}")
            return
    if not target.scored and constraint_path.exists():
        fail(
            f"{where}: {constraint_path.relative_to(DATA)} exists, but targets.py says this "
            "protein is not scored"
        )
    if not target.impact_scored and impact_path.exists():
        fail(
            f"{where}: {impact_path.relative_to(DATA)} exists, but targets.py says this "
            "gene has no impact track"
        )

    clinical_path = DATA / f"assets/clinvar/{target.slug}_clinvar.json"
    if target.clinvar_available != clinical_path.exists():
        fail(f"{where}: ClinVar asset availability differs from the table")
    mock = json.loads(mock_path.read_text())
    if target.clinvar_available and clinical_path.exists():
        check_clinvar(target, mock, json.loads(clinical_path.read_text()))

    if mock["gene"] != target.gene:
        fail(f"{where}: record is for {mock['gene']}, table says {target.gene}")

    span = mock["location"]["end"] - mock["location"]["start"] + 1
    if len(mock["sequence"]) != span:
        fail(f"{where}: {len(mock['sequence'])} bases of sequence for a {span} bp span")
    check_record(mock, where)

    if target.scored and not check_constraint(target, mock, json.loads(constraint_path.read_text())):
        return

    if target.impact_scored and not check_impact(
        target, mock, json.loads(impact_path.read_text())
    ):
        return

    nodes = glb_nodes(model_path)
    baked = {c.node for c in target.structure.chains} | (
        {"bonds"} if target.structure.bonds else set()
    )
    if nodes != baked:
        fail(f"{where}: {model_path.name} holds {sorted(nodes)}, the table says {sorted(baked)}")

    row = catalog.get(target.slug)
    if row is None:
        fail(f"{where}: no curated row")
        return
    for field, mine in (
        ("gene", target.gene),
        ("uniprot", target.uniprot),
        ("pdb", target.structure.pdb),
    ):
        if row.get(field) != mine:
            fail(f"{where}: the curated row says {field}={row.get(field)!r}, table says {mine!r}")
    if row["modelled"] != target.structure.residues:
        fail(
            f"{where}: the curated row says the model covers {row['modelled']}, "
            f"the table exports {target.structure.residues}"
        )
    if row["scored"] != target.scored:
        fail(f"{where}: the curated row says scored={row['scored']}, table says {target.scored}")
    if row["clinvar_available"] != target.clinvar_available:
        fail(f"{where}: ClinVar availability differs in the curated row and targets.py")
    if row["impact_scored"] != target.impact_scored:
        fail(
            f"{where}: the curated row says impact={row['impact_scored']}, "
            f"table says {target.impact_scored}"
        )
    if set(row["chains"]) != nodes:
        fail(
            f"{where}: the curated row paints {sorted(row['chains'])}, "
            f"{model_path.name} holds {sorted(nodes)}"
        )


def check_clinvar(target: Target, mock: dict, clinical: dict) -> None:
    """Independently rederive each allele and codon from the shipped record."""
    where = target.slug + " ClinVar"
    impact = json.loads((DATA / target.impact_asset).read_text())
    if any(clinical.get(k) != impact[k] for k in
           ("gene", "accession", "assembly", "chromosome", "start", "sequence", "runs", "complemented")):
        fail(f"{where}: mapping differs from the gene/AVI track")
    if clinical.get("protein_sequence") != mock["protein"]["translation"]:
        fail(f"{where}: protein differs")
    variants = clinical["variants"]
    if len({v["variation_id"] for v in variants}) != len(variants):
        fail(f"{where}: repeated Variation ID")
    if len(variants) + sum(clinical["excluded"].values()) != clinical["searched_records"]:
        fail(f"{where}: incomplete search")
    cds = positions(mock, mock["protein"]["segments"])
    offsets = {p: i for i, p in enumerate(cds)}
    dna = bases(mock, cds)
    complement = str.maketrans("ACGT", "TGCA")
    for v in variants:
        local = v["position"]
        if bases(mock, [local]) != v["ref"] or len(v["alt"]) != 1 or v["alt"] not in "ACGT" or v["alt"] == v["ref"]:
            fail(f"{where}: reference/allele mismatch {v['variation_id']}")
        mapped = [r["genomic"] + r["step"] * (local - r["local"]) for r in impact["runs"]
                  if r["local"] <= local < r["local"] + r["length"]]
        if mapped != [v["genomic"]]:
            fail(f"{where}: genomic mismatch {v['variation_id']}")
        for key in ("ref", "alt"):
            genomic_base = v[key].translate(complement) if impact["complemented"] else v[key]
            if genomic_base != v["genomic_" + key]:
                fail(f"{where}: strand mismatch {v['variation_id']}")
        offset = offsets.get(local)
        if offset is not None:
            codon = dna[offset // 3 * 3:offset // 3 * 3 + 3]
            ref = translate(codon)
            alt = translate(codon[:offset % 3] + v["alt"] + codon[offset % 3 + 1:])
            residue = offset // 3 + 1 if ref != "*" else None
            expected = f"p.{ref}{residue}{'=' if ref == alt else alt}" if residue else None
            if v["residue"] != residue or v["protein_change"] != expected:
                fail(f"{where}: protein mapping mismatch {v['variation_id']}")
        elif v["residue"] is not None:
            fail(f"{where}: noncoding variant has a residue")
        if not v["classification"] or not v["review_status"] or not v["accession"].startswith("VCV"):
            fail(f"{where}: missing classification provenance")
    # The identifiers ClinVar gives the conditions it names: every entry names a
    # condition some record's RCV cites, and nothing is written in its place.
    named = {n for v in variants for c in v["conditions"] for n in c["names"]}
    for name, ids in clinical.get("traits", {}).items():
        if name not in named:
            fail(f"{where}: identifiers for {name!r}, which no record names")
        if (not re.fullmatch(r"CN?\d+", ids.get("medgen", ""))
                or not re.fullmatch(r"(PS)?\d{6}", ids.get("omim", "000000"))
                or not re.fullmatch(r"MONDO:\d{7}", ids.get("mondo", "MONDO:0000000"))
                or set(ids) - {"medgen", "symbol", "omim", "mondo"}):
            fail(f"{where}: malformed identifiers for {name!r}")


def check_impact(target: Target, mock: dict, impact: dict) -> bool:
    """The per-base AVI track against the record it is filed under.

    The bake's own gates are the real ones — every score is checked against the
    reference base the Atlas returned with it. What is left for here is that the
    file on disk still belongs to this record: the same letters, the same
    coordinates, a map that covers all of them, and the biology the feature
    claims still holding.
    """
    where = target.slug
    for field, mine in (
        ("gene", target.gene),
        ("uniprot", target.uniprot),
        ("accession", target.source.accession),
        ("assembly", "GRCh38"),
        ("annotation", "GENCODE v46"),
        ("scorer", "AVI_SCORE"),
        ("score_units", "phred"),
    ):
        if impact.get(field) != mine:
            fail(f"{where}: impact track says {field}={impact.get(field)!r}, expected {mine!r}")
            return False

    # The drawn letters in increasing record position, which is how the app
    # reads them. The same as the record's own `sequence` except on a
    # minus-strand record, which stores its letters from the far end (R2.1).
    loc = mock["location"]
    if impact["sequence"] != bases(mock, list(range(loc["start"], loc["end"] + 1))):
        fail(f"{where}: the impact track's sequence is not the record's drawn letters")
        return False
    if impact["start"] != mock["location"]["start"]:
        fail(
            f"{where}: impact track starts at {impact['start']}, "
            f"the record at {mock['location']['start']}"
        )
        return False

    covered = sum(run["length"] for run in impact["runs"])
    if covered != len(mock["sequence"]):
        fail(
            f"{where}: the coordinate map covers {covered:,} of "
            f"{len(mock['sequence']):,} drawn bases"
        )
        return False

    # Every run has to stay inside the record and inside the chromosome.
    for run in impact["runs"]:
        last_local = run["local"] + run["length"] - 1
        last_genomic = run["genomic"] + run["step"] * (run["length"] - 1)
        if run["local"] < impact["start"] or last_local > mock["location"]["end"]:
            fail(f"{where}: a coordinate run leaves the record at {run['local']}")
            return False
        if min(run["genomic"], last_genomic) < 1:
            fail(f"{where}: a coordinate run leaves the chromosome at {run['genomic']}")
            return False

    start = impact["start"]
    for local, values in impact["positions"].items():
        offset = int(local) - start
        if offset < 0 or offset >= len(impact["sequence"]):
            fail(f"{where}: a score at {local} is outside the record")
            return False
        if len(values) != 3:
            fail(f"{where}: {len(values)} substitutions at {local}, expected three")
            return False

    scored = len(impact["positions"])
    if scored < len(mock["sequence"]) * 0.98:
        fail(f"{where}: only {scored:,} of {len(mock['sequence']):,} drawn bases are scored")
        return False

    # The claim the sheet makes, checked on the file rather than on the bake's
    # report of it: a splice boundary is not a quiet place and an intron
    # interior is.
    peak = {int(k): max(v) for k, v in impact["positions"].items()}
    exons = sorted((e["start"], e["end"]) for e in mock["exons"])
    exonic = [peak[p] for a, b in exons for p in range(a, b + 1) if p in peak]
    junction: list[float] = []
    interior: list[float] = []
    for (_, a), (b, _) in zip(exons, exons[1:]):
        low, high = a + 1, b - 1
        if high < low:
            continue
        edge = min(8, (high - low + 1) // 2)
        for p in range(low, high + 1):
            if p not in peak:
                continue
            (junction if p < low + edge or p > high - edge else interior).append(peak[p])
    if exonic and interior:
        if median(exonic) <= median(interior):
            fail(
                f"{where}: exons score no higher than intron interiors "
                f"({median(exonic):.1f} vs {median(interior):.1f})"
            )
            return False
        if junction and median(junction) <= median(interior):
            fail(
                f"{where}: splice boundaries score no higher than intron interiors "
                f"({median(junction):.1f} vs {median(interior):.1f})"
            )
            return False
    return True


def median(values: list[float]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def check_constraint(target: Target, mock: dict, constraint: dict) -> bool:
    """The track against the record it colours. False where nothing further is
    worth checking, because the two describe different proteins."""
    where = target.slug

    # The one that matters: the protein page is drawn from the record and
    # coloured from the track, and they are baked hours apart by different
    # tools. `AnatomyScreen` compares them too, and silently draws no colour.
    protein = mock["protein"]["translation"]
    if constraint["sequence"] != protein:
        fail(
            f"{where}: the constraint track is for a {len(constraint['sequence'])}-residue "
            f"protein, the record holds {len(protein)}. Re-run score_protein.py."
        )
        return False

    if len(constraint["positions"]) != len(protein):
        fail(f"{where}: {len(constraint['positions'])} positions for {len(protein)} residues")
    if constraint["uniprot"] != target.uniprot or constraint["gene"] != target.gene:
        fail(f"{where}: the track names {constraint['gene']}/{constraint['uniprot']}")

    expected_regions = partition(target)
    if constraint.get("regions") != expected_regions:
        fail(
            f"{where}: region table is stale. Re-run score_protein.py --metadata-only."
        )
    if constraint.get("disulfides") != [list(p) for p in target.disulfides]:
        fail(f"{where}: disulfide table is stale. Re-run --metadata-only.")
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check the baked assets against each other.")
    parser.add_argument(
        "--against", metavar="BASE_URL", default=None,
        help="also check the catalog rows a running service serves",
    )
    parser.add_argument(
        "--offline", action="store_true",
        help="skip the service check even when --against is given",
    )
    arguments = parser.parse_args()

    catalog = curated_catalog()

    if arguments.against and not arguments.offline:
        check_against(arguments.against, catalog)
        check_scenes(arguments.against)
        check_records(arguments.against)
        if problems:
            print(f"{len(problems)} problem(s):", file=sys.stderr)
            for problem in problems:
                print(f"  - {problem}", file=sys.stderr)
            raise SystemExit(1)
        raise SystemExit(0)

    extra = set(catalog) - {t.slug for t in TARGETS}
    if extra:
        fail(f"curated/catalog.json has rows with no bake: {sorted(extra)}")

    for target in TARGETS:
        check(target, catalog)
        attribution = DATA / f"assets/impact_explanations/{target.slug}.json"
        included = catalog.get(target.slug, {}).get("impact_explanations", False)
        if attribution.exists() != included:
            fail(f"{target.slug}: attribution asset and catalog availability disagree")
        elif included:
            try:
                validate_explanations((DATA / target.impact_asset).read_bytes(), json.loads(attribution.read_text()))
            except (AssertionError, KeyError, TypeError, ValueError) as error:
                fail(f"{target.slug}: invalid AVI explanations: {error}")

    # The compiled `.fsceneb` used to be checked here, out of
    # flutter_scene_generated/, because that was what shipped. Phase 5 fetches
    # it from storage instead and hook/build.dart no longer writes it, so there
    # is nothing local to open -- the same check now runs against the container
    # the service actually names, under --against. What stays offline is the
    # `.glb`: check() reads its node names and holds them to the catalog's, and
    # that is the bake's half of the contract.

    if problems:
        print(f"{len(problems)} problem(s):", file=sys.stderr)
        for problem in problems:
            print(f"  - {problem}", file=sys.stderr)
        raise SystemExit(1)

    total = sum(
        (DATA / path).stat().st_size
        for t in TARGETS
        for path in (
            [t.mock_asset]
            + ([t.constraint_asset] if t.scored else [])
            + ([t.impact_asset] if t.impact_scored else [])
            + ([f"assets/clinvar/{t.slug}_clinvar.json"] if t.clinvar_available else [])
            + ([f"assets/impact_explanations/{t.slug}.json"] if catalog[t.slug]["impact_explanations"] else [])
        )
    )
    # The models, as the bake left them. Not what a phone downloads -- the
    # `.fsceneb` compiled from each of these is about 29% of its size, and it
    # lives in storage now rather than anywhere this script can measure.
    models = sum((DATA / t.structure_asset).stat().st_size for t in TARGETS)
    scored = sum(1 for t in TARGETS if t.scored)
    tracked = sum(1 for t in TARGETS if t.impact_scored)
    print(
        f"{len(TARGETS)} targets check out, {scored} scored and "
        f"{tracked} with an impact track. "
        f"{total / 1e6:.2f} MB of records and tracks, "
        f"{models / 1e6:.2f} MB of baked models. None of it ships."
    )
