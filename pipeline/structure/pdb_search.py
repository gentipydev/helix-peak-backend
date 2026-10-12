"""Finding the PDB's experimental entries of one protein, and what each says of itself.

A protein built on demand has nobody to choose an entry for it, so the fold
page's experimental structure is chosen by rule (`pdb_policy.py`) from what
this finds:

- **Discovery.** RCSB's search, for polymer entities SIFTS maps to the UniProt
  accession, experimental results only, every page, in entity-id order. A
  sequence search can add entities SIFTS has not mapped yet; `pdb_policy` holds
  each to the protein before it counts.
- **Metadata.** RCSB's data API, a batch at a time: each entity's sequence, its
  source and host organisms, its SIFTS alignment to UniProt, and for each of its
  chains the residues the entry never placed, its helices and strands, and the
  validation measured residue by residue (RSCC and RSRZ against X-ray density,
  Q-score against an EM map, clashes); and of its entry the method, resolution,
  R-free, models, release and assembly.
- **Files.** PDBe's SIFTS-enriched mmCIF of an entry (`<id>_updated.cif`), which
  maps every residue to UniProt, read by `mmcif.py`.

What a service says is believed only whole. A page that comes back short, a
total that moves between pages, an entity the data API does not know: each is
`Incomplete`, and a service that does not answer is `Unreachable`. Both are
tried again later by the worker; neither is ever read as "the PDB holds
nothing", which only an answer that says so is (RCSB's search answers an empty
query with 204 and no body).
"""

from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from typing import Callable, Dict, FrozenSet, List, Optional, Sequence, Tuple

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
UPDATED_CIF = "https://www.ebi.ac.uk/pdbe/entry-files/download/{entry}_updated.cif"

# Rows a search page asks for, and entities a metadata query asks about. A
# batch of a hundred is about 1.4 MB of JSON, most of it validation.
PAGE = 500
BATCH = 100

# RCSB's sequence search needs this many residues; shorter, it answers nothing.
SEQUENCE_MIN = 20

# Seconds waited before each try of a request: a blip is waited out here,
# rather than spending one of a bake's tries on it.
_TRIES = (0, 2, 6)
_TIMEOUT = 120

# Get answers a URL's bytes; Post answers a JSON body's status and bytes.
Get = Callable[[str], bytes]
Post = Callable[[str, dict], Tuple[int, bytes]]


class Unreachable(RuntimeError):
    """A service did not answer, after every try. Never a verdict."""


class Incomplete(RuntimeError):
    """An answer that cannot be the whole of it. Never a verdict."""


class Missing(LookupError):
    """A service answered that it holds no such thing (a 404)."""


class TooLarge(ValueError):
    """An entry's file is larger than the bound downloads are held to."""


# ------------------------------------------------------------ transport


def _request(request: urllib.request.Request) -> Tuple[int, bytes]:
    """One request's status and body. The server's own trouble (5xx, 429) or the
    network's is asked again, twice; any other answer, a 404 among them, is the
    answer."""
    last: Optional[BaseException] = None
    for wait in _TRIES:
        time.sleep(wait)
        try:
            with urllib.request.urlopen(request, timeout=_TIMEOUT) as response:
                return response.status, response.read()
        except urllib.error.HTTPError as exc:
            if exc.code < 500 and exc.code != 429:
                return exc.code, exc.read()
            last = exc
        except OSError as exc:      # URLError, a timeout, a reset connection
            last = exc
    raise Unreachable(f"{request.full_url} did not answer: {last}")


def get(url: str) -> bytes:
    """A URL's bytes; empty for a 204 (No Content), which InterPro answers for a
    protein it has nothing on; `Missing` for a 404."""
    status, body = _request(urllib.request.Request(
        url, headers={"User-Agent": "helix-peek-pipeline"}))
    if status == 204:
        return b""
    if status == 404:
        raise Missing(url)
    if status != 200:
        raise Unreachable(f"{url} answered {status}")
    return body


def post(url: str, body: dict) -> Tuple[int, bytes]:
    return _request(urllib.request.Request(
        url, data=json.dumps(body).encode(), method="POST",
        headers={"Content-Type": "application/json", "User-Agent": "helix-peek-pipeline"}))


# ------------------------------------------------------------ discovery


@dataclass(frozen=True)
class Discovery:
    """The entities one search found: the whole of them, in entity-id order."""

    query: str
    ids: Tuple[str, ...]
    release: str            # RCSB's weekly release the search ran against
    digest: str             # sha256 of the ids, one a line

    @property
    def snapshot(self) -> dict:
        return {"query": self.query, "release": self.release, "entities": len(self.ids),
                "sha256": self.digest}


def release_of(moment: datetime) -> str:
    """The RCSB weekly release a moment falls in: the archive updates on
    Wednesdays at 00:00 UTC, so its latest Wednesday."""
    day = moment.astimezone(timezone.utc).date()
    return (day - timedelta(days=(day.weekday() - 2) % 7)).isoformat()


def _digest(ids: Sequence[str]) -> str:
    return hashlib.sha256("".join(f"{i}\n" for i in ids).encode()).hexdigest()


def _accession_query(accession: str) -> dict:
    return {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": accession}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_name",
            "operator": "exact_match", "value": "UniProt"}},
    ]}


def _sequence_query(sequence: str, identity: float) -> dict:
    return {"type": "terminal", "service": "sequence", "parameters": {
        "evalue_cutoff": 0.1, "identity_cutoff": identity, "sequence_type": "protein",
        "value": sequence}}


def _paged(query: dict, label: str, post: Post, now: Callable[[], datetime]) -> Discovery:
    """Every page of a polymer-entity search, held to its own total."""
    ids: List[str] = []
    total: Optional[int] = None
    release = release_of(now())
    while True:
        body = {"query": query, "return_type": "polymer_entity",
                "request_options": {
                    "paginate": {"start": len(ids), "rows": PAGE},
                    "results_content_type": ["experimental"],
                    "sort": [{"sort_by": "rcsb_id", "direction": "asc"}]}}
        status, raw = post(SEARCH, body)
        if status == 204:
            # RCSB's answer for "nothing matches", and only that: a page past
            # the end of a non-empty result is not asked for.
            if ids or total:
                raise Incomplete(f"{label}: the search ran dry at {len(ids)} of {total}")
            break
        if status != 200:
            raise Unreachable(f"{label}: the search answered {status}")
        try:
            said = json.loads(raw)
            found = [row["identifier"] for row in said["result_set"]]
            count = int(said["total_count"])
        except (ValueError, KeyError, TypeError) as exc:
            raise Incomplete(f"{label}: the search answered something else: {exc}") from exc
        if total is None:
            total = count
        elif count != total:
            raise Incomplete(f"{label}: the total moved from {total} to {count} between pages")
        if not found:
            raise Incomplete(f"{label}: an empty page at {len(ids)} of {total}")
        ids += found
        if len(ids) >= total:
            break
        if release_of(now()) != release:
            raise Incomplete(f"{label}: the archive was updated part-way through the search")
    if total is not None and len(ids) != total:
        raise Incomplete(f"{label}: {len(ids)} entities for a total of {total}")
    if len(set(ids)) != len(ids):
        raise Incomplete(f"{label}: an entity came back twice")
    if ids != sorted(ids):
        raise Incomplete(f"{label}: the pages are not in entity-id order")
    return Discovery(label, tuple(ids), release, _digest(ids))


def search(accession: str, post: Post = post,
           now: Callable[[], datetime] = lambda: datetime.now(timezone.utc)) -> Discovery:
    """Every experimental polymer entity SIFTS maps to `accession`."""
    return _paged(_accession_query(accession), f"UniProt {accession}", post, now)


def sequence_search(sequence: str, identity: float, post: Post = post,
                    now: Callable[[], datetime] = lambda: datetime.now(timezone.utc)
                    ) -> Discovery:
    """Every experimental polymer entity whose sequence matches `sequence` at
    `identity` or better. Empty, without asking, for a sequence too short for
    the service."""
    label = f"sequence of {len(sequence)} residues"
    if len(sequence) < SEQUENCE_MIN:
        return Discovery(label, (), release_of(now()), _digest(()))
    return _paged(_sequence_query(sequence, identity), label, post, now)


# ------------------------------------------------------------ metadata


@dataclass(frozen=True)
class Alignment:
    """An entity's residues laid on a reference sequence: each region maps
    entity residues `entity_beg..+length-1` to reference `ref_beg..`."""

    database: str
    accession: str
    isoform: Optional[str]
    regions: Tuple[Tuple[int, int, int], ...]
    provenance: Optional[str]


@dataclass(frozen=True)
class Chain:
    """One instance of an entity: a chain of the entry, by its label and
    author ids, and what the entry and its validation say residue by residue
    (numbered by label_seq_id)."""

    asym: str
    auth: str
    unobserved: FrozenSet[int]
    helices: Tuple[Tuple[int, int], ...]
    strands: Tuple[Tuple[int, int], ...]
    rscc: Dict[int, float] = field(default_factory=dict)
    rsrz: Dict[int, float] = field(default_factory=dict)
    q_score: Dict[int, float] = field(default_factory=dict)
    clashes: FrozenSet[int] = frozenset()


@dataclass(frozen=True)
class Entry:
    """What an entry says of itself, as RCSB's data API carries it."""

    entry: str
    methods: Tuple[str, ...]
    resolution: Optional[float]
    rfree: Optional[float]
    models: Optional[int]
    nmr_conformer: Optional[str]
    clashscore: Optional[float]
    ramachandran_outliers: Optional[float]     # per cent of residues
    rotamer_outliers: Optional[float]          # per cent of residues
    released: Optional[str]
    revised: Optional[str]
    revision: Optional[str]
    status: Optional[str]
    polymer_entities: Optional[int]
    nonpolymer_entities: Optional[int]
    oligomeric_count: Optional[int]


@dataclass(frozen=True)
class Entity:
    """One polymer entity, its chains, and its entry."""

    id: str                     # "1BM7_1"
    entry: Entry
    entity: str                 # "1"
    description: Optional[str]
    polymer: Optional[str]
    sequence: str               # canonical one-letter code: a modified residue as its parent
    mutations: Optional[int]
    sources: Tuple[int, ...]
    hosts: Tuple[int, ...]
    alignments: Tuple[Alignment, ...]
    chains: Tuple[Chain, ...]
    # Copies of this entity in the entry's first assembly (its chains times the
    # assembly's operators), and whether that assembly holds no other polymer.
    copies: Optional[int] = None
    alone: Optional[bool] = None
    # That assembly's point group (`C2`, `D2`, ...; `C1` is none): several
    # copies with no symmetry are molecules in one crystal, not an oligomer.
    symmetry: Optional[str] = None


_QUERY = """
query($ids: [String!]!) { polymer_entities(entity_ids: $ids) {
  rcsb_id
  entity_poly { rcsb_entity_polymer_type rcsb_mutation_count pdbx_seq_one_letter_code_can }
  rcsb_polymer_entity { pdbx_description }
  rcsb_entity_source_organism { ncbi_taxonomy_id }
  rcsb_entity_host_organism { ncbi_taxonomy_id }
  rcsb_polymer_entity_container_identifiers { entry_id entity_id }
  rcsb_polymer_entity_align { provenance_source reference_database_name
    reference_database_accession reference_database_isoform
    aligned_regions { entity_beg_seq_id ref_beg_seq_id length } }
  polymer_entity_instances {
    rcsb_polymer_entity_instance_container_identifiers { asym_id auth_asym_id }
    rcsb_polymer_instance_feature { type provenance_source
      feature_positions { beg_seq_id end_seq_id values } }
  }
  entry {
    rcsb_id
    exptl { method }
    rcsb_entry_info { resolution_combined deposited_model_count polymer_entity_count
                      nonpolymer_entity_count }
    refine { ls_R_factor_R_free }
    pdbx_nmr_representative { conformer_id }
    pdbx_vrpt_summary_geometry { clashscore percent_ramachandran_outliers
                                 percent_rotamer_outliers }
    rcsb_accession_info { initial_release_date revision_date major_revision minor_revision
                          status_code }
    assemblies { rcsb_id pdbx_struct_assembly { oligomeric_count }
                 pdbx_struct_assembly_gen { asym_id_list oper_expression }
                 rcsb_struct_symmetry { kind symbol } }
    polymer_entities { rcsb_polymer_entity_container_identifiers { asym_ids } }
  }
} }
"""


def _ranges(feature: dict) -> List[Tuple[int, int]]:
    out = []
    for position in feature.get("feature_positions") or []:
        first = position.get("beg_seq_id")
        last = position.get("end_seq_id")
        if isinstance(first, int):
            out.append((first, last if isinstance(last, int) else first))
    return out


def _values(feature: dict) -> Dict[int, float]:
    """A feature's per-residue values: each position starts at `beg_seq_id`
    and runs one residue a value."""
    out: Dict[int, float] = {}
    for position in feature.get("feature_positions") or []:
        first = position.get("beg_seq_id")
        if not isinstance(first, int):
            continue
        for offset, value in enumerate(position.get("values") or []):
            if isinstance(value, (int, float)):
                out[first + offset] = float(value)
    return out


def _chain(instance: dict) -> Chain:
    ids = instance.get("rcsb_polymer_entity_instance_container_identifiers") or {}
    unobserved, helices, strands = set(), [], []
    rscc: Dict[int, float] = {}
    rsrz: Dict[int, float] = {}
    q_score: Dict[int, float] = {}
    clashes = set()
    for feature in instance.get("rcsb_polymer_instance_feature") or []:
        kind = feature.get("type")
        if kind == "UNOBSERVED_RESIDUE_XYZ":
            for first, last in _ranges(feature):
                unobserved.update(range(first, last + 1))
        elif kind == "HELIX_P":
            helices += _ranges(feature)
        elif kind == "SHEET":
            strands += _ranges(feature)
        elif kind == "RSCC":
            rscc.update(_values(feature))
        elif kind == "RSRZ":
            rsrz.update(_values(feature))
        elif kind == "Q_SCORE":
            q_score.update(_values(feature))
        elif kind == "CLASHES":
            for first, last in _ranges(feature):
                clashes.update(range(first, last + 1))
    return Chain(
        asym=ids.get("asym_id") or "", auth=ids.get("auth_asym_id") or "",
        unobserved=frozenset(unobserved), helices=tuple(sorted(helices)),
        strands=tuple(sorted(strands)), rscc=rscc, rsrz=rsrz, q_score=q_score,
        clashes=frozenset(clashes))


def _number(value) -> Optional[float]:
    return float(value) if isinstance(value, (int, float)) else None


def _entry(raw: dict) -> Entry:
    info = raw.get("rcsb_entry_info") or {}
    accession = raw.get("rcsb_accession_info") or {}
    geometry = (raw.get("pdbx_vrpt_summary_geometry") or [{}])[0] or {}
    resolution = [r for r in info.get("resolution_combined") or [] if isinstance(r, (int, float))]
    rfree = [r.get("ls_R_factor_R_free") for r in raw.get("refine") or []
             if isinstance((r or {}).get("ls_R_factor_R_free"), (int, float))]
    nmr = raw.get("pdbx_nmr_representative") or {}
    if isinstance(nmr, list):
        nmr = nmr[0] if nmr else {}
    oligomers = None
    for assembly in raw.get("assemblies") or []:
        # The first assembly is the one the depositors (or their software)
        # name as what the molecule is in life; RCSB lists it as `-1`.
        if (assembly or {}).get("rcsb_id", "").endswith("-1"):
            count = ((assembly.get("pdbx_struct_assembly") or {}).get("oligomeric_count"))
            oligomers = count if isinstance(count, int) else None
    major, minor = accession.get("major_revision"), accession.get("minor_revision")
    return Entry(
        entry=raw.get("rcsb_id") or "",
        methods=tuple(sorted({(e or {}).get("method") or "" for e in raw.get("exptl") or []} - {""})),
        resolution=min(resolution) if resolution else None,
        rfree=min(rfree) if rfree else None,
        models=info.get("deposited_model_count"),
        nmr_conformer=nmr.get("conformer_id"),
        clashscore=_number(geometry.get("clashscore")),
        ramachandran_outliers=_number(geometry.get("percent_ramachandran_outliers")),
        rotamer_outliers=_number(geometry.get("percent_rotamer_outliers")),
        released=(accession.get("initial_release_date") or "")[:10] or None,
        revised=(accession.get("revision_date") or "")[:10] or None,
        revision=f"{major}.{minor}" if major is not None and minor is not None else None,
        status=accession.get("status_code"),
        polymer_entities=info.get("polymer_entity_count"),
        nonpolymer_entities=info.get("nonpolymer_entity_count"),
        oligomeric_count=oligomers,
    )


def operators(expression: Optional[str]) -> int:
    """How many copies an assembly's operator expression makes: `1` one, `1,2`
    two, `(1-60)` sixty, `(1-5)(6-10)` twenty-five (each group applied after
    the next)."""
    if not expression:
        return 1
    groups = re.findall(r"\(([^)]*)\)", expression) or [expression]
    total = 1
    for group in groups:
        count = 0
        for item in group.split(","):
            ends = item.strip().split("-")
            if len(ends) == 2 and ends[0].strip().isdigit() and ends[1].strip().isdigit():
                count += int(ends[1]) - int(ends[0]) + 1
            elif item.strip():
                count += 1
        total *= max(count, 1)
    return total


def _assembly(raw: dict, asyms: Sequence[str]
              ) -> Tuple[Optional[int], Optional[bool], Optional[str]]:
    """Copies of the entity with chains `asyms` in the entry's first assembly,
    whether no other polymer is in it, and its point group."""
    entry = raw.get("entry") or {}
    polymers = {a for e in entry.get("polymer_entities") or []
                for a in ((e or {}).get("rcsb_polymer_entity_container_identifiers") or {})
                .get("asym_ids") or []}
    for assembly in entry.get("assemblies") or []:
        if not (assembly or {}).get("rcsb_id", "").endswith("-1"):
            continue
        copies, others = 0, False
        for gen in assembly.get("pdbx_struct_assembly_gen") or []:
            listed = set(gen.get("asym_id_list") or [])
            times = operators(gen.get("oper_expression"))
            copies += len(listed & set(asyms)) * times
            others = others or bool((listed & polymers) - set(asyms))
        symbol = next(((sym or {}).get("symbol") for sym in assembly.get("rcsb_struct_symmetry") or []
                       if (sym or {}).get("kind") == "Global Symmetry"), None)
        return (copies or None), (not others if copies else None), symbol
    return None, None, None


def _entity(raw: dict) -> Entity:
    ids = raw.get("rcsb_polymer_entity_container_identifiers") or {}
    poly = raw.get("entity_poly") or {}
    alignments = []
    for align in raw.get("rcsb_polymer_entity_align") or []:
        regions = tuple(
            (r["entity_beg_seq_id"], r["ref_beg_seq_id"], r["length"])
            for r in align.get("aligned_regions") or []
            if all(isinstance(r.get(k), int) for k in
                   ("entity_beg_seq_id", "ref_beg_seq_id", "length")))
        alignments.append(Alignment(
            database=align.get("reference_database_name") or "",
            accession=align.get("reference_database_accession") or "",
            isoform=align.get("reference_database_isoform"),
            regions=regions, provenance=align.get("provenance_source")))
    chains = tuple(sorted((_chain(i) for i in raw.get("polymer_entity_instances") or []),
                          key=lambda c: c.asym))
    copies, alone, symmetry = _assembly(raw, [c.asym for c in chains])
    return Entity(
        id=raw["rcsb_id"], entry=_entry(raw.get("entry") or {}),
        entity=str(ids.get("entity_id") or raw["rcsb_id"].split("_")[-1]),
        description=(raw.get("rcsb_polymer_entity") or {}).get("pdbx_description"),
        polymer=poly.get("rcsb_entity_polymer_type"),
        sequence="".join((poly.get("pdbx_seq_one_letter_code_can") or "").split()),
        mutations=poly.get("rcsb_mutation_count"),
        sources=tuple(sorted({o.get("ncbi_taxonomy_id") for o in
                              raw.get("rcsb_entity_source_organism") or []
                              if isinstance((o or {}).get("ncbi_taxonomy_id"), int)})),
        hosts=tuple(sorted({o.get("ncbi_taxonomy_id") for o in
                            raw.get("rcsb_entity_host_organism") or []
                            if isinstance((o or {}).get("ncbi_taxonomy_id"), int)})),
        alignments=tuple(alignments),
        chains=chains,
        copies=copies,
        alone=alone,
        symmetry=symmetry,
    )


def metadata(ids: Sequence[str], post: Post = post, batch: int = BATCH) -> List[Entity]:
    """What RCSB's data API says of each entity, in the order asked."""
    out: List[Entity] = []
    for start in range(0, len(ids), batch):
        asked = list(ids[start:start + batch])
        status, raw = post(GRAPHQL, {"query": _QUERY, "variables": {"ids": asked}})
        if status != 200:
            raise Unreachable(f"RCSB's data API answered {status}")
        try:
            said = json.loads(raw)
        except ValueError as exc:
            raise Incomplete(f"RCSB's data API answered something else: {exc}") from exc
        if said.get("errors"):
            raise Incomplete(f"RCSB's data API: {json.dumps(said['errors'])[:300]}")
        found = {e["rcsb_id"]: e for e in ((said.get("data") or {}).get("polymer_entities") or [])
                 if isinstance(e, dict) and e.get("rcsb_id")}
        lost = [i for i in asked if i not in found]
        if lost:
            raise Incomplete(f"RCSB's data API does not know {lost[:5]} of {len(lost)}")
        out += [_entity(found[i]) for i in asked]
    return out


# ------------------------------------------------------------ files


# The largest entry file read. A whole ribosome is several times this; a
# protein's own entries are well under it.
DOWNLOAD_LIMIT = 40_000_000


def get_bounded(url: str, limit: int = DOWNLOAD_LIMIT) -> bytes:
    """`get`, refusing a body larger than `limit` (`TooLarge`) without reading it all."""
    last: Optional[BaseException] = None
    request = urllib.request.Request(url, headers={"User-Agent": "helix-peek-pipeline"})
    for wait in _TRIES:
        time.sleep(wait)
        try:
            with urllib.request.urlopen(request, timeout=_TIMEOUT) as response:
                body = response.read(limit + 1)
                if len(body) > limit:
                    raise TooLarge(f"{url} is larger than {limit:,} bytes")
                return body
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raise Missing(url) from exc
            if exc.code < 500 and exc.code != 429:
                raise Unreachable(f"{url} answered {exc.code}") from exc
            last = exc
        except OSError as exc:
            last = exc
    raise Unreachable(f"{url} did not answer: {last}")


def updated_cif(entry: str, get: Get = get_bounded) -> bytes:
    """PDBe's SIFTS-enriched mmCIF of an entry. `Missing` where PDBe has none,
    `TooLarge` past `DOWNLOAD_LIMIT`."""
    return get(UPDATED_CIF.format(entry=entry.lower()))


# ------------------------------------------------------------ UniProt


@dataclass(frozen=True)
class Domain:
    """A named stretch annotated as one folded unit, 1-based inclusive: by
    UniProt, or by Pfam (through InterPro) where UniProt names none."""

    kind: str           # Domain, Repeat, Zinc finger
    name: str
    start: int
    end: int
    source: str = "UniProt"


DOMAIN_KINDS = ("Domain", "Repeat", "Zinc finger")

INTERPRO_PFAM = ("https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{accession}/"
                 "?page_size=200")
# Pfam's entry types that are folded units; a family, a conserved site or a
# disordered region is not one.
PFAM_KINDS = {"domain": "Domain", "repeat": "Repeat"}


def domains_of(body: dict) -> List[Domain]:
    """A UniProt entry's domains, repeats and zinc fingers whose ends are exact."""
    length = len((body.get("sequence") or {}).get("value") or "")
    out = []
    for feature in body.get("features") or []:
        kind = feature.get("type")
        if kind not in DOMAIN_KINDS:
            continue
        location = feature.get("location") or {}
        start, end = location.get("start") or {}, location.get("end") or {}
        first, last = start.get("value"), end.get("value")
        if (start.get("modifier") != "EXACT" or end.get("modifier") != "EXACT"
                or not isinstance(first, int) or not isinstance(last, int)
                or not 1 <= first <= last <= length):
            continue
        out.append(Domain(kind, (feature.get("description") or kind).strip(), first, last))
    return sorted(out, key=lambda d: (d.start, d.end, d.name))


def pfam_domains(accession: str, get: Get = get) -> List[Domain]:
    """Pfam's domains and repeats on a UniProt entry, as InterPro places them.
    Each occurrence is one, from its first fragment to its last."""
    raw = get(INTERPRO_PFAM.format(accession=accession))
    if not raw:
        return []
    try:
        said = json.loads(raw)
        results = said["results"]
    except (ValueError, KeyError, TypeError) as exc:
        raise Incomplete(f"InterPro answered something else: {exc}") from exc
    if said.get("next"):
        raise Incomplete(f"InterPro has more than one page of Pfam entries for {accession}")
    out = []
    for result in results:
        meta = result.get("metadata") or {}
        kind = PFAM_KINDS.get(meta.get("type"))
        if kind is None:
            continue
        for protein in result.get("proteins") or []:
            for location in protein.get("entry_protein_locations") or []:
                fragments = [(f.get("start"), f.get("end")) for f in location.get("fragments") or []
                             if isinstance(f.get("start"), int) and isinstance(f.get("end"), int)]
                if fragments:
                    out.append(Domain(kind, meta.get("name") or meta.get("accession") or kind,
                                      min(a for a, _ in fragments), max(b for _, b in fragments),
                                      f"Pfam {meta.get('accession')}"))
    return sorted(out, key=lambda d: (d.start, d.end, d.name))
