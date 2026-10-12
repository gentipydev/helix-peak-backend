"""The fold page's words for a protein drawn from an experimental entry.

The twenty's are written by hand (`curated/catalog.json`); a protein built on
demand has its seven templated from facts (R4.5), as `alphafold.describe`
templates a predicted fold's, and nothing here knows which protein it is:

- `pdb`: the entry id. The header names it beside the span drawn
  (`PDB 1BM7 · residues 21–147 of 147`), so the sentence does not repeat it.
- `label`: what is drawn, structurally: the fold, the subunit (one chain of
  an assembly), the mature chains, or the annotated domain.
- `count`: the residues the chains drawn span, the biological span; and the
  sentence then says how many of them the entry places, which can be fewer.
- `modelled`: the span, where one chain is drawn and it is not the whole
  mature product.
- `sentence`: the method and resolution, then the placed residues or the
  bridges between chains, in two lines at most (80 characters).
- `semantics`: the same facts in full for a screen reader, and what the entry
  holds that is not drawn (partner chains, ligands).

Residues the entry does not place are "not placed", never "disordered": an
unplaced residue is a fact about the experiment, not proof of disorder.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from pipeline.structure.normalize import Examined
from pipeline.structure.pdb_policy import Candidate, CRYSTAL, EM, Inputs, NMR
from pipeline.structure.pdb_search import Entity, Entry

LINE = 80

_SHORT = {
    "X-RAY DIFFRACTION": "X-ray",
    "NEUTRON DIFFRACTION": "neutron",
    "ELECTRON CRYSTALLOGRAPHY": "electron crystallography",
    "ELECTRON MICROSCOPY": "cryo-EM",
    "SOLUTION NMR": "NMR",
    "SOLID-STATE NMR": "solid-state NMR",
}
_LONG = {
    "X-RAY DIFFRACTION": "X-ray crystallography",
    "NEUTRON DIFFRACTION": "neutron crystallography",
    "ELECTRON CRYSTALLOGRAPHY": "electron crystallography",
    "ELECTRON MICROSCOPY": "cryo-electron microscopy",
    "SOLUTION NMR": "solution NMR",
    "SOLID-STATE NMR": "solid-state NMR",
}
_WORDS = ("no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten",
          "eleven", "twelve")


def spelled(n: int) -> str:
    return _WORDS[n] if 0 <= n < len(_WORDS) else f"{n:,}"


def _capital(text: str) -> str:
    return text[:1].upper() + text[1:]


def method(methods: Sequence[str], resolution: Optional[float], long: bool = False) -> str:
    """`X-ray, 1.70 Å`, `cryo-EM, 3.1 Å`, `NMR`: the methods, and the resolution
    where the method has one."""
    names = [(_LONG if long else _SHORT).get(m, m.lower()) for m in methods]
    said = " and ".join(names)
    kinds = set(methods)
    if resolution is None or kinds <= NMR:
        return said
    places = 1 if kinds & EM and not kinds & CRYSTAL else 2
    return f"{said}, {resolution:.{places}f} Å"


def domain_phrase(name: str, kind: str) -> str:
    """`the VWFA domain`, `the protein kinase domain`, `the ANK 3 repeat`.

    UniProt writes a domain's name capitalised. Its first word is lowered only
    where it is a plain word (`Protein` kinase), never an abbreviation or a
    name (`VWFA`, `Ig-like`, `BRCT 1`, `RING-type`)."""
    words = name.split()
    first = words[0] if words else name
    if len(first) >= 3 and first.isalpha() and first[0].isupper() and first[1:].islower():
        words[0] = first.lower()
    named = " ".join(words) or name
    noun = {"Domain": "domain", "Repeat": "repeat", "Zinc finger": "zinc finger"}.get(kind, "domain")
    if named.lower().endswith(noun):
        return f"the {named}"
    return f"the {named} {noun}"


def domains_phrase(names: Sequence[str], kinds: Dict[str, str]) -> str:
    """The annotated domains a stretch holds, as the fold page's label: one by
    name, two by name, more of one kind by number, and a mixture of kinds as
    the part of the protein it is."""
    nouns = {"Domain": "domain", "Repeat": "repeat", "Zinc finger": "zinc finger"}
    if len(names) == 1:
        return domain_phrase(names[0], kinds.get(names[0], "Domain"))
    kind = {kinds.get(n, "Domain") for n in names}
    if len(kind) > 1 and len(names) > 2:
        # Repeats and domains of several kinds: a count of them says little.
        return "part of the protein"
    if len(kind) > 1:
        # Two, of two kinds: each by its own name, where the two fit.
        both = " and ".join(domain_phrase(n, kinds.get(n, "Domain")) for n in names)
        return both if len(both) <= 40 else "part of the protein"
    noun = nouns.get(kind.pop(), "domain")
    if len(names) == 2:
        bare = [domain_phrase(n, kinds.get(n, "Domain"))[4:] for n in names]
        bare = [b[: -len(noun) - 1] if b.endswith(" " + noun) else b for b in bare]
        joined = f"the {bare[0]} and {bare[1]} {noun}s"
        if len(joined) <= 40:
            return joined
    return f"{spelled(len(names))} of its {noun}s"


def tint_of(product: int, drawn: int) -> str:
    """A chain's tint. A lone chain is `mature1`, as every single-chain model of
    the twenty is. Chains drawn together take their product's place among the
    kept regions, so that the legend can name each from that same list
    (`mature1` is the first kept region)."""
    if drawn == 1:
        return "mature1"
    if product < 3:
        return f"mature{product + 1}"
    raise ValueError(f"product {product + 1} has no colour of its own")


def _fit(*forms: str) -> str:
    """The first of `forms` that fits two lines on a phone."""
    for form in forms:
        if len(form) <= LINE:
            return form
    return forms[-1]


def _bridges(examined: Examined) -> Tuple[int, int]:
    """Disulfide bridges drawn: between two chains, and inside one."""
    between = sum(1 for b in examined.bridges if b.a_node != b.b_node)
    return between, len(examined.bridges) - between


def describe(inputs: Inputs, candidate: Candidate, examined: Examined, entry: Entry,
             entity: Optional[Entity] = None) -> Tuple[dict, List[dict]]:
    """The fold page's seven words and its chains, from the entry's own numbers.

    `entity` is the one chain's entity, which says whether the entry's
    assembly is copies of this protein alone, arranged with a symmetry (a
    homo-oligomer): only then is the chain "the subunit", one of so many. Two
    molecules that merely share a crystal have none. Partner chains, which may
    be a receptor or an antibody fragment the experiment needed, are named
    only as not drawn."""
    chains = examined.chains
    placed, listed = examined.placed, examined.spanned
    how = method(examined.methods or tuple(candidate.method.split(", ")), candidate.resolution)
    how_long = method(examined.methods or tuple(candidate.method.split(", ")),
                      candidate.resolution, long=True)
    between, inside = _bridges(examined)
    copies = entity.copies if (entity is not None and entity.alone and entity.copies
                               and entity.symmetry not in (None, "C1")) else 1
    kinds = {d.name: d.kind for d in inputs.domains}

    products = len(inputs.products)
    # One of a precursor's several products, drawn alone: a hormone among the
    # hormones its precursor is cut into, which the words then name.
    one_of = candidate.scope == "product" and len(chains) == 1 and products > 1
    named = inputs.products[chains[0].product].label
    if candidate.scope == "domain":
        label = domains_phrase(list(candidate.domains), kinds)
    elif len(chains) > 1:
        label = ("the mature chains" if len(chains) == products
                 else f"{spelled(len(chains))} of its chains")
    elif one_of:
        label = "one of its products"
    elif copies > 1:
        label = "the subunit"
    else:
        label = "the fold"

    every = placed == listed
    counted = ("every residue placed" if every
               else f"{placed:,} of these {listed:,} residues placed")
    if len(chains) > 1:
        joined = (f"joined by {spelled(between)} disulfide bridge{'s' if between != 1 else ''}"
                  if between else "as the entry holds them")
        sentence = _fit(f"{how}. {_capital(spelled(len(chains)))} chains, {joined}.",
                        f"{how}. {_capital(spelled(len(chains)))} chains; {counted}.",
                        f"{_capital(spelled(len(chains)))} chains, {joined}.")
    elif candidate.scope == "domain":
        sentence = _fit(f"{how}. Only this part is drawn; {counted}.",
                        f"{how}. This part only; {counted}.",
                        f"{how}. {_capital(counted)}.")
    elif one_of:
        brief = "every residue placed" if every else f"{placed:,} of {listed:,} placed"
        sentence = _fit(f"{how}. {_capital(named)}, one of its {spelled(products)} products; "
                        f"{counted}.",
                        f"{how}. {_capital(named)}, one of {spelled(products)} products; {brief}.",
                        f"{how}. {_capital(named)}, one of {spelled(products)} products.",
                        f"{how}. One of its {spelled(products)} products; {counted}.",
                        f"{how}. {_capital(counted)}.")
    elif copies > 1:
        sentence = _fit(f"{how}. One chain of {spelled(copies)}; {counted}.",
                        f"{how}. One of {spelled(copies)} chains; {placed:,} of {listed:,} "
                        f"placed.",
                        f"{how}. {_capital(counted)}.")
    else:
        sentence = _fit(f"{how}. {_capital(counted)}.",
                        f"{how}. {placed:,} of {listed:,} residues placed.")

    # The screen reader's version: every fact, and what is not drawn.
    if len(chains) > 1:
        what = (f"{spelled(len(chains))} chains, "
                + " and ".join(f"{inputs.products[c.product].label} (residues {c.first:,} to "
                               f"{c.last:,})" for c in chains))
    else:
        chain = chains[0]
        what = (f"{label}, residues {chain.first:,} to {chain.last:,}"
                if candidate.scope == "domain"
                else f"{named}, one of its {spelled(products)} products, residues "
                     f"{chain.first:,} to {chain.last:,}" if one_of
                else f"{f'one of its {spelled(copies)} identical chains' if copies > 1 else 'its mature chain'}, "
                     f"residues {chain.first:,} to {chain.last:,}")
    parts = [f"The {inputs.display} structure in PDB {examined.entry}, by {how_long}: {what}. "
             f"{'Every residue is placed' if every else f'{placed:,} of its {listed:,} residues are placed'}."]
    if examined.bridges:
        if len(chains) > 1:
            said = []
            if between:
                said.append(f"{_capital(spelled(between))} disulfide bridge"
                            f"{'s join' if between != 1 else ' joins'} the chains")
            if inside:
                said.append(f"{spelled(inside)} more fold{'s' if inside == 1 else ''} a chain"
                            if between else
                            f"{_capital(spelled(inside))} fold{'s' if inside == 1 else ''} a chain")
            parts.append("; ".join(said) + ".")
        else:
            n = len(examined.bridges)
            parts.append("Its disulfide bridge is drawn." if n == 1
                         else f"Its {spelled(n)} disulfide bridges are drawn.")
    drawn_entities = {c.entity for c in chains}
    partners = (entry.polymer_entities or len(drawn_entities)) > len(drawn_entities)
    ligands = (entry.nonpolymer_entities or 0) > 0
    if partners and ligands:
        parts.append("The entry's other chains and its ligands are not drawn.")
    elif partners:
        parts.append("The entry's other chains are not drawn.")
    elif ligands:
        parts.append("Its ligands are not drawn.")
    parts.append("Drag to turn it.")

    first, last = (chains[0].first, chains[0].last) if len(chains) == 1 else (None, None)
    product = inputs.products[chains[0].product]
    whole = len(chains) > 1 or (first, last) == (product.start, product.end)
    chrome = {
        "pdb": examined.entry,
        "modelled": None if whole else [first, last],
        "label": label,
        "count": listed,
        "unit": "residues",
        "sentence": sentence,
        "semantics": " ".join(parts),
    }
    rows = [{"node": c.node, "tint": tint_of(c.product, len(chains))} for c in chains]
    if examined.bridges:
        rows.append({"node": "bonds", "tint": "cysteine"})
    return chrome, rows


def words(inputs: Inputs, candidate: Candidate, examined: Examined, entry: Entry,
          entity: Optional[Entity] = None) -> Dict[str, object]:
    """`describe`, as the report shows it."""
    chrome, chains = describe(inputs, candidate, examined, entry, entity)
    return {"chrome": chrome, "chains": chains}
