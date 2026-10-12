"""Reading a PDBx/mmCIF file: the format every PDB entry is archived in.

The twenty were read from legacy `.pdb` files by column (`pdb.py`). An entry
chosen for a protein built on demand is read from its mmCIF instead, because
the legacy format is neither guaranteed (large entries and extended ids have
none) nor complete: it cannot hold a chain id longer than one character, it
loses the label numbering that SIFTS maps from, and PDBe's SIFTS-enriched
file (`<id>_updated.cif`) exists only as mmCIF.

This reads the STAR syntax the archive writes and nothing more: one data
block, `loop_` tables and single items, bare and quoted values, and text
fields between lines that start with `;`. A bare `.` (inapplicable) and `?`
(unknown) are None; a quoted one is the character. It needs no package, so
it runs in every environment the pipeline has, the live worker's among them.
`test_mmcif.py` holds it to Biopython's `MMCIF2Dict` on whole entries.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

# A value on one line: a quoted string, which ends only at its quote followed
# by white space or the end of the line (so `'N'-terminal'` is one value), or
# anything else up to white space.
_VALUE = re.compile(r"""'(?:[^']|'(?=\S))*'(?=\s|$)|"(?:[^"]|"(?=\S))*"(?=\s|$)|\S+""")


class CifError(ValueError):
    """The file is not mmCIF as this reads it."""


@dataclass
class Table:
    """One category: its items' names, and one value of each per row."""

    name: str
    items: List[str] = field(default_factory=list)
    columns: Dict[str, list] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.columns[self.items[0]]) if self.items else 0

    def __contains__(self, item: str) -> bool:
        return item in self.columns

    def get(self, item: str) -> list:
        """An item's values, one a row; Nones where the category lacks it."""
        found = self.columns.get(item)
        return found if found is not None else [None] * len(self)

    def value(self, item: str) -> Optional[str]:
        """The single value of a category that has one row."""
        values = self.get(item)
        return values[0] if values else None

    def rows(self, *items: str) -> Iterator[tuple]:
        return zip(*(self.get(item) for item in items))


Block = Dict[str, Table]


def _tokens(text: str) -> Iterator[Tuple[str, Optional[str]]]:
    """(kind, value) for each token: kind is "value" or "word" (a bare token,
    which may be a keyword or a name), and a value of None is `.` or `?`."""
    lines = text.split("\n")
    i, count = 0, len(lines)
    while i < count:
        line = lines[i]
        if line.startswith(";"):
            # A text field: this line's rest, and every line to the one that
            # starts with `;` again.
            part = [line[1:]]
            i += 1
            while i < count and not lines[i].startswith(";"):
                part.append(lines[i])
                i += 1
            if i == count:
                raise CifError("a text field is never closed")
            # As written: the line break before the closing `;` belongs to
            # the delimiter, and nothing else is trimmed.
            yield "value", "\n".join(part)
            # What follows the closing `;` on its line is read as more tokens.
            line = lines[i][1:]
        for match in _VALUE.finditer(line):
            token = match.group(0)
            first = token[0]
            if first == "#":
                break                       # a comment runs to the end of the line
            if first in "'\"" and len(token) > 1 and token[-1] == first:
                yield "value", token[1:-1]
            elif token in (".", "?"):
                yield "value", None
            else:
                yield "word", token
        i += 1


def parse(text: str) -> Block:
    """The first data block of an mmCIF file, category by category.

    Categories are keyed by name without the underscore (`atom_site`), and
    their items by name without the category (`label_seq_id`).
    """
    block: Block = {}
    tokens = _tokens(text)
    started = False
    pending: List[Tuple[str, Optional[str]]] = []

    def take():
        if pending:
            return pending.pop()
        return next(tokens, None)

    while True:
        token = take()
        if token is None:
            break
        kind, value = token
        if kind == "word" and value.startswith("data_"):
            if started:
                break                       # only the first block
            started = True
            continue
        if not started:
            raise CifError("a value before the first data block")
        if kind == "word" and value == "loop_":
            names: List[str] = []
            while True:
                token = take()
                if token is None:
                    raise CifError("a loop with no items")
                if token[0] == "word" and token[1].startswith("_"):
                    names.append(token[1])
                    continue
                pending.append(token)
                break
            category, items = _category(names)
            table = block.setdefault(category, Table(category))
            if table.items:
                raise CifError(f"category {category} is written twice")
            table.items = items
            table.columns = {item: [] for item in items}
            columns = [table.columns[item] for item in items]
            width, position = len(items), 0
            while True:
                token = take()
                if token is None:
                    break
                kind, value = token
                if kind == "word" and (value.startswith("_") or value == "loop_"
                                       or value.startswith("data_")
                                       or value.startswith("save_")):
                    pending.append(token)
                    break
                columns[position].append(value)
                position = (position + 1) % width
            if position:
                raise CifError(f"loop {category} ends part-way through a row")
            continue
        if kind == "word" and value.startswith("_"):
            category, item = _split(value)
            token = take()
            if token is None:
                raise CifError(f"{value} has no value")
            if token[0] == "word" and (token[1].startswith("_") or token[1] == "loop_"):
                raise CifError(f"{value} has no value")
            table = block.setdefault(category, Table(category))
            if item in table.columns:
                raise CifError(f"{value} is written twice")
            table.items.append(item)
            table.columns[item] = [token[1]]
            continue
        raise CifError(f"an unexpected token {value!r}")
    if not started:
        raise CifError("no data block")
    return block


def _split(name: str) -> Tuple[str, str]:
    if not name.startswith("_") or "." not in name:
        raise CifError(f"{name!r} is not a category item")
    category, item = name[1:].split(".", 1)
    return category, item


def _category(names: Sequence[str]) -> Tuple[str, List[str]]:
    split = [_split(name) for name in names]
    categories = {category for category, _ in split}
    if len(categories) != 1:
        raise CifError(f"one loop holds several categories: {sorted(categories)}")
    return split[0][0], [item for _, item in split]
