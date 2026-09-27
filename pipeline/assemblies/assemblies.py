"""Assemblies: molecules made of more than one gene's chains, kept out of the catalog.

The catalog is twenty proteins, one gene each, and the walk draws only what a
gene makes (R5.1): the hemoglobin row is HBB's beta chain alone. A lab
feature about oxygen binding needs the whole tetramer, two alpha chains from
HBA1 and two beta chains from HBB, so it cannot be a catalog row. It is a row
of this table instead, `ASSEMBLIES`, which nothing that seeds or serves the
catalog reads: `targets.py` stays twenty and `/catalog` returns twenty.

An assembly is baked as a morph pair: the same molecule in two states, each
from its own PDB entry, brought into one frame so that the one can be moved
into the other.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Subunit:
    """One chain of the assembly, by the node it is drawn as."""

    node: str
    gene: str
    uniprot: str
    # The entry's own chain letter for this subunit, in every state's built
    # tetramer (see `State.biomt`).
    pdb_chain: str
    # The entries number each chain from its first residue after the initiator
    # methionine (Val1), which the precursor numbers Val2. The bake checks it:
    # the beta chains against the catalog's own HBB record.
    offset: int = 1


@dataclass(frozen=True)
class State:
    """One PDB entry of the pair, and how its tetramer is made."""

    name: str
    pdb: str
    ligand: str
    # The biological assembly to build from the entry's REMARK 350, where the
    # entry's asymmetric unit holds less than the whole molecule. The chains
    # each BIOMT operator makes are lettered in order: the first operator's
    # copy keeps the entry's letters, the next takes the following ones.
    biomt: int | None = None


@dataclass(frozen=True)
class Assembly:
    slug: str
    display: str
    subunits: tuple[Subunit, ...]
    # The first state's entry frame is the pair's; the second is superposed
    # onto it on [superpose_on], the subunits whose frame is kept.
    states: tuple[State, State]
    superpose_on: tuple[str, ...]
    # PyMOL's cartoon sampling, as the catalog's hemoglobin row has it.
    sampling: int = 4

    def subunit(self, node: str) -> Subunit:
        return next(s for s in self.subunits if s.node == node)


ASSEMBLIES: tuple[Assembly, ...] = (
    Assembly(
        # Not a catalog slug: an assembly has its own table.
        slug="hemoglobin-a",
        display="Hemoglobin A",
        subunits=(
            Subunit("alpha1", "HBA1", "P69905", "A"),
            Subunit("beta1", "HBB", "P68871", "B"),
            Subunit("alpha2", "HBA1", "P69905", "C"),
            Subunit("beta2", "HBB", "P68871", "D"),
        ),
        # Park, Yokoyama, Shibayama, Shiro and Tame (2006): the same protein at
        # 1.25 A in both states. Deoxy, T: the entry holds the tetramer. Oxy,
        # R: the entry holds one alpha and one beta chain, and its biological
        # assembly is the tetramer.
        states=(
            State("tense", "2DN2", "deoxy"),
            State("relaxed", "2DN1", "oxy", biomt=1),
        ),
        # The two crystals are different lattices, so the entries' frames do not
        # agree (the raw CA atoms are 51-59 A apart). Baldwin and Chothia's
        # convention: hold one alpha-beta dimer still, and the other turns.
        superpose_on=("alpha1", "beta1"),
    ),
)

BY_SLUG = {a.slug: a for a in ASSEMBLIES}
