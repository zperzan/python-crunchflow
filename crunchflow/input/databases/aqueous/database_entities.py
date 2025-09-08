"""Entries from the aqueous kinetics database."""

from __future__ import annotations

from dataclasses import dataclass, field

from ..common import AqueousKinetics, Reaction

__all__ = [
    "AqueousReaction",  # &Aqueous entries (stoichiometry)
    "AqueousKinetics",  # reused for &AqueousKinetics entries
]


@dataclass
class AqueousReaction:
    """
    One &Aqueous entry (reaction stoichiometry) from an aqueous kinetics database.

    Attributes
    ----------
    name : str
        Reaction label/name as it appears after '&Aqueous'.
    reaction : Reaction
        Stoichiometry map (species -> coefficient). Convention is parser-defined but
        typically negatives for reactants and positives for products.
    """

    name: str
    type: str
    logK: float
    reaction: Reaction = field(default_factory=Reaction)
