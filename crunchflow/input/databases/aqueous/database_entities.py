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

    def __str__(self) -> str:
        """
        Print as an &Aqueous block:
          &Aqueous
            name          = <name>
            type          = <type>            # if available
            stoichiometry = <coef> <Spec> ...
            keq           = <logK>            # if available
          /
        """
        KEYW = 23
        lines = ["&Aqueous"]
        lines.append(f"  {'name':<{KEYW}}= {self.name}")

        # Optional fields if present on this object
        rxn_type = getattr(self, "type", None)
        if rxn_type:
            lines.append(f"  {'type':<{KEYW}}= {rxn_type}")

        # Stoichiometry: wrap gently
        stoich = self.reaction.stoich if self.reaction else {}
        head = f"  {'stoichiometry':<{KEYW}}= "
        curr = head
        col = len(head)

        def fmt_pair(sp, coef):
            return f"{coef:g} {sp}  "

        for sp, coef in stoich.items():
            piece = fmt_pair(sp, coef)
            if col + len(piece) > 96:
                lines.append(curr.rstrip())
                curr = " " * (2 + KEYW + 2)  # indent continuation under value
                col = len(curr)
            curr += piece
            col += len(piece)
        lines.append(curr.rstrip())

        logK = getattr(self, "logK", None)
        if logK is not None:
            lines.append(f"  {'keq':<{KEYW}}= {logK:g}")

        lines.append("/")
        return "\n".join(lines)
