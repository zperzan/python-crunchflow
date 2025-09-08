"""Entries from the thermodynamic database."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Literal, Optional, Tuple

from ..common import (
    Reaction,
    ValidationError,
)

__all__ = [
    "SpeciesBase",
    "PrimarySpecies",
    "SecondarySpecies",
    "GasSpecies",
    "Minerals",
    "SurfaceComplex",
    "MineralKinetics",
    "Exchange",
    "SurfaceComplexParams",
]

# ----------------------------------------------------------------------------- #
# Species / phases
# ----------------------------------------------------------------------------- #


@dataclass
class SpeciesBase:
    """Common fields for species and phases."""

    name: str


@dataclass
class PrimarySpecies(SpeciesBase):
    """Entry from the 'primary' block"""

    size: float = 0.0
    charge: Optional[float] = None
    molar_mass: Optional[float] = None


@dataclass
class SecondarySpecies(SpeciesBase):
    """
    Secondary aqueous complex:
      - reaction: stoichiometry to primary species basis
      - logK: tabulated log10(K) aligned with the database temperature grid
    """

    reaction: Reaction = field(default_factory=Reaction)
    logK: Dict[float, float] = field(default_factory=dict)
    size: float = 0.0
    charge: Optional[float] = None
    molar_mass: Optional[float] = None


@dataclass
class GasSpecies(SpeciesBase):
    """
    Gas entry in the 'gases' block; represented by its formation/dissolution
    reaction to primary basis and a tabulated equilibrium constant.
    """

    reaction: Reaction = field(default_factory=Reaction)
    logK: Dict[float, float] = field(default_factory=dict)
    molar_volume: float = 0.0
    molar_mass: Optional[float] = None


@dataclass
class Minerals(SpeciesBase):
    """Mineral entry in the 'minerals' block; dissolution reaction and logK(T)."""

    reaction: Reaction = field(default_factory=Reaction)
    logK: Dict[float, float] = field(default_factory=dict)
    molar_volume: float = 0.0
    molar_mass: Optional[float] = None


@dataclass
class SurfaceComplex(SpeciesBase):
    """
    Surface complex from the 'surface complexation' block.
      - site: the surface site species name (e.g., '>FeOH_strong')
      - reaction: includes the site on the reactant/product side as defined
      - logK: tabulated equilibrium constant for the surface reaction
    """

    reaction: Reaction = field(default_factory=Reaction)
    logK: Dict[float, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate fields after initialization."""
        if not self.name:
            raise ValidationError("SurfaceComplex requires a non-empty name.")


@dataclass
class Exchange(SpeciesBase):
    """
    Exchange reaction from the 'Begin exchange' block.
      - name: the ion exchange species
      - reaction: includes the site on the reactant/product side as defined
      - logK: tabulated equilibrium constant for the exchange
    """

    reaction: Reaction = field(default_factory=Reaction)
    logK: Dict[float, float] = field(default_factory=dict)


@dataclass
class SurfaceComplexParams(SpeciesBase):
    """
    Surface complexation parameters from the thermodynamic database.
      - name: the ion exchange species
    """

    parameter: float = None


@dataclass
class MineralKinetics:
    """
    Mineral kinetics entry from the thermodynamic database.

    Fields mirror the database block:
      mineral: solid phase name (must match thermo section name)
      label: parallel reaction label (e.g., 'default', 'h+')
      reaction_type: 'tst' | 'monod' | 'irreversible' | 'PrecipitationOnly' | 'DissolutionOnly'
      rate25_log10: log10(rate) at 25 °C (as stored in DB)
      activation_kcal_per_mol: activation energy (kcal/mol)

    Optional fields by type:
      - TST / Irreversible / PrecipitationOnly / DissolutionOnly:
          dependence: {species -> exponent}
          affinity_dependence: (m1, m2, m3) if specified via AffinityDependence
      - Monod:
          monod_terms: [(species, Khalf)], inhibition: {species -> Kin}

    reaction / reaction_text are for optional storage of a parsed stoichiometric line,
    when provided after 'dependence :' in some entries.
    """

    mineral: str
    label: str
    reaction_type: Literal["tst", "monod", "irreversible", "PrecipitationOnly", "DissolutionOnly"]
    rate25_log10: float
    activation_kcal_per_mol: float

    # TST / Irreversible / PrecipitationOnly / DissolutionOnly
    dependence: Dict[str, float] = field(default_factory=dict)
    affinity_dependence: Optional[Tuple[float, float, float]] = None  # (m1, m2, m3)

    # Monod
    monod_terms: List[Tuple[str, float]] = field(default_factory=list)  # (species, Khalf)
    inhibition: Dict[str, float] = field(default_factory=dict)  # species -> Kin

    # Optional parsed reaction line
    reaction: Optional[Reaction] = None

    def __post_init__(self) -> None:
        """Validate and normalize fields after initialization."""
        # Normalize/validate reaction_type
        rt = self.reaction_type

        valid_reaction_types = {"tst", "monod", "irreversible", "PrecipitationOnly", "DissolutionOnly"}
        if rt not in valid_reaction_types:
            raise ValidationError(f"Unsupported reaction_type: {rt!r}")

        if not self.mineral:
            raise ValidationError("MineralKinetics.mineral must be non-empty.")
        if not self.label:
            raise ValidationError("MineralKinetics.label must be non-empty.")
