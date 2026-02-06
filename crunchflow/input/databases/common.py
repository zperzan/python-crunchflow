"""Common database entities for parsing thermodynamic and kinetics databases."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Literal, Optional, Tuple

__all__ = [
    # shared
    "ParseError",
    "ValidationError",
    "Reaction",
    "AqueousKinetics",
    # thermo
    "DebyeHuckelSet",
]

# ----------------------------------------------------------------------------- #
# Errors
# ----------------------------------------------------------------------------- #


class ParseError(ValueError):
    """Raised when a database file cannot be parsed."""


class ValidationError(ValueError):
    """Raised when parsed values fail basic validation."""


# ----------------------------------------------------------------------------- #
# Shared value objects
# ----------------------------------------------------------------------------- #


@dataclass
class Reaction:
    """
    Reaction stoichiometry in database format. Consists of one species on the
    left-hand side (lhs) and a dict of species and their coefficients on the
    right-hand side (stoich).
    """

    lhs: Optional[str] = None
    stoich: Optional[Dict[str, float]] = field(default_factory=dict)


@dataclass
class AqueousKinetics:
    """
    Single-law-first representation of an aqueous kinetic entry.

    Common, intuitive case (one law):
      - Edit fields directly on this object, e.g.:
          rate25C, orders  (for 'tst'/'irreversible')
          rate25C, monod_terms, inhibition  (for 'monod'/'MonodBiomass')

    Rare case (parallel laws):
      - Add additional laws in `parallel_laws` as dicts, each with the same
        reaction_type and the same schema as the primary law (see below).

    Schema
    ------
    label : str
        Entry name/label (often same as reaction name).
    reaction_type : {'tst','monod','irreversible','MonodBiomass'}
        Type applies to the primary law and ALL parallel laws.
    reaction : Reaction
        Stoichiometry (usually empty here; defined in &Aqueous for standalone DBs).
    rate25C : float | None
        Primary intrinsic rate (k or kmax) at 25 deg C.
    logK : float | None
        Optional equilibrium constant (some TST-style entries).
    orders : dict[str,float] | None
        Species exponents for 'tst'/'irreversible'.
    monod_terms : list[tuple[str,float]] | None
        (species, Khalf) pairs for 'monod'/'MonodBiomass'.
    inhibition : dict[str,float] | None
        (species -> Kin) for 'monod'/'MonodBiomass'.

    Monod extras (optional)
    -----------------------
    chi, biomass, bg, direction, use_metabolic_lag, lag_time, ramp_time,
    substrate_for_lag, threshold_concentration

    Parallel laws
    -------------
    parallel_laws : list[dict]
        Additional parallel rate laws (law #1, #2, ...) with the SAME
        reaction_type. Each dict may contain:
          - "rate25C": float
          - For 'tst'/'irreversible': "orders": dict[str,float]
          - For 'monod'/'MonodBiomass': "monod_terms": list[tuple[str,float]],
                                        "inhibition": dict[str,float]
          - Optional Monod extras (same keys as the primary law).
    """

    label: str
    rate25C: float
    reaction_type: Literal["tst", "monod", "irreversible", "MonodBiomass"]

    # Primary law (single-law-first)
    reaction: Optional[Reaction] = field(default_factory=Reaction)
    logK: Optional[float] = None

    # Optional tst/irreversible
    orders: Optional[Dict[str, float]] = field(default_factory=dict)
    dependence: Optional[Dict[str, float]] = field(default_factory=dict)

    # Optional Monod/biomass extras
    monod_terms: Optional[Dict[str, float]] = field(default_factory=dict)
    inhibition: Optional[Dict[str, float]] = field(default_factory=dict)
    chi: Optional[float] = None
    bq: Optional[float] = None
    direction: Optional[float] = None
    lag_time: Optional[float] = None
    ramp_time: Optional[float] = None
    threshold_concentration: Optional[float] = None
    use_metabolic_lag: Optional[bool] = None
    biomass: Optional[str] = None
    substrate_for_lag: Optional[str] = None

    # Rare case: additional parallel laws (same reaction_type)
    parallel_laws: List[Dict[str, Any]] = field(default_factory=list)

    def __post_init__(self) -> None:
        """Validate mutually exclusive fields and law consistency."""
        rt = self.reaction_type

        # Normalize simple aliases
        if rt.lower() in ("tst", "monod", "irreversible"):
            self.reaction_type = rt.lower()  # type: ignore[assignment]
        elif rt != "MonodBiomass":
            raise ValidationError(f"AqueousKinetics: unsupported reaction_type {rt!r}")

        # Avoid mixing law styles on the primary law
        if self.reaction_type in ("tst", "irreversible"):
            if self.monod_terms or self.inhibition:
                raise ValidationError("For 'tst'/'irreversible', use 'orders' (do not set monod_terms/inhibition).")
        else:  # 'monod' or 'MonodBiomass'
            if self.orders:
                raise ValidationError("For 'monod'/'MonodBiomass', use monod_terms[/inhibition] (do not set orders).")

        # Parallel laws must match reaction_type and follow the same schema
        for idx, law in enumerate(self.parallel_laws):
            if not isinstance(law, dict):
                raise ValidationError(f"parallel_laws[{idx}] must be a dict.")
            # Check mutually exclusive fields according to type
            if self.reaction_type in ("tst", "irreversible"):
                if "monod_terms" in law or "inhibition" in law:
                    raise ValidationError(f"parallel_laws[{idx}]: monod fields not allowed for '{self.reaction_type}'.")
            else:  # monod / MonodBiomass
                if "orders" in law:
                    raise ValidationError(f"parallel_laws[{idx}]: 'orders' not allowed for '{self.reaction_type}'.")

    def __str__(self) -> str:
        """
        Render as an &AqueousKinetics block (single-law-first).
        Note that parallel laws aren't printed here.
        """
        KEYW = 23

        def w(key: str, val: str) -> str:
            return f"  {key:<{KEYW}}= {val}"

        def q(name: str) -> str:
            s = str(name)
            if len(s) >= 2 and (s[0] in "'\"" and s[-1] == s[0]):
                s = s[1:-1]
            return f"'{s}'"

        out = ["&AqueousKinetics"]
        out.append(w("label", self.label or "default"))

        # Preserve MonodBiomass capitalization; 'monod' stays lowercase per examples
        rtype = self.reaction_type if self.reaction_type != "monod" else "monod"
        out.append(w("type", rtype))

        rate = 0.0 if self.rate25C is None else float(self.rate25C)
        out.append(w("rate25C", f"{rate:g}"))

        if self.reaction_type in ("monod", "MonodBiomass"):
            if self.monod_terms:
                parts = [f"{q(sp)} {float(kh):g}" for sp, kh in self.monod_terms]
                out.append(w("monod_terms", " ".join(parts)))
            if self.inhibition:
                parts = [f"{q(sp)} {float(val):g}" for sp, val in self.inhibition.items()]
                out.append(w("inhibition", " ".join(parts)))
            if self.biomass:
                out.append(w("biomass", q(self.biomass)))
            if self.bg is not None:
                out.append(w("bq", f"{float(self.bg):g}"))
            if self.chi is not None:
                out.append(w("chi", f"{float(self.chi):g}"))
            if self.direction is not None:
                out.append(w("direction", f"{float(self.direction):g}"))
            if self.use_metabolic_lag is not None:
                out.append(w("UseMetabolicLag", ".true." if self.use_metabolic_lag else ".false."))
            if self.lag_time is not None:
                out.append(w("LagTime", f"{float(self.lag_time):g}"))
            if self.ramp_time is not None:
                # Many example files use 'Ramptime' spelling
                out.append(w("Ramptime", f"{float(self.ramp_time):g}"))
            if self.threshold_concentration is not None:
                out.append(w("ThresholdConcentration", f"{float(self.threshold_concentration):g}"))
            if self.substrate_for_lag:
                out.append(w("SubstrateForLag", str(self.substrate_for_lag)))
        else:
            # tst / irreversible → species orders under 'dependence'
            if self.dependence:
                parts = [f"'{sp}'  {float(exp):g}" for sp, exp in self.dependence.items()]
                out.append(w("dependence", " ".join(parts)))

        # (Optional) hint about parallel laws without emitting them
        if self.parallel_laws:
            out.append("  !! note: additional parallel laws exist but are not printed here")

        out.append("/")
        return "\n".join(out)


# ----------------------------------------------------------------------------- #
# Thermodynamic database containers
# ----------------------------------------------------------------------------- #


@dataclass
class DebyeHuckelSet:
    """Debye–Hückel parameters (adh, bdh, bdt) tabulated on the same grid."""

    temperature_points: Tuple[float, ...]
    adh: List[float]
    bdh: List[float]
    bdt: List[float]

    def __post_init__(self) -> None:
        """Validate that all parameter arrays match the grid size."""
        n = len(self.temperature_points)
        if not (len(self.adh) == len(self.bdh) == len(self.bdt) == n):
            raise ValidationError("DebyeHuckelSet arrays must match grid size.")
