# python_crunchflow/input/databases/__init__.py
"""
Databases subpackage: typed, file-backed stores for CrunchFlow data.

Currently includes thermodynamic databases; shared value objects live in
`common.py`. Kinetics-oriented containers are also defined in `common.py`
for future use.
"""

from .aqueous.database import AqueousDatabase
from .aqueous.database_entities import AqueousReaction
from .common import (
    AqueousKinetics,
    DebyeHuckelSet,
    ParseError,
    Reaction,
    ValidationError,
)
from .thermo.database import ThermoDatabase
from .thermo.database_entities import (
    Exchange,
    GasSpecies,
    MineralKinetics,
    Minerals,
    PrimarySpecies,
    SecondarySpecies,
    SurfaceComplex,
    SurfaceComplexParams,
)

__all__ = [
    # shared/common
    "ParseError",
    "ValidationError",
    "Reaction",
    "DebyeHuckelSet",
    "AqueousKinetics",
    # thermo
    "ThermoDatabase",
    "PrimarySpecies",
    "SecondarySpecies",
    "GasSpecies",
    "Minerals",
    "SurfaceComplex",
    "MineralKinetics",
    "Exchange",
    "SurfaceComplexParams",
    # aqueous
    "AqueousDatabase",
    "AqueousReaction",
]
