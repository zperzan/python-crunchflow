"""The `crunchflow.input` subpackage provides tools for reading, writing, and modifying
CrunchFlow input files using Python classes that mirror the structure of these files.

CrunchFlow input files are divided into keyword blocks (e.g., `RUNTIME`, `FLOW`, `TRANSPORT`,
`CONDITION`), each of which controls a specific aspect of the simulation setup. This
subpackage offers Python representations for each of these blocks and supports programmatic
workflows for creating or editing `.in` files.

Some notable features include:
    - The `InputFile` class represents an entire CrunchFlow input file, with one attribute
      per block (e.g., `.runtime`, `.flow`, `.output`, etc.).
    - Each block is implemented as a dedicated class (e.g., `Flow`, `Transport`,
      `Condition`) and inherits from a shared `KeywordBlock` base class.
    - Input files can be generated from scratch, modified in memory, or parsed from existing
      files.
    - The `ThermoDatabase` class provides a way to read and represent CrunchFlow thermodynamic
      databases, including species and phases.

Most blocks in a CrunchFlow input file correspond directly to attributes of the `InputFile`
class. For example, the `FLOW` block is accessible via `InputFile.flow`, and the `RUNTIME`
block via `InputFile.runtime`. Similarly, most keywords within a block (e.g., `calculate_flow`,
`time_units`, or `permeability_x`) are represented as attributes of that block’s class. For
example, `InputFile.flow.calculate_flow` or `InputFile.runtime.time_units`.

Two notable exceptions to this structure are:
    - **Conditions**: Because multiple `CONDITION` blocks can be defined in a single input file,
      they are stored as a dictionary in `InputFile.conditions`, with each key representing the
      name of a different condition.
    - **Species in Conditions**: Within a `CONDITION` block, any number of species can be specified,
      with arbitrary names (e.g., `Na+`, `Cl-`, `HCO3-`). These are stored as key-value pairs in the
      `concentrations` dictionary of the `Condition` object, rather than as fixed attributes.
"""

# Input-file API
from .blocks import Condition, KeywordBlock

# Aqueous database API
from .databases.aqueous.database import AqueousDatabase
from .databases.aqueous.database_entities import AqueousReaction

# Thermodynamic database API
from .databases.thermo.database import ThermoDatabase
from .databases.thermo.database_entities import (
    Exchange,
    GasSpecies,
    MineralKinetics,
    Minerals,
    PrimarySpecies,
    SecondarySpecies,
    SurfaceComplex,
    SurfaceComplexParams,
)
from .inputfile import InputFile

__all__ = [
    # Input files
    "InputFile",
    "KeywordBlock",
    "Condition",
    # Thermodynamic databases
    "ThermoDatabase",
    "PrimarySpecies",
    "SecondarySpecies",
    "GasSpecies",
    "Minerals",
    "SurfaceComplex",
    "MineralKinetics",
    "Exchange",
    "SurfaceComplexParams",
    # Aqueous databases
    "AqueousDatabase",
    "AqueousReaction",
]
