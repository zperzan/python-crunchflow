# python_crunchflow/input/databases/aqueous/__init__.py
"""
Aqueous kinetics database API.

Supports standalone aqueous databases with two entry types:
  - &Aqueous          → reaction stoichiometry (AqueousReaction)
  - &AqueousKinetics  → rate laws (AqueousKinetics)

The facade `AqueousDatabase` provides read/write and simple edit helpers.
"""

from ..common import AqueousKinetics
from .database import AqueousDatabase
from .database_entities import AqueousReaction

__all__ = [
    "AqueousDatabase",
    "AqueousReaction",
    "AqueousKinetics",
]
