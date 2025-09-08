"""Aqueous database with parser/writer for aqueous databases."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .._io_utils import _fmt_float, _is_blank_or_comment, _separate_quoted_and_numeric, _strip_bang, _tokens, _unquote
from ..common import (
    AqueousKinetics,
    ParseError,
    Reaction,
)
from .database_entities import AqueousReaction

__all__ = ["AqueousDatabase"]


# ----------------------------------------------------------------------------- #
# Small local helpers
# ----------------------------------------------------------------------------- #


def _parse_key_value_pairs(line) -> tuple[str, str]:
    """Return (key, value) from 'key = value'."""
    out: Tuple[str, str]
    toks = line.split("=", 1)
    if len(toks) != 2:
        raise ParseError(f"Malformed key-value pair: '{line}'")
    key = toks[0].strip()
    value = _strip_bang(toks[1].strip()).strip()
    out = (key, value)

    return out


# ----------------------------------------------------------------------------- #
# Facade
# ----------------------------------------------------------------------------- #


@dataclass
class AqueousDatabase:
    """
    Standalone aqueous database with two entry types:

      - &Aqueous          → Reaction stoichiometry (AqueousReaction)
      - &AqueousKinetics  → Rate expressions (AqueousKinetics)

    This class includes a permissive parser/writer for the format seen in files like
    'aqueous-BhavnaSep2015.dbs'. It aims to be intuitive to edit programmatically.
    """

    reactions: Dict[str, AqueousReaction] = field(default_factory=dict)
    kinetics: Dict[str, AqueousKinetics] = field(default_factory=dict)
    source_path: Optional[str] = None

    # ------------------------------- I/O ----------------------------------- #

    @classmethod
    def from_file(cls, path: str) -> "AqueousDatabase":
        """Instantiate and read an aqueous database from a file."""
        db = cls()
        db.read(path)
        return db

    def clear(self) -> None:
        """Clear all entries and reset source path."""
        self.reactions.clear()
        self.kinetics.clear()
        self.source_path = None

    def read(self, path: str) -> None:
        """Read an aqueous kinetic database from a file."""
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()
        except Exception as exc:
            raise ParseError(f"Unable to read aqueous database '{path}': {exc}") from exc

        self.clear()
        self.source_path = path

        i, n = 0, len(lines)
        while i < n:
            line = lines[i].strip()
            low = line.lower()
            if _is_blank_or_comment(line):
                i += 1
                continue

            # --- &Aqueous block (stoichiometry) ---
            if low.startswith("&aqueous") and not low.startswith("&aqueouskinetics"):
                i = self._parse_aqueous_block(lines, i + 1)
                continue

            # --- &AqueousKinetics block (rate expressions) ---
            if low.startswith("&aqueouskinetics"):
                i = self._parse_aqueous_kinetics_block(lines, i + 1)
                continue

            # lines outside recognized blocks are ignored
            i += 1

    def write(self, path: str) -> None:
        """Write the aqueous database to a file."""
        try:
            with open(path, "w", encoding="utf-8") as f:
                self._write_aqueous_block(f)
                self._write_aqueous_kinetics_block(f)
        except Exception as exc:
            raise ParseError(f"Failed to write aqueous database '{path}': {exc}") from exc

    # ----------------------------- Parsing -------------------------------- #

    def _parse_aqueous_block(self, lines: List[str], start_idx: int) -> int:
        """
        Parse a single &Aqueous block:
            &Aqueous
            <Name>
            <Type>
            <Stoichiometry>
            <keq>
        Stops at next line that starts with '&' or EOF.
        """
        i, n = start_idx, len(lines)

        # Find name line
        while i < n and _is_blank_or_comment(lines[i]):
            i += 1
        if i >= n:
            return i
        label, name = lines[i].split("=", 1)
        name = _strip_bang(name.strip())
        if label.strip() != "name":
            raise ParseError(f"Could not identify name in &Aqueous block on line {i:d}.")
        i += 1

        # Define default values for aqueous block
        reaction_type: Optional[str] = None
        logK: float = 0.0
        reaction: Optional[Reaction] = None

        # Collect content until next '&' block or '/' or EOF
        payload: List[str] = []
        while i < n:
            s = lines[i].strip()
            if s.startswith("&") or s.startswith("/"):
                break
            if not _is_blank_or_comment(s):
                payload.append(s)
            i += 1

        # Parse payload, which should consist of key-value pairs
        j = 0
        while j < len(payload):
            line = payload[j]
            key, value = _parse_key_value_pairs(line)
            if key == "type":
                reaction_type = value
            elif key == "keq":
                try:
                    logK = float(value)
                except ValueError as exc:
                    raise ParseError(f"&Aqueous '{name}': invalid keq value.") from exc
            elif key == "stoichiometry":
                toks = _tokens(value)
                try:
                    stoich_coefs = [float(t) for t in toks[::2]]
                except ValueError as exc:
                    raise ParseError(f"Could not read '{name}' stoichiometry: {line}") from exc
                rxn_species = toks[1::2]
                stoich = {sp: coef for sp, coef in zip(rxn_species, stoich_coefs)}
                reaction = Reaction(lhs=name, stoich=stoich)

                # Check if reaction continues on next line
                while j + 1 < len(payload):
                    next_line = payload[j + 1].strip()

                    break_keys = ("&", "/", "!", "keq", "type", "name", "stoichiometry")
                    if next_line.startswith(break_keys):
                        break
                    next_toks = _tokens(next_line)
                    if not next_toks or _is_blank_or_comment(next_line):
                        j += 1
                        continue
                    try:
                        next_coefs = [float(t) for t in next_toks[::2]]
                    except ValueError as exc:
                        raise ParseError(f"Could not read '{name}' stoichiometry continuation: {next_line}") from exc
                    next_species = next_toks[1::2]
                    for sp, coef in zip(next_species, next_coefs):
                        reaction.stoich[sp] = coef
                    j += 1
            j += 1

        self.reactions[name] = AqueousReaction(name=name, logK=logK, type=reaction_type, reaction=reaction)
        return i  # next parser continues at & or EOF

    def _parse_aqueous_kinetics_block(self, lines: List[str], start_idx: int) -> int:
        """
        Parse an &AqueousKinetics block with one or more entries.
        Observed styles in the wild vary; we accept two header variants:

        A) Single header line with name, number of laws and type:
              <Name> <N_laws> <'tst'|'monod'|'irreversible'>
           followed by N_laws lines:
              - tst/irreversible: <k> <N_terms> <'Spec'> <exp> ...
              - monod:           <kmax> <N_terms> <'Spec'> <Khalf> ... ['Inhibition' <M> <'Spec'> <Kin> ...]

        B) Two header lines:
              <Name>
              <N_laws> <'tst'|'monod'|'irreversible'>
           followed by the same law lines as above.

        The block ends when the next line starts with '&' or EOF.
        """
        i, n = start_idx, len(lines)

        # Find name line
        while i < n and _is_blank_or_comment(lines[i]):
            i += 1
        if i >= n:
            return i

        label, name = lines[i].split("=", 1)
        name = _strip_bang(name.strip())
        if label.strip() != "name":
            raise ParseError(f"Could not identify name in &Aqueous block on line {i:d}.")
        i += 1

        # Instantiate aqueous kinetic object with default rate and type
        aq = AqueousKinetics(label=name, reaction_type="monod", rate25C=0.0)

        # Check to make sure default rate/type are overwritten
        found_type = False
        found_rate = False

        # Collect content until next '&' block or '/' or EOF
        while i < n:
            s = lines[i].strip()
            if s.startswith("&") or s.startswith("/"):
                break
            if _is_blank_or_comment(s):
                i += 1
                continue
            if "=" not in s:
                raise ParseError(f"Could not parse line in &AqueousKinetics block: '{s}'")

            key, value = _parse_key_value_pairs(s)

            # Assign known string attributes
            if key == "label":
                aq.label = _unquote(value)
            elif key == "type":
                found_type = True
                if value in ("tst", "monod", "irreversible", "MonodBiomass"):
                    aq.reaction_type = value  # type: ignore[assignment]
                else:
                    raise ParseError(
                        f"Unsupported reaction type '{value}' in &AqueousKinetics '{name}'.\n"
                        f"Supported types are: 'tst', 'monod', 'irreversible', 'MonodBiomass'."
                    )
            elif key == "UseMetabolicLag":
                if value == ".false.":
                    aq.use_metabolic_lag = False
                elif value == ".true.":
                    aq.use_metabolic_lag = True
            elif key == "biomass":
                aq.biomass = _unquote(value)
            elif key == "SubstrateForLag":
                aq.substrate_for_lag = value
            elif key == "monod_terms":
                toks = _tokens(value)
                quoted, nums = _separate_quoted_and_numeric(toks)
                if len(quoted) != len(nums):
                    raise ParseError(f"Mismatched Monod terms in &AqueousKinetics '{name}': {value}")
                for sp, kh in zip(quoted, nums):
                    aq.monod_terms[_unquote(sp)] = kh
            elif key == "inhibition":
                toks = _tokens(value)
                quoted, nums = _separate_quoted_and_numeric(toks)
                if len(quoted) != len(nums):
                    raise ParseError(f"Mismatched inhibition terms in &AqueousKinetics '{name}': {value}")
                for sp, kin in zip(quoted, nums):
                    aq.inhibition[sp] = kin
            elif key == "dependence":
                toks = _tokens(value)
                quoted, nums = _separate_quoted_and_numeric(toks)
                if len(quoted) != len(nums):
                    raise ParseError(f"Mismatched dependence terms in &AqueousKinetics '{name}': {value}")
                for sp, exp in zip(quoted, nums):
                    aq.dependence[sp] = exp

            # Assign known float attributes
            float_search_keys = {
                "rate25C": "rate25C",
                "chi": "chi",
                "bq": "bq",
                "LagTime": "lag_time",
                "Ramptime": "ramp_time",
                "ThresholdConcentration": "threshold_concentration",
                "direction": "direction",
            }
            for search_key, attr in float_search_keys.items():
                if key == search_key:
                    try:
                        setattr(aq, attr, float(value))
                        if key == "rate25C":
                            found_rate = True
                    except ValueError as exc:
                        raise ParseError(f"Invalid {search_key} value in &AqueousKinetics '{name}': {value}") from exc
                    break
            i += 1

        if not found_type:
            raise ParseError(f"Could not read reaction 'type' in &AqueousKinetics '{name}'.")
        if not found_rate:
            raise ParseError(f"Could not read 'rate25C' in &AqueousKinetics '{name}'.")

        self.kinetics[name] = aq

        return i

    # ----------------------------- Writing -------------------------------- #

    def _write_aqueous_block(self, f) -> None:
        if not self.reactions:
            return
        f.write("! --------------------------------------------------------------------\n")
        f.write("!\n")
        f.write("!  reaction stoichiometry\n")
        f.write("!\n")
        f.write("! --------------------------------------------------------------------\n\n")

        for name, aq_reaction in self.reactions.items():
            f.write("&Aqueous\n")
            rxn = aq_reaction.reaction
            f.write(f"  name          = {name}\n")
            f.write(f"  type          = {aq_reaction.type}\n")
            # Write stoichiometry
            stoich_line = "  stoichiometry = "
            line_len = len(stoich_line)
            for sp, coef in rxn.stoich.items():
                stoich_line += f"{_fmt_float(coef, 2)} {sp}  "
                # wrap long lines
                line_len += len(f"{_fmt_float(coef, 2)} {sp}  ")
                if line_len > 90:
                    stoich_line += "\n                  "
                    line_len = 18

            f.write(stoich_line + "\n")
            f.write(f"  keq           = {_fmt_float(aq_reaction.logK, align='left')}\n")
            f.write("/\n\n")

    def _write_aqueous_kinetics_block(self, f) -> None:
        """
        Write &AqueousKinetics entries in the standalone aqueous database format, e.g.:

          &AqueousKinetics
            name                   = AceNO3HCO3NO2
            label                  = default
            type                   = MonodBiomass
            rate25C                = 2000.0
            monod_terms            = 'tot_Acetate' 2.03E-5 'tot_NO3-' 1.06E-5 'tot_NH4+' 1.0e-6
            biomass                = 'C5H7O2NNO3(s)'
            bq                     = -0.0
            chi                    = 1
            direction              = -1
          /
        """
        if not self.kinetics:
            return
        f.write("! --------------------------------------------------------------------\n")
        f.write("!\n")
        f.write("!  rate expressions\n")
        f.write("!\n")
        f.write("! --------------------------------------------------------------------\n\n")

        def w_key(key: str, val: str) -> None:
            # fixed column for alignment (match examples with wide key column)
            KEYW = 23
            f.write(f"  {key:<{KEYW}}= {val}\n")

        for name, ak in self.kinetics.items():
            f.write("&AqueousKinetics\n")
            # Universal keys
            w_key("name", name)
            w_key("label", ak.label if ak.label else "default")
            w_key("type", ak.reaction_type)
            w_key("rate25C", f"{_fmt_float(ak.rate25C, align='left')}")

            # Reaction type-specific keys
            if ak.reaction_type in ("MonodBiomass", "monod"):
                # monod_terms
                if ak.monod_terms:
                    parts: List[str] = []
                    for sp, khalf in ak.monod_terms.items():
                        parts.append(f"'{sp}' {_fmt_float(khalf, 6, align='left')}")
                    w_key("monod_terms", "  ".join(parts))
                # inhibition
                if ak.inhibition:
                    parts = []
                    for isp, kin in ak.inhibition.items():
                        parts.append(f"'{isp}' {_fmt_float(kin, 6, align='left')}")
                    w_key("inhibition", "  ".join(parts))
                # optional biomass/monod extras
                if ak.biomass:
                    w_key("biomass", f"'{ak.biomass}'")
                if ak.chi is not None:
                    w_key("chi", f"{_fmt_float(ak.chi, 6, align='left')}")
                if ak.bq is not None:
                    w_key("bq", f"{_fmt_float(ak.bq, 6, align='left')}")
                if ak.direction is not None:
                    w_key("direction", f"{_fmt_float(ak.direction, 6, align='left')}")
                if ak.use_metabolic_lag is not None:
                    w_key("UseMetabolicLag", ".true." if ak.use_metabolic_lag else ".false.")
                if ak.lag_time is not None:
                    w_key("LagTime", f"{_fmt_float(ak.lag_time, 6, align='left')}")
                if ak.ramp_time is not None:
                    # examples often use the key 'Ramptime' (lowercase p); match that
                    w_key("Ramptime", f"{_fmt_float(ak.ramp_time, 6, align='left')}")
                if ak.threshold_concentration is not None:
                    w_key("ThresholdConcentration", f"{_fmt_float(ak.threshold_concentration, 6, align='left')}")
                if ak.substrate_for_lag:
                    w_key("SubstrateForLag", ak.substrate_for_lag)

            else:
                # tst / irreversible style → dependence (species orders)
                if ak.dependence:
                    parts = []
                    for sp, exp in ak.dependence.items():
                        parts.append(f"'{sp}'  {_fmt_float(exp, 6, align='left')}")
                    w_key("dependence", " ".join(parts))

            f.write("/\n\n")
