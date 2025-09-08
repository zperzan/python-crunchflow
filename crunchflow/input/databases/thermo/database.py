"""Module to read/write CrunchFlow-style thermodynamic databases."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .._io_utils import (
    _floats_in,
    _fmt_float,
    _is_blank_or_comment,
    _separate_quoted_and_numeric,
    _strip_bang,
    _tokens,
    _unquote,
)
from ..common import AqueousKinetics, DebyeHuckelSet, ParseError, Reaction, ValidationError
from .database_entities import (
    Exchange,
    GasSpecies,
    MineralKinetics,
    Minerals,
    PrimarySpecies,
    SecondarySpecies,
    SurfaceComplex,
    SurfaceComplexParams,
)

__all__ = ["ThermoDatabase"]


# ----------------------------------------------------------------------------- #
# Thermo database facade with inline parser/writer (CrunchFlow-style)
# ----------------------------------------------------------------------------- #


@dataclass
class ThermoDatabase:
    """
    In-memory representation of a CrunchFlow-style thermodynamic database.

    This single file contains both the data container and the lightweight
    parser/writer methods (mirroring the pattern used by input/inputfile.py).
    """

    temperature_points: Tuple[float, ...] = field(default_factory=tuple)
    activity: DebyeHuckelSet = None

    # Sections keyed by species/phase name as they appear in the file
    primary: Dict[str, PrimarySpecies] = field(default_factory=dict)
    secondary: Dict[str, SecondarySpecies] = field(default_factory=dict)
    gases: Dict[str, GasSpecies] = field(default_factory=dict)
    minerals: Dict[str, Minerals] = field(default_factory=dict)
    surface_complexes: Dict[str, SurfaceComplex] = field(default_factory=dict)
    aqueous_kinetics: Dict[str, AqueousKinetics] = field(default_factory=dict)
    mineral_kinetics: Dict[str, MineralKinetics] = field(default_factory=dict)
    exchange: Dict[str, Exchange] = field(default_factory=dict)
    surface_complex_params: Dict[str, SurfaceComplexParams] = field(default_factory=dict)

    # Optional provenance
    source_path: Optional[str] = None

    # ------------------------------- I/O ----------------------------------- #
    @classmethod
    def from_file(cls, path: str) -> "ThermoDatabase":
        """Read a database from a file and return the populated instance."""
        db = cls()
        db.read(path)
        return db

    def read(self, path: str) -> None:
        """
        Populate this instance from a thermodynamic database file.
        Overwrites any existing contents.

        Notes on expected structure (typical CrunchFlow database):
          - Temperature line (e.g., "temperature points: 0 25 60 ...")
          - Debye–Hückel arrays (adh, bdh, bdt) on the same grid
          - Sections terminated by keywords:
              End of primary
              End of secondary
              End of gases
              End of minerals
            Surface complexes bracketed by:
              Begin surface complexation
              End of surface complexation
        """
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()
        except Exception as exc:
            raise ParseError(f"Unable to read file '{path}': {exc}") from exc

        # Reset current state
        self.clear()
        self.source_path = path

        # Parse header (temperature grid + Debye–Hückel)
        idx = 0
        idx = self._parse_header(lines, idx)

        # Parse sections
        idx = self._parse_primary(lines, idx)
        idx = self._parse_secondary(lines, idx)
        idx = self._parse_gases(lines, idx)
        idx = self._parse_minerals(lines, idx)
        idx = self._parse_surface_complexation(lines, idx)
        idx = self._parse_aqueous_kinetics(lines, idx)
        idx = self._parse_mineral_kinetics(lines, idx)
        idx = self._parse_exchange(lines, idx)
        idx = self._parse_surface_complexation_params(lines, idx)

    def write(self, path: str) -> None:
        """
        Write the current database to a file in a straightforward, readable format.
        This prioritizes clarity and round-trippability for files parsed by this module.
        """
        if not self.temperature_points:
            raise ValidationError("Cannot write: header is missing.")
        try:
            with open(path, "w", encoding="utf-8") as f:
                self._write_header(f)
                self._write_primary(f)
                self._write_secondary(f)
                self._write_gases(f)
                self._write_minerals(f)
                self._write_surface_complexation(f)
                self._write_aqueous_kinetics(f)
                self._write_mineral_kinetics(f)
                self._write_exchange(f)
                self._write_surface_complexation_params(f)
        except Exception as exc:
            raise ParseError(f"Failed to write thermo database: {exc}") from exc

    # ----------------------------- Utilities ------------------------------- #

    def clear(self) -> None:
        """Clear all contents, resetting to an empty state."""
        self.temperature_points = ()
        self.primary.clear()
        self.secondary.clear()
        self.gases.clear()
        self.minerals.clear()
        self.surface_complexes.clear()
        self.aqueous_kinetics.clear()
        self.mineral_kinetics.clear()
        self.exchange.clear()
        self.surface_complex_params.clear()
        self.source_path = None

    # ----------------------------- Parsing -------------------------------- #
    def _parse_header(self, lines: List[str], start_idx: int) -> int:
        """
        Parse temperature grid and Debye–Hückel arrays.
        Returns the next line index to continue parsing sections.
        """
        n = len(lines)
        temps: Optional[List[float]] = None
        adh: Optional[List[float]] = None
        bdh: Optional[List[float]] = None
        bdt: Optional[List[float]] = None

        i = start_idx
        while i < n:
            s = lines[i].strip()
            if _is_blank_or_comment(s):
                i += 1
                continue

            low = s.lower()
            if low.startswith("temperature") or "temperature points" in low:
                vals = _floats_in(s)
                # if only label on this line, look ahead to next non-comment line
                if not vals:
                    j = i + 1
                    while j < n and _is_blank_or_comment(lines[j]):
                        j += 1
                    if j < n:
                        vals = _floats_in(lines[j])
                        i = j
                if not vals:
                    raise ParseError("Found temperature header but no numeric values.")
                temps = vals[1:]

            elif "Debye-Huckel adh" in s:
                adh = _floats_in(s)
            elif "Debye-Huckel bdh" in s:
                bdh = _floats_in(s)
            elif "Debye-Huckel bdt" in s:
                bdt = _floats_in(s)
            else:
                break

            i += 1

        if temps is None:
            raise ParseError("Temperature grid not found in header.")
        temp_values = tuple(temps)
        self.temperature_points = temp_values

        # Debye–Huckel arrays are optional in some DBs; if present, lengths must match grid.
        if not (len(adh) == len(bdh) == len(bdt) == len(temp_values)):
            raise ParseError("Debye–Huckel arrays must match temperature grid length.")
        self.activity = DebyeHuckelSet(temperature_points=temp_values, adh=list(adh), bdh=list(bdh), bdt=list(bdt))

        return i

    def _parse_primary(self, lines: List[str], start_idx: int) -> int:
        """
        Parse 'primary' section until a line containing 'End of primary'.
        Primary entries are usually given one per line, quoted name first.
        """
        i = start_idx
        n = len(lines)
        while i < n:
            s = lines[i].strip()
            if "end of primary" in s.lower():
                return i + 1
            if _is_blank_or_comment(s):
                i += 1
                continue
            toks = _tokens(s)
            if not toks:
                i += 1
                continue
            if toks and toks[0].startswith("'"):
                name = _unquote(toks[0])
                size = float(toks[1])
                charge = float(toks[2])
                molar_mass = float(toks[3])

                # Store as PrimarySpecies class
                self.primary[name] = PrimarySpecies(name=name, size=size, charge=charge, molar_mass=molar_mass)
            else:
                # If we encounter a non-quoted start, primary section likely not here.
                # Stop early and let the next section handler continue.
                return i
            i += 1
        return i

    def _parse_secondary(self, lines: List[str], start_idx: int) -> int:
        """
        Parse secondary species. Expected format is:
        <‘SpeciesName’ > <num_species> <Stoichiometric coefficient 1>  <’SpeciesName 1’>
        <Stoichiometric coefficient 2> <’SpeciesName’ 2> ... <Log K array> <Debye-Huckel size parameter>
        <Charge> <Molar mass>
        """
        grid_len = len(self.temperature_points)

        i = start_idx
        n = len(lines)
        while i < n:
            s = lines[i].strip()
            low = s.lower()
            if "end of secondary" in low:
                return i + 1
            if _is_blank_or_comment(s):
                i += 1
                continue

            toks = _tokens(s)
            if not toks:
                i += 1
                continue
            if not toks[0].startswith("'"):
                return i

            name = _unquote(toks[0])

            # Separate quoted species tokens from trailing numeric block
            quoted, nums = _separate_quoted_and_numeric(toks[1:])

            # Format for secondary species is:
            # <‘SpeciesName’ > <num_species> <Stoichiometric coefficient 1>  <’SpeciesName 1’>
            # <Stoichiometric coefficient 2> <’SpeciesName’ 2> ... <Log K array> <Debye-Huckel size parameter>
            # <Charge> <Molar mass>
            num_species = int(nums[0])

            if num_species > len(quoted):
                raise ParseError(f"Species '{name}': expected {num_species} species in reaction, found {len(quoted)}.")

            # Define j, which is the starting index of the current section within nums
            j = 1

            # Build the reaction from the quoted species
            rxn_species = quoted[:num_species]
            stoich_coefs = nums[j : num_species + j]
            rxn = Reaction(lhs=name, stoich={sp: coef for sp, coef in zip(rxn_species, stoich_coefs)})
            j += num_species

            # Get logK values
            logK = dict(zip(self.temperature_points, nums[j : grid_len + j]))
            j += grid_len

            # Finally get debye-huckel size, charge and molar mass
            dh_size = nums[j]
            charge = nums[j + 1]
            molar_mass = nums[j + 2]

            self.secondary[name] = SecondarySpecies(
                name=name,
                reaction=rxn,
                logK=logK,
                size=dh_size,
                charge=charge,
                molar_mass=molar_mass,
            )
            i += 1
        return i

    def _parse_gases(self, lines: List[str], start_idx: int) -> int:
        """
        Parse 'gases' section until 'End of gases'. Expected format is:
        <‘SpeciesName’> <Molar volume> <num_species> <Stoichiometric coefficient 1>  <’SpeciesName 1’>
        <Stoichiometric coefficient 2> <’SpeciesName’ 2> ... <Log K array> <Molar mass>
        """
        grid_len = len(self.temperature_points)

        i = start_idx
        n = len(lines)
        while i < n:
            s = lines[i].strip()
            low = s.lower()
            if "end of gases" in low:
                return i + 1
            if _is_blank_or_comment(s):
                i += 1
                continue

            toks = _tokens(s)
            if not toks:
                i += 1
                continue
            if not toks[0].startswith("'"):
                return i

            name = _unquote(toks[0])

            # Separate quoted species tokens from trailing numeric block
            quoted, nums = _separate_quoted_and_numeric(toks[1:])

            # Format for gases is:
            # <‘SpeciesName’> <Molar volume> <num_species> <Stoichiometric coefficient 1>  <’SpeciesName 1’>
            # <Stoichiometric coefficient 2> <’SpeciesName’ 2> ... <Log K array> <Molar mass>
            molar_volume = nums[0]
            num_species = int(nums[1])

            if num_species > len(quoted):
                raise ParseError(f"Gas '{name}': expected {num_species} species in reaction, found {len(quoted)}.")

            # Build the reaction from the quoted species
            # Define j, which is the starting index of the current section within nums
            j = 2
            rxn_species = quoted[:num_species]
            stoich_coefs = nums[j : num_species + j]
            rxn = Reaction(lhs=name, stoich={sp: coef for sp, coef in zip(rxn_species, stoich_coefs)})
            j += num_species

            # Get logK values
            logK = dict(zip(self.temperature_points, nums[j : grid_len + j]))
            j += grid_len

            # Finally get molar mass
            molar_mass = nums[j]

            self.gases[name] = GasSpecies(
                name=name,
                reaction=rxn,
                logK=logK,
                molar_volume=molar_volume,
                molar_mass=molar_mass,
            )
            i += 1
        return i

    def _parse_minerals(self, lines: List[str], start_idx: int) -> int:
        """
        Parse 'minerals' section until 'End of minerals'. Expected format is:
        <‘MineralName’ > <Molar Volume> <Stoichiometric coefficient 1>  <’SpeciesName 1’>
        <Stoichiometric coefficient 2> <’SpeciesName’ 2> ... <Log K array> <Molar mass>
        """
        grid_len = len(self.temperature_points)

        i = start_idx
        n = len(lines)
        while i < n:
            s = str(lines[i].strip())
            low = s.lower()
            if "end of minerals" in low:
                return i + 1
            if _is_blank_or_comment(s):
                i += 1
                continue

            toks = _tokens(s)
            if not toks:
                i += 1
                continue
            if not toks[0].startswith("'"):
                # Not in minerals yet
                return i

            name = _unquote(toks[0])

            # Separate quoted species tokens from trailing numeric block
            quoted, nums = _separate_quoted_and_numeric(toks[1:])

            # Format for minerals is:
            # <‘MineralName’> <Molar volume> <Stoichiometric coefficient 1>  <’SpeciesName 1’>
            # <Stoichiometric coefficient 2> <’SpeciesName’ 2> ... <Log K array> <Molar mass>
            molar_volume = nums[0]
            num_species = int(nums[1])

            if num_species > len(quoted):
                raise ParseError(f"Mineral '{name}': expected {num_species} species in reaction, found {len(quoted)}.")

            # Build the reaction from the quoted species
            # Define j, which is the starting index of the current section within nums
            j = 2
            rxn_species = quoted[:num_species]
            stoich_coefs = nums[j : num_species + j]
            rxn = Reaction(lhs=name, stoich={sp: coef for sp, coef in zip(rxn_species, stoich_coefs)})
            j += num_species

            # Get logK values
            logK = dict(zip(self.temperature_points, nums[j : grid_len + j]))
            j += grid_len

            # Finally get molar mass
            molar_mass = nums[j]

            self.minerals[name] = Minerals(
                name=name,
                reaction=rxn,
                logK=logK,
                molar_volume=molar_volume,
                molar_mass=molar_mass,
            )
            i += 1
        return i

    def _parse_surface_complexation(self, lines: List[str], start_idx: int) -> int:
        """
        Parse 'Begin surface complexation' ... 'End of surface complexation'.
        Each entry expected like:
          'ComplexName'  '<site_name>'  <Stoichiometric coefficient 1>  <’SpeciesName 1’>
          <Stoichiometric coefficient 2> <’SpeciesName’ 2> ... <Log K array>
        """
        grid_len = len(self.temperature_points)

        i = start_idx
        n = len(lines)

        # Seek begin marker
        while i < n and "begin surface complexation" not in lines[i].strip().lower():
            if ("end of " in lines[i].strip().lower()) or (lines[i].strip().startswith("'")):
                # allow other sections to start; surface block may be absent
                if "begin" in lines[i].strip().lower():
                    return i  # let the next parser start here
            i += 1

        if i >= n:
            # No surface complexation block
            return start_idx

        i += 1  # move past "Begin surface complexation"
        while i < n:
            s = lines[i].strip()
            low = s.lower()
            if "end of surface complexation" in low:
                return i + 1
            if _is_blank_or_comment(s):
                i += 1
                continue

            toks = _tokens(s)
            if not toks or not (toks[0].startswith("'") or toks[0].startswith('"')):
                i += 1
                continue

            name = _unquote(toks[0])

            # Separate quoted species tokens from trailing numeric block
            quoted, nums = _separate_quoted_and_numeric(toks[1:])

            if not quoted or not nums:
                raise ParseError(f"Surface complex '{name}': Could not parse reaction stoichiometry.")

            # Format for minerals is:
            # 'ComplexName'  '<site_name>'  <Stoichiometric coefficient 1>  <’SpeciesName 1’>
            # <Stoichiometric coefficient 2> <’SpeciesName’ 2> ... <Log K array>
            num_species = int(nums[0])

            if num_species > len(quoted):
                raise ParseError(f"'{name}': expected {num_species} species in reaction, found {len(quoted)}.")

            # Build the reaction from the quoted species
            # Define j, which is the starting index of the current section within nums
            j = 1
            rxn_species = quoted[:num_species]
            stoich_coefs = nums[j : num_species + j]
            rxn = Reaction(lhs=name, stoich={sp: coef for sp, coef in zip(rxn_species, stoich_coefs)})
            j += num_species

            # Get logK values
            logK = dict(zip(self.temperature_points, nums[j : grid_len + j]))

            self.surface_complexes[name] = SurfaceComplex(
                name=name,
                reaction=rxn,
                logK=logK,
            )
            i += 1

        return i

    def _parse_aqueous_kinetics(self, lines: List[str], start_idx: int) -> int:
        """
        Parse:
          Begin aqueous kinetics
          ...
          End of aqueous kinetics

        Header line per entry:
          <Label> <N_laws> <'ReactionType'> <N_species>
            <coef1> <'Spec1'> ... <coefN> <'SpecN'> [<LogK>]

        Then N_laws lines:
          - TST/irreversible: <rate> <N_terms>  <'Spec'> <exp> ...
          - monod:            <kmax> <N_terms>  <'Spec'> <Khalf> ... ['Inhibition' <M> <'Spec'> <Kin> ...]
        """
        i, n = start_idx, len(lines)

        # Seek begin marker (block may be absent)
        while i < n and "begin aqueous kinetics" not in lines[i].strip().lower():
            i += 1

        # If we get to the end of the file without finding the block, return start_idx
        if i >= n:
            return start_idx

        i += 1  # past "Begin aqueous kinetics"
        while i < n:
            s = lines[i].strip()
            low = s.lower()
            if "end of aqueous kinetics" in low:
                return i + 1
            if _is_blank_or_comment(s):
                i += 1
                continue

            toks = _tokens(s)
            if not toks:
                i += 1
                continue

            # First line in an aqueous kinetics entry consists of:
            # <Aqueous reaction label> <Number of parallel rate laws> <'ReactionType'> <Number of
            # species in reaction (integer)> <Stoichiometric coefficient> <'SpeciesName'> <Stoichiometric
            # coefficient> <'SpeciesName'> … <Log K>
            label = _unquote(toks[0])
            try:
                n_laws = int(toks[1])
                rtype = _unquote(toks[2]).strip("'\"")
                n_rxn = int(toks[3])
            except Exception as exc:
                raise ParseError(f"Aqueous kinetics '{label}': Could not read rate law header.") from exc

            # Stoichiometry pairs
            stoich: Dict[str, float] = {}
            k = 4
            for _ in range(n_rxn):
                if k + 1 >= len(toks):
                    raise ParseError(f"Aqueous kinetics '{label}': incomplete reaction stoichiometry.")
                coef = float(toks[k])
                k += 1
                sp = _unquote(toks[k])
                k += 1
                stoich[sp] = coef

            # Read logK at the end of the line
            logK_val: Optional[float] = None
            if k < len(toks):
                try:
                    logK_val = float(toks[k])
                except ValueError:
                    logK_val = None

            # There could be multiple rate laws defined for this entry. By default, store the
            # first rate law on the AqueousKinetics object itself, and any additional laws
            # into 'parallel_laws' as dicts.
            primary_rate: Optional[float] = None
            primary_orders: Optional[Dict[str, float]] = None
            primary_monod: Optional[List[Tuple[str, float]]] = None
            primary_inhib: Optional[Dict[str, float]] = None
            parallel_laws: List[dict] = []

            # Read N rate-law lines
            j, parsed = i + 1, 0
            while parsed < n_laws and j < n:
                law_line = lines[j].strip()
                if _is_blank_or_comment(law_line):
                    j += 1
                    continue
                ltoks = _tokens(law_line)
                if not ltoks:
                    j += 1
                    continue

                try:
                    rate = float(ltoks[0])
                    n_terms = int(ltoks[1])
                except Exception as exc:
                    raise ParseError(f"Aqueous kinetics '{label}': malformed rate-law line.") from exc

                p = 2
                if rtype == "monod":
                    # Monod terms
                    monod_terms: List[Tuple[str, float]] = []
                    for _ in range(n_terms):
                        sp = _unquote(ltoks[p])
                        p += 1
                        khalf = float(ltoks[p])
                        p += 1
                        monod_terms.append((sp, khalf))

                    # Read 'Inhibition' line
                    j += 1
                    if j >= n:
                        raise ParseError(f"Aqueous kinetics '{label}': expected 'Inhibition' line for monod kinetics.")
                    inhibition: Dict[str, float] = {}
                    inhib_toks = _tokens(lines[j].strip())
                    if inhib_toks and _unquote(inhib_toks[0]) == "Inhibition":
                        try:
                            num_inhibitors = int(inhib_toks[1])
                        except Exception as exc:
                            raise ParseError(
                                f"Aqueous kinetics '{label}': could not read number of inhibitors {inhib_toks[1]}"
                            ) from exc
                        q = 2  # Counter for number of inhibitors
                        for _ in range(num_inhibitors):
                            isp = _unquote(inhib_toks[q])
                            q += 1
                            kin = float(inhib_toks[q])
                            q += 1
                            inhibition[isp] = kin
                    else:
                        raise ParseError(f"Aqueous kinetics '{label}': expected 'Inhibition' line for monod kinetics.")
                    law_dict = {"rate25C": rate, "monod_terms": monod_terms, "inhibition": inhibition}

                else:
                    # TST or irreversible: orders
                    orders: Dict[str, float] = {}
                    for _ in range(n_terms):
                        sp = _unquote(ltoks[p])
                        p += 1
                        exp = float(ltoks[p])
                        p += 1
                        orders[sp] = exp
                    law_dict = {"rate25C": rate, "orders": orders}

                # first parsed law goes to primary fields; others to parallel_laws
                if parsed == 0:
                    primary_rate = law_dict.get("rate25C")
                    if "orders" in law_dict:
                        primary_orders = law_dict["orders"]  # type: ignore[assignment]
                    if "monod_terms" in law_dict:
                        primary_monod = law_dict["monod_terms"]  # type: ignore[assignment]
                        primary_inhib = law_dict.get("inhibition", {})
                else:
                    parallel_laws.append(law_dict)

                parsed += 1
                j += 1

            if parsed != n_laws:
                raise ParseError(f"Aqueous kinetics '{label}': expected {n_laws} rate-law lines, found {parsed}.")

            self.aqueous_kinetics[label] = AqueousKinetics(
                label=label,
                reaction_type=rtype,  # type: ignore[assignment]
                reaction=Reaction(stoich=stoich),
                logK=logK_val,
                rate25C=primary_rate,
                orders=primary_orders,
                monod_terms=primary_monod,
                inhibition=primary_inhib,
                parallel_laws=parallel_laws,
            )

            i = j

        return i

    def _parse_mineral_kinetics(self, lines: List[str], start_idx: int) -> int:
        """
        Parse:
          Begin mineral kinetics
          ...
          End of mineral kinetics

        Each entry is delimited by a line beginning with '+-----', then:
          <MineralName>
          label = <label>
          type  = <tst|monod|irreversible|PrecipitationOnly|DissolutionOnly>
          rate(25C) = <log10 rate at 25C>
          activation = <Ea in kcal/mole>
        Followed by, depending on type:
          - TST / Irreversible / PrecipitationOnly / DissolutionOnly:
              dependence : [<species> <exponent> ...]
              [AffinityDependence = m1] [m2] [m3]   # may span multiple lines
              [<reaction equation line with '='>]
          - Monod:
              monod_terms : [<species> <Khalf> ...]
              inhibition :  [<species> <Kin>   ...]
        """

        def parse_pairs(tokens: List[str]) -> List[Tuple[str, float]]:
            """Return [(name, value)] from alternating tokens name, number."""
            out: List[Tuple[str, float]] = []
            i = 0
            while i + 1 < len(tokens):
                name = tokens[i]
                try:
                    val = float(tokens[i + 1])
                except ValueError:
                    break
                out.append((name, val))
                i += 2
            return out

        def parse_reaction_equation(eqn: str, equal_symbol="=") -> Optional[Reaction]:
            """
            Parse a simple 'a A + b B = c C + ...' into Reaction with
            negative coeffs on LHS, positive on RHS. Best-effort only.
            """
            try:
                lhs, rhs = eqn.split(equal_symbol, 1)
            except ValueError:
                return None

            def side_to_pairs(side: str, sign: float) -> List[Tuple[str, float]]:
                parts = [p.strip() for p in side.replace("  ", " ").split(" + ")]
                pairs: List[Tuple[str, float]] = []
                for term in parts:
                    if not term:
                        continue
                    toks = term.split()
                    if not toks:
                        continue
                    # try "<coef> <name>", else assume coef = 1
                    try:
                        coef = float(toks[0])
                        name = toks[1]
                    except (ValueError, IndexError):
                        coef = 1.0
                        name = toks[0]
                    pairs.append((name, sign * coef))
                return pairs

            stoich: Dict[str, float] = {}
            for name, coef in side_to_pairs(lhs, -1.0) + side_to_pairs(rhs, +1.0):
                stoich[name] = stoich.get(name, 0.0) + coef

            return Reaction(lhs="", stoich=stoich)

        i, n = start_idx, len(lines)

        # Seek begin marker (block may be absent)
        while i < n and "begin mineral kinetics" not in lines[i].strip().lower():
            i += 1
        if i >= n:
            return start_idx

        i += 1  # move past "Begin mineral kinetics"
        while i < n:
            line = lines[i].strip()
            low = line.lower()

            if "end of mineral kinetics" in low:
                return i + 1
            elif _is_blank_or_comment(line):
                i += 1
                continue

            # Expect separator then mineral name
            elif not line.startswith("+--"):
                # If file omits '+' separators, allow mineral name directly
                mineral_name = _strip_bang(line)
            else:
                # Skip separator lines
                i += 1
                while i < n and _is_blank_or_comment(lines[i]):
                    i += 1
                if i >= n:
                    return i
                mineral_name = _strip_bang(lines[i].strip())

            # If the final mineral kinetics block is terminated by "+----", then "End of mineral kinetics"
            # is accidentally picked up as a mineral name
            if mineral_name == "End of mineral kinetics":
                return i + 1

            # Now collect the 4 standard lines (label, type, rate, activation) in any whitespace form
            label: Optional[str] = None
            rtype: Optional[str] = None
            rate25: Optional[float] = None
            Ea: Optional[float] = None

            # Advance past mineral name line
            i += 1

            # Scan forward until we've read the four header fields, then parse tail until next '+'
            # or End-of-block.
            reaction_obj: Optional[Reaction] = None
            dependence: Dict[str, float] = {}
            affinity: Optional[Tuple[float, float, float]] = None
            monod_terms: List[Tuple[str, float]] = []
            inhibition: Dict[str, float] = {}

            while i < n:
                s = lines[i].strip()
                low = s.lower()

                if low.startswith("+--"):
                    break
                if _is_blank_or_comment(s):
                    i += 1
                    continue

                s_clean = _strip_bang(s)

                # Header keys (label / type / rate / activation)
                if low.startswith("label"):
                    # label = <text>
                    if "=" in s_clean:
                        label = s_clean.split("=", 1)[1].strip()
                elif low.startswith("type"):
                    if "=" in s_clean:
                        rtype = s_clean.split("=", 1)[1].strip()
                elif low.startswith("rate(25c)") or low.startswith("rate (25c)") or low.startswith("rate"):
                    # rate(25C) = <float> [comment]
                    if "=" in s_clean:
                        try:
                            rate25 = float(s_clean.split("=", 1)[1].split()[0])
                        except ValueError:
                            raise ParseError(f"Mineral kinetics '{mineral_name}': cannot parse rate(25C).")
                elif low.startswith("activation"):
                    # activation = <float> (kcal/mole)
                    if "=" in s_clean:
                        rhs = s_clean.split("=", 1)[1].strip()
                        try:
                            Ea = float(rhs.split()[0])
                        except ValueError:
                            raise ParseError(f"Mineral kinetics '{mineral_name}': cannot parse activation energy.")

                # Type-dependent fields
                elif low.startswith("dependence"):
                    # dependence : [pairs on same line] and possibly on following lines until blank/next key
                    toks = s_clean.split(": ", 1)
                    tail = toks[1].strip() if len(toks) > 1 else ""
                    if tail:
                        for sp, exp in parse_pairs(tail.split()):
                            dependence[sp] = exp
                    # also lookahead to consume additional dependence pairs on subsequent lines
                    j = i + 1
                    while j < n:
                        ns = _strip_bang(lines[j].strip())
                        nslow = ns.lower()
                        valid_keys = (
                            "label",
                            "type",
                            "rate",
                            "activation",
                            "monod_terms",
                            "inhibition",
                            "affinitydependence",
                        )
                        if (
                            not ns
                            or nslow.startswith(valid_keys)
                            or line.startswith("+--")
                            or "end of mineral kinetics" in nslow
                        ):
                            break
                        # Split reaction described below dependence line
                        for sp, exp in parse_pairs(ns.split()):
                            dependence[sp] = exp
                        j += 1
                    i = j - 1  # -1 because outer loop will i += 1

                elif low.startswith("monod_terms"):
                    # monod_terms : <species> <Khalf> ...
                    toks = s_clean.split(": ", 1)
                    tail = toks[1].strip() if len(toks) > 1 else ""
                    monod_terms = parse_pairs(tail.split())

                elif low.startswith("inhibition"):
                    toks = s_clean.split(": ", 1)
                    tail = toks[1].strip() if len(toks) > 1 else ""
                    inhibition = {name: val for name, val in parse_pairs(tail.split())}

                elif low.startswith("affinitydependence"):
                    # AffinityDependence = m1 then next lines m2 and m3 (possibly same line)
                    nums: List[float] = []
                    if "=" in s_clean:
                        rhs = s_clean.split("=", 1)[1]
                        for tok in rhs.split():
                            try:
                                nums.append(float(tok))
                            except ValueError:
                                pass
                    # pull following lines until we have 3 numbers
                    j = i + 1
                    while len(nums) < 3 and j < n:
                        ns = _strip_bang(lines[j].strip())
                        if not ns:
                            j += 1
                            continue
                        # stop if obvious new field/entry begins
                        nslow = ns.lower()
                        valid_keys = (
                            "label",
                            "type",
                            "rate",
                            "activation",
                            "monod_terms",
                            "inhibition",
                            "affinitydependence",
                        )
                        if (
                            not ns
                            or nslow.startswith(valid_keys)
                            or line.startswith("+--")
                            or "end of mineral kinetics" in nslow
                        ):
                            break
                        for tok in ns.split():
                            try:
                                nums.append(float(tok))
                            except ValueError:
                                pass
                        j += 1
                    if len(nums) >= 3:
                        affinity = (nums[0], nums[1], nums[2])
                    i = j - 1

                elif ("=" in s_clean or "-->" in s_clean) and "affinitydependence" not in low:
                    # a stoichiometric reaction line (contains '=' and is not AffinityDependence)
                    reaction_text = s_clean
                    if "=" in reaction_text:
                        reaction_obj = parse_reaction_equation(reaction_text)
                    elif "-->" in reaction_text:
                        reaction_obj = parse_reaction_equation(reaction_text, equal_symbol="-->")

                # else: unrecognized line within entry; ignore permissively
                i += 1

            # Sanity check required header fields
            if label is None or rtype is None or rate25 is None or Ea is None:
                print(len(self.mineral_kinetics))
                raise ParseError(f"Mineral kinetics '{mineral_name}': missing one of label/type/rate/activation.")

            mk = MineralKinetics(
                mineral=mineral_name,
                label=label,
                reaction_type=rtype,  # type: ignore[assignment]
                rate25_log10=rate25,
                activation_kcal_per_mol=Ea,
                dependence=dependence,
                affinity_dependence=affinity,
                monod_terms=monod_terms,
                inhibition=inhibition,
                reaction=reaction_obj,
            )

            # Store by mineral then label
            bucket = self.mineral_kinetics.setdefault(mineral_name, {})
            bucket[label] = mk

            # If we broke out because of a separator/end, do not consume that line here.
            # Outer loop will handle it.
            if i < n and (line.startswith("+--") or "end of mineral kinetics" in lines[i].strip().lower()):
                continue

            i += 1  # advance to next line after finishing this entry

        return i

    def _parse_exchange(self, lines: List[str], start_idx: int) -> int:
        """
        Parse 'Begin exchange' ... 'End of exchange'.
        Each entry expected like:
          'ExchangerName'  '<Number of species>'  <Stoichiometric coefficient 1>  <’SpeciesName 1’>
          <Stoichiometric coefficient 2> <’SpeciesName’ 2> ... <Log K> <Log K>
        """
        i = start_idx
        n = len(lines)

        # Seek begin marker
        while i < n and "begin exchange" not in lines[i].strip().lower():
            if ("end of " in lines[i].strip().lower()) or (lines[i].strip().startswith("'")):
                # allow other sections to start if this section is missing
                if "begin" in lines[i].strip().lower():
                    return i  # let the next parser start here
            i += 1

        if i >= n:
            # No exchange block
            return start_idx

        i += 1  # move past "Begin exchange"
        while i < n:
            s = lines[i].strip()
            low = s.lower()
            if "end of exchange" in low:
                return i + 1
            if _is_blank_or_comment(s):
                i += 1
                continue

            toks = _tokens(s)
            if not toks or not toks[0].startswith("'"):
                i += 1
                continue

            name = _unquote(toks[0])

            # Separate quoted species tokens from trailing numeric block
            quoted, nums = _separate_quoted_and_numeric(toks[1:])

            if not quoted or not nums:
                raise ParseError(f"Surface complex '{name}': Could not parse reaction stoichiometry.")

            # Format for minerals is:
            num_species = int(nums[0])

            if num_species > len(quoted):
                raise ParseError(f"'{name}': expected {num_species} species in reaction, found {len(quoted)}.")

            # Build the reaction from the quoted species
            # Define j, which is the starting index of the current section within nums
            j = 1
            rxn_species = quoted[:num_species]
            stoich_coefs = nums[j : num_species + j]
            rxn = Reaction(lhs=name, stoich={sp: coef for sp, coef in zip(rxn_species, stoich_coefs)})
            j += num_species

            # Get logK values
            logK = {0: nums[j], 1: nums[j + 1]}  # exchange uses 0 and 1 as keys

            self.exchange[name] = Exchange(
                name=name,
                reaction=rxn,
                logK=logK,
            )
            i += 1

        return i

    def _parse_surface_complexation_params(self, lines: List[str], start_idx: int) -> int:
        """Parse 'Begin surface complexation parameters' ... 'End surface complexation parameters'."""
        i = start_idx
        n = len(lines)
        # Seek begin marker
        while i < n and "begin surface complexation parameters" not in lines[i].strip().lower():
            i += 1

        if i >= n:
            # No surface complexation parameters
            return start_idx

        i += 1  # move past "Begin surface complexation parameters"
        while i < n:
            s = lines[i].strip()
            low = s.lower()
            if "end surface complexation parameters" in low:
                return i + 1
            if _is_blank_or_comment(s):
                i += 1
                continue

            toks = _tokens(s)
            if not toks:
                i += 1
                continue
            if len(toks) < 2 or not toks[0].startswith(">"):
                print(toks)
                raise ValidationError("Surface complexation parameters: expected format '>Name Param ...'")

            name = toks[0]
            param = float(toks[1])

            self.surface_complex_params[name] = SurfaceComplexParams(
                name=name,
                parameter=param,
            )
            i += 1

        return i

    # ----------------------------- Writing -------------------------------- #
    def _write_header(self, f) -> None:
        assert self.temperature_points
        grid = self.temperature_points
        f.write("'temperature points' ")
        f.write(f"{len(grid):d} ")
        f.write(" " + " ".join(_fmt_float(t, 0) for t in grid) + "\n")
        dh = self.activity
        if any(v != 0.0 for v in (dh.adh + dh.bdh + dh.bdt)):
            f.write("'Debye-Huckel adh' " + " ".join(_fmt_float(v, 7) for v in dh.adh) + "\n")
            f.write("'Debye-Huckel bdh' " + " ".join(_fmt_float(v, 7) for v in dh.bdh) + "\n")
            f.write("'Debye-Huckel bdt' " + " ".join(_fmt_float(v, 7) for v in dh.bdt) + "\n")

    def _write_primary(self, f) -> None:
        # First get length of longest primary species
        max_str_len = 0
        for name in self.primary:
            if len(name) > max_str_len:
                max_str_len = len(name)
        for name in self.primary:
            size = self.primary[name].size
            charge = self.primary[name].charge
            molar_mass = self.primary[name].molar_mass
            f.write(f"'{name}'".ljust(max_str_len + 2))
            f.write(f" {_fmt_float(size, 5)} {_fmt_float(charge, 6)} {_fmt_float(molar_mass)}\n")
        f.write("'End of primary  0.0  0.0  0.0'\n")

    def _write_secondary(self, f) -> None:
        if not self.secondary:
            return
        for name, sp in self.secondary.items():
            # Write name, then number of species in reaction
            line = [f"'{name}'", f"{len(sp.reaction.stoich):d}".ljust(3)]

            # Write reaction stoichiometry
            for reactant, coefficient in sp.reaction.stoich.items():
                line.append(_fmt_float(coefficient, 4))
                line.append(f"'{reactant}'")

            # Write logK values
            if len(sp.logK) != len(self.temperature_points):
                raise ValidationError(f"Secondary '{name}': logK array length does not match temperature grid.")
            vals = list(sp.logK.values())
            line.extend(_fmt_float(v) for v in vals)

            # Write Debye–Hückel size, charge, molar mass
            line.append(_fmt_float(sp.size, 5))
            line.append(_fmt_float(sp.charge, 6))
            line.append(_fmt_float(sp.molar_mass))
            f.write(" ".join(line) + "\n")
        f.write("'End of secondary' 1 0. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.\n")

    def _write_gases(self, f) -> None:
        if not self.gases:
            return
        for name, sp in self.gases.items():
            # Write name, then molar volume, then number of species in reaction
            line = [f"'{name}'", _fmt_float(sp.molar_volume, 6), f"{len(sp.reaction.stoich):d}".ljust(3)]

            # Write reaction stoichiometry
            for reactant, coefficient in sp.reaction.stoich.items():
                line.append(_fmt_float(coefficient, 4))
                line.append(f"'{reactant}'")

            # Write logK values
            if len(sp.logK) != len(self.temperature_points):
                raise ValidationError(f"Gas '{name}': logK array length does not match temperature grid.")
            vals = list(sp.logK.values())
            line.extend(_fmt_float(v) for v in vals)

            # Write molar mass
            line.append(_fmt_float(sp.molar_mass))
            f.write(" ".join(line) + "\n")
        f.write("'End of gases' 0. 1 1. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0.\n")

    def _write_minerals(self, f) -> None:
        if not self.minerals:
            return
        for name, sp in self.minerals.items():
            # Write name, then molar volume, then number of species in reaction
            line = [f"'{name}'", _fmt_float(sp.molar_volume, 6), f"{len(sp.reaction.stoich):d}".ljust(3)]

            # Write reaction stoichiometry
            for reactant, coefficient in sp.reaction.stoich.items():
                line.append(_fmt_float(coefficient, 4))
                line.append(f"'{reactant}'")

            # Write logK values
            if len(sp.logK) != len(self.temperature_points):
                raise ValidationError(f"Mineral '{name}': logK array length does not match temperature grid.")
            vals = list(sp.logK.values())
            line.extend(_fmt_float(v) for v in vals)

            # Write molar mass
            line.append(_fmt_float(sp.molar_mass))
            f.write(" ".join(line) + "\n")
        f.write("'End of minerals' 0. 1 0. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0.\n")

    def _write_surface_complexation(self, f) -> None:
        f.write("Begin surface complexation\n")
        if self.surface_complexes:
            for name, sp in self.surface_complexes.items():
                # Write name, then number of species in reaction
                line = [f"'{name}'", f"{len(sp.reaction.stoich):d}".ljust(3)]

                # Write reaction stoichiometry
                for reactant, coefficient in sp.reaction.stoich.items():
                    line.append(_fmt_float(coefficient, 4))
                    line.append(f"'{reactant}'")

                # Write logK values
                if len(sp.logK) != len(self.temperature_points):
                    raise ValidationError(f"Mineral '{name}': logK array length does not match temperature grid.")
                vals = list(sp.logK.values())
                line.extend(_fmt_float(v) for v in vals)
                f.write(" ".join(line) + "\n")
        f.write("End of surface complexation\n")

    def _write_aqueous_kinetics(self, f) -> None:
        """Write the aqueous kinetics block."""
        f.write("Begin aqueous kinetics\n")
        if self.aqueous_kinetics:
            for label, entry in self.aqueous_kinetics.items():
                rtype = entry.reaction_type
                stoich = entry.reaction.stoich if entry.reaction else {}
                n_rxn = len(stoich)
                n_laws = 1 + len(getattr(entry, "parallel_laws", []))

                # Header line
                parts: List[str] = [f"{label}", f"{n_laws:d}", f"'{rtype}'", f"{n_rxn:d}"]
                for sp, coef in stoich.items():
                    parts.append(_fmt_float(coef))
                    parts.append(f"'{sp}'")
                if entry.logK is not None:
                    parts.append(_fmt_float(entry.logK))
                f.write(" ".join(parts) + "\n")

                # Define helper function to write a rate law
                def _write_law(
                    rate: float,
                    orders: Optional[Dict[str, float]] = None,
                    monod_terms: Optional[List[Tuple[str, float]]] = None,
                    inhibition: Optional[Dict[str, float]] = None,
                ) -> None:
                    if rtype == "monod":
                        mt = monod_terms or []
                        line = [_fmt_float(rate), f"{len(mt):d}"]
                        for sp, khalf in mt:
                            line.append(f"'{sp}' {_fmt_float(khalf)}")
                        inhib = inhibition or {}
                        if inhib:
                            line.append(f"\n   'Inhibition'  {len(inhib):d}")
                            for isp, kin in inhib.items():
                                line.append(f"'{isp}' {_fmt_float(kin)}")
                        f.write(" ".join(line) + "\n")
                    else:
                        od = orders or {}
                        line = [_fmt_float(rate), f"{len(od):d}"]
                        for sp, exp in od.items():
                            line.append(f"'{sp}' {_fmt_float(exp)} ")
                        f.write(" ".join(line) + "\n")

                # Write primary rate law
                if rtype == "monod":
                    _write_law(entry.rate25C, monod_terms=entry.monod_terms, inhibition=entry.inhibition)
                else:
                    _write_law(entry.rate25C, orders=entry.orders)

                # Write any parallel rate laws
                for law in getattr(entry, "parallel_laws", []):
                    if rtype == "monod":
                        _write_law(
                            law["rate25C"],
                            monod_terms=law.get("monod_terms") or [],
                            inhibition=law.get("inhibition") or {},
                        )
                    else:
                        _write_law(law["rate25C"], orders=law.get("orders") or {})
        f.write("End of aqueous kinetics\n")

    def _write_mineral_kinetics(self, f) -> None:
        r"""
        Write the 'Begin mineral kinetics' block using MineralKinetics records.

        Typical format:
          +----------------------------------------------------
          <MineralName>
            label = <label>
            type  = <tst|monod|irreversible|PrecipitationOnly|DissolutionOnly>
            rate(25C) = <log10 rate at 25C>
            activation = <Ea>  (kcal/mole)
            [dependence : <species> <exp> ...]             # for TST/Irrev/PrecOnly/DissOnly
            [AffinityDependence = m1 \n  m2 \n  m3]        # optional, those types only
            [<stoichiometric reaction line with '='>]      # optional
            -- OR (for monod) --
            monod_terms : [<species> <Khalf> ...]
            inhibition :  [<species> <Kin>   ...]
          +----------------------------------------------------
          ...
          End of mineral kinetics
        """
        if not getattr(self, "mineral_kinetics", None):
            return

        f.write("Begin mineral kinetics\n")
        sep = "+---------------------------------------------------"

        # Stable order: mineral then label
        for mineral in sorted(self.mineral_kinetics.keys()):
            entries = self.mineral_kinetics[mineral]
            for label in sorted(entries.keys()):
                mk = entries[label]

                f.write(sep + "\n")
                f.write(f"{mk.mineral}\n")
                f.write(f"  label = {mk.label}\n")
                f.write(f"  type  = {mk.reaction_type}\n")
                f.write(f"  rate(25C) = {_fmt_float(mk.rate25_log10, 1)}\n")
                f.write(f"  activation = {_fmt_float(mk.activation_kcal_per_mol, 1)}  (kcal/mole)\n")

                if mk.reaction_type == "monod":
                    # Monod terms
                    f.write("  monod_terms :")
                    if mk.monod_terms:
                        f.write("  " + "  ".join(f"{sp} {khalf:g}" for sp, khalf in mk.monod_terms) + "\n")
                    else:
                        f.write("\n")

                    # Inhibition terms
                    f.write("  inhibition :")
                    if mk.inhibition:
                        f.write("  " + "  ".join(f"{sp} {kin:g}" for sp, kin in mk.inhibition.items()) + "\n")
                    else:
                        f.write("\n")

                    # Write inhibition reaction
                    if mk.reaction and mk.reaction.stoich:
                        # Reconstruct: negatives on LHS, positives on RHS
                        lhs_terms = []
                        rhs_terms = []
                        for sp, coef in mk.reaction.stoich.items():
                            if coef < 0:
                                lhs_terms.append(f"{_fmt_float(abs(coef), 3)} {sp}")
                            elif coef > 0:
                                rhs_terms.append(f"{_fmt_float(coef, 3)} {sp}")
                        left = " + ".join(lhs_terms) if lhs_terms else "0"
                        right = " + ".join(rhs_terms) if rhs_terms else "0"
                        f.write(f"  {left} -->  {right}\n")
                else:
                    # Dependence (species orders)
                    f.write("  dependence :")
                    if mk.dependence:
                        dependence_parts = "  ".join(f"{sp} {_fmt_float(exp, 3)}" for sp, exp in mk.dependence.items())
                        f.write("  " + dependence_parts + "\n")
                    else:
                        f.write("\n")

                    # Optional AffinityDependence (m1, m2, m3), printed over 3 lines like manual examples
                    if mk.affinity_dependence:
                        m1, m2, m3 = mk.affinity_dependence
                        f.write(f"  AffinityDependence = {m1:g}\n")
                        f.write(f"   {m2:g}\n")
                        f.write(f"   {m3:g}\n")

                    # Write dependence reaction
                    if mk.reaction and mk.reaction.stoich:
                        # Reconstruct: negatives on LHS, positives on RHS
                        lhs_terms = []
                        rhs_terms = []
                        for sp, coef in mk.reaction.stoich.items():
                            if coef < 0:
                                lhs_terms.append(f"{_fmt_float(abs(coef), 3)} {sp}")
                            elif coef > 0:
                                rhs_terms.append(f"{_fmt_float(coef, 3)} {sp}")
                        left = " + ".join(lhs_terms) if lhs_terms else "0"
                        right = " + ".join(rhs_terms) if rhs_terms else "0"
                        f.write(f"  {left} =  {right}\n")

        f.write(sep + "\n")
        f.write("End of mineral kinetics\n")

    def _write_exchange(self, f) -> None:
        if not self.exchange:
            return
        f.write("Begin exchange\n")
        # First get length of longest primary species
        max_str_len = 0
        for name in self.exchange:
            if len(name) > max_str_len:
                max_str_len = len(name)

        for name, sp in self.exchange.items():
            # Write name, then number of species in reaction
            line = [f"'{name}'".ljust(max_str_len + 2), f"{len(sp.reaction.stoich):d}".ljust(3)]

            # Write reaction stoichiometry
            for reactant, coefficient in sp.reaction.stoich.items():
                line.append(_fmt_float(coefficient, 4))
                line.append(f"'{reactant}'")

            # Write logK values
            vals = list(sp.logK.values())
            line.extend(_fmt_float(v, 6) for v in vals)

            f.write(" ".join(line) + "\n")
        f.write("'End of exchange\n")

    def _write_surface_complexation_params(self, f) -> None:
        if not self.exchange:
            return
        f.write("Begin surface complexation parameters\n")
        # First get length of longest primary species
        max_str_len = 0
        for name in self.surface_complex_params:
            if len(name) > max_str_len:
                max_str_len = len(name)

        for name, sp in self.surface_complex_params.items():
            # Write name, then number of species in reaction
            line = f"{name}".ljust(max_str_len + 2)
            line += f" {_fmt_float(sp.parameter, 5)}\n"
            f.write(line)
        f.write("'End surface complexation parameters\n")
