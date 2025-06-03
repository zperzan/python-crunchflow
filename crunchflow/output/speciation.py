"""
Module for parsing speciation blocks from CrunchFlow output files.

Provides the Speciation class, which extracts geochemical species data and
associated properties from each condition included within an input file.
"""

import pandas as pd


def get_float(line):
    """
    Extract a float value from a line formatted like 'Label = Value'.

    Parameters
    ----------
    line : str
        The line to parse.

    Returns
    -------
    float or None
        The extracted float, or None if conversion fails.
    """
    parts = line.split("=")
    if len(parts) > 1:
        try:
            return float(parts[1].strip())
        except ValueError:
            return None
    return None


class Speciation:
    """
    A class for parsing and storing the results from a single speciation block
    within a CrunchFlow output file.

    This class processes the full text of a speciation block and extracts:
    - Primary species and their total concentrations.
    - Secondary species with molality, activity, and activity coefficient.
    - Gas species and their partial pressures.
    - Mineral phases and their saturation indices.
    - Exchange sites, exchangers, and surface complexation data.
    - General geochemical properties such as temperature, porosity, pH, etc.

    Parameters
    ----------
    block_text : str
        The full text of a speciation block from a CrunchFlow output file.

    Attributes
    ----------
    totcon : pandas.DataFrame
        Total molality of primary species. Index is species name.

    conc : pandas.DataFrame
        Molality, activity, and activity coefficient of secondary species.

    gas_conc : pandas.DataFrame
        Partial pressure of gas species, in bars.

    saturation : pandas.DataFrame
        Saturation indices of minerals, in log format.

    exchangers : pandas.DataFrame
        Equivalent concentrations of exchange sites.

    exchange_site_conc : pandas.DataFrame
        Concentrations of primary species in exchange sites.

    surface_complex : pandas.DataFrame
        Concentrations of surface complex species.

    temperature, porosity, pH, pe, eh, total_charge : float
        Scalar geochemical properties extracted from the block.

    primary_species, secondary_species, minerals, gases : list of str
        Lists of species in each category.

    Methods
    -------
    __str__()
        Returns a summary string with primary and secondary species counts.
    """

    def __init__(self, block_text):
        """
        Initialize a Speciation object from a block of CrunchFlow output text.

        Parameters
        ----------
        block_text : str
            The full text content of a single speciation block.
        """
        self.block_text = block_text
        self.primary_species = []
        self.secondary_species = []
        self.minerals = []
        self.gases = []
        self.totcon = pd.DataFrame(columns=["molality"])
        self.conc = pd.DataFrame(columns=["molality", "activity", "activity_coefficient"])
        self.gas_conc = pd.DataFrame(columns=["partial_pressure"])
        self.saturation = pd.DataFrame(columns=["saturation_index"])
        self.exchangers = pd.DataFrame(columns=["equiv/kgw", "equiv/g solid", "equiv/m^3 bulk"])
        self.exchange_site_conc = pd.DataFrame(columns=["mol/kgw", "mol/g solid", "equiv/g solid"])
        self.surface_complex = pd.DataFrame(columns=["Sites/kgw", "Moles/g solid", "Moles/m^3 bulk"])

        self.temperature = None
        self.porosity = None
        self.liquid_saturation = None
        self.liquid_density = None
        self.solid_density = None
        self.solid_solution_ratio = None
        self.ionic_strength = None
        self.pH = None
        self.pe = None
        self.eh = None
        self.total_charge = None
        self.conversion = None

        self._parse_block()

    def _parse_block(self):
        """
        Parse all lines in the block text and dispatch them to the appropriate section
        handlers. This identifies the current section of the block and triggers parsing
        of temperature, concentrations, etc.
        """
        lines = [line.strip() for line in self.block_text.splitlines() if line.strip()]
        section = "initial"

        for line in lines:
            if "Total Aqueous Concentrations of Primary Species" in line:
                section = "primary"
            elif "Concentrations of Individual Species, Exchangers, and Surface" in line:
                section = "secondary"
            elif "Partial pressure of gases" in line:
                section = "gas"
            elif "Saturation state of minerals" in line:
                section = "mineral"
            elif line.strip().startswith("Exchangers"):
                section = "exchangers"
            elif "Total Concentrations in Exchange Sites" in line:
                section = "exchange sites"
            elif line.strip().startswith("Surface complex"):
                section = "surface_complex"

            # For each section, capture data within that section
            elif section == "initial":
                # Capture initial metadata
                if "Temperature (C)" in line:
                    self.temperature = get_float(line)
                elif "Porosity" in line:
                    self.porosity = get_float(line)
                elif "Liquid Saturation" in line:
                    self.liquid_saturation = get_float(line)
                elif "Liquid Density" in line:
                    self.liquid_density = get_float(line)
                elif "Solid Density" in line:
                    self.solid_density = get_float(line)
                elif "Solid:Solution Ratio" in line:
                    self.solid_solution_ratio = get_float(line)
                elif "Ionic Strength" in line:
                    self.ionic_strength = get_float(line)
                elif "Solution pH" in line:
                    self.pH = get_float(line)
                elif "Solution pe" in line:
                    self.pe = get_float(line)
                elif "Solution Eh" in line:
                    self.eh = get_float(line)
                elif "Total Charge" in line:
                    self.total_charge = get_float(line)
                elif "Conversion (M->m)" in line:
                    self.conversion = get_float(line)
            elif section == "primary":
                if line.strip("-") and not line.startswith("Species"):
                    parts = line.split()
                    self.totcon.loc[parts[0]] = float(parts[1])
            elif section == "secondary":
                if line.strip("-") and not line.startswith("Species") and not line.startswith("Log"):
                    parts = line.split()
                    # Determine whether this is an exchanger or an aqueous species
                    if parts[-1] == "Exchange":
                        self.conc.loc[parts[0], "molality"] = parts[3]
                        self.conc.loc[parts[0], "activity"] = parts[4]
                    elif parts[-1] == "Aqueous":
                        self.conc.loc[parts[0]] = [float(parts[i]) for i in [3, 4, 5]]
            elif section == "gas":
                if line.strip(""):
                    parts = line.split()
                    self.gas_conc.loc[parts[0]] = float(parts[1])
            elif section == "mineral":
                if (
                    line.strip("")
                    and not line.startswith("SPECIATION OF")
                    and not line.startswith("INITIAL AND BOUNDARY CONDITIONS")
                ):
                    parts = line.split()
                    self.saturation.loc[parts[0]] = float(parts[1])
            elif section == "exchangers":
                if line.strip("-") and not line.strip().startswith("Exchangers"):
                    parts = line.split()
                    self.exchangers.loc[parts[0]] = [float(parts[i]) for i in [1, 2, 3]]
            elif section == "exchange sites":
                if line.strip("-") and not line.strip().startswith("Primary Species"):
                    parts = line.split()
                    self.exchange_site_conc.loc[parts[0]] = [float(parts[i]) for i in [1, 2, 3]]
            elif section == "surface_complex":
                if line.strip("-") and "Total Concentrations on Surface" not in line and "Primary Species" not in line:
                    parts = line.split()
                    self.surface_complex.loc[parts[0]] = [float(parts[i]) for i in [1, 2, 3]]

        # Get primary and secondary species
        self.primary_species = self.totcon.index.tolist()
        self.secondary_species = self.conc.index.tolist()
        self.minerals = self.saturation.index.tolist()
        self.gases = self.gas_conc.index.tolist()

    def __str__(self):
        """
        Return a short summary string for the speciation block.

        Includes counts of primary and secondary species parsed.
        """
        output_str = f"SpeciationBlock with {len(self.primary_species)} primary "
        output_str += f"species and {len(self.secondary_species)} secondary species."

        return output_str
