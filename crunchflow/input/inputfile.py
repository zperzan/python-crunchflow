"""A module for creating and saving CrunchFlow input files."""

import os
from importlib.metadata import version

from crunchflow.input.blocks import (
    AqueousKinetics,
    BoundaryConditions,
    Condition,
    DatabaseBlock,
    Discretization,
    Erosion,
    Flow,
    Gases,
    InitialConditions,
    IonExchange,
    Isotopes,
    Minerals,
    Output,
    Pest,
    Porosity,
    PrimarySpecies,
    Runtime,
    SecondarySpecies,
    SurfaceComplexation,
    Temperature,
    Title,
    Transport,
)


def format_class_name(name):
    """Convert a class name to the appropriate format for writing a block
    keyword in an input file. If the name is in camel case, an underscore is
    inserted between words. The name is then converted to uppercase.

    Parameters
    ----------
    name : str
        The name of the class to convert.

    Returns
    -------
    str
        The class name converted to uppercase with underscores between words.
    """
    import re

    # DatabaseBlock is a special case, because we want to distinguish it from the
    # separate Database() class for manipulating databases
    if name == "DatabaseBlock":
        return "DATABASE"
    else:
        # Regex patter recognizes lower case followed by upper case, so
        # "InitialConditions" becomes "INITIAL_CONDITIONS"
        return re.sub(r"(?<!^)(?=[A-Z])", "_", name).upper()


class InputFile:
    """The main class for creating and saving CrunchFlow input files."""

    def __init__(self):
        self.title = Title()
        self.runtime = Runtime()
        self.database_block = DatabaseBlock()
        self.output = Output()
        self.discretization = Discretization()
        self.flow = Flow()
        self.transport = Transport()
        self.primary_species = PrimarySpecies()
        self.secondary_species = SecondarySpecies()
        self.minerals = Minerals()
        self.gases = Gases()
        self.initial_conditions = InitialConditions()
        self.boundary_conditions = BoundaryConditions()
        self.ion_exchange = IonExchange()
        self.surface_complexation = SurfaceComplexation()
        self.aqueous_kinetics = AqueousKinetics()
        self.conditions = {}
        self.porosity = Porosity()
        self.temperature = Temperature()
        self.pest = Pest()
        self.erosion = Erosion()
        self.isotopes = Isotopes()
        # Initialize other blocks as needed

    def set_block(self, block_name, parameters):
        """Set the parameters of a block in an InputFile instance.

        Parameters
        ----------
        block_name : str
            The name of the block to set.
        parameters : dict
            A dictionary of parameters to set for the block.
        """
        if block_name == "condition":
            name = parameters.pop("name", None)
            if name:
                condition = Condition(name)
                condition.set_parameters(parameters)
                self.conditions[name] = condition
        else:
            block = getattr(self, block_name)
            block.set_parameters(parameters)

    def __str__(self):
        """Return a string representation of the CrunchFlow input file."""
        result = []
        for attr, value in self.__dict__.items():
            # Conditions are handled as dict since there can be multiple
            # condition blocks, each with a different name
            if isinstance(value, dict):
                for condition in value.values():
                    if any(val for val in condition.__dict__.values() if val):
                        # Note that the "CONDITION  <name>" is included in
                        # Condition.__str__ method, so it is omitted here
                        result.append(f"{str(condition)}\nEND")

            # Otherwise, check if the block has any attributes set
            # and if so, format the block name then print it
            elif any(val for val in value.__dict__.values() if val):
                block_name = format_class_name(value.__class__.__name__)
                result.append(f"{block_name}\n{str(value)}\nEND")
        return "\n\n".join(result)

    def save(self, filename, path=".", update_pestcontrol=False):
        """Write this CrunchFlow run to an input file.

        Parameters
        ----------
        filename : str
            The name of the output file to save.
        path : str, optional
            The path to the output file. Default is the current directory.
        update_pestcontrol : bool, optional
            Whether to update PestControl.ant with the name of the
            CrunchFlow input file. Default is False.

        Returns
        -------
        None
            The CrunchFlow input file is saved to disk.
        """
        full_path = os.path.join(path, filename)
        with open(full_path, "w") as file:
            file.write("! CrunchFlow input file\n")
            cf_version = version("crunchflow")
            file.write("! Generated automatically by python-crunchflow v%s\n" % cf_version)
            file.write(str(self) + "\n")

        # Update PestControl.ant if requested
        if update_pestcontrol:
            folder = os.path.dirname(full_path)
            pestfile = os.path.join(folder, "PestControl.ant")
            with open(pestfile, "w") as file:
                file.write("%s \n" % os.path.basename(full_path))

    @classmethod
    def load(cls, filename, path=".", warnings=True):
        """Read a CrunchFlow input file and create a Run instance.

        Parameters
        ----------
        filename : str
            The name of the input file to read.
        path : str, optional
            The path to the input file. Default is the current directory.
        warnings : bool, optional
            Whether to print warnings for unrecognized blocks or attributes.
            Default is True.

        Returns
        -------
        InputFile
            A Run instance with the data read from the input file.
        """
        instance = cls()
        full_path = os.path.join(path, filename)

        with open(full_path, "r") as file:
            lines = file.readlines()

        current_block = None
        current_block_name = None

        # Define attributes and blocks to be handled as special cases
        # SpeciesBlock and KineticsBlock instances are handled differently below
        species_blocks = ["primary_species", "secondary_species", "gases"]
        kinetics_blocks = ["minerals", "aqueous_kinetics"]

        # Some attributes can be set multiple times within a single block
        multiply_defined = [
            "time_series",
            "pressure",
            "mineral",
            "primary",
            "D_25",
            "tortuosityMP",
            "permeability_x",
            "permeability_y",
            "permeability_z",
        ]

        # List of condition attributes that are not species
        condition_attributes = ["units", "equilibrate_surface", "temperature", "set_porosity", "set_saturation"]

        for line in lines:
            line = line.strip()

            # Once stripped, each line should be one of four categories:
            # (1) empty, commented
            # (2) "END"
            # (3) a block name
            # (4) otherwise, we assume we are within a keyword block

            # Category (1): empty or commented lines
            if line.startswith("!") or line.startswith("#") or not line:
                continue

            # Category (2): "END" line
            if line.endswith("END"):
                current_block = None
                current_block_name = None
                continue

            # Category (3): block name
            if not current_block:
                parts = line.split(maxsplit=2)
                block_name = parts[0].lower()
                if block_name == "condition" and len(parts) > 1:
                    condition_name = parts[1]
                    condition = Condition(condition_name)
                    instance.conditions[condition_name] = condition
                    current_block = condition
                    current_block_name = "condition"
                elif block_name == "title":
                    current_block = instance.title
                    current_block_name = block_name
                elif block_name == "database":
                    current_block = instance.database_block
                    current_block_name = block_name
                elif hasattr(instance, block_name):
                    current_block = getattr(instance, block_name)
                    current_block_name = block_name
                else:
                    if warnings:
                        print(f"\tWarning: Unrecognized block name '{block_name}'")
                continue

            # Category (4): within a keyword block
            if current_block and current_block_name:
                # First, take care of the special cases: title, initial_conditions,
                # database, species_blocks (see above) and kinetics_blocks (see above)
                if current_block_name in ["title", "database"]:
                    current_block.set_parameters(line)

                elif current_block_name == "initial_conditions":
                    instance.initial_conditions.conditions.append(line.strip())

                elif current_block_name in species_blocks:
                    parts = line.split()
                    current_block.species.append(parts[0])

                elif current_block_name == "surface_complexation":
                    parts = line.split()
                    species = parts[0]
                    details = " ".join(parts[1:]) if len(parts) > 1 else ""
                    current_block.species.append(species)
                    current_block.species_dict[species] = details

                elif current_block_name in kinetics_blocks:
                    parts = line.split()
                    species = parts[0]
                    details = {}

                    # Assume the default label
                    # This will be updated below if it is included in details
                    label = "default"

                    # If there are details associated with the species_dict, then parse them
                    if len(parts) > 1:
                        # Loop through kinetic options (e.g., -activation, -rate, etc.)
                        # and store this information in a dictionary
                        for i in range(1, len(parts), 2):
                            key = parts[i].lstrip("-")
                            value = parts[i + 1] if (i + 1) < len(parts) else None
                            if key == "label":
                                label = value
                            else:
                                details[key] = value
                    current_block.set_parameters({species: {**details, "label": label}})

                # All other blocks, split on whitespace
                else:
                    parts = line.split()
                    if len(parts) > 1:
                        attribute = parts[0]
                        value = " ".join(parts[1:])

                        if current_block_name != "condition":
                            # In everything but conditions block, replace hyphens
                            # with underscores
                            attribute = attribute.replace("-", "_")

                        if current_block_name == "condition" and attribute not in condition_attributes:
                            current_block.concentrations[attribute] = value
                            current_block.species.append(attribute)
                        elif attribute in multiply_defined:
                            cur_val = getattr(current_block, attribute)
                            cur_val.append(value)
                        elif hasattr(current_block, attribute):
                            setattr(current_block, attribute, value)
                        else:
                            # If the attribute isn't found, check if there's a capitalization issue
                            # and try to set the attribute with the correct capitalization
                            lower_attrs = [attr.lower() for attr in current_block.__dict__.keys()]
                            if attribute.lower() in lower_attrs:
                                correct_attribute = [
                                    attr for attr in current_block.__dict__.keys() if attr.lower() == attribute.lower()
                                ][0]
                                setattr(current_block, correct_attribute, value)
                            # If it still can't be found, issue a warning and save it in the
                            # "other" attribute of the KeywordBlock
                            else:
                                if warnings:
                                    print(
                                        f"\tWarning: Unrecognized attribute '{attribute}' "
                                        f"in block '{current_block_name}'"
                                    )
                                current_block.other[attribute] = value
                    # If time_series_print is passed on its own without a list of species
                    # then set the time_series_print attribute to True, so that all species
                    # are printed to file
                    elif parts[0] in ["time_series_print"]:
                        current_block.time_series_print = True
                    else:
                        if warnings:
                            print(
                                f"\tWarning: Attribute '{line}' in block '{current_block_name}' "
                                f"does not have a value associated with it"
                            )

        return instance
