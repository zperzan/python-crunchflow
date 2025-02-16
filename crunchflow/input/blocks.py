"""A submodule for defining classes that represent blocks in a CrunchFlow input file."""


class KeywordBlock:
    """Base class for all keyword blocks in a CrunchFlow input file"""

    def __init__(self):
        self.other = {}

    def set_parameters(self, parameters):
        """Set the parameters for the block

        Parameters
        ----------
        parameters : dict
            Dictionary of parameters to set
        """
        for key, value in parameters.items():
            key = key.lower()
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                self.other[key] = value

    def __str__(self):
        """Return a string representation of the block"""
        result = []
        for attr, value in self.__dict__.items():
            # Only print values that are not empty strings, empty lists or None
            # This allows users to set false values as either False or 'false'
            if attr != "other" and (value or isinstance(value, bool)):
                if isinstance(value, list):
                    result.append(f"{attr:<22} {' '.join(map(str, value))}")
                elif isinstance(value, bool):
                    result.append(f"{attr:<22} {str(value).lower()}")
                else:
                    result.append(f"{attr:<22} {value}")
        for key, value in self.other.items():
            result.append(f"{key:<22} {value}")
        return "\n".join(result)


class Title(KeywordBlock):
    """Title block for the input file"""

    def __init__(self):
        super().__init__()
        self.title = ""

    def set_parameters(self, parameters):
        """Set parameters for the Title block"""
        self.title = parameters.strip()

    def __str__(self):
        """Return the title as a string"""
        return self.title


class DatabaseBlock(KeywordBlock):
    """Typically, the database is specified in the Runtime() block, but
    some input files use a separate DATABASE block. This class is not to be
    confused with the Database() class for manipulating geochemical databases
    """

    def __init__(self):
        super().__init__()
        self.database = ""

    def set_parameters(self, parameters):
        """Set the database for the block"""
        self.database = parameters.strip()

    def __str__(self):
        """Return the database block as a string"""
        return self.database


class Runtime(KeywordBlock):
    """Runtime block for the input file

    This block contains parameters that control the runtime of the simulation
    """

    def __init__(self):
        super().__init__()
        self.aqueousdatabase = ""
        self.Benchmark = ""
        self.ChromeTopes = ""
        self.coordinate = ""
        self.coordinates = ""
        self.correction_max = ""
        self.courant_number = ""
        self.database = ""
        self.database_sweep = ""
        self.debye_huckel = ""
        self.density_module = ""
        self.dissolution_max = ""
        self.Duan = ""
        self.fix_saturation = ""
        self.generic_rates = ""
        self.gimrt = ""
        self.gimrt_pc = ""
        self.gimrt_pclevel = ""
        self.gimrt_solver = ""
        self.graphics = ""
        self.hindmarsh = ""
        self.JennyDruhan = ""
        self.kinetic_database = ""
        self.lag_activity = ""
        self.later_inputfiles = ""
        self.master = ""
        self.master_variable = ""
        self.OvershootTolerance = ""
        self.pc = ""
        self.pclevel = ""
        self.preconditioner = ""
        self.precondition_level = ""
        self.Qingyun = ""
        self.reaction_path = ""
        self.read_saturationfile = ""
        self.ResidualTolerance = ""
        self.restart = ""
        self.save_restart = ""
        self.screen_output = ""
        self.SetSurfaceAreaConstant = ""
        self.set_saturation = ""
        self.solver = ""
        self.speciate_only = ""
        self.steady_state = ""
        self.timestep_max = ""
        self.timestep_init = ""
        self.time_tolerance = ""
        self.time_units = ""

    # Redefine __str__ to ensure certain underscores are replaced with
    # hyphens in the output (e.g., debye_huckel -> debye-huckel)
    def __str__(self):
        """Return the runtime block as a string"""
        result = []
        for attr, value in self.__dict__.items():
            if attr == "debye_huckel" and value:
                result.append(f"debye-huckel           {value}")
            elif attr != "other" and (value or isinstance(value, bool)):
                if isinstance(value, list):
                    result.append(f"{attr:<22} {' '.join(map(str, value))}")
                elif isinstance(value, bool):
                    result.append(f"{attr:<22} {str(value).lower()}")
                else:
                    result.append(f"{attr:<22} {value}")
        for key, value in self.other.items():
            result.append(f"{key:<22} {value}")
        return "\n".join(result)


class Output(KeywordBlock):
    """Output block for the input file

    This block contains parameters that control the output of the simulation
    """

    def __init__(self):
        super().__init__()
        self.FluxWeightedConcentrationSpecies = ""
        self.ime_series_interval = ""
        self.MakeMovie = ""
        self.spatial_profile = ""
        self.time_units = ""
        self.time_series = []
        self.time_series_at_node = ""
        self.time_series_interval = ""
        self.time_series_output = ""
        self.time_series_print = ""
        self.time_series_units = ""
        self.WriteFluxWeightedConcentration = ""

    def set_parameters(self, parameters):
        """Set the parameters for the Output block"""
        for key, value in parameters.items():
            key = key.lower()
            if key == "time_series":
                self.time_series.append(value)
            elif hasattr(self, key):
                setattr(self, key, value)
            else:
                self.other[key] = value

    def __str__(self):
        """Return the output block as a string, as it would be represented
        in a CrunchFlow input file
        """
        result = []
        for attr, value in self.__dict__.items():
            exceptions = ["other", "time_series", "time_series_print"]
            if attr not in exceptions and (value or isinstance(value, bool)):
                if isinstance(value, list):
                    result.append(f"{attr:<22} {' '.join(map(str, value))}")
                elif isinstance(value, bool):
                    result.append(f"{attr:<22} {str(value).lower()}")
                else:
                    result.append(f"{attr:<22} {value}")
        if self.time_series_print:
            # If time_series_print == True, then print it on its own line
            # to ensure all species are printed to file
            if isinstance(self.time_series_print, bool) and self.time_series_print:
                result.append("time_series_print")
            else:
                result.append(f"time_series_print      {self.time_series_print}")
        for ts in self.time_series:
            result.append(f"time_series            {ts}")
        for key, value in self.other.items():
            result.append(f"{key:<22} {value}")
        return "\n".join(result)


class Discretization(KeywordBlock):
    """Discretization block for the input file"""

    def __init__(self):
        super().__init__()
        self.distance_units = ""
        self.xzones = ""
        self.yzones = ""
        self.zzones = ""


class IonExchange(KeywordBlock):
    """Ion exchange block for the input file"""

    def __init__(self):
        super().__init__()
        self.exchange = ""
        self.convention = ""


class Condition(KeywordBlock):
    """Condition block for the input file"""

    def __init__(self, name):
        super().__init__()
        self.concentrations = {}
        self.constraints = {}
        self.equilibrate_surface = ""
        self.name = name
        self.set_porosity = ""
        self.set_saturation = ""
        self.species = []
        self.temperature = ""
        self.units = ""

    def set_parameters(self, parameters):
        """Set the parameters for the Condition block"""
        for key, value in parameters.items():
            condition_attributes = ["units", "equilibrate_surface", "temperature", "set_porosity", "set_saturation"]
            if hasattr(self, key):
                setattr(self, key, value)
            elif key in condition_attributes:
                setattr(self, key, value)
            else:
                self.concentrations[key] = value
        self.species = list(self.concentrations.keys())

    def __str__(self):
        """Return the condition block as a string"""
        result = [f"CONDITION              {self.name}"]
        if self.units:
            result.append(f"units                  {self.units}")
        if self.equilibrate_surface:
            result.append(f"equilibrate_surface    {self.equilibrate_surface}")
        if self.temperature:
            result.append(f"temperature            {self.temperature}")
        if self.set_porosity:
            result.append(f"set_porosity           {self.set_porosity}")
        if self.set_saturation:
            result.append(f"set_saturation         {self.set_saturation}")
        for key, value in self.concentrations.items():
            result.append(f"{key:<22} {value}")
        return "\n".join(result)

    def __repr__(self):
        """Return the name of the condition"""
        return self.name


class Transport(KeywordBlock):
    """Transport block for the input file"""

    def __init__(self):
        super().__init__()
        self.anisotropy_ratioY = ""
        self.anisotropy_ratioZ = ""
        self.calculate_diffusion = ""
        self.cementation_exponent = ""
        self.constant_tortuosity = ""
        self.diffusion_activation = ""
        self.dispersion = ""
        self.dispersivity = ""
        self.distance_units = ""
        self.D_25 = []
        self.D_MP = ""
        self.fix_diffusion = ""
        self.formation_factor = ""
        self.gas_diffusion = ""
        self.set_porosity = ""
        self.tortuosity = ""
        self.read_TortuosityFile = ""
        self.threshold_porosity = ""
        self.time_units = ""
        self.tortuosity_above = ""
        self.tortuosity_below = ""
        self.tortuosityMP = []

    def __str__(self):
        """Return the transport block as a string"""
        result = []
        for attr, value in self.__dict__.items():
            exceptions = ["other", "D_25", "tortuosityMP"]
            if attr not in exceptions and (value or isinstance(value, bool)):
                if isinstance(value, list):
                    result.append(f"{attr:<22} {' '.join(map(str, value))}")
                elif isinstance(value, bool):
                    result.append(f"{attr:<22} {str(value).lower()}")
                else:
                    result.append(f"{attr:<22} {value}")
        for d in self.D_25:
            result.append(f"D_25                   {d}")
        for t in self.tortuosityMP:
            result.append(f"tortuosityMP           {t}")
        for key, value in self.other.items():
            result.append(f"{key:<22} {value}")
        return "\n".join(result)


class Flow(KeywordBlock):
    """Flow block for the input file"""

    def __init__(self):
        super().__init__()
        self.calculate_flow = ""
        self.constant_flow = ""
        self.constant_gasflow = ""
        self.distance_units = ""
        self.gaspump = ""
        self.gravity = ""
        self.infiltration = ""
        self.initialize_hydrostatic = ""
        self.permeability_x = []
        self.permeability_y = []
        self.permeability_z = []
        self.porosity_update = ""
        self.pressure = []
        self.pump = ""
        self.read_GasVelocityFile = ""
        self.read_PermeabilityFile = ""
        self.read_VelocityFile = ""
        self.space_units = ""
        self.time_units = ""

    def set_parameters(self, parameters):
        """Set the parameters for the Flow block"""
        for key, value in parameters.items():
            key = key.lower()
            if key in ["pressure", "permeability_x", "permeability_y", "permeability_z"]:
                self.pressure.append(value)
            elif hasattr(self, key):
                setattr(self, key, value)
            else:
                self.other[key] = value

    def __str__(self):
        """Return the flow block as a string"""
        exceptions = ["other", "pressure", "permeability_x", "permeability_y", "permeability_z"]
        result = []
        for attr, value in self.__dict__.items():
            if attr not in exceptions and (value or isinstance(value, bool)):
                if isinstance(value, list):
                    result.append(f"{attr:<22} {' '.join(map(str, value))}")
                elif isinstance(value, bool):
                    result.append(f"{attr:<22} {str(value).lower()}")
                else:
                    result.append(f"{attr:<22} {value}")

        # Print pressure, since it can be defined multiple times
        for press in self.pressure:
            line = "pressure               "
            for p in press.split():
                line += f"{p:<8} "
            result.append(line.strip())

        # Print permeability, since it can be defined multiple times
        for perm in self.permeability_x:
            result.append(f"permeability_x         {perm}")
        for perm in self.permeability_y:
            result.append(f"permeability_y         {perm}")
        for perm in self.permeability_z:
            result.append(f"permeability_z         {perm}")

        # Finally, the other category
        for key, value in self.other.items():
            result.append(f"{key:<22} {value}")
        return "\n".join(result)


class Temperature(KeywordBlock):
    """Temperature block for the input file"""

    def __init__(self):
        super().__init__()
        self.read_TemperatureFile = ""
        self.set_temperature = ""
        self.temperature_gradient = ""


class Porosity(KeywordBlock):
    """Porosity block for the input file"""

    def __init__(self):
        super().__init__()
        self.fix_microporosity = ""
        self.fix_porosity = ""
        self.mineral_update = ""
        self.minimum_porosity = ""
        self.MultiplyPorosityTortuosity = ""
        self.porosity_exponent = ""
        self.porosity_threshold = ""
        self.porosity_update = ""
        self.read_PorosityFile = ""
        self.set_porosity = ""
        self.UpdateDDL = ""
        self.update_porosity = ""


class Pest(KeywordBlock):
    """Pest block for the input file"""

    def __init__(self):
        super().__init__()
        self.CreatePestInstructionFile = ""
        self.CreatePestExchangeFile = ""
        self.exchange = ""


class Erosion(KeywordBlock):
    """Erosion block for the input file"""

    def __init__(self):
        super().__init__()
        self.read_BurialFile = ""


class Isotopes(KeywordBlock):
    """Isotopes block for the input file"""

    def __init__(self):
        super().__init__()
        self.isotope_time_series = ""
        self.mineral = []
        self.primary = []

    def set_parameters(self, parameters):
        """Set the parameters for the Isotopes block"""
        for key, value in parameters.items():
            key = key.lower()
            if key == "mineral":
                self.mineral.append(value)
            elif key == "primary":
                self.primary.append(value)
            elif hasattr(self, key):
                setattr(self, key, value)
            else:
                self.other[key] = value

    def __str__(self):
        """Return the isotopes block as a string"""
        result = []
        for attr, value in self.__dict__.items():
            if attr not in ["other", "mineral", "primary"] and (value or isinstance(value, bool)):
                if isinstance(value, list):
                    result.append(f"{attr:<22} {' '.join(map(str, value))}")
                elif isinstance(value, bool):
                    result.append(f"{attr:<22} {str(value).lower()}")
                else:
                    result.append(f"{attr:<22} {value}")
        for p in self.primary:
            result.append(f"primary                {p}")
        for m in self.mineral:
            result.append(f"mineral                {m}")
        for key, value in self.other.items():
            result.append(f"{key:<22} {value}")
        return "\n".join(result)


class InitialConditions(KeywordBlock):
    """Initial conditions block for the input file"""

    def __init__(self):
        super().__init__()
        self.conditions = []

    def set_parameters(self, parameters):
        """Set the parameters for the InitialConditions block"""
        self.conditions = parameters.get("conditions", [])

    def __str__(self):
        """Return the initial conditions block as a string"""
        result = []
        for cond in self.conditions:
            name, *coords = cond.split()
            coords_fmt = " ".join(c.ljust(8) for c in coords)
            line = f"{name:<22} {coords_fmt}"
            result.append(line.strip())
        return "\n".join(result)


class BoundaryConditions(KeywordBlock):
    """Boundary conditions block for the input file"""

    def __init__(self):
        super().__init__()
        self.x_begin = ""
        self.x_end = ""
        self.y_begin = ""
        self.y_end = ""
        self.z_begin = ""
        self.z_end = ""


class SpeciesBlock(KeywordBlock):
    """Base class for all species blocks in a CrunchFlow input file"""

    def __init__(self):
        super().__init__()
        self.species = []

    def set_parameters(self, parameters):
        """Set the parameters for the species block"""
        self.species = parameters.get("species", [])

    def __str__(self):
        """Return the species block as a string"""
        return "\n".join(self.species)


class PrimarySpecies(SpeciesBlock):
    """Primary species block for the input file"""

    pass


class SecondarySpecies(SpeciesBlock):
    """Secondary species block for the input file"""

    pass


class Gases(SpeciesBlock):
    """Gases block for the input file"""

    pass


class SurfaceComplexation(KeywordBlock):
    """Surface complexation block for the input file"""

    def __init__(self):
        super().__init__()
        self.species = []
        self.species_dict = {}

    def set_parameters(self, parameters):
        """Set the parameters for the SurfaceComplexation block"""
        for species, details in parameters.items():
            self.species.append(species)
            self.species_dict[species] = details

    def __str__(self):
        """Return the surface complexation block as a string"""
        result = []
        for species, details in self.species_dict.items():
            result.append(f"{species:<22} {details}")
        return "\n".join(result)


class KineticsBlock(KeywordBlock):
    """Base class for all kinetics blocks in a CrunchFlow input file"""

    def __init__(self):
        super().__init__()
        self.species_dict = {}
        self.species = []

    def set_parameters(self, parameters):
        """Set the parameters for the kinetics block"""
        for species, details in parameters.items():
            label = details.pop("label", None)
            if species not in self.species_dict:
                self.species_dict[species] = {}
                self.species.append(species)
            if label:
                self.species_dict[species][label] = details
            else:
                self.species_dict[species]["default"] = details

    def update_parameters(self, species, label, new_details):
        """Update the parameters for a species in the kinetics block"""
        if species in self.species_dict and label in self.species_dict[species]:
            self.species_dict[species][label].update(new_details)
        else:
            print(f"Species '{species}' with label '{label}' not found.")

    def __str__(self):
        """Return the kinetics block as a string"""
        result = []
        for species, labels in self.species_dict.items():
            if "default" in labels and not labels["default"]:
                result.append(f"{species:<20}")
            else:
                for label, details in labels.items():
                    if label == "default" and not details:
                        result.append(f"{species:<20}")
                    else:
                        details_str = " ".join(f"-{k} {v}" for k, v in details.items())
                        label_str = f"-label {label} " if label != "default" else ""
                        result.append(f"{species:<20} {label_str}{details_str}")
        return "\n".join(result)


class Minerals(KineticsBlock):
    """Minerals block for the input file"""

    pass


class AqueousKinetics(KineticsBlock):
    """Aqueous kinetics block for the input file"""

    pass
