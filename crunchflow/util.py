"""Various utilities for working with CrunchFlow files."""

import os
import re

# Prefixes of the spatial_profile files CrunchFlow writes. Each is followed by
# a time-step index and an output suffix, e.g. 'volume740.tec'.
SPATIAL_PROFILE_PREFIXES = (
    "Aq_conc",
    "Aq_totconc",
    "AqRate",
    "AqSat",
    "area",
    "conc",
    "DarcyFlux",
    "DelGBiomass",
    "divergence",
    "exchange",
    "fTBiomass",
    "gasdiffflux",
    "gases",
    "gases_conc",
    "liquid_saturation",
    "Mineralarea",
    "Mineralconc",
    "MineralPercent",
    "Mineralrate",
    "Mineralsaturation",
    "Mineralvolumefraction",
    "MineralVolfraction",
    "permeability",
    "pH",
    "porosity",
    "porositychange",
    "pressure",
    "pressure_head",
    "rate",
    "saturation",
    "speciation",
    "surface",
    "temperature",
    "tortuosity",
    "totcon",
    "totexchange",
    "TotMineral",
    "totsurface",
    "velocity",
    "velocityx",
    "velocityy",
    "velocityz",
    "volume",
    "water_content",
    "water_potential",
    "WeightPercent",
)

# Suffixes used for spatial_profile output by different versions of CrunchFlow
SPATIAL_PROFILE_SUFFIXES = (".out", ".tec", ".dat")

# Output files CrunchFlow writes under a fixed name, without a time-step index
FIXED_OUTPUT_FILES = (
    "CrunchJunk2.out",
    "fort.123",
    "initial_condition_Richards.out",
    "initial_condition_Richards.tec",
)


def correct_exponent(filename, folder=".", verbose="med"):
    """Correct triple digit exponents within a file. CrunchFlow
    has trouble outputting triple-digit exponents and omits
    the 'E'. For example, '2.5582E-180' prints as '2.5582-180'.

    Parameters
    ----------
    filename : str
        name of the file to be processed
    folder : str
        folder containing the file, either relative or absolute path
    verbose : {'med', 'high', 'low'}
        Print each correction as it's performed ('high'), print total
        number of corrections ('med'), or print nothing ('low'). The
        default is 'med'

    Returns
    -------
        None. Modifies the file in place.
    """
    crunch_file = os.path.join(folder, filename)
    n_repl = 0  # Count number of replacements

    # Pre-compile regex to save time
    neg_search = re.compile(r"([0-9][0-9])\-([0-9][0-9][0-9])")
    pos_search = re.compile(r"([0-9][0-9])\+([0-9][0-9][0-9])")

    # Read in crunch input file
    with open(crunch_file, "r") as fin:
        cf = fin.readlines()

    with open(crunch_file, "w") as fout:
        for line in cf:
            # before any changes, store line for printing
            tmp = line

            if re.search(neg_search, line):
                line = re.sub(neg_search, r"\1E-\2", line)
                n_repl += 1

            if re.search(pos_search, line):
                line = re.sub(pos_search, r"\1E+\2", line)
                n_repl += 1

            if verbose == "high":
                print("Original: \n\t " + tmp)
                print("New: \n\t " + line)

            fout.write(line)

    if verbose == "med":
        print("Made {} replacements in {}".format(n_repl, crunch_file))


def clear_output(folder=".", suffixes=SPATIAL_PROFILE_SUFFIXES, dry_run=False, verbose=True):
    """Delete the output files that CrunchFlow writes into a run folder, so
    that a simulation can be re-run from a clean directory.

    Removes the spatial_profile files that CrunchFlow names by appending a
    time-step index to a fixed prefix (for example ``volume740.tec`` or
    ``pH12.out``), along with the few output files it writes under a fixed
    name. If a ``PestControl.ant`` file is present, the main output log named
    within it is removed as well.

    Files are matched by name only: a file is deleted just when its name is a
    known CrunchFlow output prefix followed by one or more digits and one of
    `suffixes`. Only `folder` itself is searched, so subdirectories are left
    untouched.

    Note that ``time_series`` output files are not removed, because CrunchFlow
    takes their names from the input file rather than generating them.

    Parameters
    ----------
    folder : str
        folder to clear, either a relative or absolute path. The default is
        the current working directory
    suffixes : sequence of str
        file suffixes to treat as spatial_profile output. The default is
        ('.out', '.tec', '.dat')
    dry_run : bool
        if True, report the files that would be deleted without deleting
        them. The default is False
    verbose : bool
        print the number of files deleted. The default is True

    Returns
    -------
    list of str
        paths of the deleted files, sorted by name. If `dry_run` is True,
        the paths that would have been deleted.
    """
    if not os.path.isdir(folder):
        raise NotADirectoryError("No such folder: {}".format(folder))

    # Accept a single suffix passed as a bare string
    if isinstance(suffixes, str):
        suffixes = (suffixes,)

    targets = set(FIXED_OUTPUT_FILES)

    if suffixes:
        pattern = re.compile(
            r"^({})[0-9]+({})$".format(
                "|".join(re.escape(prefix) for prefix in SPATIAL_PROFILE_PREFIXES),
                "|".join(re.escape(suffix) for suffix in suffixes),
            )
        )
        targets.update(name for name in os.listdir(folder) if pattern.match(name))

    # PestControl.ant holds the name of the input file, next to which
    # CrunchFlow writes the main output log as run_name.out
    pest_control = os.path.join(folder, "PestControl.ant")
    if os.path.isfile(pest_control):
        with open(pest_control, "r") as fin:
            infile = fin.read().strip()
        if infile:
            targets.add(os.path.splitext(infile)[0] + ".out")

    deleted = []
    for name in sorted(targets):
        path = os.path.join(folder, name)
        if os.path.isfile(path):
            if not dry_run:
                os.remove(path)
            deleted.append(path)

    if verbose:
        verb = "Would delete" if dry_run else "Deleted"
        print("{} {} output file(s) in {}".format(verb, len(deleted), folder))

    return deleted
