"""Various utilities for working with CrunchFlow files."""

import os
import re


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
