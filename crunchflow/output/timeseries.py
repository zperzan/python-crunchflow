"""A module for loading and plotting CrunchFlow time series output files."""

import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def get_ts_coords(tsfile):
    """Given a CrunchFlow time series ouput file, return the coordinate
    at which that time series was output.

    Parameters
    ----------
    tsfile : str
        filename containing timeseries output

    Returns
    -------
    coords : tuple of int
        Coordinates of the form (x, y, z)
    """
    # Open the file and read in the first line
    with open(tsfile) as f:
        firstline = f.readline()
        if "Flux weighted" in firstline:
            # First split on colon
            coord_str = firstline.split(":")[-1]

            # Then split on both space and hyphen. Sometimes CrunchFlow outputs
            # these coordinates as "1-100" and sometimes as "1- 100"
            raw_coord_list = re.split(r"[-\s]", coord_str)

            # Remove empty strings
            coord_list = [x for x in raw_coord_list if x]

            if len(coord_list) == 6:
                coords = (
                    "%s-%s" % (coord_list[0], coord_list[1]),
                    "%s-%s" % (coord_list[2], coord_list[3]),
                    "%s-%s" % (coord_list[4], coord_list[5]),
                )
            else:
                coords = coord_str

        else:
            # Final 3 fields are x, y, and z
            fields = firstline.split()[-3:]
            x = int(fields[0])
            y = int(fields[1])
            z = int(fields[2])

            # Return tuple
            coords = (x, y, z)

    return coords


def get_ts_duplicates(tsfile):
    """Find duplicate columns in a time series file. Useful when a user
    specifies `time_series_print all` in CrunchFlow; the time series file
    includes primary species printed twice.

    Parameters
    ----------
    tsfile : str
        path to the time series file

    Returns
    -------
    columns : list
        list of column headings without duplicates
    nondup_indices : list
        list of column indices without duplicates
    """
    # Open the file and read in the second line
    with open(tsfile) as f:
        for i, line in enumerate(f):
            if i == 1:
                columns = line.split()

    # Get list of duplicates and how many times each is seen
    seen = {}  # Count occurrences of each item
    duplicates = []  # List of duplicate items
    dup_indices = []  # List of first occurrences of each duplicate item

    for col in columns:
        if col not in seen:
            seen[col] = 1
        else:
            if seen[col] == 1:
                duplicates.append(col)
                idx = columns.index(col)
                dup_indices.append(idx)
            seen[col] += 1  # Do not go through if seen[col] == 1 again

    # Total number of columns in the file, including duplicates
    totcols = len(columns)
    nondup_indices = list(range(totcols))  # List of non-duplicated indices

    # Delete each duplicate, but reverse sort to avoid throwing off indices
    # after deleting earlier elements
    for idx in sorted(dup_indices, reverse=True):
        del columns[idx]
        del nondup_indices[idx]

    return columns, nondup_indices


class timeseries:
    """The deprecated timeseries class for working
    with CrunchFlow time series output files.
    """

    def __init__(self, tsfile, folder="."):
        raise DeprecationWarning(
            "The crunchflow.output.timeseries class has been deprecated. Use crunchflow.output.TimeSeries instead."
        )


class TimeSeries:
    """The timeseries class for working with CrunchFlow time
    series output files.

    Attributes
    ----------
    coords : tuple of int
        x, y and z coordinates of the time series
    timeunit : str
        time unit used in the CrunchFlow input
    unit : str
        Concentration units included in the file. Automatically set to the
        default CrunchFlow concentration units (mol/kgw)
    species : list of str
        list of aqueous species in the file
    data : ndarray of float
        Numpy array of all data. First col is the time step and remaining
        cols are species in the same order as self.species list
    df : dataframe of float
        Pandas dataframe of all data. Index is the time step and columns
        #are the aqueous species

    Methods
    -------
    convert_mgL(database='datacom.dbs', folder='.')
        Convert time series concentrations from mol/kgw to mg/L (ppm).
    plot(species, units='mg/L', **kwargs)
        Plot the time series of one or more species.

    Examples
    --------
    >>> ts = TimeSeries('Well1-1.txt')
    >>> ts.convert_mgL()
    >>> calcium = ts.df['Ca++']
    >>> ts.plot('Ca++')
    """

    def __init__(self, tsfile, folder="."):
        """Read in and get basic info about the timeseries file `tsfile`.

        Parameters
        ----------
        tsfile : str
            Name of the CrunchFlow time series file
        folder : str
            Path to the CrunchFlow time series file
        """
        tsfilepath = os.path.join(folder, tsfile)

        # Get coordinates at which time series was output
        self.coords = get_ts_coords(tsfilepath)

        # Get list of duplicates and their indices
        self.columns, indices = get_ts_duplicates(tsfilepath)

        # Assume that, if there are duplicates, it's because user
        # specified `time_series_print all` in CrunchFlow, in which
        # case non-duplicate columns are printed in log format
        ncols = len(self.columns)
        lastcol = max(indices) + 1

        # If idx of the last column is greater than the # cols,
        # duplicates were deleted, so set logformat = True
        if lastcol > ncols:
            logformat = True
        else:
            logformat = False

        # List of species (columns without the time column)
        self.species = self.columns[1:]

        # Set time and concentration units
        t = self.columns[0]
        self.timeunit = t[t.find("(") + 1 : t.find(")")]
        self.unit = "mol/L"

        # Load data into numpy array
        self.data = np.genfromtxt(
            tsfilepath, skip_header=2, usecols=indices, missing_values=["Infinity", "NaN", "-Infinity"]
        )

        # If necessary, convert from log to real format
        if logformat:
            if self.columns[1] != "pH":
                self.data[:, 1:] = 10 ** self.data[:, 1:]
            else:
                # Skip first two cols, which are time and pH
                self.data[:, 2:] = 10 ** self.data[:, 2:]

        # Load into a pandas dataframe as well
        self.df = pd.DataFrame(data=self.data[:, 1:], index=self.data[:, 0], columns=self.species)
        self.df.index.name = "time"

    def convert_mgL(self, database="datacom.dbs", folder=".", warnings=True):
        """Convert time series concentrations from mol/kgw to mg/L (ppm).
        Note that this assumes that 1 kg water = 1 L water.

        Parameters
        ----------
        database : str
            name of the CrunchFlow database. The default is 'datacom.dbs'
        folder : str
            path to the database. The default is current directory.
        warnings : bool
            whether to print warnings for species not found in the database.

        Returns
        -------
            None. Modifies timeseries object in place.
        """
        databasepath = os.path.join(folder, database)

        # If units are already mg/L, no need to do anything
        if self.unit == "mg/L":
            return

        # Check if database exists
        if not os.path.exists(databasepath):
            raise OSError("Could not find " + databasepath)

        molar_mass = {}

        # Open the database and get the molar mass of each species
        with open(databasepath) as db:
            for line in db:
                # Skip blank lines
                if not line.strip():
                    continue

                for spec in self.species:
                    # Database format is, e.g., "'Ca++' 6.0  2.0    40.0780",
                    # where the last value is the molar mass
                    if line.split()[0] == "'{}'".format(spec):
                        molar_mass[spec] = float(line.split()[-1])

        # Delete keys with molar masses of 0 (e.g., tracers)
        # and do not convert them to mg/L
        del_keys = []
        for key, value in molar_mass.items():
            if value == 0:
                del_keys.append(key)  # Cannot delete key within loop, otherwise
                # it changes size on each iteration
        for key in del_keys:
            del molar_mass[key]

        for spec in self.species:
            if spec not in molar_mass.keys():
                if warnings:
                    print("Warning -- Did not convert {} to mg/L".format(spec))
            else:
                idx = self.columns.index(spec)
                # Only need to convert .data since .data and .df are linked
                self.data[:, idx] = self.data[:, idx] * molar_mass[spec] * 1000

        # Update the unit attribute
        self.unit = "mg/L"

    def plot(self, species, units="mg/L", **kwargs):
        """Plot the time series of one or more species.

        Parameters
        ----------
        species : str or list of str
            Either single species or list of species to be plotted
        units : str
            Concentration units to use for plotting. The default is 'mg/L'
        **kwargs : dict
            keyword arguments passed to plt.subplots (e.g., figsize)

        Returns
        -------
        fig : pyplot object
            figure handle for current plot
        ax : pyplot object
            axis handle for current plot
        """
        if units == "mg/L" and self.unit != "mg/L":
            # Raise error if cannot find datacom.dbs
            if not os.path.exists("./datacom.dbs"):
                raise OSError(
                    "Could not find default database. \
                Plot with other units or convert to mg/L first using the \
                convert_mgL method. See convert_mgL.__doc__ for more info."
                )

            self.convert_mgL()

        # Accept both str and list input, so if str, convert to list
        if isinstance(species, str):
            species = [species]

        fig, ax = plt.subplots(**kwargs)

        for spec in species:
            ax.plot(self.df.index, self.df[spec], label=spec)

        ax.set(xlabel="Time ({})".format(self.timeunit), ylabel="Concentration ({})".format(units))
        ax.legend()

        return fig, ax


if __name__ == "__main__":
    print(TimeSeries.__doc__)
