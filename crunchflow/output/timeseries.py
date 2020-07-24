import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_ts_coords(tsfile):
    """Given a CrunchFlow time series ouput file, return the coordinate 
    at which that time series was output.
    
    params:
        tsfile [str]: filename containing timeseries output
    
    returns:
        coords [tuple(int)]: Coordinate of the form (x, y, z)"""
    
    # Open the file and read in the first line
    with open(tsfile) as f:
        for i, line in enumerate(f):
            if i == 0:
                fields = line.split()[-3:]
                
                # Final 3 fields are x, y, and z
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
    
    params:
        tsfile [str]: path to the time series file
    
    returns:
       columns [list]: list of column headings without duplicates
       nondup_indices [list]: list of column indices without duplicates
    """
    
    # Open the file and read in the second line
    with open(tsfile) as f:
        for i, line in enumerate(f):
            if i == 1:
                columns = line.split()
               
    # Get list of duplicates and how many times each is seen
    seen = {} # Count occurrences of each item
    duplicates = [] # List of duplicate items
    dup_indices = [] # List of first occurrences of each duplicate item
    
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
    nondup_indices = list(range(totcols)) # List of non-duplicated indices

    # Delete each duplicate, but reverse sort to avoid throwing off indices 
    # after deleting earlier elements
    for idx in sorted(dup_indices, reverse=True):
        del columns[idx]
        del nondup_indices[idx]

    return columns, nondup_indices

    
class timeseries:
    """This is the timeseries class for working with CrunchFlow time 
    series output files.

    Available methods include:
        __init__(tsfile, folder='.')
        convert_mgL(database='datacom.dbs', folder='.')
        plot(species, units='mg/L', **kwargs)

    Example usage:
        >>> ts = timeseries('Well1-1.txt')
        >>> ts.convert_mgL()
        >>> calcium = ts.df['Ca++']
        >>> ts.plot('Ca++')
    """
    
    def __init__(self, tsfile, folder='.'):
        """Read in and get basic info about the timeseries file `tsfile`.
        
        The __init__ method sets the following attributes:
            coords [tuple(int)]: x, y and z coordinates of the time series
            timeunit [str]: time unit used in the CrunchFlow input
            unit [str]: default CrunchFlow concentration unit is mol/kgw
            species [list(str)]: list of aqueous species in the file
            data [array(float)]: Numpy array of all data. First col is the 
                time step and remaining cols are species in the same order 
                as self.species list
            df [dataframe(float)]: Pandas dataframe of all data. Index is 
                the time step and columns are the aqueous species
                                   
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
        self.timeunit = t[t.find("(")+1:t.find(")")]
        self.unit = 'mol/L'
        
        # Load data into numpy array
        self.data = np.genfromtxt(tsfilepath, skip_header=2, usecols=indices, 
                                  missing_values=['Infinity', 'NaN', '-Infinity'])
        
        # If necessary, convert from log to real format
        if logformat:
            if self.columns[1] != 'pH':
                self.data[:, 1:] = 10**self.data[:, 1:]
            else:
                # Skip first two cols, which are time and pH
                self.data[:, 2:] = 10**self.data[:, 2:]
        
        # Load into a pandas dataframe as well
        self.df = pd.DataFrame(data=self.data[:, 1:], 
                               index=self.data[:, 0], 
                               columns=self.species)
        self.df.index.name = 'time'
    
    
    def convert_mgL(self, database='datacom.dbs', folder='.'):
        """Convert time series concentrations from mol/kgw to mg/L (ppm). 
        Note that this assumes that 1 kg water = 1 L water.
        
        params:
            database [str]: name of the CrunchFlow database. Default: 'datacom.dbs'
            folder [str]: path to the database. Default: current directory.
        """

        databasepath = os.path.join(folder, database)

        # If units are already mg/L, no need to do anything
        if self.unit == 'mg/L':
            return
        
        # Check if database exists
        if not os.path.exists(databasepath):
            raise OSError('Could not find ' + databasepath)
        
        molar_mass = {}

        # Open the database and get the molar mass of each species
        with open(databasepath) as db:
            for line in db:
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
                del_keys.append(key) # Cannot delete key within loop, otherwise 
                                     # it changes size on each iteration
        for key in del_keys:
            del molar_mass[key]

        for spec in self.species:
            if spec not in molar_mass.keys():
                print('Warning -- Did not convert {} to mg/L'.format(spec))
            else:
                idx = self.columns.index(spec)
                # Only need to convert .data since .data and .df are linked
                self.data[:, idx] = self.data[:, idx]*molar_mass[spec]*1000

        # Update the unit attribute
        self.unit = 'mg/L'
        
        return
        
        
    def plot(self, species, units='mg/L', **kwargs):
        """Plot the time series of species.
        
        params:
            species [str|list(str)]: either str or list of species to be plotted
            units [str]: units to use for plotting. Default: 'mg/L'
            **kwargs: keyword arguments passed to plt.subplots (e.g., figsize)

        returns:
            fig [matplotlib figure]: figure handle for current plot
            ax [matplotlib axis]: axis handle for current plot
        """
        
        if units == 'mg/L' and self.unit != 'mg/L':
            # Raise error if cannot find datacom.dbs
            if not os.path.exists('./datacom.dbs'):
                raise OSError('Could not find default database. \
                Plot with other units or convert to mg/L first using the \
                convert_mgL method. See convert_mgL.__doc__ for more info.')

            self.convert_mgL()
        
        # Accept both str and list input, so if str, convert to list
        if isinstance(species, str):
            species = [species]
        
        fig, ax = plt.subplots(**kwargs)
        
        for spec in species:
            ax.plot(self.df.index, self.df[spec], label=spec)
        
        ax.set(xlabel='Time ({})'.format(self.timeunit), 
               ylabel='Concentration ({})'.format(units))
        ax.legend();
        
        return fig, ax


if __name__ == '__main__':
    print(timeseries.__doc__)
