import os
import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def get_tec_metadata(file):
    """Given a crunch .tec output file, read it in and return a list of the 
    variables (e.g., 'X-Perm', 'Y-Perm', etc) included in the output file.
    
    Parameters
    ----------
    file : str
        Filename to read in
        
    Returns
    -------
    columns : list, str
        ordered list of variables included in filename 
    title : str
        Value included in this file, as output by CrunchFlow
    """

    # Instantiate empty list of columns
    columns = []
    
    # Open file and read line by line
    with open(file) as f:
        for line in f:
            if "TITLE" in line:
                title = line.split('"')[1]
            
            if "VARIABLES" in line:
                columns = line.split('"')
    
    # Remove x, y, z
    columns.remove('VARIABLES = ')
    columns.remove('X')
    columns.remove('Y')
    columns.remove('Z')
    
    # Remove whitespace
    remove_fields = [col for col in columns if col.isspace()]

    for col in remove_fields:
        columns.remove(col)
    
    return columns, title


def get_tec_output_times(crunch_log, folder='.'):
    """Given a log of CrunchFlow terminal output, get the time steps 
    associated with .tec file numbers. For example, during a model run, 
    CrunchFlow will print progress to the screen, including blocks that 
    look like:
    
        >>> WRITING OUTPUT FILES
        >>> Time (days) =  9.000E+00
        >>> File number  =            2
    
    This function combs through that terminal output to create a list of 
    times at which .tec files were written.
    
    Parameters
    ----------
    crunch_log : str
        Filename of the crunchflow terminal output
    folder : str
        Folder in which 'crunch_log' is located.The default is the
        current directory
        
    Returns
    -------
    output_times : list of float 
        Ordered list of the output times (in CrunchFlow time units) 
        corresponding to each .tec file number.

    """
    
    output_times = {}
    
    filepath = os.path.join(folder, crunch_log)
    with open(filepath, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if 'WRITING OUTPUT FILES' in line:
                time = float(lines[i+1].split()[3])
                fileno = int(lines[i+2].split()[3])
                
                output_times[fileno] = time
                
    return output_times


class tec:
    """This is the tec class for working with CrunchFlow .tec files.

    Attributes
    ----------
    times : list of int 
        relative times at which .tec files were output, sorted
    files : list of str
        ordered list of all files matching fileprefix
    columns : list of str
        list of all vars included within each .tec file
    nx : int
        # of grid cells in the x direction
    ny : int
        # of grid cells in the y direction
    nz : int
        # of grid cells in the z direction
    coords : ndarray of float 
        (nx*ny*nz, 2) array of x, y coordinates of each node
    griddedX : ndarray of float 
        (ny, nx) array of x coordinates; used for plotting
    griddedY : ndarray of float 
        (ny, nx) array of y coordinates; used for plotting
        
    Methods
    -------
    plot(variable, time)
        Plot .tec output from a single time step.
    plot_series(variable)
        Plot .tec output for a series of time steps.
    extract(variable, time)
        Return the spatial profile of var at the given time.
    outline(variable, value, time)
        Get line segments that outline all the regions equal to 
        the provided value.

    Examples
    --------
    >>> vol = tec('volume')
    >>> print(vol.columns)
    >>> vol.plot('Calcite')
    >>> calcite = vol.extract('Calcite', time=2)
    """
    
    def __init__(self, fileprefix, folder='.', output_times=None):
        """Read in and get basic info about all .tec files matching `fileprefix`.
        For example, `tec('volume')` will read in all files matching 
        'volume[0-9]+.tec'
        
        Parameters
        ----------
        fileprefix : str
            file prefix of the tec file to read in, without the timestep 
            number or ".tec" file ending.
        folder : str
            folder in which the .tec files are located. The default is the 
            current directory
        output_times : list of float 
            list of the actual output times at which the .tec files were 
            output, in CrunchFlow time units
        """
             
        # Glob doesn't support regex, so match manually using re and os
        search_str = fileprefix + '[0-9]+.tec'
        unsort_files = [f for f in os.listdir(folder) if re.search(search_str, f)]
        
        # Get the output time in sequence for each file
        st = len(fileprefix)
        self.times = [int(f[st:-4]) for f in unsort_files]
        self.times.sort()
        
        # Set output times if provided
        if output_times is None:
            self.output_times = None
        else:
            self.output_times = output_times

        # Create sored list of files using the times (which are sorted numerically)
        # This solves the problem that sorted str print as ['vol1', 'vol10', 'vol2', etc.]
        self.files = [fileprefix + str(t) + '.tec' for t in self.times]
        self.files = [os.path.join(folder, f) for f in self.files]

        # Get the columns (variables) available in each .tec file
        self.columns, self.title = get_tec_metadata(self.files[0])
        
        # Get the grid size
        with open(self.files[0]) as f:
            for i, line in enumerate(f):
                # Go to third line
                if i == 2:
                    # Split on comma
                    fields = line.split(',') 
                    
                    # Workaround for CrunchFlow typo in AqRate.tec files
                    if 'Aqueous Rate' not in self.title:
                        # Extract numeric digits from each field
                        self.nx = int(re.sub('\D', '', fields[0]))
                        self.ny = int(re.sub('\D', '', fields[1]))
                        self.nz = int(re.sub('\D', '', fields[2]))
                        
                    else:
                        self.nx = int(re.sub('\D', '', fields[1]))
                        self.ny = int(re.sub('\D', '', fields[2]))
                        self.nz = np.nan


        # Get gridded X and Y for plotting
        self.coords = np.loadtxt(self.files[0], skiprows=3, usecols=[0, 1])
        self.griddedX = self.coords[:, 0].reshape(self.ny, self.nx)
        self.griddedY = self.coords[:, 1].reshape(self.ny, self.nx)


    def plot(self, var, time=None, plot_type='image', figsize=(12,3), 
             **kwargs):
        """Plot .tec output from a single time step.
        
        Parameters
        ----------
        var : str
            variable to plot (e.g., 'H+')
        time : int
            which time slice to show
        figsize : tuple of int
            figure size; default (12, 3)
        plot_type : str
            whether to display an image ('image') or a color-filled 
            contour plot ('contourf')
        **kwargs : dict
            args passed on to matplotlib plot function (either contourf 
            or imshow)
        
        Returns
        -------
        fig : pyplot object
            matplotlib fig handle for the plot
        ax : pyplot object
            matplotlib ax handle for the plot
        """
        
        # If time is not given, plot the first time slice
        if time is None:
            time = self.times[0]
        
        if plot_type not in ['image', 'contour']:
            raise ValueError("plot_type must be either 'image' or 'contour'")
        
        # Get column number of var to plot
        col_no = self.columns.index(var) + 3 # Add 3 after deleting X, Y and Z
        itime = self.times.index(time) # Index of time, if self.times[0] != 0
        
        # Read in gridded data to numpy array
        z = np.loadtxt(self.files[itime], skiprows=3, usecols=col_no)
        z = z.reshape(self.ny, self.nx)
        z[np.isinf(z)] = np.nan # Mark inf values as nan
        
        fig, ax = plt.subplots(figsize=figsize)

        # Plot data as color-filled contour
        if plot_type == 'contour':
            # Define colorbar steps
            steps = np.linspace(np.nanmin(z), np.nanmax(z), 12)

            # Ensure that colorbar steps are strictly increasing
            if np.nanmin(z) == np.nanmax(z):
                if np.nanmin(z) == 0:
                    steps = np.arange(12)
                else:
                    factors = np.arange(12)/100 + 1
                    steps = steps * factors

            CS = ax.contourf(self.griddedX, self.griddedY, z, steps, **kwargs)
        
        # Plot data as image
        elif plot_type == 'image':
            extent = [self.griddedX[0,0], self.griddedX[-1,-1], 
                    self.griddedY[0,0], self.griddedY[-1,-1]]
            CS = ax.imshow(z, origin='lower', extent=extent, **kwargs)
        
        ax.set(aspect='equal')
        
        # Add colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.2)
        plt.colorbar(CS, cax=cax)
        
        # Build title
        title=self.title
        
        # Add the output time if it was provided
        if self.output_times is not None:
            itime = self.times.index(time)
            title += " @ t=" + str(self.output_times[itime])
        
        ax.set(title=title+' -- '+var)
        
        return fig, ax
    
    
    def plot_series(self, var, times=None, plot_type='image', 
                    figsize=None, **kwargs):
        """Plot .tec output for a series of time steps.
        
        Parameters
        ----------
        var : str
            variable to plot (e.g., 'H+')
        times : list
            time steps to plot, must be either components of self.times or 
            self.output_times; defaults to all time steps
        plot_type : str
            whether to display an image ('image') or a contour-filled contour 
            plot ('contourf')
        figsize : tuple of int
            figure size; default (12, 1.5 * # of timesteps)
        **kwargs : dict
            args passed on to matplotlib's imshow
            
        Returns
        -------
        fig : pyplot object
            matplotlib fig handle for the plot
        axes : list of pyplot objects
            list of matplotlib axes handles for each subplot
        """
        
        # Check input
        if plot_type not in ['image', 'contour']:
            raise ValueError("plot_type must be either 'image' or 'contour'")

        if times is None:
            times = self.times

        # Set up figsize if not provided
        if figsize is None:
            figsize = (12,1.8*len(times))

        # Interpret whether times are in output_times or self.times
        in_self_times=True
        for time in times:
            # Assume first that they are in self.times
            if time not in self.times:
                in_self_times=False
                if self.output_times is None:
                    raise ValueError("Could not find {} in self.times".format(time))
                else:
                    if time not in self.output_times:
                        raise ValueError("Could not find {} in self.output_times".format(time))
                    
        # If times were given as components of self.output_times, convert 
        # to self.times equivalent
        if not in_self_times:
            times = [self.times[self.output_times.index(t)] for t in times]
    
        # Get column number of var to plot
        col_no = self.columns.index(var) + 3 
              
        # Set up minimum/maximum values observed across all time slices
        # Used for min/max of color range
        minz = np.Inf
        maxz = np.NINF
        
        # z is a dict of numpy arrays
        z = {}
        for time in times:
            # Get index of time within self.times
            itime = self.times.index(time) 
            z[time] = np.loadtxt(self.files[itime], skiprows=3, usecols=col_no)
            z[time] = z[time].reshape(self.ny, self.nx)
            z[time][np.isinf(z[time])] = np.nan # Mark inf values as nan
            
            # Update min/max for color range
            if np.nanmin(z[time]) < minz:
                minz = np.nanmin(z[time])
            if np.nanmax(z[time]) > maxz:
                maxz = np.nanmax(z[time])

        # If it's a contour plot, setup the color intervals/steps
        if plot_type == 'contour':
            # Define colorbar steps
            steps = np.linspace(minz, maxz, 12)

            # Ensure that colors are strictly increasing for the colorbar
            if minz == maxz:
                if minz == 0:
                    steps = np.arange(12)
                else:
                    factors = np.arange(12)/100 + 1
                    steps = steps * factors

        fig, axes = plt.subplots(len(times), sharex=True, figsize=figsize)
        
        # Loop through each time step and add to plot
        for i, time in enumerate(times):
            # itime is the index of time within self.times
            # i is the index of time within times
            itime = self.times.index(time) 
            
            # In case plot_series is called when there's only one time step
            if len(times) > 1:
                ax = axes[i]
            else:
                ax = axes
            
            # Add a title to the first plot 
            # (and to subsequent plots if output_times was provided)
            if i == 0:
                title = self.title
                
                # Add the output time if it was provided when tec was init
                if self.output_times is not None:
                    itime = self.times.index(time)
                    title += " @ t=" + str(self.output_times[itime])
                ax.set(title=title+' -- '+var)
            else:
                if self.output_times is not None:
                    title = "t = " + str(self.output_times[itime])
                    ax.set(title=title)
            
            
            if plot_type == 'contour':
                CS = ax.contourf(self.griddedX, self.griddedY, z[time], 
                                 steps, **kwargs)

            elif plot_type == 'image':                
                extent = [self.griddedX[0,0], self.griddedX[-1,-1], 
                        self.griddedY[0,0], self.griddedY[-1,-1]]
                CS = ax.imshow(z[time], origin='lower', extent=extent, 
                            vmin=minz, vmax=maxz, **kwargs)

            ax.set(aspect='equal')
        
            # Add colorbar
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="2%", pad=0.2)
            plt.colorbar(CS, cax=cax)

        
        return fig, axes
        
        
    def extract(self, var, time=1):
        """Return the spatial profile of var at the given time.
        
        Parameters
        ----------
        var : str
            variable to return (e.g., 'Porosity')
        time : int
            time slice at which to return the profile
        
        Returns
        -------
        data : ndarray of float
            (ny, nx) numpy array of var
        
        """
        if var not in self.columns:
            raise ValueError('{} not found'.format(var))
        
        # Get column number of var to retrieve
        col_no = self.columns.index(var) + 3 # Add 3 because we deleted X, Y and Z cols        
        data = np.loadtxt(self.files[time - 1], skiprows=3, usecols=col_no)
        
        # Reshape to (ny, nx)
        data = data.reshape(self.ny, self.nx)
        
        return data


    def outline(self, var, value=None, time=1):
        """For a given .tec file, get line segments that outline all the regions 
        equal to the provided value. Useful for generating outlines of a areas 
        of a model domain sharing a single attribute (e.g., permeability) and 
        later adding the outlines to a plot of a different variable.
        
        Parameters
        ----------
        var : str
            variable to outline (e.g., 'Porosity')
        value : float 
            value to outline. The default is least-frequent value in array
        time : int
            time slice at which to create an outline
        
        Returns
        -------
        segments : ndarray of float
            array of (x,y) coordinate pairs for each line segment that comprises 
            the outline

        Examples
        --------
        >>> # Get stratigraphy from permeability map
        >>> perm = tec('permeability')
        >>> segments = perm.outline('X-Perm')
        >>> 
        >>> # Plot stratigraphy outlines on O2 map
        >>> conc = tec('conc')
        >>> fig, ax = conc.plot('O2(aq)')
        >>> ax.plot(segments[:,0], segments[:,1])
        
        """
        
        # Extract the data for this var and time slice
        data = self.extract(var, time=time)
        
        # Assume least-frequent value if it isn't provided explicitly
        if value is None:
            vals, counts = np.unique(data, return_counts=True)
            value = vals[np.argmin(counts)]
            
        # Check that value is in data
        if np.count_nonzero(data == value) == 0:
            raise ValueError("{} not found in array".format(value))

        # Mask pixels to outline
        mask = (data == value)

        # Get coordinates of where to draw horizontal and vertical segments
        # I.e., where adjacent pixels are not equal to each other
        ver_seg = np.where(mask[:,1:] != mask[:, :-1])
        hor_seg = np.where(mask[1:,:] != mask[:-1, :])

        # Create list of each line segment, 
        # separated by NaN [(start_coord), (end_coord), (nan nan), ... ]
        l = []
        for p in zip(*hor_seg):
            l.append((p[1], p[0]+1))
            l.append((p[1]+1, p[0]+1))
            l.append((np.nan,np.nan))

        # and the same for vertical segments
        for p in zip(*ver_seg):
            l.append((p[1]+1, p[0]))
            l.append((p[1]+1, p[0]+1))
            l.append((np.nan, np.nan))
        
        # Convert list to array
        segments = np.array(l) 

        # Before rescaling, get image extent
        x0 = self.griddedX[0, 0]
        x1 = self.griddedX[-1, -1] 
        y0 = self.griddedY[0, 0]
        y1 = self.griddedY[-1, -1] 

        # Rescale points according to the extent
        segments[:,0] = x0 + (x1-x0) * segments[:,0] / data.shape[1]
        segments[:,1] = y0 + (y1-y0) * segments[:,1] / data.shape[0]

        return segments


if __name__ == '__main__':
    print(tec.__doc__)
