python-crunchflow -- A Python Library for Working with the CrunchFlow Reactive Transport Code
=============================================================================================

``python-crunchflow`` is an open source python package designed for working with the the reactive
transport model `CrunchFlow <https://github.com/CISteefel/CrunchTope>`_. This library provides a set
of tools for managing and manipulating CrunchFlow input files, as well as quickly reading and
plotting model output.

If you use the ``crunchflow`` package in published work, please cite the paper for which it was originally developed:

Perzan, Z., Babey, T., Caers, J., Bargar, J.R. and Maher, K., 2021, Local and global sensitivity analysis of a reactive transport model simulating floodplain redox cycling, *Water Resources Research*, doi: `10.1029/2021WR029723 <https://dx.doi.org/10.1029/2021WR029723>`_

This library consists of three main subpackages:

- ``crunchflow.output`` contains the ``TimeSeries`` and ``SpatialProfile`` classes for working with CrunchFlow ``time_series`` and ``spatial_profile`` output, respectively. These include concentration histories at a given grid cell as well as spatially-variable geochemical conditions and sediment properties at a given snapshot in time.
- ``crunchflow.input`` contains the ``InputFile`` and ``KeywordBlock`` classes for working with CrunchFlow input files.
- ``crunchflow.cli`` contains command-line tools for managing CrunchFlow runs, invoked as ``python -m crunchflow.cli``.

Contents
========

.. toctree::
   :maxdepth: 3
   :hidden:

   Home page <self>
   Tutorials <tutorials>

- **Tutorials**. Most of the functionality of this package is documented in the `Tutorials <tutorials.html>`_ section. There is one tutorial each for the ``crunchflow.input`` and ``crunchflow.output`` subpackages.
- **Documentation**. Some documentation of individual classes and functions is also included in the `API Reference <autoapi/index.html>`_.

Command-line tools
==================

The ``crunchflow.cli`` subpackage is run as ``python -m crunchflow.cli``. Its available
subcommands are listed by ``python -m crunchflow.cli --help``.

``clear-output`` deletes the output files that CrunchFlow writes into a run folder, so that a
simulation can be re-run from a clean directory:

.. code-block:: bash

   # Clear the current directory
   python -m crunchflow.cli clear-output

   # Check what would be removed from another folder first
   python -m crunchflow.cli clear-output --folder runs/case1 --dry-run

Files are matched by name, so a file is deleted just when its name is a known CrunchFlow output
prefix followed by a time-step index and a recognized suffix (``.out``, ``.tec``, or ``.dat`` by
default). Subdirectories are left untouched. Note that ``time_series`` output files are not
removed, because CrunchFlow takes their names from the input file rather than generating them.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
