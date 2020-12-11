python-crunchflow -- A Python Library for Processing CrunchFlow Output
======================================================================

python-crunchflow is an open source python package designed for processing and visualizing output from the reactive transport model `CrunchFlow <https://bitbucket.org/crunchflow/>`_. This library provides a set of tools for quickly plotting time series and geochemical output produced by CrunchFlow. 

This library consists of three main modules: 

- ``output.tec``. Class for working with CrunchFlow ``spatial_profile`` output, which include spatially-variable geochemical conditions and sediment properties at a given point in time.
- ``output.timeseries``. Class for working with CrunchFlow ``time_series`` output, which consists of temporally-variable aqueous concentrations at a specific node.
- ``util``. General utilities for parsing CrunchFlow input files and manipulating output.

Contents:

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Spatial (tec) Output
====================
.. automodule:: output.tec
    :members:

Time Series Output
==================
.. automodule:: output.timeseries
    :members:

Utilities
=========
.. automodule:: util
    :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
