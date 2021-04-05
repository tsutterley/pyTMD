==============
pyTMD Overview
==============

``pyTMD`` is a Python-based tidal prediction software that reads OTIS, GOT and FES
formatted tidal solutions for predicting ocean and load tides and can use IERS
conventions for calculating radial pole tide displacements.
It was developed with the goal of supporting science applications for
airborne and satellite altimetry.
``pyTMD`` provides data access utilities for ascii, netCDF4, HDF5, and geotiff
formats.
``pyTMD`` also provides some very high-level plotting programs through the
use of `Jupyter Notebooks <./Examples.html>`_.

The program provides the interface between spatial and temporal coordinates and
the output ocean, load and pole tide variables as shown in the flowcharts below:

.. graphviz::
    :caption: Ocean and Load Tide Framework
    :align: center

    digraph {
        S [label="Spatial Coordinates" shape=box style="filled" color="darkorchid"]
        E [label="Temporal Values" shape=box style="filled" color="darkorchid"]
        M [label="Tide Model" shape=box style="filled" color="darkorchid"]
        T [label="pyTMD" shape=box style="filled" color="gray"]
        P [label="Tide Predictions" shape=box style="filled" color="mediumseagreen"]
        H [label="Tide Heights" shape=box style="filled" color="mediumseagreen"]
        C [label="Average Tidal Currents" shape=box style="filled" color="mediumseagreen"]
        S -> T
        E -> T
        M -> T
        T -> P
        T -> H
        T -> C
    }

.. graphviz::
    :caption: Pole Tide Framework
    :align: center

    digraph {
        S [label="Spatial Coordinates" shape=box style="filled" color="darkorchid"]
        E [label="Temporal Values" shape=box style="filled" color="darkorchid"]
        T [label="pyTMD" shape=box style="filled" color="gray"]
        P [label="Pole Tide Displacements" shape=box style="filled" color="mediumseagreen"]
        S -> T
        E -> T
        T -> P
    }
