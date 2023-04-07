==============
pyTMD Overview
==============

``pyTMD`` is a Python-based tidal prediction software.
It was developed with the goal of supporting science applications for
airborne and satellite altimetry originally as part of a
NASA Postdoctoral Program (NPP) Fellowship.
``pyTMD`` provides data access utilities for ascii, netCDF4, HDF5, and geotiff
formats.
``pyTMD`` also provides some very high-level plotting programs through the
use of `Jupyter Notebooks <../user_guide/Examples.html>`_.

The program provides the interface between spatial and temporal coordinates and
the output ocean, load, solid Earth and pole tide variables as shown in the flowcharts below:

.. graphviz::
    :caption: Ocean and Load Tide Framework
    :align: center

    digraph {
        S [label="Spatial Coordinates" shape=box style="filled" color="#7570b3"]
        E [label="Temporal Values" shape=box style="filled" color="#7570b3"]
        M [label="Tide Model" shape=box style="filled" color="#7570b3"]
        T [label="pyTMD" shape=box style="filled" color="gray"]
        P [label="Tide Predictions" shape=box style="filled" color="#1b9e77"]
        H [label="Tide Heights" shape=box style="filled" color="#1b9e77"]
        C [label="Average Tidal Currents" shape=box style="filled" color="#1b9e77"]
        S -> T
        E -> T
        M -> T
        T -> P
        T -> H
        T -> C
    }

.. graphviz::
    :caption: Solid Earth and Pole Tide Framework
    :align: center

    digraph {
        S [label="Spatial Coordinates" shape=box style="filled" color="#7570b3"]
        E [label="Temporal Values" shape=box style="filled" color="#7570b3"]
        T [label="pyTMD" shape=box style="filled" color="gray"]
        D [label="Solid Earth Displacements" shape=box style="filled" color="#1b9e77"]
        P [label="Pole Tide Displacements" shape=box style="filled" color="#1b9e77"]
        S -> T
        E -> T
        T -> D
        T -> P
    }
