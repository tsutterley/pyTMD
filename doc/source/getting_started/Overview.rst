==============
pyTMD Overview
==============

``pyTMD`` is a Python-based tidal prediction software.
It was developed with the goal of supporting science applications for
airborne and satellite altimetry originally as part of a
NASA Postdoctoral Program (NPP) Fellowship.
``pyTMD`` provides data access utilities for ascii, netCDF4, HDF5, parquet,
and geotiff formats.
High-level plotting programs are also provided through the
use of `Jupyter Notebooks <../user_guide/Examples.html>`_.

This software provides an interface between spatial and temporal coordinates and
the output ocean, load, solid Earth and pole tide variables as shown in the flowcharts below:

.. graphviz::
    :caption: Ocean and Load Tide Framework
    :align: center

    digraph {
        S [label="Spatial Coordinates"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        E [label="Temporal Values"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        M [label="Tide Model"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        T [label="pyTMD"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        P [label="Tide Predictions"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        H [label="Tide Heights"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        C [label="Average Tidal Currents"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        S -> T [arrowsize=0.8]
        E -> T [arrowsize=0.8]
        M -> T [arrowsize=0.8]
        T -> P [arrowsize=0.8]
        T -> H [arrowsize=0.8]
        T -> C [arrowsize=0.8]
    }

.. graphviz::
    :caption: Solid Earth and Pole Tide Framework
    :align: center

    digraph {
        S [label="Spatial Coordinates"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        E [label="Temporal Values"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#7570b3"]
        T [label="pyTMD"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="gray"]
        D [label="Solid Earth Displacements"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        P [label="Pole Tide Displacements"
            fontname="Lato"
            fontsize=11
            shape=box
            style="filled"
            color="#1b9e77"]
        S -> T [arrowsize=0.8]
        E -> T [arrowsize=0.8]
        T -> D [arrowsize=0.8]
        T -> P [arrowsize=0.8]
    }
