#!/usr/bin/env python
u"""
convert_crs.py
Written by Tyler Sutterley (12/2023)
Converts points to and from Coordinates Reference Systems (CRS)

CALLING SEQUENCE:
    x, y = convert_crs(lon, lat, PROJ, 'F')
    lon, lat = convert_crs(x, y, PROJ, 'B')

INPUTS:
    i1: longitude ('F') or projection easting x ('B')
    i2: latitude ('F') or projection northing y ('B')
    PROJ: spatial reference system code for coordinate transformations
    BF: backwards ('B') or forward ('F') translations

OPTIONS:
    EPSG: spatial reference system code for input (F) and output (B) coordinates

OUTPUTS:
    o1: projection easting x ('F') or longitude ('B')
    o2: projection northing y ('F') or latitude ('B')

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
        https://pyproj4.github.io/pyproj/

UPDATE HISTORY:
    Updated 03/2023: deprecated in favor of crs class
    Updated 03/2023: add basic variable typing to function inputs
        renamed coordinate reference system conversion functions
    Updated 02/2023: use named exception before passing to custom
    Updated 11/2022: place some imports within try/except statements
        use f-strings for formatting verbose or ascii output
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 09/2021: added function for using custom projections
    Updated 06/2021: added 3413 for new 1km Greenland model from ESR
    Updated 08/2020: using conversion protocols following pyproj-2 updates
        https://pyproj4.github.io/pyproj/stable/gotchas.html
    Updated 07/2020: added function docstrings. changed function name
    Updated 03/2020: remove commented coordinate conversion functions
    Updated 11/2019: using pyproj for coordinate conversions
    Written 09/2017
"""
from __future__ import annotations

import warnings
import numpy as np
import pyTMD.crs

def convert_crs(
        i1: np.ndarray,
        i2: np.ndarray,
        PROJ: str,
        BF: str,
        EPSG: int | str = 4326
    ):
    """
    Converts points to and from Coordinates Reference Systems (CRS)

    Parameters
    ----------
    i1: np.ndarray
        Longitude (``'F'``) or projected x-coordinates (``'B'``)
    i2: np.ndarray
        Latitude (``'F'``) or projected y-coordinates (``'B'``)
    PROJ: str
        Spatial reference system code for coordinate transformations
    BF: str
        Direction of translation

            - ``'B'``: backwards
            - ``'F'``: forwards
    EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
        EPSG code for input (``'F'``) or output (``'B'``) coordinate system

    Returns
    -------
    o1: np.ndarray
        Projected x-coordinates (``'F'``) or longitude (``'B'``)
    o2: np.ndarray
        Projected y-coordinates (``'F'``) or latitude (``'B``')
    """
    warnings.warn("Deprecated. Please use pyTMD.crs().convert instead",
        DeprecationWarning)
    # call updated function to not break current workflows
    return pyTMD.crs().convert(i1, i2, PROJ, BF, EPSG=EPSG)

# PURPOSE: try to get the projection information
def crs_from_input(PROJECTION: int | str):
    """
    Attempt to get the Coordinate Reference System for an input code

    Parameters
    ----------
    PROJECTION: int or str
        Coordinate Reference System code
    """
    warnings.warn("Deprecated. Please use pyTMD.crs().from_input instead",
        DeprecationWarning)
    # call updated function to not break current workflows
    return pyTMD.crs().from_input(PROJECTION)
