#!/usr/bin/env python
u"""
convert_ll_xy.py
Written by Tyler Sutterley (03/2023)
Wrapper function to convert points to and from Coordinates Reference Systems

CALLING SEQUENCE:
    x,y = convert_ll_xy(lon,lat,PROJ,'F')
    lon,lat = convert_ll_xy(x,y,PROJ,'B')

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
    Updated 03/2023: deprecated in favor of convert_crs.py
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
import warnings
import pyTMD.convert_crs

def convert_ll_xy(*args, **kwargs):
    """
    Converts lat/lon points to and from projected coordinates

    Parameters
    ----------
    i1: float
        Longitude (``'F'``) or projected x-coordinates (``'B'``)
    i2: float
        Latitude (``'F'``) or projected y-coordinates (``'B'``)
    PROJ: str
        Spatial reference system code for coordinate transformations
    BF: str
        Direction of translation

            - ``'B'``: backwards
            - ``'F'``: forwards
    EPSG: int, default 4326 (WGS84 Latitude/Longitude)
        EPSG code for input (``'F'``) or output (``'B'``) coordinate system

    Returns
    -------
    o1: float
        Projected x-coordinates (``'F'``) or longitude (``'B'``)
    o2: float
        Projected y-coordinates (``'F'``) or latitude (``'B``')
    """
    warnings.filterwarnings("module")
    warnings.warn("Deprecated. Please use pyTMD.convert_crs instead",
        DeprecationWarning)
    warnings.filterwarnings("ignore")
    # call updated function to not break current workflows
    return pyTMD.convert_crs(*args, **kwargs)
