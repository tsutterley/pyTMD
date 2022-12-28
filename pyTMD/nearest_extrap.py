#!/usr/bin/env python
u"""
nearest_extrap.py
Written by Tyler Sutterley (04/2022)
Uses kd-trees for nearest-neighbor extrapolation of tide model data

CALLING SEQUENCE:
    data = nearest_extrap(ilon,ilat,idata,lon,lat)

INPUTS:
    x: x-coordinates of tidal model
    y: y-coordinates of tidal model
    data: tide model data
    XI: output x-coordinates
    YI: output y-coordinates

OPTIONS:
    fill_value: invalid value
    dtype: output data type
    cutoff: return only neighbors within distance [km]
    EPSG: projection of tide model data

OUTPUT:
    DATA: extrapolated data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/

PROGRAM DEPENDENCIES:
    spatial.py: utilities for reading and writing spatial data

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
        use valid tide model points when creating kd-trees
    Updated 02/2022: fix equirectangular case for cutoffs near poles
    Updated 05/2021: set ellipsoidal major axis to WGS84 in kilometers
    Updated 03/2021: add checks to prevent runtime exception
        where there are no valid points within the input bounds
        or no points to be extrapolated
        replaced numpy bool/int to prevent deprecation warnings
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Written 12/2020
"""
import warnings
import pyTMD.interpolate

# PURPOSE: Nearest-neighbor extrapolation of valid data to output data
def nearest_extrap(*args, **kwargs):
    """
    Nearest-neighbor extrapolation of valid model data

    Parameters
    ----------
    x: float
        x-coordinates of tidal model
    y: float
        y-coordinates of tidal model
    data: float
        Tide model data
    XI: float
        Output x-coordinates
    YI: float
        Output y-coordinates
    fill_value: float, default np.nan
        Invalid value
    dtype: obj, default np.float64
        Output data type
    cutoff: float, default np.inf
        return only neighbors within distance [km]
    EPSG: str, default '4326'
        projection of tide model data

    Returns
    -------
    DATA: float
        interpolated data
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.interpolate instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.interpolate.extrapolate(*args, **kwargs)

# PURPOSE: calculate Euclidean distances between points
def distance_matrix(c1, c2):
    """
    Calculate Euclidean distances between points

    Parameters
    ----------
    c1: float
        first set of coordinates
    c2: float
        second set of coordinates

    Returns
    -------
    c: float
        Euclidean distance
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.interpolate instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.interpolate._distance(c1, c2)
