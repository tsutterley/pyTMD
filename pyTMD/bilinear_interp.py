#!/usr/bin/env python
u"""
bilinear_interp.py
Written by Tyler Sutterley (07/2022)
Bilinear interpolation of input data to output coordinates

CALLING SEQUENCE:
    data = bilinear_interp(ilon,ilat,idata,lon,lat)

INPUTS:
    ilon: longitude of tidal model
    ilat: latitude of tidal model
    idata: tide model data
    lat: output latitude
    lon: output longitude

OPTIONS:
    fill_value: invalid value
    dtype: output data type

OUTPUT:
    data: interpolated data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 07/2022: verify that input data is a masked array
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 03/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 12/2020: using numpy isclose to check corner points
    Updated 08/2020: check that output coordinates are within bounds
        allow small extrapolations if individual grid cells are invalid
    Updated 07/2020: split into separate function
    Updated 06/2020: use argmin and argmax in bilinear interpolation
    Updated 09/2017: Rewritten in Python
"""
import warnings
import pyTMD.interpolate

# PURPOSE: bilinear interpolation of input data to output data
def bilinear_interp(*args, **kwargs):
    """
    Bilinear interpolation of input data to output coordinates

    Parameters
    ----------
    ilon: float
        longitude of tidal model
    ilat: float
        latitude of tidal model
    idata: float
        tide model data
    lat: float
        output latitude
    lon: float
        output longitude
    fill_value: float, default np.nan
        invalid value
    dtype: obj, default np.float64
        output data type

    Returns
    -------
    data: float
        interpolated data
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.interpolate instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.interpolate.bilinear(*args, **kwargs)