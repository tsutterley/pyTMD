#!/usr/bin/env python
u"""
read_netcdf_model.py
Written by Tyler Sutterley (12/2022)

Reads files for a tidal model and makes initial calculations to run tide program
Includes functions to extract tidal harmonic constants from OTIS tide models for
    given locations
netCDF4 files can be been compressed using gzip

Reads netCDF4 ATLAS tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/

INPUTS:
    ilon: longitude to interpolate
    ilat: latitude to interpolate
    grid_file: grid file for model (can be gzipped)
    model_files: list of model files for each constituent (can be gzipped)

OPTIONS:
    type: tidal variable to run
        z: heights
        u: horizontal transport velocities
        U: horizontal depth-averaged transport
        v: vertical transport velocities
        V: vertical depth-averaged transport
    method: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, nearest: scipy regular grid interpolations
    extrapoalte: extrapolate model using nearest-neighbors
    cutoff: extrapolation cutoff in kilometers
        set to np.inf to extrapolate for all points
    compressed: input netCDF4 files are gzip compressed
    scale: scaling factor for converting to output units

OUTPUTS:
    amplitude: amplitudes of tidal constituents
    phase: phases of tidal constituents
    D: bathymetry of tide model
    constituents: list of model constituents

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html

PROGRAM DEPENDENCIES:
    bilinear_interp.py: bilinear interpolation of data to coordinates
    nearest_extrap.py: nearest-neighbor extrapolation of data to coordinates

UPDATE HISTORY:
    Updated 12/2022: refactored tide read programs under io
    Updated 11/2022: place some imports within try/except statements
        use f-strings for formatting verbose or ascii output
    Updated 07/2022: fix setting of masked array data to NaN
    Updated 05/2022: reformat arguments to extract_netcdf_constants definition
        changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 12/2021: adjust longitude convention based on model longitude
    Updated 09/2021: fix cases where there is no mask on constituent files
    Updated 07/2021: added check that tide model files are accessible
    Updated 06/2021: add warning for tide models being entered as string
    Updated 05/2021: added option for extrapolation cutoff in kilometers
    Updated 03/2021: add extrapolation check where there are no invalid points
        prevent ComplexWarning for fill values when calculating amplitudes
        simplified inputs to be similar to binary OTIS read program
    Updated 02/2021: set invalid values to nan in extrapolation
        replaced numpy bool to prevent deprecation warning
    Updated 12/2020: added valid data extrapolation with nearest_extrap
        replace tostring with tobytes to fix DeprecationWarning
    Updated 11/2020: create function to read bathymetry and spatial coordinates
    Updated 09/2020: set bounds error to false for regular grid interpolations
        adjust dimensions of input coordinates to be iterable
        reduce number of interpolations by copying bathymetry mask to variables
    Updated 08/2020: replaced griddata with scipy regular grid interpolators
    Updated 07/2020: added function docstrings. separate bilinear interpolation
        changed TYPE variable to keyword argument. update griddata interpolation
    Updated 06/2020: use argmin and argmax in bilinear interpolation
    Written 09/2019
"""
import os
import copy
import warnings
import numpy as np
import pyTMD.io

# PURPOSE: extract tidal harmonic constants from tide models at coordinates
def extract_netcdf_constants(ilon, ilat,
    grid_file=None,
    model_files=None,
    **kwargs):
    """
    Reads files for ATLAS netCDF4 tidal models

    Makes initial calculations to run the tide program

    Spatially interpolates tidal constituents to input coordinates

    Parameters
    ----------
    ilon: float
        longitude to interpolate
    ilat: float
        latitude to interpolate
    grid_file: str or NoneType, default None
        grid file for model
    model_files: list or NoneType, default None
        list of model files for each constituent
    type: str, default 'z'
        Tidal variable to read

            - ``'z'``: heights
            - ``'u'``: horizontal transport velocities
            - ``'U'``: horizontal depth-averaged transport
            - ``'v'``: vertical transport velocities
            - ``'V'``: vertical depth-averaged transport
    method: str, default 'spline'
        Interpolation method

            - ``'bilinear'``: quick bilinear interpolation
            - ``'spline'``: scipy bivariate spline interpolation
            - ``'linear'``, ``'nearest'``: scipy regular grid interpolations
    extrapolate: bool, default False
        Extrapolate model using nearest-neighbors
    cutoff: float, default 10.0
        Extrapolation cutoff in kilometers

        Set to np.inf to extrapolate for all points
    compressed: bool, default False
        Input files are gzip compressed
    scale: float, default 1.0
        Scaling factor for converting to output units

    Returns
    -------
    amplitude: float
        amplitudes of tidal constituents
    phase: float
        phases of tidal constituents
    D: float
        bathymetry of tide model
    constituents: list
        list of model constituents
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)

    # set default keyword arguments
    kwargs.setdefault('type', 'z')
    kwargs.setdefault('method', 'spline')
    kwargs.setdefault('extrapolate', False)
    kwargs.setdefault('cutoff', 10.0)
    kwargs.setdefault('compressed', True)
    kwargs.setdefault('scale', 1.0)
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(TYPE='type', METHOD='method',
        EXTRAPOLATE='extrapolate', CUTOFF='cutoff',
        GZIP='compressed', SCALE='scale')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn(f"""Deprecated keyword argument {old}.
                Changed to '{new}'""", DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])

    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.ATLAS.extract_constants(ilon, ilat,
        grid_file=grid_file, model_files=model_files, **kwargs)

# PURPOSE: wrapper function to extend an array
def extend_array(input_array, step_size):
    """
    Wrapper function to extend an array

    Parameters
    ----------
    input_array: float
        array to extend
    step_size: float
        step size between elements of array

    Returns
    -------
    temp: float
        extended array
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.ATLAS.extend_array(input_array, step_size)

# PURPOSE: wrapper function to extend a matrix
def extend_matrix(input_matrix):
    """
    Wrapper function to extend a matrix

    Parameters
    ----------
    input_matrix: float
        matrix to extend

    Returns
    -------
    temp: float
        extended matrix
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.ATLAS.extend_matrix(input_matrix)

# PURPOSE: read grid file
def read_netcdf_grid(input_file, variable, **kwargs):
    """
    Read grid file to extract model coordinates and bathymetry

    Parameters
    ----------
    input_file: str
        input grid file
    variable: str
        Tidal variable to read

            - ``'z'``: heights
            - ``'u'``: horizontal transport velocities
            - ``'U'``: horizontal depth-averaged transport
            - ``'v'``: vertical transport velocities
            - ``'V'``: vertical depth-averaged transport

    compressed: bool, default False
        Input file is gzip compressed

    Returns
    -------
    lon: float
        longitudinal coordinates of input grid
    lat: float
        latitudinal coordinates of input grid
    bathymetry: float
        model bathymetry
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.ATLAS.read_netcdf_grid(input_file, variable, **kwargs)

# PURPOSE: read elevation file to extract real and imaginary components for
# constituent
def read_elevation_file(input_file, **kwargs):
    """
    Read elevation file to extract real and imaginary components for constituent

    Parameters
    ----------
    input_file: str
        input elevation file
    compressed: bool, default False
        Input file is gzip compressed

    Returns
    -------
    h: float
        tidal elevation
    con: str
        tidal constituent ID
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.ATLAS.read_netcdf_elevation(input_file, **kwargs)

# PURPOSE: read transport file to extract real and imaginary components for
# constituent
def read_transport_file(input_file, variable, **kwargs):
    """
    Read transport file to extract real and imaginary components for constituent

    Parameters
    ----------
    input_file: str
        input transport file
    variable: str
        Tidal variable to read

            - ``'u'``: horizontal transport velocities
            - ``'U'``: horizontal depth-averaged transport
            - ``'v'``: vertical transport velocities
            - ``'V'``: vertical depth-averaged transport

    compressed: bool, default False
        Input file is gzip compressed

    Returns
    -------
    tr: float
        tidal transport
    con: str
        tidal constituent ID
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.ATLAS.read_netcdf_transport(input_file, variable, **kwargs)
