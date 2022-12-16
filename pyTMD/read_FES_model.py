#!/usr/bin/env python
u"""
read_FES_model.py
Written by Tyler Sutterley (12/2022)

Reads files for a tidal model and makes initial calculations to run tide program
Includes functions to extract tidal harmonic constants from the
    FES (Finite Element Solution) tide models for given locations
ascii and netCDF4 files can be been compressed using gzip

Reads ascii and netCDF4 FES tidal solutions provided by AVISO
    https://www.aviso.altimetry.fr/data/products/auxiliary-products/
        global-tide-fes.html

INPUTS:
    ilon: longitude to interpolate
    ilat: latitude to interpolate
    model_files: list of model files for each constituent

OPTIONS:
    type: tidal variable to run
        z: heights
        u: horizontal transport velocities
        v: vertical transport velocities
    version: model version to run
        FES1999
        FES2004
        FES2012
        FES2014
        EOT20
    method: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, nearest: scipy regular grid interpolations
    extrapolate: extrapolate model using nearest-neighbors
    cutoff: extrapolation cutoff in kilometers
        set to np.inf to extrapolate for all points
    compressed: input files are gzip compressed
    scale: scaling factor for converting to output units

OUTPUTS:
    amplitude: amplitudes of tidal constituents
    phase: phases of tidal constituents

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
    Updated 05/2022: reformat arguments to extract_FES_constants definition
        changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
        fix netCDF4 masks for nan values
    Updated 01/2022: added global Empirical Ocean Tide model (EOT20)
    Updated 12/2021: adjust longitude convention based on model longitude
    Updated 07/2021: added check that tide model files are accessible
    Updated 06/2021: add warning for tide models being entered as string
    Updated 05/2021: added option for extrapolation cutoff in kilometers
    Updated 03/2021: add extrapolation check where there are no invalid points
        prevent ComplexWarning for fill values when calculating amplitudes
        simplified inputs to be similar to binary OTIS read program
        replaced numpy bool/int to prevent deprecation warnings
        use uuid for reading from gzipped netCDF4 files
    Updated 02/2021: set invalid values to nan in extrapolation
        replaced numpy bool to prevent deprecation warning
    Updated 12/2020: added nearest-neighbor data extrapolation
    Updated 09/2020: set bounds error to false for regular grid interpolations
        adjust dimensions of input coordinates to be iterable
    Updated 08/2020: replaced griddata with scipy regular grid interpolators
    Written 07/2020
"""
import os
import copy
import warnings
import numpy as np
import pyTMD.io

# PURPOSE: extract tidal harmonic constants from tide models at coordinates
def extract_FES_constants(ilon, ilat, model_files=None, **kwargs):
    """
    Reads files for an ascii or netCDF4 tidal model

    Makes initial calculations to run the tide program

    Spatially interpolates tidal constituents to input coordinates

    Parameters
    ----------
    ilon: float
        longitude to interpolate
    ilat: float
        latitude to interpolate
    model_files: list or NoneType, default None
        list of model files for each constituent
    type: str, default 'z'
        Tidal variable to read

            - ``'z'``: heights
            - ``'u'``: horizontal transport velocities
            - ``'v'``: vertical transport velocities
    version: str or NoneType, default None
        Model version to read

            - ``'FES1999'``
            - ``'FES2004'``
            - ``'FES2012'``
            - ``'FES2014'``
            - ``'EOT20'``
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
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)

    # set default keyword arguments
    kwargs.setdefault('type', 'z')
    kwargs.setdefault('version', None)
    kwargs.setdefault('method', 'spline')
    kwargs.setdefault('extrapolate', False)
    kwargs.setdefault('cutoff', 10.0)
    kwargs.setdefault('compressed', False)
    kwargs.setdefault('scale', 1.0)
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(TYPE='type',VERSION='version',
        METHOD='method',EXTRAPOLATE='extrapolate',CUTOFF='cutoff',
        GZIP='compressed',SCALE='scale')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn(f"""Deprecated keyword argument {old}.
                Changed to '{new}'""", DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])
    # call renamed version to not break workflows
    return pyTMD.io.FES.extract_constants(ilon, ilat,
        model_files=model_files, **kwargs)

# PURPOSE: read FES ascii tide model grid files
def read_ascii_file(input_file, **kwargs):
    """
    Read FES (Finite Element Solution) tide model file

    Parameters
    ----------
    input_file: str
        model file
    compressed: bool, default False
        Input file is gzip compressed

    Returns
    -------
    hc: complex form of tidal constituent oscillation
    lon: longitude of tidal model
    lat: latitude of tidal model
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.FES.read_ascii_file(input_file, **kwargs)

# PURPOSE: read FES netCDF4 tide model files
def read_netcdf_file(input_file, **kwargs):
    """
    Read FES (Finite Element Solution) tide model netCDF4 file

    Parameters
    ----------
    input_file: str
        model file
    type: str or NoneType, default None
        Tidal variable to read

            - ``'z'``: heights
            - ``'u'``: horizontal transport velocities
            - ``'v'``: vertical transport velocities
    version: str or NoneType, default None
        FES model version
    compressed: bool, default False
        Input file is gzip compressed

    Returns
    -------
    hc: complex
        complex form of tidal constituent oscillation
    lon: float
        longitude of tidal model
    lat: float
        latitude of tidal model
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.FES.read_netcdf_file(input_file, **kwargs)
