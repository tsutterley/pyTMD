#!/usr/bin/env python
u"""
read_GOT_model.py
Written by Tyler Sutterley (12/2022)

Reads files for Richard Ray's Global Ocean Tide (GOT) models and makes initial
    calculations to run the tide program
Includes functions to extract tidal harmonic constants out of a tidal model for
    given locations

INPUTS:
    ilon: longitude to interpolate
    ilat: latitude to interpolate
    model_files: list of model files for each constituent

OPTIONS:
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
    constituents: list of model constituents

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/

PROGRAM DEPENDENCIES:
    bilinear_interp.py: bilinear interpolation of data to coordinates
    nearest_extrap.py: nearest-neighbor extrapolation of data to coordinates

UPDATE HISTORY:
    Updated 12/2022: refactored tide read programs under io
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: reformat arguments to extract_GOT_constants definition
        changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
    Updated 12/2021: adjust longitude convention based on model longitude
    Updated 07/2021: added check that tide model files are accessible
    Updated 06/2021: add warning for tide models being entered as string
    Updated 05/2021: added option for extrapolation cutoff in kilometers
    Updated 03/2021: add extrapolation check where there are no invalid points
        prevent ComplexWarning for fill values when calculating amplitudes
        simplified inputs to be similar to binary OTIS read program
        replaced numpy bool/int to prevent deprecation warnings
    Updated 02/2021: set invalid values to nan in extrapolation
        replaced numpy bool to prevent deprecation warning
    Updated 12/2020: added valid data extrapolation with nearest_extrap
    Updated 09/2020: set bounds error to false for regular grid interpolations
        adjust dimensions of input coordinates to be iterable
    Updated 08/2020: replaced griddata with scipy regular grid interpolators
    Updated 07/2020: added function docstrings. separate bilinear interpolation
        update griddata interpolation. add option for compression
    Updated 06/2020: use argmin and argmax in bilinear interpolation
    Updated 11/2019: find invalid mask points for each constituent
    Updated 09/2019: output as numpy masked arrays instead of nan-filled arrays
    Updated 07/2019: interpolate fill value mask with bivariate splines
    Updated 12/2018: python3 compatibility updates for division and zip
    Updated 10/2018: added scale as load tides are in mm and ocean are in cm
    Updated 08/2018: added multivariate spline interpolation option
    Written 07/2018
"""
from __future__ import division

import os
import re
import copy
import warnings
import numpy as np
import pyTMD.io

# PURPOSE: extract tidal harmonic constants out of GOT model at coordinates
def extract_GOT_constants(ilon, ilat, model_files=None, **kwargs):
    """
    Reads files for Richard Ray's Global Ocean Tide (GOT) models

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
    constituents: list
        list of model constituents
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)

    # set default keyword arguments
    kwargs.setdefault('method', 'spline')
    kwargs.setdefault('extrapolate', False)
    kwargs.setdefault('cutoff', 10.0)
    kwargs.setdefault('compressed', False)
    kwargs.setdefault('scale', 1.0)
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(METHOD='method',
        EXTRAPOLATE='extrapolate',CUTOFF='cutoff',
        GZIP='compressed',SCALE='scale')
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
    return pyTMD.io.GOT.extract_constants(ilon, ilat,
        model_files=model_files, **kwargs)

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
    return pyTMD.io.GOT.extend_array(input_array, step_size)

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
    return pyTMD.io.GOT.extend_matrix(input_matrix)

# PURPOSE: read GOT model grid files
def read_GOT_grid(input_file, **kwargs):
    """
    Read Richard Ray's Global Ocean Tide (GOT) model file

    Parameters
    ----------
    input_file: str
        Model file
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
    cons: str
        tidal constituent ID
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.io.GOT.read_ascii_file(input_file, **kwargs)