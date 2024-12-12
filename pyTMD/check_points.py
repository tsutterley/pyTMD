#!/usr/bin/env python
u"""
check_points.py
Written by Tyler Sutterley (12/2024)
Check if points are within a tide model domain

OTIS format tidal solutions provided by Oregon State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/
Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
or Finite Element Solution (FES) models provided by AVISO

INPUTS:
    x: x-coordinates in projection EPSG
    y: y-coordinates in projection EPSG

OPTIONS:
    DIRECTORY: working data directory for tide models
    MODEL: Tide model to use
    GZIP: Tide model files are gzip compressed
    DEFINITION_FILE: Tide model definition file for use
    EPSG: input coordinate system
        default: 3031 Polar Stereographic South, WGS84
    METHOD: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, nearest: scipy regular grid interpolations

OUTPUTS:
    valid: array describing if input coordinate is within model domain

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    crs.py: Coordinate Reference System (CRS) routines
    io/model.py: retrieves tide model parameters for named tide models
    io/OTIS.py: extract tidal harmonic constants from OTIS tide models
    io/ATLAS.py: extract tidal harmonic constants from ATLAS netcdf models
    io/GOT.py: extract tidal harmonic constants from GSFC GOT models
    io/FES.py: extract tidal harmonic constants from FES tide models
    interpolate.py: interpolation routines for spatial data

UPDATE HISTORY:
    Updated 12/2024: deprecated in favor of compute.tide_masks
    Updated 09/2024: use JSON database for known model parameters
        drop support for the ascii definition file format
    Updated 07/2024: renamed format for ATLAS to ATLAS-compact
        renamed format for netcdf to ATLAS-netcdf
        renamed format for FES to FES-netcdf and added FES-ascii
        renamed format for GOT to GOT-ascii and added GOT-netcdf
    Updated 04/2024: use wrapper to importlib for optional dependencies
    Updated 12/2023: use new crs class for coordinate reprojection
    Updated 08/2023: changed ESR netCDF4 format to TMD3 format
    Updated 04/2023: using pathlib to define and expand paths
    Updated 03/2023: add basic variable typing to function inputs
    Updated 12/2022: refactored tide read programs under io
        refactored bilinear interpolation routine
    Updated 11/2022: place some imports within try/except statements
        use f-strings for formatting verbose or ascii output
    Updated 05/2022: added ESR netCDF4 formats to list of model types
        updated keyword arguments to read tide model programs
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 09/2021: refactor to use model class for files and attributes
    Updated 07/2021: added check that tide model directory is accessible
    Updated 06/2021: add try/except for input projection strings
    Written 05/2021
"""
from __future__ import print_function, annotations

import pathlib
import warnings
import pyTMD.compute

# PURPOSE: compute tides at points and times using tide model algorithms
def check_points(*args, **kwargs):
    """
    Check if points are within a tide model domain

    Parameters
    ----------
    x: np.ndarray
        x-coordinates in projection EPSG
    y: np.ndarray
        y-coordinates in projection EPSG
    DIRECTORY: str or NoneType, default None
        working data directory for tide models
    MODEL: str or NoneType, default None
        Tide model to use
    GZIP: bool, default False
        Tide model files are gzip compressed
    DEFINITION_FILE: str or NoneType, default None
        Tide model definition file for use
    EPSG: str or int, default: 3031 (Polar Stereographic South, WGS84)
        Input coordinate system
    METHOD: str, default 'spline'
        interpolation method

            - ```bilinear```: quick bilinear interpolation
            - ```spline```: scipy bivariate spline interpolation
            - ```linear```, ```nearest```: scipy regular grid interpolations

    Returns
    -------
    valid: bool
        array describing if input coordinate is within model domain
    """
    warnings.warn("Deprecated. Please use pyTMD.compute instead",
        DeprecationWarning)
    return pyTMD.compute.tide_masks(*args, **kwargs)
