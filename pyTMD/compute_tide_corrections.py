#!/usr/bin/env python
u"""
compute_tide_corrections.py
Written by Tyler Sutterley (01/2024)
Calculates tidal elevations for correcting elevation or imagery data

Ocean and Load Tides
Uses OTIS format tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/
Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
or Finite Element Solution (FES) models provided by AVISO

Long-Period Equilibrium Tides (LPET)
Calculates long-period equilibrium tidal elevations for correcting
elevation or imagery data from the summation of fifteen spectral lines
    https://doi.org/10.1111/j.1365-246X.1973.tb03420.x

Load Pole Tides (LPT)
Calculates radial load pole tide displacements following IERS Convention
(2010) guidelines for correcting elevation or imagery data
    https://iers-conventions.obspm.fr/chapter7.php

Ocean Pole Tides (OPT)
Calculates radial ocean pole load tide displacements following IERS Convention
(2010) guidelines for correcting elevation or imagery data
    https://iers-conventions.obspm.fr/chapter7.php

Ocean Pole Tides (SET)
Calculates radial Solid Earth tide displacements following IERS Convention
(2010) guidelines for correcting elevation or imagery data
    https://iers-conventions.obspm.fr/chapter7.php

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
    time.py: utilities for calculating time operations
    spatial: utilities for reading, writing and operating on spatial data
    utilities.py: download and management utilities for syncing files
    arguments.py: load the nodal corrections for tidal constituents
    astro.py: computes the basic astronomical mean longitudes
    crs.py: Coordinate Reference System (CRS) routines
    predict.py: predict tide values using harmonic constants
    io/model.py: retrieves tide model parameters for named tide models
    io/OTIS.py: extract tidal harmonic constants from OTIS tide models
    io/ATLAS.py: extract tidal harmonic constants from netcdf models
    io/GOT.py: extract tidal harmonic constants from GSFC GOT models
    io/FES.py: extract tidal harmonic constants from FES tide models
    interpolate.py: interpolation routines for spatial data

UPDATE HISTORY:
    Updated 01/2024: made the inferrence of minor constituents an option
        refactored lunisolar ephemerides functions
        deprecated in favor of refactored compute.py module
    Updated 12/2023: use new crs class for coordinate reprojection
    Updated 08/2023: changed ESR netCDF4 format to TMD3 format
    Updated 05/2023: use timescale class for time conversion operations
        use defaults from eop module for pole tide and EOP files
        add option for using higher resolution ephemerides from JPL
    Updated 04/2023: added function for radial solid earth tides
        using pathlib to define and expand paths
    Updated 03/2023: add basic variable typing to function inputs
        added function for long-period equilibrium tides
        added function for radial load pole tides
        added function for radial ocean pole tides
    Updated 12/2022: refactored tide read and prediction programs
    Updated 11/2022: place some imports within try/except statements
        use f-strings for formatting verbose or ascii output
    Updated 05/2022: added ESR netCDF4 formats to list of model types
        updated keyword arguments to read tide model programs
        added option to apply flexure to heights for applicable models
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 12/2021: added function to calculate a tidal time series
        verify coordinate dimensions for each input data type
        added option for converting from LORAN times to UTC
    Updated 09/2021: refactor to use model class for files and attributes
    Updated 07/2021: can use numpy datetime arrays as input time variable
        added function for determining the input spatial variable type
        added check that tide model directory is accessible
    Updated 06/2021: added new Gr1km-v2 1km Greenland model from ESR
        add try/except for input projection strings
    Updated 05/2021: added option for extrapolation cutoff in kilometers
    Updated 03/2021: added TPXO9-atlas-v4 in binary OTIS format
        simplified netcdf inputs to be similar to binary OTIS read program
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 12/2020: added valid data extrapolation with nearest_extrap
    Updated 11/2020: added model constituents from TPXO9-atlas-v3
    Updated 08/2020: using builtin time operations.
        calculate difference in leap seconds from start of epoch
        using conversion protocols following pyproj-2 updates
    Updated 07/2020: added function docstrings, FES2014 and TPXO9-atlas-v2
        use merged delta time files combining biannual, monthly and daily files
    Updated 03/2020: added TYPE, TIME, FILL_VALUE and METHOD options
    Written 03/2020
"""
from __future__ import print_function, annotations

import warnings
import numpy as np
import pyTMD.compute

# PURPOSE: wrapper function for computing corrections
def compute_corrections(
        x: np.ndarray, y: np.ndarray, delta_time: np.ndarray,
        CORRECTION: str = 'ocean',
        **kwargs
    ):
    """
    Wrapper function to compute tide corrections at points and times

    Parameters
    ----------
    x: np.ndarray
        x-coordinates in projection EPSG
    y: np.ndarray
        y-coordinates in projection EPSG
    delta_time: np.ndarray
        seconds since EPOCH or datetime array
    CORRECTION: str, default 'ocean'
        Correction type to compute

            - ``'ocean'``: ocean tide from model constituents
            - ``'load'``: load tide from model constituents
            - ``'LPET'``: long-period equilibrium tide
            - ``'LPT'``: solid earth load pole tide
            - ``'OPT'``: ocean pole tide
            - ``'SET'``: solid earth tide
    **kwargs: dict
        keyword arguments for correction functions

    Returns
    -------
    correction: np.ndarray
        tidal correction at coordinates and time in meters
    """
    warnings.warn("Deprecated. Please use pyTMD.compute instead",
        DeprecationWarning)
    if CORRECTION.lower() in ('ocean', 'load'):
        return pyTMD.compute.tide_elevations(x, y, delta_time, **kwargs)
    elif (CORRECTION.upper() == 'LPET'):
        return pyTMD.compute.LPET_elevations(x, y, delta_time, **kwargs)
    elif (CORRECTION.upper() == 'LPT'):
        return pyTMD.compute.LPT_displacements(x, y, delta_time, **kwargs)
    elif (CORRECTION.upper() == 'OPT'):
        return pyTMD.compute.OPT_displacements(x, y, delta_time, **kwargs)
    elif (CORRECTION.upper() == 'SET'):
        return pyTMD.compute.SET_displacements(x, y, delta_time, **kwargs)
    else:
        raise ValueError(f'Unrecognized correction type: {CORRECTION}')

def compute_tide_corrections(*args, **kwargs):
    return pyTMD.compute.tide_elevations(*args, **kwargs)

def compute_LPET_corrections(*args, **kwargs):
    return pyTMD.compute.LPET_elevations(*args, **kwargs)

def compute_LPT_corrections(*args, **kwargs):
    return pyTMD.compute.LPT_displacements(*args, **kwargs)

def compute_OPT_corrections(*args, **kwargs):
    return pyTMD.compute.OPT_displacements(*args, **kwargs)

def compute_SET_corrections(*args, **kwargs):
    return pyTMD.compute.SET_displacements(*args, **kwargs)
