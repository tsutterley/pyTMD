#!/usr/bin/env python
u"""
load_nodal_corrections.py
Written by Tyler Sutterley (04/2023)
Calculates the nodal corrections for tidal constituents
Modification of ARGUMENTS fortran subroutine by Richard Ray 03/1999

CALLING SEQUENCE:
    pu,pf,G = load_nodal_corrections(MJD, constituents)

INPUTS:
    MJD: Modified Julian Day of input date
    constituents: tidal constituent IDs

OUTPUTS:
    pu,pf: nodal corrections for the constituents
    G: phase correction in degrees

OPTIONS:
    deltat: time correction for converting to Ephemeris Time (days)
    corrections: use nodal corrections from OTIS/ATLAS or GOT models

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    astro.py: computes the basic astronomical mean longitudes

REFERENCES:
    A. T. Doodson and H. Warburg, "Admiralty Manual of Tides", HMSO, (1941).
    P. Schureman, "Manual of Harmonic Analysis and Prediction of Tides"
        US Coast and Geodetic Survey, Special Publication, 98, (1958).
    M. G. G. Foreman and R. F. Henry, "The harmonic analysis of tidal model
        time series", Advances in Water Resources, 12, (1989).
    G. D. Egbert and S. Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).

UPDATE HISTORY:
    Updated 04/2023: using renamed astro mean_longitudes function
        deprecated in favor of pyTMD.arguments function
    Updated 03/2023: add basic variable typing to function inputs
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: added ESR netCDF4 formats to list of model types
        changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 12/2020: fix k1 for FES models
    Updated 08/2020: change time variable names to not overwrite functions
        update nodal corrections for FES models
    Updated 07/2020: added function docstrings.  add shallow water constituents
    Updated 09/2019: added netcdf option to corrections option
    Updated 08/2018: added correction option ATLAS for localized OTIS solutions
    Updated 07/2018: added option to use GSFC GOT nodal corrections
    Updated 09/2017: Rewritten in Python
    Rewritten in Matlab by Lana Erofeeva 01/2003
    Written by Richard Ray 03/1999
"""
from __future__ import annotations

import copy
import warnings
import pyTMD.arguments

def load_nodal_corrections(*args, **kwargs):
    """
    Calculates the nodal corrections for tidal constituents
    [1]_ [2]_ [3]_ [4]_

    Parameters
    ----------
    MJD: np.ndarray
        modified julian day of input date
    constituents: list
        tidal constituent IDs
    deltat: float or np.ndarray, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS/ATLAS or GOT models

    Returns
    -------
    pu: np.ndarray
        nodal correction for the constituent amplitude
    pf: np.ndarray
        nodal correction for the constituent phase
    G: np.ndarray
        phase correction in degrees

    References
    ----------
    .. [1] A. T. Doodson and H. D. Warburg, "Admiralty Manual of Tides",
        HMSO, London, (1941).
    .. [2] P. Schureman, "Manual of Harmonic Analysis and Prediction of Tides,"
        *US Coast and Geodetic Survey*, Special Publication, 98, (1958).
    .. [3] M. G. G. Foreman and R. F. Henry, "The harmonic analysis of tidal
        model time series," *Advances in Water Resources*, 12(3), 109--120,
        (1989). `doi: 10.1016/0309-1708(89)90017-1
        <https://doi.org/10.1016/0309-1708(89)90017-1>`_
    .. [4] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of
        Barotropic Ocean Tides," *Journal of Atmospheric and Oceanic
        Technology*, 19(2), 183--204, (2002).
        `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

    .. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
    """
    # raise warning for deprecated function call
    warnings.filterwarnings("module")
    warnings.warn("Deprecated. Please use pyTMD.arguments instead",
        DeprecationWarning)
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(DELTAT='deltat',CORRECTIONS='corrections')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn(f"""Deprecated keyword argument {old}.
                Changed to '{new}'""", DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])
    warnings.filterwarnings("ignore")
    # call updated function to not break current workflows
    return pyTMD.arguments(*args, **kwargs)
