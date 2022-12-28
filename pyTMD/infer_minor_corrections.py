#!/usr/bin/env python
u"""
infer_minor_corrections.py (11/2022)
Return correction for minor constituents based on Richard Ray's PERTH3 code
    PERTH: PREdict Tidal Heights

CALLING SEQUENCE:
    dh = infer_minor_corrections(t,zmajor,constituents)

INPUTS:
    t: days relative to Jan 1, 1992 (48622 MJD)
    zmajor: Complex HC for given constituents/points
    constituents: tidal constituent IDs

OUTPUT:
    dh: height from minor constituents

OPTIONS:
    deltat: time correction for converting to Ephemeris Time (days)
    corrections: use nodal corrections from OTIS/ATLAS or GOT models

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    calc_astrol_longitudes.py: computes the basic astronomical mean longitudes

REFERENCES:
    A. T. Doodson and H. Warburg, "Admiralty Manual of Tides", HMSO, (1941).
    P. Schureman, "Manual of Harmonic Analysis and Prediction of Tides"
        US Coast and Geodetic Survey, Special Publication, 98, (1958).
    M. G. G. Foreman and R. F. Henry, "The harmonic analysis of tidal model
        time series", Advances in Water Resources, 12, (1989).

UPDATE HISTORY:
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 08/2020: change time variable names to not overwrite functions
        update nodal corrections for FES models
    Updated 07/2020: added function docstrings
        reduce list of minor constituents if in list of major values
    Updated 11/2019: output as numpy masked arrays instead of nan-filled arrays
    Updated 08/2018: added correction option ATLAS for localized OTIS solutions
    Updated 07/2018: added option to use GSFC GOT nodal corrections
        use the number of dates if calculating a tidal time series at a point
    Updated 09/2017: Rewritten in Python
"""
import copy
import warnings
import pyTMD.predict

def infer_minor_corrections(t, zmajor, constituents, **kwargs):
    """
    Calculate the tidal corrections for minor constituents inferred using
    major constituents

    Parameters
    ----------
    t: float
        days relative to 1992-01-01T00:00:00
    zmajor: complex
        Complex HC for given constituents/points
    constituents: list
        tidal constituent IDs
    deltat: float, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS/ATLAS or GOT models

    Returns
    -------
    dh: float
        Height from minor constituents

    References
    ----------
    .. [1] A. T. Doodson and H. Warburg, "Admiralty Manual of Tides", HMSO, (1941).
    .. [2] P. Schureman, "Manual of Harmonic Analysis and Prediction of Tides"
        *US Coast and Geodetic Survey*, Special Publication, 98, (1958).
    .. [3] Foreman and Henry, "The harmonic analysis of tidal model time
        series," *Advances in Water Resources*, 12(3), 109--120, (1989).
        `doi: 10.1016/0309-1708(89)90017-1
        <https://doi.org/10.1016/0309-1708(89)90017-1>`_
    """
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.predict instead",DeprecationWarning)
    # set default keyword arguments
    kwargs.setdefault('deltat', 0.0)
    kwargs.setdefault('corrections', 'OTIS')
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(DELTAT='deltat',CORRECTIONS='corrections')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn(f"""Deprecated keyword argument {old}.
                Changed to '{new}'""", DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])
    # call renamed version to not break workflows
    return pyTMD.predict.infer_minor(t, zmajor, constituents, **kwargs)
