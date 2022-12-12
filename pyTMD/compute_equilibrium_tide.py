#!/usr/bin/env python
u"""
compute_equilibrium_tide.py (04/2022)
Calculates the long-period equilibrium ocean tides
Uses the summation of fifteen tidal spectral lines from the
    Cartwright-Tayler-Edden tables to compute the long-period tides

INPUTS:
    t: days relative to Jan 1, 1992 (48622 MJD)
    lat: latitudes in degrees

OUTPUTS:
    lpet: long-period equilibrium tide in meters

REFERENCES:
    Cartwright & Tayler, Geophys. J. R.A.S., 23, 45, 1971.
    Cartwright & Edden, Geophys. J. R.A.S., 33, 253, 1973.

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2020: check sizes of inputs to check if using grid or drift data
        drift: drift buoys or satellite/airborne altimetry (time per data point)
        grid: spatial grids or images (single time for all data points)
    Updated 08/2020: use Load love number value of 0.693 for gamma_2
    Written 08/2020
"""
import warnings
import pyTMD.predict

def compute_equilibrium_tide(t, lat):
    """
    Parameters
    ----------
    t: float
        time (days relative to January 1, 1992)
    lat: float
        latitudes (degrees)

    Returns
    -------
    lpet: float
        long-period equilibrium tide in meters

    References
    ----------
    .. [1] Cartwright & Tayler, Geophys. J. R.A.S., 23, 45, 1971.
    .. [2] Cartwright & Edden, Geophys. J. R.A.S., 33, 253, 1973.
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.predict instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.predict.equilibrium_tide(t, lat)
