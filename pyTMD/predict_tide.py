#!/usr/bin/env python
u"""
predict_tide.py (05/2022)
Predict tides at a single time using harmonic constants

CALLING SEQUENCE:
    ht = predict_tide(t,hc,con)

INPUTS:
    t: days relative to Jan 1, 1992 (48622mjd)
    hc: harmonic constant vector (complex)
    constituents: tidal constituent IDs

OUTPUT:
    ht: tide values reconstructed using the nodal corrections

OPTIONS:
    deltat: time correction for converting to Ephemeris Time (days)
    corrections: use nodal corrections from OTIS/ATLAS or GOT models

REFERENCES:
    G. D. Egbert and S. Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    load_constituent.py: loads parameters for a given tidal constituent
    load_nodal_corrections.py: loads nodal corrections for tidal constituents

UPDATE HISTORY:
    Updated 05/2022: added ESR netCDF4 formats to list of model types
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 09/2020: append output mask over each constituent
    Updated 08/2020: change time variable names to not overwrite functions
    Updated 07/2020: added function docstrings
    Updated 11/2019: can output an array of heights with a single time stamp
        such as for estimating tide height maps from imagery
    Updated 09/2019: added netcdf option to CORRECTIONS option
    Updated 08/2018: added correction option ATLAS for localized OTIS solutions
    Updated 07/2018: added option to use GSFC GOT nodal corrections
    Updated 09/2017: Rewritten in Python
"""
import warnings
import pyTMD.predict

def predict_tide(*args, **kwargs):
    """
    Predict tides at a single time using harmonic constants

    Parameters
    ----------
    t: float
        days relative to 1992-01-01T00:00:00
    hc: complex
        harmonic constant vector
    constituents: list
        tidal constituent IDs
    deltat: float, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS/ATLAS or GOT models

    Returns
    -------
    ht: float
        tide values reconstructed using the nodal corrections

    References
    ----------
    .. [1] Egbert and Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides," *Journal of Atmospheric and Oceanic Technology*,
        19(2), 183--204, (2002).
        `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

    .. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.predict instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.predict.map(*args, **kwargs)
