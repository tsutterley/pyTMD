#!/usr/bin/env python
u"""
predict_tidal_ts.py (05/2022)
Predict tidal time series at a location using harmonic constants

CALLING SEQUENCE:
    ht = predict_tidal_ts(t,hc,con)

INPUTS:
    t: days relative to Jan 1, 1992 (48622mjd)
    hc: harmonic constant vector (complex)
    constituents: tidal constituent IDs

OUTPUT:
    ht: tidal time series reconstructed using the nodal corrections

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
    Updated 11/2019: output as numpy masked arrays instead of nan-filled arrays
    Updated 09/2019: added netcdf option to CORRECTIONS option
    Updated 08/2018: added correction option ATLAS for localized OTIS solutions
    Updated 07/2018: added option to use GSFC GOT nodal corrections
    Updated 09/2017: Rewritten in Python
"""
import numpy as np
from pyTMD.load_constituent import load_constituent
from pyTMD.load_nodal_corrections import load_nodal_corrections

def predict_tidal_ts(t, hc, constituents, deltat=0.0, corrections='OTIS'):
    """
    Predict tidal time series at a single location using harmonic constants

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
    corrections: str, default ''
        use nodal corrections from OTIS/ATLAS or GOT models

    Returns
    -------
    ht: float
        tidal time series reconstructed using the nodal corrections

    References
    ----------
    .. [1] Egbert and Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).
    """

    nt = len(t)
    # load the nodal corrections
    # convert time to Modified Julian Days (MJD)
    pu,pf,G = load_nodal_corrections(t + 48622.0, constituents,
        deltat=deltat, corrections=corrections)
    # allocate for output time series
    ht = np.ma.zeros((nt))
    ht.mask = np.zeros((nt),dtype=bool)
    # for each constituent
    for k,c in enumerate(constituents):
        if corrections in ('OTIS','ATLAS','ESR','netcdf'):
            # load parameters for each constituent
            amp,ph,omega,alpha,species = load_constituent(c)
            # add component for constituent to output tidal time series
            th = omega*t*86400.0 + ph + pu[:,k]
        elif corrections in ('GOT','FES'):
            th = G[:,k]*np.pi/180.0 + pu[:,k]
        # sum over all tides at location
        ht.data[:] += pf[:,k]*hc.real[0,k]*np.cos(th) - \
            pf[:,k]*hc.imag[0,k]*np.sin(th)
        ht.mask[:] |= np.any(hc.real.mask[0,k] | hc.imag.mask[0,k])
    # return the tidal time series
    return ht
