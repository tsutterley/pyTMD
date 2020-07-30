#!/usr/bin/env python
u"""
predict_tide.py (07/2020)
Predict tides at a single time using harmonic constants

CALLING SEQUENCE:
    ht = predict_tide(time,hc,con)

INPUTS:
    time: days relative to Jan 1, 1992 (48622mjd)
    hc: harmonic constant vector (complex)
    constituents: tidal constituent IDs

OUTPUT:
    ht: tide values reconstructed using the nodal corrections

OPTIONS:
    DELTAT: time correction for converting to Ephemeris Time (days)
    CORRECTIONS: use nodal corrections from OTIS/ATLAS or GOT models

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    load_constituent.py: loads parameters for a given tidal constituent
    load_nodal_corrections.py: loads nodal corrections for tidal constituents

UPDATE HISTORY:
    Updated 07/2020: added function docstrings
    Updated 11/2019: can output an array of heights with a single time stamp
        such as for estimating tide height maps from imagery
    Updated 09/2019: added netcdf option to CORRECTIONS option
    Updated 08/2018: added correction option ATLAS for localized OTIS solutions
    Updated 07/2018: added option to use GSFC GOT nodal corrections
    Updated 09/2017: Rewritten in Python
"""
import numpy as np
from pyTMD.load_constituent import load_constituent
from pyTMD.load_nodal_corrections import load_nodal_corrections

def predict_tide(time,hc,constituents,DELTAT=0.0,CORRECTIONS='OTIS'):
    """
    Predict tides at a single time using harmonic constants

    Arguments
    ---------
    time: days relative to 1992-01-01T00:00:00
    hc: harmonic constant vector (complex)
    constituents: tidal constituent IDs

    Keyword arguments
    -----------------
    DELTAT: time correction for converting to Ephemeris Time (days)
    CORRECTIONS: use nodal corrections from OTIS/ATLAS or GOT models

    Returns
    -------
    ht: tide values reconstructed using the nodal corrections
    """

    #-- number of points and number of constituents
    npts,nc = np.shape(hc)
    #-- load the nodal corrections
    pu,pf,G = load_nodal_corrections(time + 48622.0, constituents,
        DELTAT=DELTAT, CORRECTIONS=CORRECTIONS)
    #-- allocate for output tidal elevation
    ht = np.ma.zeros((npts))
    #-- for each constituent
    for k,c in enumerate(constituents):
        if CORRECTIONS in ('OTIS','ATLAS','netcdf'):
            #-- load parameters for each constituent
            amp,ph,omega,alpha,species = load_constituent(c)
            #-- add component for constituent to output tidal elevation
            th = omega*time*86400.0 + ph + pu[0,k]
        elif CORRECTIONS in ('GOT','FES'):
            th = G[0,k]*np.pi/180.0 + pu[0,k]
        #-- sum over all tides
        ht += pf[0,k]*hc.real[:,k]*np.cos(th) - pf[0,k]*hc.imag[:,k]*np.sin(th)
    #-- return the tidal elevation after removing singleton dimensions
    return np.squeeze(ht)
