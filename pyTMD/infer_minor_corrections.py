#!/usr/bin/env python
u"""
infer_minor_corrections.py (07/2020)
Return correction for minor constituents based on Richard Ray's PERTH3 code
    PERTH: PREdict Tidal Heights

CALLING SEQUENCE:
    dh = infer_minor_corrections(time,zmajor,constituents)

INPUTS:
    constituents: tidal constituent IDs
    zmajor: Complex HC for GIVEN constituents/points
    time: days relative to Jan 1, 1992 (48622mjd)

OUTPUT:
    dh: height from minor constituents

OPTIONS:
    DELTAT: time correction for converting to Ephemeris Time (days)
    CORRECTIONS: use nodal corrections from OTIS/ATLAS or GOT models

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
    Updated 07/2020: added function docstrings
        reduce list of minor constituents if in list of major values
    Updated 11/2019: output as numpy masked arrays instead of nan-filled arrays
    Updated 08/2018: added correction option ATLAS for localized OTIS solutions
    Updated 07/2018: added option to use GSFC GOT nodal corrections
        use the number of dates if calculating a tidal time series at a point
    Updated 09/2017: Rewritten in Python
"""
import numpy as np
from pyTMD.calc_astrol_longitudes import calc_astrol_longitudes

def infer_minor_corrections(time,zmajor,constituents,DELTAT=0.0,CORRECTIONS=''):
    """
    Calculate the tidal corrections for minor constituents inferred using
    major constituents

    Arguments
    ---------
    constituents: tidal constituent IDs
    zmajor: Complex HC for GIVEN constituents/points
    time: days relative to 1992-01-01T00:00:00

    Keyword arguments
    -----------------
    DELTAT: time correction for converting to Ephemeris Time (days)
    CORRECTIONS: use nodal corrections from OTIS/ATLAS or GOT models

    Returns
    -------
    dh: height from minor constituents
    """
    #-- degrees to radians
    dtr = np.pi/180.0
    #-- number of constituents
    npts,nc = np.shape(zmajor)
    nt = 1 if (np.ndim(time) == 0) else len(time)
    #-- number of data points to calculate if running time series/drift/map
    n = nt if ((npts == 1) & (nt > 1)) else npts
    #-- allocate for output elevation correction
    dh = np.ma.zeros((n))
    #-- convert time from days relative to Jan 1, 1992 to modified Julian days
    time_mjd = 48622.0 + time
    cindex = ['q1','o1','p1','k1','n2','m2','s2','k2']
    #-- re-order major tides to correspond to order of cindex
    z8 = np.ma.zeros((n,8),dtype=np.complex64)
    ni = 0
    for i,c in enumerate(cindex):
        j = [j for j,val in enumerate(constituents) if val == c]
        if j:
            j1, = j
            z8[:,i] = zmajor[:,j1]
            ni += 1

    if (ni < 6):
        raise Exception('Not enough constituents for inference')

    #-- list of minor constituents
    minor = ['2q1','sigma1','rho1','m1','m1','chi1','pi1','phi1','theta1','j1',
        'oo1','2n2','mu2','nu2','lambda2','l2','l2','t2']
    #-- only add minor consituents that are not on the list of major values
    minor_indices = [i for i,m in enumerate(minor) if m not in constituents]

    #-- relationship between major and minor constituent amplitude and phase
    zmin = np.zeros((n,18),dtype=np.complex64)
    zmin[:,0] = 0.263*z8[:,0] - 0.0252*z8[:,1]#-- 2Q1
    zmin[:,1] = 0.297*z8[:,0] - 0.0264*z8[:,1]#-- sigma1
    zmin[:,2] = 0.164*z8[:,0] + 0.0048*z8[:,1]#-- rho1
    zmin[:,3] = 0.0140*z8[:,1] + 0.0101*z8[:,3]#-- M1
    zmin[:,4] = 0.0389*z8[:,1] + 0.0282*z8[:,3]#-- M1
    zmin[:,5] = 0.0064*z8[:,1] + 0.0060*z8[:,3]#-- chi1
    zmin[:,6] = 0.0030*z8[:,1] + 0.0171*z8[:,3]#-- pi1
    zmin[:,7] = -0.0015*z8[:,1] + 0.0152*z8[:,3]#-- phi1
    zmin[:,8] = -0.0065*z8[:,1] + 0.0155*z8[:,3]#-- theta1
    zmin[:,9] = -0.0389*z8[:,1] + 0.0836*z8[:,3]#-- J1
    zmin[:,10] = -0.0431*z8[:,1] + 0.0613*z8[:,3]#-- OO1
    zmin[:,11] = 0.264*z8[:,4] - 0.0253*z8[:,5]#-- 2N2
    zmin[:,12] = 0.298*z8[:,4] - 0.0264*z8[:,5]#-- mu2
    zmin[:,13] = 0.165*z8[:,4] + 0.00487*z8[:,5]#-- nu2
    zmin[:,14] = 0.0040*z8[:,5] + 0.0074*z8[:,6]#-- lambda2
    zmin[:,15] = 0.0131*z8[:,5] + 0.0326*z8[:,6]#-- L2
    zmin[:,16] = 0.0033*z8[:,5] + 0.0082*z8[:,6]#-- L2
    zmin[:,17] = 0.0585*z8[:,6]#-- t2

    hour = (time % 1)*24.0
    t1 = 15.0*hour
    t2 = 30.0*hour
    #-- set function for astrological longitudes
    ASTRO5 = True if CORRECTIONS in ('GOT','FES') else False
    #-- convert from Modified Julian Dates into Ephemeris Time
    S,H,P,omega,pp = calc_astrol_longitudes(time_mjd+DELTAT, ASTRO5=ASTRO5)

    #-- determine equilibrium tidal arguments
    arg = np.zeros((n,18))
    arg[:,0] = t1 - 4.0*S + H + 2.0*P - 90.0#-- 2Q1
    arg[:,1] = t1 - 4.0*S + 3.0*H - 90.0#-- sigma1
    arg[:,2] = t1 - 3.0*S + 3.0*H - P - 90.0#-- rho1
    arg[:,3] = t1 - S + H - P + 90.0#-- M1
    arg[:,4] = t1 - S + H + P + 90.0#-- M1
    arg[:,5] = t1 - S + 3.0*H - P + 90.0#-- chi1
    arg[:,6] = t1 - 2.0*H + pp - 90.0#-- pi1
    arg[:,7] = t1 + 3.0*H + 90.0#-- phi1
    arg[:,8] = t1 + S - H + P + 90.0#-- theta1
    arg[:,9] = t1 + S + H - P + 90.0#-- J1
    arg[:,10] = t1 + 2.0*S + H + 90.0#-- OO1
    arg[:,11] = t2 - 4.0*S + 2.0*H + 2.0*P#-- 2N2
    arg[:,12] = t2 - 4.0*S + 4.0*H#-- mu2
    arg[:,13] = t2 - 3.0*S + 4.0*H - P#-- nu2
    arg[:,14] = t2 - S + P + 180.0#-- lambda2
    arg[:,15] = t2 - S + 2.0*H - P + 180.0#-- L2
    arg[:,16] = t2 - S + 2.0*H + P#-- L2
    arg[:,17] = t2 - H + pp#-- t2

    #-- determine nodal corrections f and u
    sinn = np.sin(omega*dtr)
    cosn = np.cos(omega*dtr)
    sin2n = np.sin(2.0*omega*dtr)
    cos2n = np.cos(2.0*omega*dtr)

    f = np.ones((n,18))
    f[:,0] = np.sqrt((1.0 + 0.189*cosn - 0.0058*cos2n)**2 +
        (0.189*sinn - 0.0058*sin2n)**2)#-- 2Q1
    f[:,1] = f[:,0]#-- sigma1
    f[:,2] = f[:,0]#-- rho1
    f[:,3] = np.sqrt((1.0 + 0.185*cosn)**2 + (0.185*sinn)**2)#-- M1
    f[:,4] = np.sqrt((1.0 + 0.201*cosn)**2 + (0.201*sinn)**2)#-- M1
    f[:,5] = np.sqrt((1.0 + 0.221*cosn)**2 + (0.221*sinn)**2)#-- chi1
    f[:,9] = np.sqrt((1.0 + 0.198*cosn)**2 + (0.198*sinn)**2)#-- J1
    f[:,10] = np.sqrt((1.0 + 0.640*cosn + 0.134*cos2n)**2 +
        (0.640*sinn + 0.134*sin2n)**2)#-- OO1
    f[:,11] = np.sqrt((1.0 - 0.0373*cosn)**2 + (0.0373*sinn)**2)#-- 2N2
    f[:,12] = f[:,11]#-- mu2
    f[:,13] = f[:,11]#-- nu2
    f[:,15] = f[:,11]#-- L2
    f[:,16] = np.sqrt((1.0 + 0.441*cosn)**2 + (0.441*sinn)**2)#-- L2

    u = np.zeros((n,18))
    u[:,0] = np.arctan2(0.189*sinn - 0.0058*sin2n,
        1.0 + 0.189*cosn - 0.0058*sin2n)/dtr#-- 2Q1
    u[:,1] = u[:,0]#-- sigma1
    u[:,2] = u[:,0]#-- rho1
    u[:,3] = np.arctan2( 0.185*sinn, 1.0 + 0.185*cosn)/dtr#-- M1
    u[:,4] = np.arctan2(-0.201*sinn, 1.0 + 0.201*cosn)/dtr#-- M1
    u[:,5] = np.arctan2(-0.221*sinn, 1.0 + 0.221*cosn)/dtr#-- chi1
    u[:,9] = np.arctan2(-0.198*sinn, 1.0 + 0.198*cosn)/dtr#-- J1
    u[:,10] = np.arctan2(-0.640*sinn - 0.134*sin2n,
        1.0 + 0.640*cosn + 0.134*cos2n)/dtr#-- OO1
    u[:,11] = np.arctan2(-0.0373*sinn, 1.0 - 0.0373*cosn)/dtr#-- 2N2
    u[:,12] = u[:,11]#-- mu2
    u[:,13] = u[:,11]#-- nu2
    u[:,15] = u[:,11]#-- L2
    u[:,16] = np.arctan2(-0.441*sinn, 1.0 + 0.441*cosn)/dtr#-- L2

    #-- sum over the minor tidal constituents of interest
    for k in minor_indices:
        th = (arg[:,k] + u[:,k])*dtr
        dh += zmin.real[:,k]*f[:,k]*np.cos(th)-zmin.imag[:,k]*f[:,k]*np.sin(th)
    #-- return the inferred elevation
    return dh
