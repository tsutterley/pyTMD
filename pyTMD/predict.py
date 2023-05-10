#!/usr/bin/env python
u"""
predict.py
Written by Tyler Sutterley (04/2023)
Prediction routines for ocean, load, equilibrium and solid earth tides

REFERENCES:
    G. D. Egbert and S. Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    arguments.py: loads nodal corrections for tidal constituents
    astro.py: computes the basic astronomical mean longitudes
    constants.py: calculate reference parameters for common ellipsoids
    load_constituent.py: loads parameters for a given tidal constituent
    spatial.py: utilities for working with geospatial data

UPDATE HISTORY:
    Updated 04/2023: using renamed astro mean_longitudes function
        using renamed arguments function for nodal corrections
        adding prediction routine for solid earth tides
        output solid earth tide corrections as combined XYZ components
    Updated 03/2023: add basic variable typing to function inputs
    Updated 12/2022: merged prediction functions into a single module
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
from __future__ import annotations

import numpy as np
import pyTMD.astro
from pyTMD.arguments import arguments
from pyTMD.constants import constants
from pyTMD.load_constituent import load_constituent

# PURPOSE: Predict tides at single times
def map(t: float | np.ndarray,
        hc: np.ndarray,
        constituents: list | np.ndarray,
        deltat: float | np.ndarray = 0.0,
        corrections: str = 'OTIS'
    ):
    """
    Predict tides at a single time using harmonic constants [1]_

    Parameters
    ----------
    t: float or np.ndarray
        days relative to 1992-01-01T00:00:00
    hc: np.ndarray
        harmonic constant vector
    constituents: list or np.ndarray
        tidal constituent IDs
    deltat: float or np.ndarray, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS/ATLAS or GOT/FES models

    Returns
    -------
    ht: np.ndarray
        tide values reconstructed using the nodal corrections

    References
    ----------
    .. [1] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of
        Barotropic Ocean Tides," *Journal of Atmospheric and Oceanic
        Technology*, 19(2), 183--204, (2002).
        `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

    .. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
    """
    # number of points and number of constituents
    npts,nc = np.shape(hc)
    # load the nodal corrections
    # convert time to Modified Julian Days (MJD)
    pu,pf,G = arguments(t + 48622.0, constituents,
        deltat=deltat, corrections=corrections)
    # allocate for output tidal elevation
    ht = np.ma.zeros((npts))
    ht.mask = np.zeros((npts),dtype=bool)
    # for each constituent
    for k,c in enumerate(constituents):
        if corrections in ('OTIS','ATLAS','ESR','netcdf'):
            # load parameters for each constituent
            amp, ph, omega, alpha, species = load_constituent(c)
            # add component for constituent to output tidal elevation
            th = omega*t*86400.0 + ph + pu[0,k]
        elif corrections in ('GOT','FES'):
            th = G[0,k]*np.pi/180.0 + pu[0,k]
        # sum over all tides
        ht.data[:] += pf[0,k]*hc.real[:,k]*np.cos(th) - \
            pf[0,k]*hc.imag[:,k]*np.sin(th)
        ht.mask[:] |= (hc.real.mask[:,k] | hc.imag.mask[:,k])
    # return the tidal elevation after removing singleton dimensions
    return np.squeeze(ht)

# PURPOSE: Predict tides at drift bouys or altimetry points
def drift(t: float | np.ndarray,
        hc: np.ndarray,
        constituents: list | np.ndarray,
        deltat: float | np.ndarray = 0.0,
        corrections: str = 'OTIS'
    ):
    """
    Predict tides at multiple times and locations using harmonic
    constants [1]_

    Parameters
    ----------
    t: float or np.ndarray
        days relative to 1992-01-01T00:00:00
    hc: np.ndarray
        harmonic constant vector
    constituents: list or np.ndarray
        tidal constituent IDs
    deltat: float or np.ndarray, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS/ATLAS or GOT/FES models

    Returns
    -------
    ht: np.ndarray
        tidal time series reconstructed using the nodal corrections

    References
    ----------
    .. [1] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of
        Barotropic Ocean Tides," *Journal of Atmospheric and Oceanic
        Technology*, 19(2), 183--204, (2002).
        `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

    .. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
    """
    nt = len(t)
    # load the nodal corrections
    # convert time to Modified Julian Days (MJD)
    pu,pf,G = arguments(t + 48622.0, constituents,
        deltat=deltat, corrections=corrections)
    # allocate for output time series
    ht = np.ma.zeros((nt))
    ht.mask = np.zeros((nt),dtype=bool)
    # for each constituent
    for k,c in enumerate(constituents):
        if corrections in ('OTIS','ATLAS','ESR','netcdf'):
            # load parameters for each constituent
            amp, ph, omega, alpha, species = load_constituent(c)
            # add component for constituent to output tidal elevation
            th = omega*t*86400.0 + ph + pu[:,k]
        elif corrections in ('GOT','FES'):
            th = G[:,k]*np.pi/180.0 + pu[:,k]
        # sum over all tides
        ht.data[:] += pf[:,k]*hc.real[:,k]*np.cos(th) - \
            pf[:,k]*hc.imag[:,k]*np.sin(th)
        ht.mask[:] |= (hc.real.mask[:,k] | hc.imag.mask[:,k])
    # return tides
    return ht

# PURPOSE: Predict a tidal time series at a location
def time_series(t: float | np.ndarray,
        hc: np.ndarray,
        constituents: list | np.ndarray,
        deltat: float | np.ndarray = 0.0,
        corrections: str = 'OTIS'
    ):
    """
    Predict tidal time series at a single location using harmonic
    constants [1]_

    Parameters
    ----------
    t: float or np.ndarray
        days relative to 1992-01-01T00:00:00
    hc: np.ndarray
        harmonic constant vector
    constituents: list or np.ndarray
        tidal constituent IDs
    deltat: float or np.ndarray, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS/ATLAS or GOT/FES models

    Returns
    -------
    ht: np.ndarray
        tidal time series reconstructed using the nodal corrections

    References
    ----------
    .. [1] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of
        Barotropic Ocean Tides," *Journal of Atmospheric and Oceanic
        Technology*, 19(2), 183--204, (2002).
        `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

    .. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
    """
    nt = len(t)
    # load the nodal corrections
    # convert time to Modified Julian Days (MJD)
    pu,pf,G = arguments(t + 48622.0, constituents,
        deltat=deltat, corrections=corrections)
    # allocate for output time series
    ht = np.ma.zeros((nt))
    ht.mask = np.zeros((nt),dtype=bool)
    # for each constituent
    for k,c in enumerate(constituents):
        if corrections in ('OTIS','ATLAS','ESR','netcdf'):
            # load parameters for each constituent
            amp, ph, omega, alpha, species = load_constituent(c)
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

# PURPOSE: infer the minor corrections from the major constituents
def infer_minor(
        t: float | np.ndarray,
        zmajor: np.ndarray,
        constituents: list | np.ndarray,
        **kwargs
    ):
    """
    Calculate the tidal corrections for minor constituents inferred using
    major constituents [1]_ [2]_ [3]_ [4]_

    Parameters
    ----------
    t: float or np.ndarray
        days relative to 1992-01-01T00:00:00
    zmajor: np.ndarray
        Complex HC for given constituents/points
    constituents: list
        tidal constituent IDs
    deltat: float or np.ndarray, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS/ATLAS or GOT/FES models

    Returns
    -------
    dh: np.ndarray
        Height from minor constituents

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
    # set default keyword arguments
    kwargs.setdefault('deltat', 0.0)
    kwargs.setdefault('corrections', 'OTIS')

    # degrees to radians
    dtr = np.pi/180.0
    # number of constituents
    npts, nc = np.shape(zmajor)
    nt = len(np.atleast_1d(t))
    # number of data points to calculate if running time series/drift/map
    n = nt if ((npts == 1) & (nt > 1)) else npts
    # allocate for output elevation correction
    dh = np.ma.zeros((n))
    # convert time from days relative to Jan 1, 1992 to Modified Julian Days
    MJD = 48622.0 + t
    # major constituents used for inferring minor tides
    cindex = ['q1', 'o1', 'p1', 'k1', 'n2', 'm2', 's2', 'k2', '2n2']
    # re-order major tides to correspond to order of cindex
    z = np.ma.zeros((n,9),dtype=np.complex64)
    nz = 0
    for i,c in enumerate(cindex):
        j = [j for j,val in enumerate(constituents) if (val == c)]
        if j:
            j1, = j
            z[:,i] = zmajor[:,j1]
            nz += 1

    if (nz < 6):
        raise Exception('Not enough constituents for inference')

    # list of minor constituents
    minor = ['2q1','sigma1','rho1','m12','m11','chi1','pi1','phi1','theta1',
        'j1','oo1','2n2','mu2','nu2','lambda2','l2','l2','t2','eps2','eta2']
    # only add minor constituents that are not on the list of major values
    minor_indices = [i for i,m in enumerate(minor) if m not in constituents]

    # relationship between major and minor constituent amplitude and phase
    zmin = np.zeros((n,20),dtype=np.complex64)
    zmin[:,0] = 0.263*z[:,0] - 0.0252*z[:,1]# 2Q1
    zmin[:,1] = 0.297*z[:,0] - 0.0264*z[:,1]# sigma1
    zmin[:,2] = 0.164*z[:,0] + 0.0048*z[:,1]# rho1
    zmin[:,3] = 0.0140*z[:,1] + 0.0101*z[:,3]# M12
    zmin[:,4] = 0.0389*z[:,1] + 0.0282*z[:,3]# M11
    zmin[:,5] = 0.0064*z[:,1] + 0.0060*z[:,3]# chi1
    zmin[:,6] = 0.0030*z[:,1] + 0.0171*z[:,3]# pi1
    zmin[:,7] = -0.0015*z[:,1] + 0.0152*z[:,3]# phi1
    zmin[:,8] = -0.0065*z[:,1] + 0.0155*z[:,3]# theta1
    zmin[:,9] = -0.0389*z[:,1] + 0.0836*z[:,3]# J1
    zmin[:,10] = -0.0431*z[:,1] + 0.0613*z[:,3]# OO1
    zmin[:,11] = 0.264*z[:,4] - 0.0253*z[:,5]# 2N2
    zmin[:,12] = 0.298*z[:,4] - 0.0264*z[:,5]# mu2
    zmin[:,13] = 0.165*z[:,4] + 0.00487*z[:,5]# nu2
    zmin[:,14] = 0.0040*z[:,5] + 0.0074*z[:,6]# lambda2
    zmin[:,15] = 0.0131*z[:,5] + 0.0326*z[:,6]# L2
    zmin[:,16] = 0.0033*z[:,5] + 0.0082*z[:,6]# L2
    zmin[:,17] = 0.0585*z[:,6]# t2
    # additional coefficients for FES models
    if kwargs['corrections'] in ('FES',):
        # spline coefficients for admittances
        mu2 = [0.069439968323, 0.351535557706, -0.046278307672]
        nu2 = [-0.006104695053, 0.156878802427, 0.006755704028]
        l2 = [0.077137765667, -0.051653455134, 0.027869916824]
        t2 = [0.180480173707, -0.020101177502, 0.008331518844]
        lda2 = [0.016503557465, -0.013307812292, 0.007753383202]
        zmin[:,12] = mu2[0]*z[:,7] + mu2[1]*z[:,4] + mu2[2]*z[:,5]# mu2
        zmin[:,13] = nu2[0]*z[:,7] + nu2[1]*z[:,4] + nu2[2]*z[:,5]# nu2
        zmin[:,14] = lda2[0]*z[:,7] + lda2[1]*z[:,4] + lda2[2]*z[:,5]# lambda2
        zmin[:,16] = l2[0]*z[:,7] + l2[1]*z[:,4] + l2[2]*z[:,5]# L2
        zmin[:,17] = t2[0]*z[:,7] + t2[1]*z[:,4] + t2[2]*z[:,5]# t2
        zmin[:,18] = 0.53285*z[:,8] - 0.03304*z[:,4]# eps2
        zmin[:,19] = -0.0034925*z[:,5] + 0.0831707*z[:,7]# eta2

    hour = (t % 1)*24.0
    t1 = 15.0*hour
    t2 = 30.0*hour
    # set function for astronomical longitudes
    ASTRO5 = True if kwargs['corrections'] in ('GOT','FES') else False
    # convert from Modified Julian Dates into Ephemeris Time
    S,H,P,omega,pp = pyTMD.astro.mean_longitudes(MJD + kwargs['deltat'],
        ASTRO5=ASTRO5)

    # determine equilibrium tidal arguments
    arg = np.zeros((n,20))
    arg[:,0] = t1 - 4.0*S + H + 2.0*P - 90.0# 2Q1
    arg[:,1] = t1 - 4.0*S + 3.0*H - 90.0# sigma1
    arg[:,2] = t1 - 3.0*S + 3.0*H - P - 90.0# rho1
    arg[:,3] = t1 - S + H - P + 90.0# M12
    arg[:,4] = t1 - S + H + P + 90.0# M11
    arg[:,5] = t1 - S + 3.0*H - P + 90.0# chi1
    arg[:,6] = t1 - 2.0*H + pp - 90.0# pi1
    arg[:,7] = t1 + 3.0*H + 90.0# phi1
    arg[:,8] = t1 + S - H + P + 90.0# theta1
    arg[:,9] = t1 + S + H - P + 90.0# J1
    arg[:,10] = t1 + 2.0*S + H + 90.0# OO1
    arg[:,11] = t2 - 4.0*S + 2.0*H + 2.0*P# 2N2
    arg[:,12] = t2 - 4.0*S + 4.0*H# mu2
    arg[:,13] = t2 - 3.0*S + 4.0*H - P# nu2
    arg[:,14] = t2 - S + P + 180.0# lambda2
    arg[:,15] = t2 - S + 2.0*H - P + 180.0# L2
    arg[:,16] = t2 - S + 2.0*H + P# L2
    arg[:,17] = t2 - H + pp# t2
    arg[:,18] = t2 - 5.0*S + 4.0*H + P # eps2
    arg[:,19] = t2 + S + 2.0*H - pp # eta2

    # determine nodal corrections f and u
    sinn = np.sin(omega*dtr)
    cosn = np.cos(omega*dtr)
    sin2n = np.sin(2.0*omega*dtr)
    cos2n = np.cos(2.0*omega*dtr)

    f = np.ones((n,20))
    f[:,0] = np.sqrt((1.0 + 0.189*cosn - 0.0058*cos2n)**2 +
        (0.189*sinn - 0.0058*sin2n)**2)# 2Q1
    f[:,1] = f[:,0]# sigma1
    f[:,2] = f[:,0]# rho1
    f[:,3] = np.sqrt((1.0 + 0.185*cosn)**2 + (0.185*sinn)**2)# M12
    f[:,4] = np.sqrt((1.0 + 0.201*cosn)**2 + (0.201*sinn)**2)# M11
    f[:,5] = np.sqrt((1.0 + 0.221*cosn)**2 + (0.221*sinn)**2)# chi1
    f[:,9] = np.sqrt((1.0 + 0.198*cosn)**2 + (0.198*sinn)**2)# J1
    f[:,10] = np.sqrt((1.0 + 0.640*cosn + 0.134*cos2n)**2 +
        (0.640*sinn + 0.134*sin2n)**2)# OO1
    f[:,11] = np.sqrt((1.0 - 0.0373*cosn)**2 + (0.0373*sinn)**2)# 2N2
    f[:,12] = f[:,11]# mu2
    f[:,13] = f[:,11]# nu2
    f[:,15] = f[:,11]# L2
    f[:,16] = np.sqrt((1.0 + 0.441*cosn)**2 + (0.441*sinn)**2)# L2

    u = np.zeros((n,20))
    u[:,0] = np.arctan2(0.189*sinn - 0.0058*sin2n,
        1.0 + 0.189*cosn - 0.0058*sin2n)/dtr# 2Q1
    u[:,1] = u[:,0]# sigma1
    u[:,2] = u[:,0]# rho1
    u[:,3] = np.arctan2( 0.185*sinn, 1.0 + 0.185*cosn)/dtr# M12
    u[:,4] = np.arctan2(-0.201*sinn, 1.0 + 0.201*cosn)/dtr# M11
    u[:,5] = np.arctan2(-0.221*sinn, 1.0 + 0.221*cosn)/dtr# chi1
    u[:,9] = np.arctan2(-0.198*sinn, 1.0 + 0.198*cosn)/dtr# J1
    u[:,10] = np.arctan2(-0.640*sinn - 0.134*sin2n,
        1.0 + 0.640*cosn + 0.134*cos2n)/dtr# OO1
    u[:,11] = np.arctan2(-0.0373*sinn, 1.0 - 0.0373*cosn)/dtr# 2N2
    u[:,12] = u[:,11]# mu2
    u[:,13] = u[:,11]# nu2
    u[:,15] = u[:,11]# L2
    u[:,16] = np.arctan2(-0.441*sinn, 1.0 + 0.441*cosn)/dtr# L2

    if kwargs['corrections'] in ('FES',):
        # additional astronomical terms for FES models
        II = np.arccos(0.913694997 - 0.035692561*np.cos(omega*dtr))
        at1 = np.arctan(1.01883*np.tan(omega*dtr/2.0))
        at2 = np.arctan(0.64412*np.tan(omega*dtr/2.0))
        xi = -at1 - at2 + omega*dtr
        xi[xi > np.pi] -= 2.0*np.pi
        nu = at1 - at2
        I2 = np.tan(II/2.0)
        Ra1 = np.sqrt(1.0 - 12.0*(I2**2)*np.cos(2.0*(P - xi)) + 36.0*(I2**4))
        P2 = np.sin(2.0*(P - xi))
        Q2 = 1.0/(6.0*(I2**2)) - np.cos(2.0*(P - xi))
        R = np.arctan(P2/Q2)

        f[:,0] = np.sin(II)*(np.cos(II/2.0)**2)/0.38 # 2Q1
        f[:,1] = f[:,0] # sigma1
        f[:,2] = f[:,0] # rho1
        f[:,3] = f[:,0] # M12
        f[:,4] = np.sin(2.0*II)/0.7214 # M11
        f[:,5] = f[:,4] # chi1
        f[:,9] = f[:,5] # J1
        f[:,10] = np.sin(II)*np.power(np.sin(II/2.0),2.0)/0.01640 # OO1
        f[:,11] = np.power(np.cos(II/2.0),4.0)/0.9154 # 2N2
        f[:,12] = f[:,11] # mu2
        f[:,13] = f[:,11] # nu2
        f[:,14] = f[:,11] # lambda2
        f[:,15] = f[:,11]*Ra1 # L2
        f[:,18] = f[:,11] # eps2
        f[:,19] = np.power(np.sin(II),2.0)/0.1565 # eta2

        u[:,0] = (2.0*xi - nu)/dtr # 2Q1
        u[:,1] = u[:,0] # sigma1
        u[:,2] = u[:,0] # rho1
        u[:,3] = u[:,0] # M12
        u[:,4] = -nu/dtr # M11
        u[:,5] = u[:,4] # chi1
        u[:,9] = u[:,4] # J1
        u[:,10] = (-2.0*xi - nu)/dtr # OO1
        u[:,11] = (2.0*xi - 2.0*nu)/dtr # 2N2
        u[:,12] = u[:,11] # mu2
        u[:,13] = u[:,11] # nu2
        u[:,14] = (2.0*xi - 2.0*nu)/dtr # lambda2
        u[:,15] = (2.0*xi - 2.0*nu - R)/dtr# L2
        u[:,18] = u[:,12] # eps2
        u[:,19] = -2.0*nu/dtr # eta2

    # sum over the minor tidal constituents of interest
    for k in minor_indices:
        th = (arg[:,k] + u[:,k])*dtr
        dh += zmin.real[:,k]*f[:,k]*np.cos(th)-zmin.imag[:,k]*f[:,k]*np.sin(th)
    # return the inferred elevation
    return dh

# PURPOSE: estimate long-period equilibrium tides
def equilibrium_tide(t: np.ndarray, lat: np.ndarray):
    """
    Compute the long-period equilibrium tides the summation of fifteen
    tidal spectral lines from Cartwright-Tayler-Edden tables [1]_ [2]_

    Parameters
    ----------
    t: np.ndarray
        time (days relative to January 1, 1992)
    lat: np.ndarray
        latitudes (degrees north)

    Returns
    -------
    lpet: np.ndarray
        long-period equilibrium tide in meters

    References
    ----------
    .. [1] D. E. Cartwright and R. J. Tayler,
        "New Computations of the Tide-generating Potential,"
        *Geophysical Journal of the Royal Astronomical Society*,
        23(1), 45--73. (1971). `doi: 10.1111/j.1365-246X.1971.tb01803.x
        <https://doi.org/10.1111/j.1365-246X.1971.tb01803.x>`_
    .. [2] D. E. Cartwright and A. C. Edden,
        "Corrected Tables of Tidal Harmonics,"
        *Geophysical Journal of the Royal Astronomical Society*,
        33(3), 253--264, (1973). `doi: 10.1111/j.1365-246X.1973.tb03420.x
        <https://doi.org/10.1111/j.1365-246X.1973.tb03420.x>`_
    """
    # longitude of moon
    # longitude of sun
    # longitude of lunar perigee
    # longitude of ascending lunar node
    PHC = np.array([290.21,280.12,274.35,343.51])
    DPD = np.array([13.1763965,0.9856473,0.1114041,0.0529539])

    # number of input points
    nt = len(np.atleast_1d(t))
    nlat = len(np.atleast_1d(lat))
    # compute 4 principal mean longitudes in radians at delta time (SHPN)
    SHPN = np.zeros((4,nt))
    for N in range(4):
        # convert time from days relative to 1992-01-01 to 1987-01-01
        ANGLE = PHC[N] + (t + 1826.0)*DPD[N]
        SHPN[N,:] = np.pi*np.mod(ANGLE, 360.0)/180.0

    # assemble long-period tide potential from 15 CTE terms greater than 1 mm
    # nodal term is included but not the constant term.
    PH = np.zeros((nt))
    Z = np.zeros((nt))
    Z += 2.79*np.cos(SHPN[3,:]) - 0.49*np.cos(SHPN[1,:] - \
        283.0*np.pi/180.0) - 3.10*np.cos(2.0*SHPN[1,:])
    PH += SHPN[0,:]
    Z += -0.67*np.cos(PH - 2.0*SHPN[1,:] + SHPN[2,:]) - \
        (3.52 - 0.46*np.cos(SHPN[3,:]))*np.cos(PH - SHPN[2,:])
    PH += SHPN[0,:]
    Z += - 6.66*np.cos(PH) - 2.76*np.cos(PH + SHPN[3,:]) - \
        0.26 * np.cos(PH + 2.*SHPN[3,:]) - 0.58 * np.cos(PH - 2.*SHPN[1,:]) - \
        0.29 * np.cos(PH - 2.*SHPN[2,:])
    PH += SHPN[0,:]
    Z += - 1.27*np.cos(PH - SHPN[2,:]) - \
        0.53*np.cos(PH - SHPN[2,:] + SHPN[3,:]) - \
        0.24*np.cos(PH - 2.0*SHPN[1,:] + SHPN[2,:])

    # Multiply by gamma_2 * normalization * P20(lat)
    gamma_2 = 0.693
    P20 = 0.5*(3.0*np.sin(lat*np.pi/180.0)**2 - 1.0)
    # calculate long-period equilibrium tide and convert to meters
    if (nlat != nt):
        lpet = gamma_2*np.sqrt((4.0+1.0)/(4.0*np.pi))*np.outer(P20,Z/100.0)
    else:
        lpet = gamma_2*np.sqrt((4.0+1.0)/(4.0*np.pi))*P20*(Z/100.0)
    # return the long-period equilibrium tides
    return lpet

# get IERS parameters
_iers = constants(ellipsoid='IERS', units='MKS')

# PURPOSE: estimate solid Earth tides due to gravitational attraction
def solid_earth_tide(
        t: np.ndarray,
        XYZ: np.ndarray,
        SXYZ: np.ndarray,
        LXYZ: np.ndarray,
        a_axis: float = _iers.a_axis,
        tide_system: str = 'tide_free',
        **kwargs
    ):
    """
    Compute the solid Earth tides due to the gravitational attraction
    of the moon and sun [1]_ [2]_ [3]_ [4]_

    Parameters
    ----------
    t: np.ndarray
        Time (days relative to January 1, 1992)
    XYZ: np.ndarray
        Cartesian coordinates of the station (meters)
    SXYZ: np.ndarray
        Earth-centered Earth-fixed coordinates of the sun (meters)
    LXYZ: np.ndarray
        Earth-centered Earth-fixed coordinates of the moon (meters)
    a_axis: float, default 6378136.3
        Semi-major axis of the Earth (meters)
    tide_system: str, default 'tide_free'
        Permanent tide system for the output solid Earth tide

            - ``'tide_free'``: no permanent direct and indirect tidal potentials
            - ``'mean_tide'``: permanent tidal potentials (direct and indirect)

    Returns
    -------
    dxt: np.ndarray
        Solid Earth tide in meters in Cartesian coordinates

    References
    ----------
    .. [1] P. M. Mathews, B. A. Buffett, T. A. Herring and I. I Shapiro,
        "Forced nutations of the Earth: Influence of inner core dynamics:
        1. Theory", *Journal of Geophysical Research: Solid Earth*,
        96(B5), 8219--8242, (1991). `doi: 10.1029/90JB01955
        <https://doi.org/10.1029/90JB01955>`_
    .. [2] P. M. Mathews, V. Dehant and J. M. Gipson,
        "Tidal station displacements", *Journal of Geophysical
        Research: Solid Earth*, 102(B9), 20469--20477, (1997).
        `doi: 10.1029/97JB01515 <https://doi.org/10.1029/97JB01515>`_
    .. [3] J. C. Ries, R. J. Eanes, C. K. Shum and M. M. Watkins,
        "Progress in the determination of the gravitational
        coefficient of the Earth", *Geophysical Research Letters*,
        19(6), 529--531, (1992). `doi: 10.1029/92GL00259
        <https://doi.org/10.1029/92GL00259>`_
    .. [4] J. M. Wahr, "Body tides on an elliptical, rotating, elastic
        and oceanless Earth", *Geophysical Journal of the Royal
        Astronomical Society*, 64(3), 677--703, (1981).
        `doi: 10.1111/j.1365-246X.1981.tb02690.x
        <https://doi.org/10.1111/j.1365-246X.1981.tb02690.x>`_
    """
    # set default keyword arguments
    # nominal Love and Shida numbers
    kwargs.setdefault('h2', 0.6078)
    kwargs.setdefault('l2', 0.0847)
    kwargs.setdefault('h3', 0.292)
    kwargs.setdefault('l3', 0.015)
    # mass ratios between earth and sun/moon
    kwargs.setdefault('mass_ratio_solar', 332946.0482)
    kwargs.setdefault('mass_ratio_lunar', 0.0123000371)
    # validate output tide system
    assert tide_system in ('tide_free', 'mean_tide')
    # number of input coordinates
    nt = len(np.atleast_1d(t))
    # convert time to Modified Julian Days (MJD)
    MJD = t + 48622.0
    # scalar product of input coordinates with sun/moon vectors
    radius = np.sqrt(XYZ[:,0]**2 + XYZ[:,1]**2 + XYZ[:,2]**2)
    solar_radius = np.sqrt(SXYZ[:,0]**2 + SXYZ[:,1]**2 + SXYZ[:,2]**2)
    lunar_radius = np.sqrt(LXYZ[:,0]**2 + LXYZ[:,1]**2 + LXYZ[:,2]**2)
    solar_scalar = (XYZ[:,0]*SXYZ[:,0] + XYZ[:,1]*SXYZ[:,1] +
        XYZ[:,2]*SXYZ[:,2])/(radius*solar_radius)
    lunar_scalar = (XYZ[:,0]*LXYZ[:,0] + XYZ[:,1]*LXYZ[:,1] +
        XYZ[:,2]*LXYZ[:,2])/(radius*lunar_radius)
    # compute new h2 and l2 (Mathews et al., 1997)
    cosphi = np.sqrt(XYZ[:,0]**2 + XYZ[:,1]**2)/radius
    h2 = kwargs['h2'] - 0.0006*(1.0 - 3.0/2.0*cosphi**2)
    l2 = kwargs['l2'] + 0.0002*(1.0 - 3.0/2.0*cosphi**2)
    # compute P2 terms
    P2_solar = 3.0*(h2/2.0 - l2)*solar_scalar**2 - h2/2.0
    P2_lunar = 3.0*(h2/2.0 - l2)*lunar_scalar**2 - h2/2.0
    # compute P3 terms
    P3_solar = 5.0/2.0*(kwargs['h3'] - 3.0*kwargs['l3'])*solar_scalar**3 + \
        3.0/2.0*(kwargs['l3'] - kwargs['h3'])*solar_scalar
    P3_lunar = 5.0/2.0*(kwargs['h3'] - 3.0*kwargs['l3'])*lunar_scalar**3 + \
        3.0/2.0*(kwargs['l3'] - kwargs['h3'])*lunar_scalar
    # compute terms in direction of sun/moon vectors
    X2_solar = 3.0*l2*solar_scalar
    X2_lunar = 3.0*l2*lunar_scalar
    X3_solar = 3.0*kwargs['l3']/2.0*(5.0*solar_scalar**2 - 1.0)
    X3_lunar = 3.0*kwargs['l3']/2.0*(5.0*lunar_scalar**2 - 1.0)
    # factors for sun and moon using IAU estimates of mass ratios
    F2_solar = kwargs['mass_ratio_solar']*a_axis*(a_axis/solar_radius)**3
    F2_lunar = kwargs['mass_ratio_lunar']*a_axis*(a_axis/lunar_radius)**3
    F3_solar = kwargs['mass_ratio_solar']*a_axis*(a_axis/solar_radius)**4
    F3_lunar = kwargs['mass_ratio_lunar']*a_axis*(a_axis/lunar_radius)**4
    # compute total displacement (Mathews et al. 1997)
    dxt = np.zeros((nt, 3))
    for i in range(3):
        S2 = F2_solar*(X2_solar*SXYZ[:,i]/solar_radius+P2_solar*XYZ[:,i]/radius)
        L2 = F2_lunar*(X2_lunar*LXYZ[:,i]/lunar_radius+P2_lunar*XYZ[:,i]/radius)
        S3 = F3_solar*(X3_solar*SXYZ[:,i]/solar_radius+P3_solar*XYZ[:,i]/radius)
        L3 = F3_lunar*(X3_lunar*LXYZ[:,i]/lunar_radius+P3_lunar*XYZ[:,i]/radius)
        dxt[:,i] = S2 + L2 + S3 + L3
    # corrections for out-of-phase portions of the Love and Shida numbers
    dxt += _out_of_phase_diurnal(XYZ, SXYZ, LXYZ, F2_solar, F2_lunar)
    dxt += _out_of_phase_semidiurnal(XYZ, SXYZ, LXYZ, F2_solar, F2_lunar)
    # corrections for the latitudinal dependence
    dxt += _latitude_dependence(XYZ, SXYZ, LXYZ, F2_solar, F2_lunar)
    # corrections for the frequency dependence
    dxt += _frequency_dependence_diurnal(XYZ, MJD)
    dxt += _frequency_dependence_long_period(XYZ, MJD)
    # convert the permanent tide system if specified
    if (tide_system.lower() == 'mean_tide'):
        dxt += _free_to_mean(XYZ, h2, l2)
    # return the solid earth tide
    return dxt

def _out_of_phase_diurnal(
        XYZ: np.ndarray,
        SXYZ: np.ndarray,
        LXYZ: np.ndarray,
        F2_solar: np.ndarray,
        F2_lunar: np.ndarray
    ):
    """
    Computes the out-of-phase corrections induced by mantle
    anelasticity in the diurnal band

    Parameters
    ----------
    XYZ: np.ndarray
        Cartesian coordinates of the station (meters)
    SXYZ: np.ndarray
        Earth-centered Earth-fixed coordinates of the sun (meters)
    LXYZ: np.ndarray
        Earth-centered Earth-fixed coordinates of the moon (meters)
    F2_solar: np.ndarray
        Factors for the sun
    F2_lunar: np.ndarray
        Factors for the moon
    """
    # Love and Shida number corrections
    dhi = -0.0025
    dli = -0.0007
    # Compute the normalized position vector of coordinates
    radius = np.sqrt(np.sum(XYZ**2, axis=1))
    sinphi = XYZ[:,2]/radius
    cosphi = np.sqrt(XYZ[:,0]**2 + XYZ[:,1]**2)/radius
    cos2phi = cosphi**2 - sinphi**2
    sinla = XYZ[:,1]/cosphi/radius
    cosla = XYZ[:,0]/cosphi/radius
    # Compute the normalized position vector of the Sun/Moon
    solar_radius = np.sqrt(np.sum(SXYZ**2, axis=1))
    lunar_radius = np.sqrt(np.sum(LXYZ**2, axis=1))
    # calculate offsets
    dr_solar = -3.0*dhi*sinphi*cosphi*F2_solar*SXYZ[:,2]* \
        (SXYZ[:,0]*sinla-SXYZ[:,1]*cosla)/solar_radius**2
    dr_lunar = -3.0*dhi*sinphi*cosphi*F2_lunar*LXYZ[:,2]* \
        (LXYZ[:,0]*sinla-LXYZ[:,1]*cosla)/lunar_radius**2
    dn_solar = -3.0*dli*cos2phi*F2_solar*SXYZ[:,2]* \
        (SXYZ[:,0]*sinla-SXYZ[:,1]*cosla)/solar_radius**2
    dn_lunar = -3.0*dli*cos2phi*F2_lunar*LXYZ[:,2]* \
        (LXYZ[:,0]*sinla-LXYZ[:,1]*cosla)/lunar_radius**2
    de_solar = -3.0*dli*sinphi*F2_solar*SXYZ[:,2]* \
        (SXYZ[:,0]*cosla+SXYZ[:,1]*sinla)/solar_radius**2
    de_lunar = -3.0*dli*sinphi*F2_lunar*LXYZ[:,2]* \
        (LXYZ[:,0]*cosla+LXYZ[:,1]*sinla)/lunar_radius**2
    # add solar and lunar offsets
    DR = dr_solar + dr_lunar
    DN = dn_solar + dn_lunar
    DE = de_solar + de_lunar
    # compute corrections
    DX = DR*cosla*cosphi - DE*sinla - DN*cosla*sinphi
    DY = DR*sinla*cosphi + DE*cosla - DN*sinla*sinphi
    DZ = DR*sinphi + DN*cosphi
    # return the corrections
    return np.c_[DX, DY, DZ]

def _out_of_phase_semidiurnal(
        XYZ: np.ndarray,
        SXYZ: np.ndarray,
        LXYZ: np.ndarray,
        F2_solar: np.ndarray,
        F2_lunar: np.ndarray
    ):
    """
    Computes the out-of-phase corrections induced by mantle
    anelasticity in the semi-diurnal band

    Parameters
    ----------
    XYZ: np.ndarray
        Cartesian coordinates of the station (meters)
    SXYZ: np.ndarray
        Earth-centered Earth-fixed coordinates of the sun (meters)
    LXYZ: np.ndarray
        Earth-centered Earth-fixed coordinates of the moon (meters)
    F2_solar: np.ndarray
        Factors for the sun
    F2_lunar: np.ndarray
        Factors for the moon
    """
    # Love and Shida number corrections
    dhi = -0.0022
    dli = -0.0007
    # Compute the normalized position vector of coordinates
    radius = np.sqrt(np.sum(XYZ**2, axis=1))
    sinphi = XYZ[:,2]/radius
    cosphi = np.sqrt(XYZ[:,0]**2 + XYZ[:,1]**2)/radius
    sinla = XYZ[:,1]/cosphi/radius
    cosla = XYZ[:,0]/cosphi/radius
    cos2la = cosla**2 - sinla**2
    sin2la = 2.0*cosla*sinla
    # Compute the normalized position vector of the Sun/Moon
    solar_radius = np.sqrt(np.sum(SXYZ**2, axis=1))
    lunar_radius = np.sqrt(np.sum(LXYZ**2, axis=1))
    # calculate offsets
    dr_solar = -3.0/4.0*dhi*cosphi**2*F2_solar * \
        ((SXYZ[:,0]**2-SXYZ[:,1]**2)*sin2la-2.0*SXYZ[:,0]*SXYZ[:,1]*cos2la) / \
        solar_radius**2
    dr_lunar = -3.0/4.0*dhi*cosphi**2*F2_lunar * \
        ((LXYZ[:,0]**2-LXYZ[:,1]**2)*sin2la-2.0*LXYZ[:,0]*LXYZ[:,1]*cos2la) / \
        lunar_radius**2
    dn_solar = 3.0/2.0*dli*sinphi*cosphi*F2_solar * \
        ((SXYZ[:,0]**2-SXYZ[:,1]**2)*sin2la-2.0*SXYZ[:,0]*SXYZ[:,1]*cos2la) / \
        solar_radius**2
    dn_lunar = 3.0/2.0*dli*sinphi*cosphi*F2_lunar * \
        ((LXYZ[:,0]**2-LXYZ[:,1]**2)*sin2la-2.0*LXYZ[:,0]*LXYZ[:,1]*cos2la) / \
        lunar_radius**2
    de_solar = -3.0/2.0*dli*cosphi*F2_solar * \
        ((SXYZ[:,0]**2-SXYZ[:,1]**2)*cos2la+2.0*SXYZ[:,0]*SXYZ[:,1]*sin2la) / \
        solar_radius**2
    de_lunar = -3.0/2.0*dli*cosphi*F2_lunar * \
        ((LXYZ[:,0]**2-LXYZ[:,1]**2)*cos2la+2.0*LXYZ[:,0]*LXYZ[:,1]*sin2la) / \
        lunar_radius**2
    # add solar and lunar offsets
    DR = dr_solar + dr_lunar
    DN = dn_solar + dn_lunar
    DE = de_solar + de_lunar
    # compute corrections
    DX = DR*cosla*cosphi - DE*sinla - DN*cosla*sinphi
    DY = DR*sinla*cosphi + DE*cosla - DN*sinla*sinphi
    DZ = DR*sinphi + DN*cosphi
    # return the corrections
    return np.c_[DX, DY, DZ]

def _latitude_dependence(
        XYZ: np.ndarray,
        SXYZ: np.ndarray,
        LXYZ: np.ndarray,
        F2_solar: np.ndarray,
        F2_lunar: np.ndarray
    ):
    """
    Computes the corrections induced by the latitude of the
    dependence given by L\ :sup:`1`

    Parameters
    ----------
    XYZ: np.ndarray
        Cartesian coordinates of the station (meters)
    SXYZ: np.ndarray
        Earth-centered Earth-fixed coordinates of the sun (meters)
    LXYZ: np.ndarray
        Earth-centered Earth-fixed coordinates of the moon (meters)
    F2_solar: np.ndarray
        Factors for the sun
    F2_lunar: np.ndarray
        Factors for the moon
    """
    # Love/Shida number corrections (diurnal and semi-diurnal)
    l1d = 0.0012
    l1sd = 0.0024
    # Compute the normalized position vector of coordinates
    radius = np.sqrt(np.sum(XYZ**2, axis=1))
    sinphi = XYZ[:,2]/radius
    cosphi = np.sqrt(XYZ[:,0]**2 + XYZ[:,1]**2)/radius
    sinla = XYZ[:,1]/cosphi/radius
    cosla = XYZ[:,0]/cosphi/radius
    cos2la = cosla**2 - sinla**2
    sin2la = 2.0*cosla*sinla
    # Compute the normalized position vector of the Sun/Moon
    solar_radius = np.sqrt(np.sum(SXYZ**2, axis=1))
    lunar_radius = np.sqrt(np.sum(LXYZ**2, axis=1))
    # calculate offsets for the diurnal band
    dn_d_solar = -l1d*sinphi**2*F2_solar*SXYZ[:,2] * \
        (SXYZ[:,0]*cosla+SXYZ[:,1]*sinla)/solar_radius**2
    dn_d_lunar = -l1d*sinphi**2*F2_lunar*LXYZ[:,2] * \
        (LXYZ[:,0]*cosla+LXYZ[:,1]*sinla)/lunar_radius**2
    de_d_solar = l1d*sinphi*(cosphi**2-sinphi**2)*F2_solar*SXYZ[:,2] * \
        (SXYZ[:,0]*sinla-SXYZ[:,1]*cosla)/solar_radius**2
    de_d_lunar = l1d*sinphi*(cosphi**2-sinphi**2)*F2_lunar*LXYZ[:,2] * \
        (LXYZ[:,0]*sinla-LXYZ[:,1]*cosla)/lunar_radius**2
    # calculate offsets for the semi-diurnal band
    dn_s_solar = -l1sd/2.0*sinphi*cosphi*F2_solar * \
        ((SXYZ[:,0]**2-SXYZ[:,1]**2)*cos2la+2.0*SXYZ[:,0]*SXYZ[:,1]*sin2la) / \
        solar_radius**2
    dn_s_lunar =-l1sd/2.0*sinphi*cosphi*F2_lunar * \
        ((LXYZ[:,0]**2-LXYZ[:,1]**2)*cos2la+2.0*LXYZ[:,0]*LXYZ[:,1]*sin2la) / \
        lunar_radius**2
    de_s_solar =-l1sd/2.0*sinphi**2*cosphi*F2_solar * \
        ((SXYZ[:,0]**2-SXYZ[:,1]**2)*sin2la-2.0*SXYZ[:,0]*SXYZ[:,1]*cos2la) / \
        solar_radius**2
    de_s_lunar =-l1sd/2.0*sinphi**2*cosphi*F2_lunar * \
        ((LXYZ[:,0]**2-LXYZ[:,1]**2)*sin2la-2.0*LXYZ[:,0]*LXYZ[:,1]*cos2la) / \
        lunar_radius**2
    # add solar and lunar offsets (diurnal and semi-diurnal)
    DN = 3.0*(dn_d_solar + dn_d_lunar + dn_s_solar + dn_s_lunar)
    DE = 3.0*(de_d_solar + de_d_lunar + de_s_solar + de_s_lunar)
    # compute combined diurnal and semi-diurnal corrections
    DX = -DE*sinla - DN*cosla*sinphi
    DY = DE*cosla - DN*sinla*sinphi
    DZ = DN*cosphi
    # return the corrections
    return np.c_[DX, DY, DZ]

def _frequency_dependence_diurnal(
        XYZ: np.ndarray,
        MJD: np.ndarray
    ):
    """
    Computes the in-phase and out-of-phase corrections induced by mantle
    anelasticity in the diurnal band

    Parameters
    ----------
    XYZ: np.ndarray
        Cartesian coordinates of the station (meters)
    MJD: np.ndarray
        Modified Julian Day (MJD)
    """
    # number of time steps
    nt = len(np.atleast_1d(MJD))
    # Corrections to Diurnal Tides for Frequency Dependence
    # of Love and Shida Number Parameters
    # table 7.3a of IERS conventions
    table = np.array([
        [-3.0, 0.0, 2.0, 0.0, 0.0, -0.01, 0.0, 0.0, 0.0],
        [-3.0, 2.0, 0.0, 0.0, 0.0, -0.01, 0.0, 0.0, 0.0],
        [-2.0, 0.0, 1.0, -1.0, 0.0, -0.02, 0.0, 0.0, 0.0],
        [-2.0, 0.0, 1.0, 0.0, 0.0, -0.08, 0.0, -0.01, 0.01],
        [-2.0, 2.0, -1.0, 0.0, 0.0, -0.02, 0.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0,-1.0, 0.0, -0.10, 0.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0, 0.0, 0.0, -0.51, 0.0, -0.02, 0.03],
        [-1.0, 2.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0],
        [0.0, -2.0, 1.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0],
        [0.0, 0.0, -1.0, 0.0, 0.0, 0.02, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.06, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 1.0, 0.0, 0.01, 0.0, 0.0, 0.0],
        [0.0, 2.0, -1.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0],
        [1.0, -3.0, 0.0, 0.0, 1.0, -0.06, 0.0, 0.0, 0.0],
        [1.0, -2.0, 0.0, -1.0, 0.0, 0.01, 0.0, 0.0, 0.0],
        [1.0, -2.0, 0.0, 0.0, 0.0, -1.23, -0.07, 0.06, 0.01],
        [1.0, -1.0, 0.0, 0.0,-1.0, 0.02, 0.0, 0.0, 0.0],
        [1.0, -1.0, 0.0, 0.0, 1.0, 0.04, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, -1.0, 0.0, -0.22, 0.01, 0.01, 0.0],
        [1.0, 0.0, 0.0, 0.0, 0.0, 12.00, -0.80, -0.67, -0.03],
        [1.0, 0.0, 0.0, 1.0, 0.0, 1.73, -0.12, -0.10, 0.0],
        [1.0, 0.0, 0.0, 2.0, 0.0, -0.04, 0.0, 0.0, 0.0],
        [1.0, 1.0, 0.0, 0.0, -1.0, -0.50, -0.01, 0.03, 0.0],
        [1.0, 1.0, 0.0, 0.0, 1.0, 0.01, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 1.0, -1.0, -0.01, 0.0, 0.0, 0.0],
        [1.0, 2.0, -2.0, 0.0, 0.0, -0.01, 0.0, 0.0, 0.0],
        [1.0, 2.0, 0.0, 0.0, 0.0, -0.11, 0.01, 0.01, 0.0],
        [2.0, -2.0, 1.0, 0.0, 0.0, -0.01, 0.0, 0.0, 0.0],
        [2.0, 0.0,-1.0, 0.0, 0.0, -0.02, 0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    ])
    # get phase angles
    S, H, P, TAU, ZNS, PS = pyTMD.astro.phase_angles(MJD)
    # Compute the normalized position vector of coordinates
    radius = np.sqrt(np.sum(XYZ**2, axis=1))
    sinphi = XYZ[:,2]/radius
    cosphi = np.sqrt(XYZ[:,0]**2 + XYZ[:,1]**2)/radius
    sinla = XYZ[:,1]/cosphi/radius
    cosla = XYZ[:,0]/cosphi/radius
    zla = np.arctan2(XYZ[:,1], XYZ[:,0])
    # compute corrections (Mathews et al. 1997)
    DX = np.zeros((nt))
    DY = np.zeros((nt))
    DZ = np.zeros((nt))
    # iterate over rows in the table
    for i, row in enumerate(table):
        thetaf = TAU + S*row[0] + H*row[1] + P*row[2] + \
            ZNS*row[3] + PS*row[4]
        dr = 2.0*row[5]*sinphi*cosphi*np.sin(thetaf + zla) + \
            2.0*row[6]*sinphi*cosphi*np.cos(thetaf + zla)
        dn = row[7]*(cosphi**2 - sinphi**2)*np.sin(thetaf + zla) + \
            row[8]*(cosphi**2 - sinphi**2)*np.cos(thetaf + zla)
        de = row[7]*sinphi*np.cos(thetaf + zla) - \
            row[8]*sinphi*np.sin(thetaf + zla)
        DX += 1e-3*(dr*cosla*cosphi - de*sinla - dn*cosla*sinphi)
        DY += 1e-3*(dr*sinla*cosphi + de*cosla - dn*sinla*sinphi)
        DZ += 1e-3*(dr*sinphi + dn*cosphi)
    # return the corrections
    return np.c_[DX, DY, DZ]

def _frequency_dependence_long_period(
        XYZ: np.ndarray,
        MJD: np.ndarray
    ):
    """
    Computes the in-phase and out-of-phase corrections induced by mantle
    anelasticity in the long-period band

    Parameters
    ----------
    XYZ: np.ndarray
        Cartesian coordinates of the station (meters)
    MJD: np.ndarray
        Modified Julian Day (MJD)
    """
    # number of time steps
    nt = len(np.atleast_1d(MJD))
    # Corrections to Long-Peroid Tides for Frequency Dependence
    # of Love and Shida Number Parameters
    # table 7.3b of IERS conventions
    table = np.array([
        [0.0, 0.0, 0.0, 1.0, 0.0, 0.47, 0.23, 0.16, 0.07],
        [0.0, 2.0, 0.0, 0.0, 0.0, -0.20, -0.12, -0.11, -0.05],
        [1.0, 0.0, -1.0, 0.0, 0.0, -0.11, -0.08, -0.09, -0.04],
        [2.0, 0.0, 0.0, 0.0, 0.0, -0.13, -0.11, -0.15, -0.07],
        [2.0, 0.0, 0.0, 1.0, 0.0, -0.05, -0.05, -0.06, -0.03]
    ])
    # get phase angles
    S, H, P, TAU, ZNS, PS = pyTMD.astro.phase_angles(MJD)
    # Compute the normalized position vector of coordinates
    radius = np.sqrt(np.sum(XYZ**2, axis=1))
    sinphi = XYZ[:,2]/radius
    cosphi = np.sqrt(XYZ[:,0]**2 + XYZ[:,1]**2)/radius
    sinla = XYZ[:,1]/cosphi/radius
    cosla = XYZ[:,0]/cosphi/radius
    # compute corrections (Mathews et al. 1997)
    DX = np.zeros((nt))
    DY = np.zeros((nt))
    DZ = np.zeros((nt))
    # iterate over rows in the table
    for i, row in enumerate(table):
        thetaf = S*row[0] + H*row[1] + P*row[2] + ZNS*row[3] + PS*row[4]
        dr = row[5]*(3.0*sinphi**2 - 1.0)*np.cos(thetaf)/2.0 + \
            row[7]*(3.0*sinphi**2 - 1.0)*np.sin(thetaf)/2.0
        dn = row[6]*(2.0*cosphi*sinphi)*np.cos(thetaf) + \
            row[8]*(2.0*cosphi*sinphi)*np.sin(thetaf)
        de = 0.0
        DX += 1e-3*(dr*cosla*cosphi - de*sinla - dn*cosla*sinphi)
        DY += 1e-3*(dr*sinla*cosphi + de*cosla - dn*sinla*sinphi)
        DZ += 1e-3*(dr*sinphi + dn*cosphi)
    # return the corrections
    return np.c_[DX, DY, DZ]

def _free_to_mean(
        XYZ: np.ndarray,
        h2: float | np.ndarray,
        l2: float | np.ndarray
    ):
    """
    Calculate offsets for converting the permanent tide from
    a tide-free to a mean-tide state

    Parameters
    ----------
    XYZ: np.ndarray
        Cartesian coordinates of the station (meters)
    h2: float or np.ndarray
        Degree-2 Love number of vertical displacement
    l2: float or np.ndarray
        Degree-2 Love (Shida) number of horizontal displacement
    """
    # Compute the normalized position vector of coordinates
    radius = np.sqrt(np.sum(XYZ**2, axis=1))
    sinphi = XYZ[:,2]/radius
    cosphi = np.sqrt(XYZ[:,0]**2 + XYZ[:,1]**2)/radius
    sinla = XYZ[:,1]/cosphi/radius
    cosla = XYZ[:,0]/cosphi/radius
    # time-independent constituent of amplitude (Mathews et al. 1997)
    H0 = -0.31460
    # in Mathews et al. (1997): dR0=-0.1196 m with h2=0.6026
    dR0 = np.sqrt(5.0/(4.0*np.pi))*h2*H0
    # in Mathews et al. (1997): dN0=-0.0247 m with l2=0.0831
    dN0 = np.sqrt(45.0/(16.0*np.pi))*l2*H0
    # use double angle formula for sin(2*phi)
    dr = dR0*(3.0/2.0*sinphi**2 - 1.0/2.0)
    dn = 2.0*dN0*cosphi*sinphi
    # compute as an additive correction (Mathews et al. 1997)
    DX = dr*cosla*cosphi - dn*cosla*sinphi
    DY = dr*sinla*cosphi - dn*sinla*sinphi
    DZ = dr*sinphi + dn*cosphi
    # return the corrections
    return np.c_[DX, DY, DZ]
