#!/usr/bin/env python
u"""
arguments.py
Written by Tyler Sutterley (08/2024)
Calculates the nodal corrections for tidal constituents
Modification of ARGUMENTS fortran subroutine by Richard Ray 03/1999

CALLING SEQUENCE:
    pu, pf, G = arguments(MJD, constituents)

INPUTS:
    MJD: Modified Julian Day of input date
    constituents: tidal constituent IDs

OUTPUTS:
    pu, pf: nodal corrections for the constituents
    G: phase correction in degrees

OPTIONS:
    deltat: time correction for converting to Ephemeris Time (days)
    corrections: use nodal corrections from OTIS, FES or GOT models

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
    Updated 08/2024: add support for constituents in PERTH5 tables
        add back nodal arguments from PERTH3 for backwards compatibility
    Updated 01/2024: add function to create arguments coefficients table
        add function to calculate the arguments for minor constituents
        include multiples of 90 degrees as variable following Ray 2017
        add function to calculate the Doodson numbers for constituents
        add option to return None and not raise error for Doodson numbers
        moved constituent parameters function from predict to arguments
        added more constituent parameters for OTIS/ATLAS predictions
    Updated 12/2023: made keyword argument for selecting M1 coefficients
    Updated 08/2023: changed ESR netCDF4 format to TMD3 format
    Updated 04/2023: using renamed astro mean_longitudes function
        function renamed from original load_nodal_corrections
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

import numpy as np
import pyTMD.astro

def arguments(
        MJD: np.ndarray,
        constituents: list | np.ndarray,
        **kwargs
    ):
    """
    Calculates the nodal corrections for tidal constituents
    [1]_ [2]_ [3]_ [4]_

    Parameters
    ----------
    MJD: np.ndarray
        modified Julian day of input date
    constituents: list
        tidal constituent IDs
    deltat: float or np.ndarray, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS, FES or GOT models
    M1: str, default 'perth5'
        coefficients to use for M1 tides

                - ``'Doodson'``
                - ``'Ray'``
                - ``'perth5'``

    Returns
    -------
    pu: np.ndarray
        nodal angle correction
    pf: np.ndarray
        nodal factor correction
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
    # set default keyword arguments
    kwargs.setdefault('deltat', 0.0)
    kwargs.setdefault('corrections', 'OTIS')
    kwargs.setdefault('M1', 'perth5')

    # set function for astronomical longitudes
    # use ASTRO5 routines if not using an OTIS type model
    ASTRO5 = kwargs['corrections'] not in ('OTIS','ATLAS','TMD3','netcdf')
    # convert from Modified Julian Dates into Ephemeris Time
    s, h, p, n, pp = pyTMD.astro.mean_longitudes(MJD + kwargs['deltat'],
        ASTRO5=ASTRO5)

    # number of temporal values
    nt = len(np.atleast_1d(MJD))
    # initial time conversions
    hour = 24.0*np.mod(MJD, 1)
    # convert from hours solar time into mean lunar time in degrees
    tau = 15.0*hour - s + h
    # variable for multiples of 90 degrees (Ray technical note 2017)
    k = 90.0 + np.zeros((nt))

    # determine equilibrium arguments
    fargs = np.c_[tau, s, h, p, n, pp, k]
    G = np.dot(fargs, coefficients_table(constituents, **kwargs))

    # set nodal corrections
    # determine nodal corrections f and u for each model type
    pu, pf = nodal(n, p, constituents, **kwargs)

    # return values as tuple
    return (pu, pf, G)

def minor_arguments(
        MJD: np.ndarray,
        **kwargs
    ):
    """
    Calculates the nodal corrections for minor tidal constituents
    in order to infer their values [1]_ [2]_ [3]_ [4]_

    Parameters
    ----------
    MJD: np.ndarray
        modified Julian day of input date
    deltat: float or np.ndarray, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS, FES or GOT models

    Returns
    -------
    pu: np.ndarray
        nodal angle correction
    pf: np.ndarray
        nodal factor correction
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
    # set default keyword arguments
    kwargs.setdefault('deltat', 0.0)
    kwargs.setdefault('corrections', 'OTIS')

    # degrees to radians
    dtr = np.pi/180.0
    # set function for astronomical longitudes
    # use ASTRO5 routines if not using an OTIS type model
    ASTRO5 = kwargs['corrections'] not in ('OTIS','ATLAS','TMD3','netcdf')
    # convert from Modified Julian Dates into Ephemeris Time
    s, h, p, n, pp = pyTMD.astro.mean_longitudes(MJD + kwargs['deltat'],
        ASTRO5=ASTRO5)

    # number of temporal values
    nt = len(np.atleast_1d(MJD))
    # initial time conversions
    hour = 24.0*np.mod(MJD, 1)
    # convert from hours solar time into mean lunar time in degrees
    tau = 15.0*hour - s + h
    # variable for multiples of 90 degrees (Ray technical note 2017)
    k = 90.0 + np.zeros((nt))

    # determine equilibrium arguments
    fargs = np.c_[tau, s, h, p, n, pp, k]
    arg = np.dot(fargs, _minor_table())

    # determine nodal corrections f and u
    sinn = np.sin(n*dtr)
    cosn = np.cos(n*dtr)
    sin2n = np.sin(2.0*n*dtr)
    cos2n = np.cos(2.0*n*dtr)

    # nodal factor corrections for minor constituents
    f = np.ones((nt, 20))
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

    # nodal angle corrections for minor constituents
    u = np.zeros((nt, 20))
    u[:,0] = np.arctan2(0.189*sinn - 0.0058*sin2n,
        1.0 + 0.189*cosn - 0.0058*sin2n)# 2Q1
    u[:,1] = u[:,0]# sigma1
    u[:,2] = u[:,0]# rho1
    u[:,3] = np.arctan2( 0.185*sinn, 1.0 + 0.185*cosn)# M12
    u[:,4] = np.arctan2(-0.201*sinn, 1.0 + 0.201*cosn)# M11
    u[:,5] = np.arctan2(-0.221*sinn, 1.0 + 0.221*cosn)# chi1
    u[:,9] = np.arctan2(-0.198*sinn, 1.0 + 0.198*cosn)# J1
    u[:,10] = np.arctan2(-0.640*sinn - 0.134*sin2n,
        1.0 + 0.640*cosn + 0.134*cos2n)# OO1
    u[:,11] = np.arctan2(-0.0373*sinn, 1.0 - 0.0373*cosn)# 2N2
    u[:,12] = u[:,11]# mu2
    u[:,13] = u[:,11]# nu2
    u[:,15] = u[:,11]# L2
    u[:,16] = np.arctan2(-0.441*sinn, 1.0 + 0.441*cosn)# L2

    if kwargs['corrections'] in ('FES',):
        # additional astronomical terms for FES models
        II = np.arccos(0.913694997 - 0.035692561*np.cos(n*dtr))
        at1 = np.arctan(1.01883*np.tan(n*dtr/2.0))
        at2 = np.arctan(0.64412*np.tan(n*dtr/2.0))
        xi = -at1 - at2 + n*dtr
        xi = np.arctan2(np.sin(xi), np.cos(xi))
        nu = at1 - at2
        I2 = np.tan(II/2.0)
        Ra1 = np.sqrt(1.0 - 12.0*(I2**2)*np.cos(2.0*(p - xi)) + 36.0*(I2**4))
        P2 = np.sin(2.0*(p - xi))
        Q2 = 1.0/(6.0*(I2**2)) - np.cos(2.0*(p - xi))
        R = np.arctan(P2/Q2)

        # nodal factor corrections for minor constituents
        f[:,0] = np.sin(II)*(np.cos(II/2.0)**2)/0.38 # 2Q1
        f[:,1] = f[:,0] # sigma1
        f[:,2] = f[:,0] # rho1
        f[:,3] = f[:,0] # M12
        f[:,4] = np.sin(2.0*II)/0.7214 # M11
        f[:,5] = f[:,4] # chi1
        f[:,9] = f[:,4] # J1
        f[:,10] = np.sin(II)*np.power(np.sin(II/2.0),2.0)/0.01640 # OO1
        f[:,11] = np.power(np.cos(II/2.0),4.0)/0.9154 # 2N2
        f[:,12] = f[:,11] # mu2
        f[:,13] = f[:,11] # nu2
        f[:,14] = f[:,11] # lambda2
        f[:,15] = f[:,11]*Ra1 # L2
        f[:,18] = f[:,11] # eps2
        f[:,19] = np.power(np.sin(II),2.0)/0.1565 # eta2

        # nodal angle corrections for minor constituents
        u[:,0] = (2.0*xi - nu) # 2Q1
        u[:,1] = u[:,0] # sigma1
        u[:,2] = u[:,0] # rho1
        u[:,3] = u[:,0] # M12
        u[:,4] = -nu # M11
        u[:,5] = u[:,4] # chi1
        u[:,9] = u[:,4] # J1
        u[:,10] = (-2.0*xi - nu) # OO1
        u[:,11] = (2.0*xi - 2.0*nu) # 2N2
        u[:,12] = u[:,11] # mu2
        u[:,13] = u[:,11] # nu2
        u[:,14] = (2.0*xi - 2.0*nu) # lambda2
        u[:,15] = (2.0*xi - 2.0*nu - R)# L2
        u[:,18] = u[:,12] # eps2
        u[:,19] = -2.0*nu # eta2
    elif kwargs['corrections'] in ('GOT',):
        f[:,18] = f[:,11] # eps2
        f[:,19] = np.sqrt((1.0 + 0.436*cosn)**2 + (0.436*sinn)**2) # eta2
        u[:,18] = u[:,11] # eps2
        u[:,19] = np.arctan(-0.436*sinn/(1.0 + 0.436*cosn)) # eta2

    # return values as tuple
    return (u, f, arg)

def coefficients_table(
        constituents: list | tuple | np.ndarray | str,
        **kwargs
    ):
    """
    Doodson table coefficients for tidal constituents [1]_ [2]_

    Parameters
    ----------
    constituents: list, tuple, np.ndarray or str
        tidal constituent IDs
    corrections: str, default 'OTIS'
        use coefficients from OTIS, FES or GOT models

    Returns
    -------
    coef: np.ndarray
        Doodson coefficients (Cartwright numbers) for each constituent

    References
    ----------
    .. [1] A. T. Doodson and H. Lamb, "The harmonic development of
        the tide-generating potential", *Proceedings of the Royal Society
        of London. Series A, Containing Papers of a Mathematical and
        Physical Character*, 100(704), 305--329, (1921).
        `doi: 10.1098/rspa.1921.0088 <https://doi.org/10.1098/rspa.1921.0088>`_
    .. [2] A. T. Doodson and H. D. Warburg, "Admiralty Manual of Tides",
        HMSO, London, (1941).
    """
    # set default keyword arguments
    kwargs.setdefault('corrections', 'OTIS')

    # modified Doodson coefficients for constituents
    # using 7 index variables: tau, s, h, p, n, pp, k
    # tau: mean lunar time
    # s: mean longitude of moon
    # h: mean longitude of sun
    # p: mean longitude of lunar perigee
    # n: mean longitude of ascending lunar node
    # pp: mean longitude of solar perigee
    # k: 90-degree phase
    coefficients = {}
    coefficients['z0'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['node'] = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 2.0]
    coefficients['omega0'] = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 2.0]
    # With p'
    coefficients['sa'] = [0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0]
    # # Without p'
    # coefficients['sa'] = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['ssa'] = [0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    # With p'
    coefficients['sta'] = [0.0, 0.0, 3.0, 0.0, 0.0, -1.0, 0.0]
    # # Without p'
    # coefficients['sta'] = [0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['st'] = [0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['msm'] = [0.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['mm'] = [0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0]
    # annual sideline
    coefficients['msfa'] = [0.0, 2.0, -3.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['msf'] = [0.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    # annual sideline
    coefficients['msfb'] = [0.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['mfa'] = [0.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['mf-'] = [0.0, 2.0, 0.0, 0.0, -1.0, 0.0, 0.0]
    coefficients['mf'] = [0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # nodal line
    coefficients['mf+'] = [0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    # nodal line
    coefficients['mfn'] = [0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    coefficients['mfb'] = [0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['sn0'] = [0.0, 3.0, -2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['mst'] = [0.0, 3.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['mt'] = [0.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['mtm'] = [0.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['msqm'] = [0.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['mq'] = [0.0, 4.0, 0.0, -2.0, 0.0, 0.0, 0.0]
    coefficients['2smn0'] = [0.0, 5.0, -4.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['msp'] = [0.0, 5.0, -2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['mp'] = [0.0, 5.0, 0.0, -3.0, 0.0, 0.0, 0.0]
    coefficients['2qj1'] = [1.0, -6.0, 0.0, 3.0, 0.0, 0.0, 1.0]
    coefficients['2qk1'] = [1.0, -5.0, 0.0, 2.0, 0.0, 0.0, 1.0]
    coefficients['2qs1'] = [1.0, -5.0, 1.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['2qp1'] = [1.0, -5.0, 2.0, 2.0, 0.0, 0.0, -1.0]
    coefficients['2oj1'] = [1.0, -4.0, 0.0, 1.0, 0.0, 0.0, 1.0]
    coefficients['alpha1'] = [1.0, -4.0, 2.0, 1.0, 0.0, 0.0, -1.0]
    coefficients['2ok1'] = [1.0, -3.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    # 3rd degree terms
    coefficients["2q1'"] = [1.0, -3.0, 0.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['2q1'] = [1.0, -3.0, 0.0, 2.0, 0.0, 0.0, -1.0]
    coefficients['sigma1'] = [1.0, -3.0, 2.0, 0.0, 0.0, 0.0, -1.0]
    # 3rd degree terms
    coefficients["q1'"] = [1.0, -2.0, 0.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['q1'] = [1.0, -2.0, 0.0, 1.0, 0.0, 0.0, -1.0]
    coefficients['rho1'] = [1.0, -2.0, 2.0, -1.0, 0.0, 0.0, -1.0]
    coefficients['np1'] = [1.0, -2.0, 2.0, 1.0, 0.0, 0.0, 1.0]
    coefficients['opk1'] = [1.0, -1.0, -2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['oa1'] = [1.0, -1.0, -1.0, 0.0, 0.0, 0.0, -1.0]
    # O1 nodal line
    coefficients['o1n'] = [1.0, -1.0, 0.0, 0.0, -1.0, 0.0, -1.0]
    # O1 nodal line
    coefficients['o1-'] = [1.0, -1.0, 0.0, 0.0, -1.0, 0.0, -1.0]
    coefficients['o1'] = [1.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0]
    # conjugate to nodal line
    coefficients['o1+'] = [1.0, -1.0, 0.0, 0.0, 1.0, 0.0, -1.0]
    # 3rd degree terms
    coefficients["o1'"] = [1.0, -1.0, 0.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['ob1'] = [1.0, -1.0, 1.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['ms1'] = [1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['mp1'] = [1.0, -1.0, 2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['tau1'] = [1.0, -1.0, 2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['beta1'] = [1.0, 0.0, -2.0, 1.0, 0.0, 0.0, 1.0]
    coefficients['pqo1'] = [1.0, 0.0, -2.0, 1.0, 0.0, 0.0, -1.0]
    coefficients['2oq1'] = [1.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0]
    # 3rd degree terms
    coefficients["m1'"] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['m1'] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['m1a'] = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
    coefficients['m1b'] = [1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0]
    coefficients['no1'] = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
    coefficients['chi1'] = [1.0, 0.0, 2.0, -1.0, 0.0, 0.0, 1.0]
    coefficients['2pk1'] = [1.0, 1.0, -4.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['pi1'] = [1.0, 1.0, -3.0, 0.0, 0.0, 1.0, -1.0]
    coefficients['tk1'] = [1.0, 1.0, -3.0, 0.0, 0.0, 1.0, -1.0]
    coefficients['s1-1'] = [1.0, 1.0, -2.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['p1'] = [1.0, 1.0, -2.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['s1-'] = [1.0, 1.0, -1.0, 0.0, 0.0, -1.0, 1.0]
    if kwargs['corrections'] in ('OTIS','ATLAS','TMD3','netcdf'):
        coefficients['s1'] = [1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 1.0]
    else:
        # Doodson's phase
        coefficients['s1'] = [1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['s1+'] = [1.0, 1.0, -1.0, 0.0, 0.0, 1.0, 1.0]
    coefficients['ojm1'] = [1.0, 1.0, 0.0, -2.0, 0.0, 0.0, -1.0]
    # 3rd degree terms
    coefficients["k1'"] = [1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['k1-'] = [1.0, 1.0, 0.0, 0.0, -1.0, 0.0, -1.0]
    coefficients['k1'] = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['s1+1'] = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['k1+'] = [1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0]
    coefficients['k1n'] = [1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0]
    coefficients['k1++'] = [1.0, 1.0, 0.0, 0.0, 2.0, 0.0, -1.0]
    coefficients['psi1'] = [1.0, 1.0, 1.0, 0.0, 0.0, -1.0, 1.0]
    coefficients['rp1'] = [1.0, 1.0, 1.0, 0.0, 0.0, -1.0, -1.0]
    coefficients['phi1'] = [1.0, 1.0, 2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['kp1'] = [1.0, 1.0, 2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['opq1'] = [1.0, 2.0, -2.0, -1.0, 0.0, 0.0, -1.0]
    coefficients['the1'] = [1.0, 2.0, -2.0, 1.0, 0.0, 0.0, 1.0]
    coefficients['theta1'] = [1.0, 2.0, -2.0, 1.0, 0.0, 0.0, 1.0]
    coefficients['mq1'] = [1.0, 2.0, 0.0, -1.0, 0.0, 0.0, 1.0]
    coefficients['j1'] = [1.0, 2.0, 0.0, -1.0, 0.0, 0.0, 1.0]
    # 3rd degree terms
    coefficients["j1'"] = [1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['2po1'] = [1.0, 3.0, -4.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['so1'] = [1.0, 3.0, -2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['oo1'] = [1.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['kpq1'] = [1.0, 4.0, -2.0, -1.0, 0.0, 0.0, 1.0]
    coefficients['ups1'] = [1.0, 4.0, 0.0, -1.0, 0.0, 0.0, 1.0]
    coefficients['kq1'] = [1.0, 4.0, 0.0, -1.0, 0.0, 0.0, 1.0]
    coefficients['2jo1'] = [1.0, 5.0, 0.0, -2.0, 0.0, 0.0, -1.0]
    coefficients['kjq1'] = [1.0, 5.0, 0.0, -2.0, 0.0, 0.0, -1.0]
    coefficients['2ook1'] = [1.0, 5.0, 0.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['2oop1'] = [1.0, 5.0, 2.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['2jq1'] = [1.0, 6.0, 0.0, -3.0, 0.0, 0.0, -1.0]
    coefficients['2mn2s2'] = [2.0, -5.0, 4.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['2nk2'] = [2.0, -4.0, 0.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['3m(sk)2'] = [2.0, -4.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['3mks2'] = [2.0, -4.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2ns2'] = [2.0, -4.0, 2.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['3m2s2'] = [2.0, -4.0, 4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['oq2'] = [2.0, -3.0, 0.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['mnk2'] = [2.0, -3.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['3n2'] = [2.0, -3.0, 0.0, 3.0, 0.0, 0.0, 0.0]
    coefficients['eps2'] = [2.0, -3.0, 2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['mns2'] = [2.0, -3.0, 2.0, 1.0, 0.0, 0.0, 0.0]
    # coefficients['mns2'] = [2.0, -3.0, 4.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['mnus2'] = [2.0, -3.0, 4.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2ml2s2'] = [2.0, -3.0, 4.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['mnk2s2'] = [2.0, -3.0, 4.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['2ms2k2'] = [2.0, -2.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2mk2'] = [2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['o2'] = [2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 2.0]
    # 3rd degree terms
    coefficients["2n2'"] = [2.0, -2.0, 0.0, 1.0, 0.0, 0.0, 1.0]
    coefficients['2n2'] = [2.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['2nm2'] = [2.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['2ms2'] = [2.0, -2.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['mu2'] = [2.0, -2.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2mk2s2'] = [2.0, -2.0, 4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['omg2'] = [2.0, -2.0, 3.0, 0.0, 0.0, -1.0, 0.0]
    coefficients['nsk2'] = [2.0, -1.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['snk2'] = [2.0, -1.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['na2'] = [2.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0]
    # 3rd degree terms
    coefficients["n2'"] = [2.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['n2'] = [2.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['2ml2'] = [2.0, -1.0, 0.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['nb2'] = [2.0, -1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['nu2'] = [2.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['mmun2'] = [2.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['3m(sn)2'] = [2.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['nks2'] = [2.0, -1.0, 2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['nkp2'] = [2.0, -1.0, 2.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['mkl2s2'] = [2.0, -1.0, 4.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['2kn2s2'] = [2.0, -1.0, 4.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['msk2'] = [2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['op2'] = [2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['m2-2'] = [2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['gamma2'] = [2.0, 0.0, -2.0, 2.0, 0.0, 0.0, 2.0]
    coefficients['gam2'] = [2.0, 0.0, -2.0, 2.0, 0.0, 0.0, 2.0]
    coefficients['alp2'] = [2.0, 0.0, -1.0, 0.0, 0.0, 1.0, 2.0]
    coefficients['alpha2'] = [2.0, 0.0, -1.0, 0.0, 0.0, 1.0, 2.0]
    coefficients['m2a'] = [2.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0]
    coefficients['m2-1'] = [2.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['ma2'] = [2.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    # M2 nodal line
    coefficients['m2-'] = [2.0, 0.0, 0.0, 0.0, -1.0, 0.0, 2.0]
    coefficients['m2'] = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['ko2'] = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # conjugate to nodal
    coefficients['m2+'] = [2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 2.0]
    # 3rd degree terms
    coefficients["m2'"] = [2.0, 0.0, 0.0, 1.0, 0.0, 0.0, -1.0]
    coefficients['mb2'] = [2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['m2+1'] = [2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['m2b'] = [2.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0]
    coefficients['beta2'] = [2.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0]
    coefficients['delta2'] = [2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['m2+2'] = [2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['mks2'] = [2.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['m2(ks)2'] = [2.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2snmk2'] = [2.0, 1.0, -4.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['2sn(mk)2'] = [2.0, 1.0, -4.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['snm2'] = [2.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['lambda2'] = [2.0, 1.0, -2.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['2mn2'] = [2.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['l2'] = [2.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['l2a'] = [2.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0]
    # 3rd degree terms
    coefficients["l2'"] = [2.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['nkm2'] = [2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['l2b'] = [2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['lb2'] = [2.0, 1.0, 1.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['2sk2'] = [2.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2t2'] = [2.0, 2.0, -4.0, 0.0, 0.0, 2.0, 0.0]
    coefficients['s2-2'] = [2.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['s2-1'] = [2.0, 2.0, -3.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['t2'] = [2.0, 2.0, -3.0, 0.0, 0.0, 1.0, 0.0]
    coefficients['kp2'] = [2.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['s2r'] = [2.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['s2'] = [2.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['r2'] = [2.0, 2.0, -1.0, 0.0, 0.0, -1.0, 2.0]
    coefficients['s2+1'] = [2.0, 2.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['s2+2'] = [2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['k2'] = [2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['kb2'] = [2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2ks2'] = [2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2pmn2'] = [2.0, 3.0, -4.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['msnu2'] = [2.0, 3.0, -4.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['msn2'] = [2.0, 3.0, -2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['zeta2'] = [2.0, 3.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['eta2'] = [2.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['mkn2'] = [2.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['kj2'] = [2.0, 3.0, 0.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['2kmsn2'] = [2.0, 3.0, 2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2km(sn)2'] = [2.0, 3.0, 2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2sm2'] = [2.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2ms2n2'] = [2.0, 4.0, -2.0, -2.0, 0.0, 0.0, 0.0]
    coefficients['skm2'] = [2.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2j2'] = [2.0, 4.0, 0.0, -2.0, 0.0, 0.0, 2.0]
    coefficients['2k2'] = [2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2snu2'] = [2.0, 5.0, -6.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['3(sm)n2'] = [2.0, 5.0, -6.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['2sn2'] = [2.0, 5.0, -4.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['skn2'] = [2.0, 5.0, -2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2kn2'] = [2.0, 5.0, 0.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['3s2m2'] = [2.0, 6.0, -6.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2sk2m2'] = [2.0, 6.0, -4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2oq3'] = [3.0, -4.0, 0.0, 1.0, 0.0, 0.0, 1.0]
    # compound 3 O1
    coefficients['o3'] = [3.0, -3.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['nq3'] = [3.0, -3.0, 0.0, 2.0, 0.0, 0.0, -1.0]
    coefficients['muo3'] = [3.0, -3.0, 2.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['mq3'] = [3.0, -2.0, 0.0, 1.0, 0.0, 0.0, -1.0]
    coefficients['no3'] = [3.0, -2.0, 0.0, 1.0, 0.0, 0.0, -1.0]
    coefficients['mnp3'] = [3.0, -2.0, 2.0, 1.0, 0.0, 0.0, 1.0]
    coefficients['2op3'] = [3.0, -1.0, -2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['2os3'] = [3.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    # Q1+M1+S1
    coefficients['qms3'] = [3.0, -1.0, -1.0, 2.0, 0.0, 0.0, 2.0]
    coefficients['mo3-'] = [3.0, -1.0, 0.0, 0.0, -1.0, 0.0, -1.0]
    coefficients['mo3'] = [3.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['mo3+'] = [3.0, -1.0, 0.0, 0.0, 1.0, 0.0, -1.0]
    coefficients['2mk3'] = [3.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['e3n'] = [3.0, -1.0, 0.0, 1.0, -1.0, 0.0, 0.0]
    # 3rd degree terms
    coefficients['e3'] = [3.0, -1.0, 0.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['2no3'] = [3.0, -1.0, 0.0, 2.0, 0.0, 0.0, 1.0]
    coefficients['2nkm3'] = [3.0, -1.0, 0.0, 2.0, 0.0, 0.0, 1.0]
    coefficients['2ms3'] = [3.0, -1.0, 1.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['2mp3'] = [3.0, -1.0, 2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['ns3'] = [3.0, 0.0, -1.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['2oj3'] = [3.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0]
    # 2M2 - M1
    coefficients['2mm3'] = [3.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0]
    coefficients['m3'] = [3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # 3rd degree terms
    coefficients["m3'"] = [3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['nk3'] = [3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['mp3'] = [3.0, 1.0, -2.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['so3'] = [3.0, 1.0, -2.0, 0.0, 0.0, 0.0, -1.0]
    # 3rd degree terms
    coefficients['lambda3'] = [3.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['ms3'] = [3.0, 1.0, -1.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['pjrho3'] = [3.0, 1.0, 0.0, -2.0, 0.0, 0.0, -1.0]
    # 3rd degree terms
    coefficients['l3'] = [3.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['mk3-'] = [3.0, 1.0, 0.0, 0.0, -1.0, 0.0, 1.0]
    coefficients['mk3'] = [3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['mk3+'] = [3.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0]
    # 3rd degree terms
    coefficients['l3b'] = [3.0, 1.0, 0.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['mks3'] = [3.0, 1.0, 1.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['mkp3'] = [3.0, 1.0, 2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['nso3'] = [3.0, 2.0, -2.0, 1.0, 0.0, 0.0, 1.0]
    coefficients['2mq3'] = [3.0, 2.0, 0.0, -1.0, 0.0, 0.0, 1.0]
    # 3rd degree terms
    coefficients['f3'] = [3.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0]
    # 3rd degree terms
    coefficients['j3'] = [3.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0]
    # s3 perturbation
    coefficients['2t3'] = [3.0, 3.0, -5.0, 0.0, 0.0, 0.0, 2.0]
    # s3 perturbation
    coefficients['t3'] = [3.0, 3.0, -4.0, 0.0, 0.0, 0.0, 2.0]
    # = 2SK3
    coefficients['sp3'] = [3.0, 3.0, -4.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['s3'] = [3.0, 3.0, -3.0, 0.0, 0.0, 0.0, 0.0]
    # 3rd degree terms
    coefficients["s3'"] = [3.0, 3.0, -3.0, 0.0, 0.0, 0.0, 2.0]
    # s3 perturbation
    coefficients['r3'] = [3.0, 3.0, -2.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['sk3'] = [3.0, 3.0, -2.0, 0.0, 0.0, 0.0, 1.0]
    # s3 perturbation
    coefficients['2r3'] = [3.0, 3.0, -1.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['k3'] = [3.0, 3.0, 0.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['2so3'] = [3.0, 5.0, -4.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['2jp3'] = [3.0, 5.0, -2.0, -2.0, 0.0, 0.0, 1.0]
    coefficients['kso3'] = [3.0, 5.0, -2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['2jk3'] = [3.0, 5.0, -1.0, -2.0, 0.0, 0.0, 0.0]
    coefficients['2ko3'] = [3.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['2sq3'] = [3.0, 6.0, -4.0, -1.0, 0.0, 0.0, 1.0]
    coefficients['o4'] = [4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2qm4'] = [4.0, -4.0, 0.0, 2.0, 0.0, 0.0, 2.0]
    coefficients['4ms4'] = [4.0, -4.0, 4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['4m2s4'] = [4.0, -4.0, 4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2mnk4'] = [4.0, -3.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['moq4'] = [4.0, -3.0, 0.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['2mns4'] = [4.0, -3.0, 2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['2mns4'] = [4.0, -3.0, 4.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2mnus4'] = [4.0, -3.0, 4.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['3mk4'] = [4.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2om4'] = [4.0, -2.0, 0.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['n4'] = [4.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['3ms4'] = [4.0, -2.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['msnk4'] = [4.0, -1.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['mpq4'] = [4.0, -1.0, -2.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['mn4'] = [4.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['mnu4'] = [4.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2mls4'] = [4.0, -1.0, 2.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['mnks4'] = [4.0, -1.0, 2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['2msk4'] = [4.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['mop4'] = [4.0, 0.0, -2.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['ma4'] = [4.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['m4'] = [4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['mb4'] = [4.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2mks4'] = [4.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['sn4'] = [4.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['3mn4'] = [4.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['ml4'] = [4.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['nk4'] = [4.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['2smk4'] = [4.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2pm4'] = [4.0, 2.0, -4.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['mt4'] = [4.0, 2.0, -3.0, 0.0, 0.0, 1.0, 0.0]
    coefficients['ms4'] = [4.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['mr4'] = [4.0, 2.0, -1.0, 0.0, 0.0, -1.0, 2.0]
    coefficients['mk4'] = [4.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2snm4'] = [4.0, 3.0, -4.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['2msn4'] = [4.0, 3.0, -2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['sl4'] = [4.0, 3.0, -2.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['2mkn4'] = [4.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['mkj4'] = [4.0, 3.0, 0.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['2t4'] = [4.0, 4.0, -6.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['t4'] = [4.0, 4.0, -5.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['st4'] = [4.0, 4.0, -5.0, 0.0, 0.0, 1.0, 0.0]
    coefficients['s4'] = [4.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['r4'] = [4.0, 4.0, -3.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2r4'] = [4.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['sk4'] = [4.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['3ks4'] = [4.0, 4.0, -1.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['k4'] = [4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['3sm4'] = [4.0, 6.0, -6.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2(kj)4'] = [4.0, 6.0, 0.0, -2.0, 0.0, 0.0, 2.0]
    coefficients['2no5'] = [5.0, -3.0, 0.0, 2.0, 0.0, 0.0, -1.0]
    coefficients['mno5'] = [5.0, -2.0, 0.0, 1.0, 0.0, 0.0, -1.0]
    coefficients['2mq5'] = [5.0, -2.0, 0.0, 1.0, 0.0, 0.0, -1.0]
    coefficients['2nkms5'] = [5.0, -2.0, 2.0, 2.0, 0.0, 0.0, 1.0]
    coefficients['2mo5'] = [5.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['mnm5'] = [5.0, -1.0, 0.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['2nk5'] = [5.0, -1.0, 0.0, 2.0, 0.0, 0.0, 1.0]
    coefficients['3ms5'] = [5.0, -1.0, 1.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['3mp5'] = [5.0, -1.0, 2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['nso5'] = [5.0, 0.0, -2.0, 1.0, 0.0, 0.0, -1.0]
    coefficients['m5'] = [5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['mnk5'] = [5.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
    coefficients['mb5'] = [5.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['mso5'] = [5.0, 1.0, -2.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['2mp5'] = [5.0, 1.0, -2.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['2ms5'] = [5.0, 1.0, -1.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['2mk5'] = [5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    # N2 + K2 + S1
    coefficients['nks5'] = [5.0, 2.0, -1.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['nsk5'] = [5.0, 2.0, -2.0, -1.0, 0.0, 0.0, 1.0]
    coefficients['msm5'] = [5.0, 2.0, -2.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['snk5'] = [5.0, 2.0, -2.0, 1.0, 0.0, 0.0, -1.0]
    coefficients['3mq5'] = [5.0, 2.0, 0.0, -1.0, 0.0, 0.0, 1.0]
    coefficients['msp5'] = [5.0, 3.0, -4.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['msk5'] = [5.0, 3.0, -2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['3km5'] = [5.0, 3.0, 0.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['2sp5'] = [5.0, 5.0, -6.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['t5'] = [5.0, 5.0, -6.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['s5'] = [5.0, 5.0, -5.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['r5'] = [5.0, 5.0, -4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2sk5'] = [5.0, 5.0, -4.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['k5'] = [5.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['o6'] = [6.0, -6.0, 0.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['2(mn)k6'] = [6.0, -4.0, 0.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['5mks6'] = [6.0, -4.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2(mn)s6'] = [6.0, -4.0, 2.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['5m2s6'] = [6.0, -4.0, 4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['3mnk6'] = [6.0, -3.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['n6'] = [6.0, -3.0, 0.0, 3.0, 0.0, 0.0, 0.0]
    coefficients['3mns6'] = [6.0, -3.0, 2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['2nmls6'] = [6.0, -3.0, 2.0, 1.0, 0.0, 0.0, 2.0]
    coefficients['3nks6'] = [6.0, -3.0, 2.0, 3.0, 0.0, 0.0, 0.0]
    coefficients['3mnus6'] = [6.0, -3.0, 4.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['4mk6'] = [6.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2nm6'] = [6.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['m2n6'] = [6.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['4ms6'] = [6.0, -2.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2mlns6'] = [6.0, -2.0, 2.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['2mnls6'] = [6.0, -2.0, 2.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['2nmks6'] = [6.0, -2.0, 2.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['2msnk6'] = [6.0, -1.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['2mn6'] = [6.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['2mnu6'] = [6.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2mn6'] = [6.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['3mls6'] = [6.0, -1.0, 2.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['2mno6'] = [6.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2mnks6'] = [6.0, -1.0, 2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['3msk6'] = [6.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['ma6'] = [6.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['m6'] = [6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2nk6'] = [6.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['mb6'] = [6.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['3mks6'] = [6.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['msn6'] = [6.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['4mn6'] = [6.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2ml6'] = [6.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['mnk6'] = [6.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['mkn6'] = [6.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['mknu6'] = [6.0, 1.0, 2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2(ms)k6'] = [6.0, 2.0, -4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2mt6'] = [6.0, 2.0, -3.0, 0.0, 0.0, 1.0, 0.0]
    coefficients['2ms6'] = [6.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2mr6'] = [6.0, 2.0, -1.0, 0.0, 0.0, -1.0, 2.0]
    coefficients['2mk6'] = [6.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2sn6'] = [6.0, 3.0, -4.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['3msn6'] = [6.0, 3.0, -2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['msl6'] = [6.0, 3.0, -2.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['snk6'] = [6.0, 3.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['mkl6'] = [6.0, 3.0, 0.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['3mkn6'] = [6.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2sm6'] = [6.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['msk6'] = [6.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2km6'] = [6.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2(ms)n6'] = [6.0, 5.0, -4.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2mskn6'] = [6.0, 5.0, -2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2t6'] = [6.0, 6.0, -8.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['t6'] = [6.0, 6.0, -7.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['s6'] = [6.0, 6.0, -6.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['r6'] = [6.0, 6.0, -5.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2r6'] = [6.0, 6.0, -4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2sk6'] = [6.0, 6.0, -4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['k6'] = [6.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['4sm6'] = [6.0, 8.0, -8.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2mno7'] = [7.0, -2.0, 0.0, 1.0, 0.0, 0.0, -1.0]
    coefficients['4mk7'] = [7.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['3mo7'] = [7.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['2nmk7'] = [7.0, -1.0, 0.0, 2.0, 0.0, 0.0, 1.0]
    # = 2MNK7 = MNKO7
    coefficients['m7'] = [7.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
    coefficients['3mp7'] = [7.0, 1.0, -2.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['2mso7'] = [7.0, 1.0, -2.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['3mk7'] = [7.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    # = MSKO7 (noaa)
    coefficients['2msk7'] = [7.0, 3.0, -2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['msko7'] = [7.0, 3.0, -2.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['2t7'] = [7.0, 7.0, -9.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['3sp7'] = [7.0, 7.0, -8.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['t7'] = [7.0, 7.0, -8.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['s7'] = [7.0, 7.0, -7.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['r7'] = [7.0, 7.0, -6.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['3sk7'] = [7.0, 7.0, -6.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['2r7'] = [7.0, 7.0, -5.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['3r7'] = [7.0, 7.0, -4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['4r7'] = [7.0, 7.0, -3.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['k7'] = [7.0, 7.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['3m2ns8'] = [8.0, -4.0, 2.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['4mns8'] = [8.0, -3.0, 2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['5mk8'] = [8.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2mn8'] = [8.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['2(mn)8'] = [8.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
    coefficients['5ms8'] = [8.0, -2.0, 2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['3mn8'] = [8.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['3mnu8'] = [8.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['3mn8'] = [8.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['4mls8'] = [8.0, -1.0, 2.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['ma8'] = [8.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['m8'] = [8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['mb8'] = [8.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2msn8'] = [8.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['3ml8'] = [8.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['2mnk8'] = [8.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['3mt8'] = [8.0, 2.0, -3.0, 0.0, 0.0, 1.0, 0.0]
    coefficients['3ms8'] = [8.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['3mk8'] = [8.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2smn8'] = [8.0, 3.0, -4.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['2msl8'] = [8.0, 3.0, -2.0, -1.0, 0.0, 0.0, 2.0]
    coefficients['msnk8'] = [8.0, 3.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['4msn8'] = [8.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0]
    coefficients['2(ms)8'] = [8.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2msk8'] = [8.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['2(mk)8'] = [8.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['s8'] = [8.0, 8.0, -8.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['k8'] = [8.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['4mo9'] = [9.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['2(mn)k9'] = [9.0, -1.0, 0.0, 2.0, 0.0, 0.0, 1.0]
    coefficients['2m2nk9'] = [9.0, -1.0, 0.0, 2.0, 0.0, 0.0, 1.0]
    coefficients['3mnk9'] = [9.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
    coefficients['4mp9'] = [9.0, 1.0, -2.0, 0.0, 0.0, 0.0, -1.0]
    coefficients['4ms9'] = [9.0, 1.0, -1.0, 0.0, 0.0, 0.0, 2.0]
    coefficients['4mk9'] = [9.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['3msk9'] = [9.0, 3.0, -2.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['s9'] = [9.0, 9.0, -9.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['k9'] = [9.0, 9.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    coefficients['4mn10'] = [10.0, -1, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['m10'] = [10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['3msn10'] = [10.0, 1.0, -2, 1.0, 0.0, 0.0, 0.0]
    coefficients['4ms10'] = [10.0, 2.0, -2, 0.0, 0.0, 0.0, 0.0]
    coefficients['4mk10'] = [10.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['3msl10'] = [10.0, 3.0, -2, -1.0, 0.0, 0.0, 2.0]
    coefficients['3m2s10'] = [10.0, 4.0, -4, 0.0, 0.0, 0.0, 0.0]
    coefficients['3msk10'] = [10.0, 4.0, -2, 0.0, 0.0, 0.0, 0.0]
    coefficients['s10'] = [10.0, 10.0, -10.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['4msk11'] = [11.0, 3.0, -2, 0.0, 0.0, 0.0, 1.0]
    coefficients['s11'] = [11.0, 11.0, -11.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['5mn12'] = [12.0, -1, 0.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['m12'] = [12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['4msn12'] = [12.0, 1.0, -2, 1.0, 0.0, 0.0, 0.0]
    coefficients['4mns12'] = [12.0, 1.0, -2, 1.0, 0.0, 0.0, 0.0]
    coefficients['5ms12'] = [12.0, 2.0, -2, 0.0, 0.0, 0.0, 0.0]
    coefficients['4msl12'] = [12.0, 3.0, -2, -1.0, 0.0, 0.0, 2.0]
    coefficients['4m2s12'] = [12.0, 4.0, -4, 0.0, 0.0, 0.0, 0.0]
    coefficients['s12'] = [12.0, 12.0, -12.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['m14'] = [14.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['5msn14'] = [14.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0]
    coefficients['6ms14'] = [14.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]
    coefficients['5m2s14'] = [14.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0]

    # set constituents to be iterable
    if isinstance(constituents, str):
        constituents = [constituents]
    # allocate for output coefficients
    nc = len(constituents)
    coef = np.zeros((7, nc))
    # for each constituent of interest
    for i, c in enumerate(constituents):
        try:
            coef[:,i] = coefficients[c]
        except KeyError:
            raise ValueError(f'Unsupported constituent: {c}')

    # return Doodson coefficients for constituents
    return coef

def doodson_number(
        constituents: str | list | np.ndarray,
        **kwargs
    ):
    """
    Calculates the Doodson or Cartwright number for
    tidal constituents [1]_

    Parameters
    ----------
    constituents: str, list or np.ndarray
        tidal constituent ID(s)
    corrections: str, default 'OTIS'
        use arguments from OTIS, FES or GOT models
    formalism: str, default 'Doodson'
        constituent identifier formalism

            - ``'Cartwright'``
            - ``'Doodson'``
    raise_error: bool, default True
        Raise exception if constituent is unsupported

    Returns
    -------
    numbers: float, np.ndarray or dict
        Doodson or Cartwright number for each constituent

    References
    ----------
    .. [1] A. T. Doodson and H. Lamb, "The harmonic development of
        the tide-generating potential", *Proceedings of the Royal Society
        of London. Series A, Containing Papers of a Mathematical and
        Physical Character*, 100(704), 305--329, (1921).
        `doi: 10.1098/rspa.1921.0088 <https://doi.org/10.1098/rspa.1921.0088>`_
    """
    # set default keyword arguments
    kwargs.setdefault('corrections', 'OTIS')
    kwargs.setdefault('formalism', 'Doodson')
    kwargs.setdefault('raise_error', True)
    # validate inputs
    assert kwargs['formalism'].title() in ('Cartwright', 'Doodson'), \
        f'Unknown formalism {kwargs["formalism"]}'
    # get the coefficients of coefficients
    if isinstance(constituents, str):
        # try to get the Doodson coefficients for constituent
        try:
            coefficients = coefficients_table(constituents.lower(), **kwargs)
        except ValueError as exc:
            if kwargs['raise_error']:
                raise ValueError(f'Unsupported constituent {constituents}')
            else:
                return None
        # extract identifier in formalism
        if (kwargs['formalism'] == 'Cartwright'):
            # extract Cartwright number
            numbers = np.array(coefficients[:6,0])
        elif (kwargs['formalism'] == 'Doodson'):
            # convert from coefficients to Doodson number
            numbers = _to_doodson_number(coefficients[:,0], **kwargs)
    else:
        # output dictionary with Doodson numbers
        numbers = {}
        # for each input constituent
        for i,c in enumerate(constituents):
            # try to get the Doodson coefficients for constituent
            try:
                coefficients = coefficients_table(c.lower(), **kwargs)
            except ValueError as exc:
                if kwargs['raise_error']:
                    raise ValueError(f'Unsupported constituent {c}')
                else:
                    numbers[c] = None
                    continue
            # convert from coefficients to Doodson number
            if (kwargs['formalism'] == 'Cartwright'):
                # extract Cartwright number
                numbers[c] = np.array(coefficients[:6,0])
            elif (kwargs['formalism'] == 'Doodson'):
                # convert from coefficients to Doodson number
                numbers[c] = _to_doodson_number(coefficients[:,0], **kwargs)
    # return the Doodson or Cartwright number
    return numbers

# PURPOSE: compute the nodal corrections
def nodal(
        n: np.ndarray,
        p: np.ndarray,
        constituents: list | tuple | np.ndarray | str,
        **kwargs
    ):
    """
    Calculates the nodal corrections for tidal constituents
    [1]_ [2]_ [3]_ [4]_

    Calculates factors for compound tides using recursion

    Parameters
    ----------
    n: np.ndarray
        mean longitude of ascending lunar node (degrees)
    p: np.ndarray
        mean longitude of lunar perigee (degrees)
    constituents: list, tuple, np.ndarray or str
        tidal constituent IDs
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS, FES or GOT models
    M1: str, default 'perth5'
        coefficients to use for M1 tides

                - ``'Doodson'``
                - ``'Ray'``
                - ``'perth5'``

    Returns
    -------
    f: np.ndarray
        nodal factor correction
    u: np.ndarray
        nodal angle correction

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
    kwargs.setdefault('corrections', 'OTIS')
    kwargs.setdefault('M1', 'perth5')
    # set correction type
    OTIS_TYPE = kwargs['corrections'] in ('OTIS','ATLAS','TMD3','netcdf')
    FES_TYPE = kwargs['corrections'] in ('FES',)
    PERTH3_TYPE = kwargs['corrections'] in ('perth3',)

    # degrees to radians
    dtr = np.pi/180.0
    # trigonometric factors for nodal corrections
    sinn = np.sin(n*dtr)
    cosn = np.cos(n*dtr)
    sin2n = np.sin(2.0*n*dtr)
    cos2n = np.cos(2.0*n*dtr)
    sin3n = np.sin(3.0*n*dtr)
    sinp  = np.sin(p*dtr)
    cosp  = np.cos(p*dtr)
    sin2p = np.sin(2.0*p*dtr)
    cos2p = np.cos(2.0*p*dtr)

    # set constituents to be iterable
    if isinstance(constituents, str):
        constituents = [constituents]

    # set nodal corrections
    nt = len(np.atleast_1d(n))
    nc = len(constituents)
    # nodal factor correction
    f = np.zeros((nt, nc))
    # nodal angle correction
    u = np.zeros((nt, nc))

    # additional astronomical terms for FES models
    II = np.arccos(0.913694997 - 0.035692561*np.cos(n*dtr))
    at1 = np.arctan(1.01883*np.tan(n*dtr/2.0))
    at2 = np.arctan(0.64412*np.tan(n*dtr/2.0))
    xi = -at1 - at2 + n*dtr
    xi = np.arctan2(np.sin(xi), np.cos(xi))
    nu = at1 - at2
    I2 = np.tan(II/2.0)
    Ra1 = np.sqrt(1.0 - 12.0*(I2**2)*np.cos(2.0*(p - xi)) + 36.0*(I2**4))
    P2 = np.sin(2.0*(p - xi))
    Q2 = 1.0/(6.0*(I2**2)) - np.cos(2.0*(p - xi))
    R = np.arctan(P2/Q2)
    P_prime = np.sin(2.0*II)*np.sin(nu)
    Q_prime = np.sin(2.0*II)*np.cos(nu) + 0.3347
    nu_prime = np.arctan(P_prime/Q_prime)
    P_sec = (np.sin(II)**2)*np.sin(2.0*nu)
    Q_sec = (np.sin(II)**2)*np.cos(2.0*nu) + 0.0727
    nu_sec = 0.5*np.arctan(P_sec/Q_sec)

    # compute standard nodal corrections f and u
    for i, c in enumerate(constituents):
        if c in ('msf','tau1','p1','theta1','lambda2','s2') and OTIS_TYPE:
            term1 = 0.0
            term2 = 1.0
        elif c in ('p1','s2') and (FES_TYPE or PERTH3_TYPE):
            term1 = 0.0
            term2 = 1.0
        elif c in ('mm','msm') and OTIS_TYPE:
            term1 = 0.0
            term2 = 1.0 - 0.130*cosn
        elif c in ('mm','msm') and FES_TYPE:
            term1 = 0.0
            term2 = (2.0/3.0 - np.power(np.sin(II),2.0))/0.5021
        elif c in ('mm','msm'):
            term1 = -0.0534*sin2p - 0.0219*np.sin((2.0*p-n)*dtr)
            term2 = 1.0 - 0.1308*cosn - 0.0534*cos2p - 0.0219*np.cos((2.0*p-n)*dtr)
        elif c in ('mf','msqm','msp','mq','mtm') and OTIS_TYPE:
            f[:,i] = 1.043 + 0.414*cosn
            u[:,i] = dtr*(-23.7*sinn + 2.7*sin2n - 0.4*sin3n)
            continue
        elif c in ('mf','msqm','msp','mq','mt','mtm') and FES_TYPE:
            f[:,i] = np.power(np.sin(II),2.0)/0.1578
            u[:,i] = -2.0*xi
            continue
        elif c in ('mf','msqm','msp','mq'):
            term1 = -0.04324*sin2p - 0.41465*sinn - 0.03873*sin2n
            term2 = 1.0 + 0.04324*cos2p + 0.41465*cosn + 0.03873*cos2n
        elif c in ('mt',) and OTIS_TYPE:
            term1 = -0.203*sinn - 0.040*sin2n
            term2 = 1.0 + 0.203*cosn + 0.040*cos2n
        elif c in ('mt','mtm',):
            term1 = -0.018*sin2p - 0.4145*sinn - 0.040*sin2n
            term2 = 1.0 + 0.018*cos2p + 0.4145*cosn + 0.040*cos2n
        elif c in ('msf',) and FES_TYPE:
            f[:,i] = 1.0
            u[:,i] = (2.0*xi - 2.0*nu)
            continue
        elif c in ('msf',):
            # linear tide and not compound
            term1 = 0.137*sinn
            term2 = 1.0
        elif c in ('mst',):
            term1 = -0.380*sin2p - 0.413*sinn - 0.037*sin2n
            term2 = 1.0 + 0.380*cos2p + 0.413*cosn + 0.037*cos2n
        elif c in ('o1','so3','op2') and OTIS_TYPE:
            term1 = 0.189*sinn - 0.0058*sin2n
            term2 = 1.0 + 0.189*cosn - 0.0058*cos2n
            f[:,i] = np.sqrt(term1**2 + term2**2) # O1
            u[:,i] = dtr*(10.8*sinn - 1.3*sin2n + 0.2*sin3n)
            continue
        elif c in ('o1','so3','op2','2q1','q1','rho1','sigma1') and FES_TYPE:
            f[:,i] = np.sin(II)*(np.cos(II/2.0)**2)/0.38
            u[:,i] = (2.0*xi - nu)
            continue
        elif c in ('q1','o1') and PERTH3_TYPE:
            f[:,i] = 1.009 + 0.187*cosn - 0.015*cos2n
            u[:,i] = dtr*(10.8*sinn - 1.3*sin2n)
            continue
        elif c in ('o1','so3','op2'):
            term1 = 0.1886*sinn - 0.0058*sin2n - 0.0065*sin2p
            term2 = 1.0 + 0.1886*cosn - 0.0058*cos2n - 0.0065*cos2p
        elif c in ('2q1','q1','rho1','sigma1') and OTIS_TYPE:
            f[:,i] = np.sqrt((1.0 + 0.188*cosn)**2 + (0.188*sinn)**2)
            u[:,i] = np.arctan(0.189*sinn/(1.0 + 0.189*cosn))
            continue
        elif c in ('2q1','q1','rho1','sigma1'):
            term1 = 0.1886*sinn
            term2 = 1.0 + 0.1886*cosn
        elif c in ('tau1',):
            term1 = 0.219*sinn
            term2 = 1.0 - 0.219*cosn
        elif c in ('beta1',):
            term1 = 0.226*sinn
            term2 = 1.0 + 0.226*cosn
        elif c in ('m1',) and (kwargs['M1'] == 'Doodson'):
            # A. T. Doodson's coefficients for M1 tides
            term1 = sinp + 0.2*np.sin((p-n)*dtr)
            term2 = 2.0*cosp + 0.4*np.cos((p-n)*dtr)
        elif c in ('m1',) and (kwargs['M1'] == 'Ray'):
            # R. Ray's coefficients for M1 tides (perth3)
            term1 = 0.64*sinp + 0.135*np.sin((p-n)*dtr)
            term2 = 1.36*cosp + 0.267*np.cos((p-n)*dtr)
        elif c in ('m1',) and (kwargs['M1'] == 'perth5'):
            # assumes M1 argument includes p
            term1 = -0.2294*sinn - 0.3594*sin2p - 0.0664*np.sin((2.0*p-n)*dtr)
            term2 = 1.0 + 0.1722*cosn + 0.3594*cos2p + 0.0664*np.cos((2.0*p-n)*dtr)
        elif c in ('chi1',) and OTIS_TYPE:
            term1 = -0.221*sinn
            term2 = 1.0 + 0.221*cosn
        elif c in ('chi1','theta1','j1') and FES_TYPE:
            f[:,i] = np.sin(2.0*II) / 0.7214
            u[:,i] = -nu
            continue
        elif c in ('chi1',):
            term1 = -0.250*sinn
            term2 = 1.0 + 0.193*cosn
        elif c in ('p1',):
            term1 = -0.0112*sinn
            term2 = 1.0 - 0.0112*cosn
        elif c in ('k1','sk3','2sk5') and OTIS_TYPE:
            term1 = -0.1554*sinn + 0.0029*sin2n
            term2 = 1.0 + 0.1158*cosn - 0.0029*cos2n
        elif c in ('k1','sk3','2sk5') and FES_TYPE:
            temp1 = 0.8965*np.power(np.sin(2.0*II),2.0)
            temp2 = 0.6001*np.sin(2.0*II)*np.cos(nu)
            f[:,i] = np.sqrt(temp1 + temp2 + 0.1006)
            u[:,i] = -nu_prime
            continue
        elif c in ('k1',) and PERTH3_TYPE:
            f[:,i] = 1.006 + 0.115*cosn - 0.009*cos2n
            u[:,i] = dtr*(-8.9*sinn + 0.7*sin2n)
            continue
        elif c in ('k1','sk3','2sk5'):
            term1 = -0.1554*sinn + 0.0031*sin2n
            term2 = 1.0 + 0.1158*cosn - 0.0028*cos2n
        elif c in ('j1','theta1'):
            term1 = -0.227*sinn
            term2 = 1.0 + 0.169*cosn
        elif c in ('oo1','ups1') and OTIS_TYPE:
            term1 = -0.640*sinn - 0.134*sin2n
            term2 = 1.0 + 0.640*cosn + 0.134*cos2n
        elif c in ('oo1','ups1') and FES_TYPE:
            f[:,i] = np.sin(II)*np.power(np.sin(II/2.0),2.0)/0.01640
            u[:,i] = -2.0*xi - nu
            continue
        elif c in ('oo1','ups1'):
            term1 = -0.640*sinn - 0.134*sin2n - 0.150*sin2p
            term2 = 1.0 + 0.640*cosn + 0.134*cos2n + 0.150*cos2p
        elif c in ('m2','2n2','mu2','n2','nu2','lambda2','ms4','eps2','2sm6',
                '2sn6','mp1','mp3','sn4') and FES_TYPE:
            f[:,i] = np.power(np.cos(II/2.0),4.0)/0.9154
            u[:,i] = 2.0*xi - 2.0*nu
            continue
        elif c in ('m2','n2') and PERTH3_TYPE:
            f[:,i] = 1.000 - 0.037*cosn
            u[:,i] = dtr*(-2.1*sinn)
            continue
        elif c in ('m2','2n2','mu2','n2','nu2','lambda2','ms4','eps2','2sm6',
                '2sn6','mp1','mp3','sn4'):
            term1 = -0.03731*sinn + 0.00052*sin2n
            term2 = 1.0 - 0.03731*cosn + 0.00052*cos2n
        elif c in ('l2','sl4') and OTIS_TYPE:
            term1 = -0.25*sin2p - 0.11*np.sin((2.0*p-n)*dtr) - 0.04*sinn
            term2 = 1.0 - 0.25*cos2p - 0.11*np.cos((2.0*p - n)*dtr) - 0.04*cosn
        elif c in ('l2','sl4') and FES_TYPE:
            f[:,i] = Ra1*np.power(np.cos(II/2.0),4.0)/0.9154
            u[:,i] = 2.0*xi - 2.0*nu - R
            continue
        elif c in ('l2','sl4'):
            term1 = -0.25*sin2p - 0.11*np.sin((2.0*p-n)*dtr) - 0.037*sinn
            term2 = 1.0 - 0.25*cos2p - 0.11*np.cos((2.0*p-n)*dtr) - 0.037*cosn
        elif c in ('l2b',):
            # for when l2 is split into two constituents
            term1 = 0.441*sinn
            term2 = 1.0 + 0.441*cosn
        elif c in ('k2','sk4','2sk6','kp1') and OTIS_TYPE:
            term1 = -0.3108*sinn - 0.0324*sin2n
            term2 = 1.0 + 0.2852*cosn + 0.0324*cos2n
        elif c in ('k2','sk4','2sk6','kp1') and FES_TYPE:
            term1 = 19.0444 * np.power(np.sin(II),4.0)
            term2 = 2.7702 * np.power(np.sin(II),2.0) * np.cos(2.0*nu)
            f[:,i] = np.sqrt(term1 + term2 + 0.0981)
            u[:,i] = -2.0*nu_sec
            continue
        elif c in ('k2',) and PERTH3_TYPE:
            f[:,i] = 1.024 + 0.286*cosn + 0.008*cos2n
            u[:,i] = dtr*(-17.7*sinn + 0.7*sin2n)
            continue
        elif c in ('k2','sk4','2sk6','kp1'):
            term1 = -0.3108*sinn - 0.0324*sin2n
            term2 = 1.0 + 0.2853*cosn + 0.0324*cos2n
        elif c in ('gamma2',):
            term1 = 0.147*np.sin(2.0*(n-p)*dtr)
            term2 = 1.0 + 0.147*np.cos(2.0*(n-p)*dtr)
        elif c in ('delta2',):
            term1 = 0.505*sin2p + 0.505*sinn - 0.165*sin2n
            term2 = 1.0 - 0.505*cos2p - 0.505*cosn + 0.165*cos2n
        elif c in ('eta2','zeta2') and FES_TYPE:
            f[:,i] = np.power(np.sin(II),2.0)/0.1565
            u[:,i] = -2.0*nu
            continue
        elif c in ('eta2','zeta2'):
            term1 = -0.436*sinn
            term2 = 1.0 + 0.436*cosn
        elif c in ('s2',):
            term1 = 0.00225*sinn
            term2 = 1.0 + 0.00225*cosn
        elif c in ("m1'",):
            # Linear 3rd degree terms
            term1 = -0.01815*sinn
            term2 = 1.0 - 0.27837*cosn
        elif c in ("q1'",):
            # Linear 3rd degree terms
            term1 = 0.3915*sinn + 0.033*sin2n + 0.061*sin2p
            term2 = 1.0 + 0.3915*cosn + 0.033*cos2n + 0.06*cos2p
        elif c in ("j1'",):
            # Linear 3rd degree terms
            term1 = -0.438*sinn - 0.033*sin2n
            term2 = 1.0 + 0.372*cosn + 0.033*cos2n
        elif c in ("2n2'",):
            # Linear 3rd degree terms
            term1 = 0.166*sinn
            term2 = 1.0 + 0.166*cosn
        elif c in ("n2'",):
            # Linear 3rd degree terms
            term1 = 0.1705*sinn - 0.0035*sin2n - 0.0176*sin2p
            term2 = 1.0 + 0.1705*cosn - 0.0035*cos2n - 0.0176*cos2p
        elif c in ("l2'"):
            # Linear 3rd degree terms
            term1 = -0.2495*sinn
            term2 = 1.0 + 0.1315*cosn
        elif c in ('m3',) and FES_TYPE:
            f[:,i] = np.power(np.cos(II/2.0), 6.0) / 0.8758
            u[:,i] = (3.0*xi - 3.0*nu)
            continue
        elif c in ('m3','e3'):
            # Linear 3rd degree terms
            term1 = -0.05644*sinn
            term2 = 1.0 - 0.05644*cosn
        elif c in ('j3','f3'):
            term1 = -0.464*sinn - 0.052*sin2n
            term2 = 1.0 + 0.387*cosn + 0.052*cos2n
        elif c in ('l3',):
            term1 = -0.373*sin2p - 0.164*np.sin((2.0*p-n)*dtr)
            term2 = 1.0 - 0.373*cos2p - 0.164*np.cos((2.0*p-n)*dtr)
        elif c in ('mfdw',):
            # special test of Doodson-Warburg formula
            f[:,i] = 1.043 + 0.414*cosn
            u[:,i] = dtr*(-23.7*sinn + 2.7*sin2n - 0.4*sin3n)
            continue
        elif c in ('so1','2so3','2po1'):
            # compound tides calculated using recursion
            parents = ['o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]
            u[:,i] = -utmp[:,0]
            continue
        elif c in ('o3',):
            # compound tides calculated using recursion
            parents = ['o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3
            u[:,i] = 3.0*utmp[:,0]
            continue
        elif c in ('2k2'):
            # compound tides calculated using recursion
            parents = ['k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2
            u[:,i] = 2.0*utmp[:,0]
            continue
        elif c in ('tk1'):
            # compound tides calculated using recursion
            parents = ['k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]
            u[:,i] = -utmp[:,0]
            continue
        elif c in ('2oop1'):
            # compound tides calculated using recursion
            parents = ['oo1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2
            u[:,i] = 2.0*utmp[:,0]
            continue
        elif c in ('oq2'):
            # compound tides calculated using recursion
            parents = ['o1','q1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('2oq1'):
            # compound tides calculated using recursion
            parents = ['o1','q1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('ko2'):
            # compound tides calculated using recursion
            parents = ['o1','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('opk1',):
            # compound tides calculated using recursion
            parents = ['o1','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] - utmp[:,1]
            continue
        elif c in ('2ook1',):
            # compound tides calculated using recursion
            parents = ['oo1','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('kj2',):
            # compound tides calculated using recursion
            parents = ['k1','j1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('kjq1'):
            # compound tides calculated using recursion
            parents = ['k1','j1','q1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1] * ftmp[:,2]
            u[:,i] = utmp[:,0] + utmp[:,1] - utmp[:,2]
            continue
        elif c in ('k3',):
            # compound tides calculated using recursion
            parents = ['k1','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in('m4','mn4','mns2','2ms2','mnus2','mmus2','2ns2','n4','mnu4',
                'mmu4','2mt6','2ms6','msn6','mns6','2mr6','msmu6','2mp3','2ms3',
                '2mp5','2msp7','2(ms)8','2ms8'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2
            u[:,i] = 2.0*utmp[:,0]
            continue
        elif c in ('msn2','snm2','nsm2'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2
            u[:,i] = 0.0
            continue
        elif c in ('mmun2','2mn2'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3
            u[:,i] = utmp[:,0]
            continue
        elif c in ('2sm2',):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]
            u[:,i] = -utmp[:,0]
            continue
        elif c in ('m6','2mn6','2mnu6','2mmu6','2nm6','mnnu6','mnmu6','3ms8',
                '3mp7','2msn8','3ms5','3mp5','3ms4','3m2s2','3m2s10','2mn2s2'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3
            u[:,i] = 3.0*utmp[:,0]
            continue
        elif c in ('m8','ma8','3mn8','3mnu8','3mmu8','2mn8','2(mn):8','3msn10',
                '4ms10','2(mn)S10','4m2s12'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**4
            u[:,i] = 4.0*utmp[:,0]
            continue
        elif c in ('m10','4mn10','5ms12','4msn12','4mns12'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**5
            u[:,i] = 5.0*utmp[:,0]
            continue
        elif c in ('m12','5mn12','6ms14','5msn14'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**6
            u[:,i] = 6.0*utmp[:,0]
            continue
        elif c in ('m14',):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**7
            u[:,i] = 7.0*utmp[:,0]
            continue
        elif c in ('mo3','no3','mso5'):
            # compound tides calculated using recursion
            parents = ['m2','o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('no1','nso3'):
            # compound tides calculated using recursion
            parents = ['m2','o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] - utmp[:,1]
            continue
        elif c in ('mq3','nq3'):
            # compound tides calculated using recursion
            parents = ['m2','q1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('2mq3',):
            # compound tides calculated using recursion
            parents = ['m2','q1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('2no3',):
            # compound tides calculated using recursion
            parents = ['m2','o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('2mo5','2no5','mno5','2mso7','2(ms):o9'):
            # compound tides calculated using recursion
            parents = ['m2','o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('2mno7','3mo7'):
            # compound tides calculated using recursion
            parents = ['m2','o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3 * ftmp[:,1]
            u[:,i] = 3.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('mk3','nk3','msk5','nsk5'):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('mnk5','2mk5','2nk5','2msk7'):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('2mk3',):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('3mk7','2mnk7','2nmk7','3nk7','3msk9'):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3 * ftmp[:,1]
            u[:,i] = 3.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('3msk7',):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3 * ftmp[:,1]
            u[:,i] = 3.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('4mk9','3mnk9','2m2nk9','2(mn):k9','3nmk9','4msk11'):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**4 * ftmp[:,1]
            u[:,i] = 4.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('3km5',):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]**3
            u[:,i] = utmp[:,0] + 3.0*utmp[:,1]
            continue
        elif c in ('mk4','nk4','mks2'):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('msk2','2smk4','msk6','snk6'):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] - utmp[:,1]
            continue
        elif c in ('mnk6','2mk6','2msk8','msnk8'):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('mnk2','2mk2'):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('mkn2','nkm2'):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = utmp[:,1]
            continue
        elif c in ('skm2',):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = -utmp[:,0] + utmp[:,1]
            continue
        elif c in ('3mk8','2mnk8'):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3 * ftmp[:,1]
            u[:,i] = 3.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('m2(ks):2',):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]**2
            u[:,i] = utmp[:,0] + 2.0*utmp[:,1]
            continue
        elif c in ('2ms2k2',):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]**2
            u[:,i] = 2.0*utmp[:,0] - 2.0*utmp[:,1]
            continue
        elif c in ('mko5','msko7'):
            # compound tides calculated using recursion
            parents = ['m2','k2','o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1] * ftmp[:,2]
            u[:,i] = utmp[:,0] + utmp[:,1] + utmp[:,2]
            continue
        elif c in ('ml4','msl6'):
            # compound tides calculated using recursion
            parents = ['m2','l2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('2ml2',):
            # compound tides calculated using recursion
            parents = ['m2','l2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('2ml6','2ml2s2','2mls4','2msl8'):
            # compound tides calculated using recursion
            parents = ['m2','l2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('2nmls6','3mls6','2mnls6','3ml8','2mnl8','3msl10'):
            # compound tides calculated using recursion
            parents = ['m2','l2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3 * ftmp[:,1]
            u[:,i] = 3.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('4msl12',):
            # compound tides calculated using recursion
            parents = ['m2','l2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**4 * ftmp[:,1]
            u[:,i] = 4.0*utmp[:,0] + utmp[:,1]
            continue
        else:
            # default for linear tides
            term1 = 0.0
            term2 = 1.0

        # calculate factors for linear tides
        # and parent waves in compound tides
        f[:,i] = np.sqrt(term1**2 + term2**2)
        u[:,i] = np.arctan2(term1, term2)

    # return corrections for constituents
    return (u, f)

def _arguments_table(**kwargs):
    """
    Arguments table for tidal constituents [1]_ [2]_

    Parameters
    ----------
    corrections: str, default 'OTIS'
        use arguments from OTIS, FES or GOT models

    Returns
    -------
    coef: np.ndarray
        Doodson coefficients (Cartwright numbers) for each constituent

    References
    ----------
    .. [1] A. T. Doodson and H. Lamb, "The harmonic development of
        the tide-generating potential", *Proceedings of the Royal Society
        of London. Series A, Containing Papers of a Mathematical and
        Physical Character*, 100(704), 305--329, (1921).
        `doi: 10.1098/rspa.1921.0088 <https://doi.org/10.1098/rspa.1921.0088>`_
    .. [2] A. T. Doodson and H. D. Warburg, "Admiralty Manual of Tides",
        HMSO, London, (1941).
    """
    # set default keyword arguments
    kwargs.setdefault('corrections', 'OTIS')

    # constituents array (not all are included in tidal program)
    cindex = ['sa', 'ssa', 'mm', 'msf', 'mf', 'mt', 'alpha1', '2q1', 'sigma1',
        'q1', 'rho1', 'o1', 'tau1', 'm1', 'chi1', 'pi1', 'p1', 's1', 'k1',
        'psi1', 'phi1', 'theta1', 'j1', 'oo1', '2n2', 'mu2', 'n2', 'nu2', 'm2a',
        'm2', 'm2b', 'lambda2', 'l2', 't2', 's2', 'r2', 'k2', 'eta2', 'mns2',
        '2sm2', 'm3', 'mk3', 's3', 'mn4', 'm4', 'ms4', 'mk4', 's4', 's5', 'm6',
        's6', 's7', 's8', 'm8', 'mks2', 'msqm', 'mtm', 'n4', 'eps2', 'z0']
    # modified Doodson coefficients for constituents
    # using 7 index variables: tau, s, h, p, n, pp, k
    # tau: mean lunar time
    # s: mean longitude of moon
    # h: mean longitude of sun
    # p: mean longitude of lunar perigee
    # n: mean longitude of ascending lunar node
    # pp: mean longitude of solar perigee
    # k: 90-degree phase
    coef = coefficients_table(cindex, **kwargs)
    # return the coefficient table
    return coef

def _minor_table(**kwargs):
    """
    Arguments table for minor tidal constituents [1]_ [2]_

    Returns
    -------
    coef: np.ndarray
        Doodson coefficients (Cartwright numbers) for each constituent

    References
    ----------
    .. [1] A. T. Doodson and H. Lamb, "The harmonic development of
        the tide-generating potential", *Proceedings of the Royal Society
        of London. Series A, Containing Papers of a Mathematical and
        Physical Character*, 100(704), 305--329, (1921).
        `doi: 10.1098/rspa.1921.0088 <https://doi.org/10.1098/rspa.1921.0088>`_
    .. [2] A. T. Doodson and H. D. Warburg, "Admiralty Manual of Tides",
        HMSO, London, (1941).
    """
    # modified Doodson coefficients for constituents
    # using 7 index variables: tau, s, h, p, n, pp, k
    # tau: mean lunar time
    # s: mean longitude of moon
    # h: mean longitude of sun
    # p: mean longitude of lunar perigee
    # n: mean longitude of ascending lunar node
    # pp: mean longitude of solar perigee
    # k: 90-degree phase
    minor = ['2q1', 'sigma1', 'rho1', 'm1b', 'm1a', 'chi1', 'pi1',
        'phi1', 'theta1', 'j1', 'oo1', '2n2', 'mu2', 'nu2', 'lambda2',
        'l2a', 'l2b', 't2', 'eps2', 'eta2']
    coef = coefficients_table(minor, **kwargs)
    # return the coefficient table
    return coef

def _constituent_parameters(c: str, **kwargs):
    """
    Loads parameters for a given tidal constituent

    Parameters
    ----------
    c: str
        tidal constituent ID
    raise_error: bool, default False
        Raise exception if constituent is unsupported

    Returns
    -------
    amplitude: float
        amplitude of equilibrium tide for tidal constituent (meters)
    phase: float
        phase of tidal constituent (radians)
    omega: float
        angular frequency of constituent (radians)
    alpha: float
        load love number of tidal constituent
    species: float
        spherical harmonic dependence of quadrupole potential

    References
    ----------
    .. [1] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of
        Barotropic Ocean Tides," *Journal of Atmospheric and Oceanic
        Technology*, 19(2), 183--204, (2002).
        `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

    .. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
    """
    # default keyword arguments
    kwargs.setdefault('raise_error', False)
    # constituents array that are included in tidal program
    cindex = ['m2', 's2', 'k1', 'o1', 'n2', 'p1', 'k2', 'q1', '2n2', 'mu2',
        'nu2', 'l2', 't2', 'j1', 'm1', 'oo1', 'rho1', 'mf', 'mm', 'ssa',
        'm4', 'ms4', 'mn4', 'm6', 'm8', 'mk3', 's6', '2sm2', '2mk3',
        'msf', 'sa', 'mt', '2q1']
    # species type (spherical harmonic dependence of quadrupole potential)
    _species = np.array([2, 2, 1, 1, 2, 1, 2, 1, 2, 2, 2, 2, 2, 1, 1,
        1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    # Load Love numbers
    # alpha = correction factor for first order load tides
    _alpha = np.array([0.693, 0.693, 0.736, 0.695, 0.693, 0.706, 0.693,
        0.695, 0.693, 0.693, 0.693, 0.693, 0.693, 0.695, 0.695, 0.695, 0.695,
        0.693, 0.693, 0.693, 0.693, 0.693, 0.693, 0.693, 0.693, 0.693, 0.693,
        0.693, 0.693, 0.693, 0.693, 0.693, 0.693])
    # omega: angular frequency of constituent, in radians
    _omega = np.array([1.405189e-04, 1.454441e-04, 7.292117e-05, 6.759774e-05,
        1.378797e-04, 7.252295e-05, 1.458423e-04, 6.495854e-05, 1.352405e-04,
        1.355937e-04, 1.382329e-04, 1.431581e-04, 1.452450e-04, 7.556036e-05,
        7.028195e-05, 7.824458e-05, 6.531174e-05, 0.053234e-04, 0.026392e-04,
        0.003982e-04, 2.810377e-04, 2.859630e-04, 2.783984e-04, 4.215566e-04,
        5.620755e-04, 2.134402e-04, 4.363323e-04, 1.503693e-04, 2.081166e-04,
        4.925200e-06, 1.990970e-07, 7.962619e-06, 6.231934e-05])
    # Astronomical arguments (relative to t0 = 1 Jan 0:00 1992)
    # phases for each constituent are referred to the time when the phase of
    # the forcing for that constituent is zero on the Greenwich meridian
    _phase = np.array([1.731557546, 0.000000000, 0.173003674, 1.558553872,
        6.050721243, 6.110181633, 3.487600001, 5.877717569, 4.086699633,
        3.463115091, 5.427136701, 0.553986502, 0.052841931, 2.137025284,
        2.436575100, 1.929046130, 5.254133027, 1.756042456, 1.964021610,
        3.487600001, 3.463115091, 1.731557546, 1.499093481, 5.194672637,
        6.926230184, 1.904561220, 0.000000000, 4.551627762, 3.809122439,
        4.551627762, 6.232786837, 3.720064066, 3.91369596])
    # amplitudes of equilibrium tide in meters
    # _amplitude = np.array([0.242334,0.112743,0.141565,0.100661,0.046397,
    _amplitude = np.array([0.2441, 0.112743, 0.141565, 0.100661, 0.046397,
        0.046848, 0.030684, 0.019273, 0.006141, 0.007408, 0.008811, 0.006931,
        0.006608, 0.007915, 0.007915, 0.004338, 0.003661, 0.042041, 0.022191,
        0.019567, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.003681, 0.003104,
        0.008044, 0.002565])

    # map between input constituent and cindex
    j = [j for j,val in enumerate(cindex) if (val == c.lower())]
    # set the values for the constituent
    if j:
        amplitude, = _amplitude[j]
        phase, = _phase[j]
        omega, = _omega[j]
        alpha, = _alpha[j]
        species, = _species[j]
    elif kwargs['raise_error']:
        raise ValueError(f'Unsupported constituent {c}')
    else:
        amplitude = 0.0
        phase = 0.0
        omega = 0.0
        alpha = 0.0
        species = 0
    # return the values for the constituent
    return (amplitude, phase, omega, alpha, species)

def _to_doodson_number(coef: list | np.ndarray, **kwargs):
    """
    Converts Cartwright numbers into a Doodson number

    Parameters
    ----------
    coef: list or np.ndarray
        Doodson coefficients (Cartwright numbers) for constituent
    raise_error: bool, default True
        Raise exception if constituent is unsupported

    Returns
    -------
    DO: float
        Doodson number for constituent
    """
    # default keyword arguments
    kwargs.setdefault('raise_error', True)
    # assert length and verify array
    coef = np.array(coef[:6])
    # add 5 to values following Doodson convention (prevent negatives)
    coef[1:] += 5
    # check for unsupported constituents
    if (np.any(coef < 0) or np.any(coef > 10)) and kwargs['raise_error']:
        raise ValueError('Unsupported constituent')
    elif (np.any(coef < 0) or np.any(coef > 10)):
        return None
    else:
        # convert to single number and round off floating point errors
        DO = np.sum([v*10**(2-o) for o,v in enumerate(coef)])
        return np.round(DO, decimals=3)

def _from_doodson_number(DO: float | np.ndarray):
    """
    Converts Doodson numbers into Cartwright numbers

    Parameters
    ----------
    DO: float or np.ndarray
        Doodson number for constituent

    Returns
    -------
    coef: np.ndarray
        Doodson coefficients (Cartwright numbers) for constituent
    """
    # convert from Doodson number to Cartwright numbers
    # multiply by 1000 to prevent floating point errors
    coef = np.array([np.mod(1e3*DO, 10**(6-o))//10**(5-o) for o in range(6)])
    # remove 5 from values following Doodson convention
    coef[1:] -= 5
    return coef
