#!/usr/bin/env python
u"""
arguments.py
Written by Tyler Sutterley (01/2024)
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
        use nodal corrections from OTIS/ATLAS or GOT models
    M1: str, default 'Ray'
        coefficients to use for M1 tides

                - ``'Doodson'``
                - ``'Ray'``

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
    kwargs.setdefault('M1', 'Ray')

    # constituents array (not all are included in tidal program)
    cindex = ['sa', 'ssa', 'mm', 'msf', 'mf', 'mt', 'alpha1', '2q1', 'sigma1',
        'q1', 'rho1', 'o1', 'tau1', 'm1', 'chi1', 'pi1', 'p1', 's1', 'k1',
        'psi1', 'phi1', 'theta1', 'j1', 'oo1', '2n2', 'mu2', 'n2', 'nu2', 'm2a',
        'm2', 'm2b', 'lambda2', 'l2', 't2', 's2', 'r2', 'k2', 'eta2', 'mns2',
        '2sm2', 'm3', 'mk3', 's3', 'mn4', 'm4', 'ms4', 'mk4', 's4', 's5', 'm6',
        's6', 's7', 's8', 'm8', 'mks2', 'msqm', 'mtm', 'n4', 'eps2', 'z0']

    # set function for astronomical longitudes
    ASTRO5 = True if kwargs['corrections'] in ('GOT', 'FES') else False
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
    # degrees to radians
    dtr = np.pi/180.0

    # determine equilibrium arguments
    fargs = np.c_[tau, s, h, p, n, pp, k]
    arg = np.dot(fargs, _arguments_table(**kwargs))

    # determine nodal corrections f and u
    sinn = np.sin(n*dtr)
    cosn = np.cos(n*dtr)
    sin2n = np.sin(2.0*n*dtr)
    cos2n = np.cos(2.0*n*dtr)
    sin3n = np.sin(3.0*n*dtr)

    # set nodal corrections
    # nodal factor correction
    f = np.zeros((nt, 60))
    # nodal angle correction
    u = np.zeros((nt, 60))
    # determine nodal corrections f and u for each model type
    if kwargs['corrections'] in ('OTIS', 'ATLAS', 'TMD3', 'netcdf'):
        # nodal factors
        f[:,0] = 1.0 # Sa
        f[:,1] = 1.0 # Ssa
        f[:,2] = 1.0 - 0.130*cosn # Mm
        f[:,3] = 1.0 # MSf
        f[:,4] = 1.043 + 0.414*cosn # Mf
        temp1 = (1.0 + 0.203*cosn + 0.040*cos2n)**2
        temp2 = (0.203*sinn + 0.040*sin2n)**2
        f[:,5] = np.sqrt(temp1 + temp2) # Mt
        f[:,6] = 1.0 # alpha1
        f[:,7] = np.sqrt((1.0 + 0.188*cosn)**2 + (0.188*sinn)**2) # 2Q1
        f[:,8] = f[:,7] # sigma1
        f[:,9] = f[:,7] # q1
        f[:,10] = f[:,7] # rho1
        temp1 = (1.0 + 0.189*cosn - 0.0058*cos2n)**2
        temp2 = (0.189*sinn - 0.0058*sin2n)**2
        f[:,11] = np.sqrt(temp1 + temp2) # O1
        f[:,12] = 1.0 # tau1
        if (kwargs['M1'] == 'Doodson'):
            # A. T. Doodson's coefficients for M1 tides
            Mtmp1 = 2.0*np.cos(p*dtr) + 0.4*np.cos((p-n)*dtr)
            Mtmp2 = np.sin(p*dtr) + 0.2*np.sin((p-n)*dtr)
        elif (kwargs['M1'] == 'Ray'):
            # R. Ray's coefficients for M1 tides
            Mtmp1 = 1.36*np.cos(p*dtr) + 0.267*np.cos((p-n)*dtr)
            Mtmp2 = 0.64*np.sin(p*dtr) + 0.135*np.sin((p-n)*dtr)
        f[:,13] = np.sqrt(Mtmp1**2 + Mtmp2**2) # M1
        f[:,14] = np.sqrt((1.0+0.221*cosn)**2+(0.221*sinn)**2) # chi1
        f[:,15] = 1.0 # pi1
        f[:,16] = 1.0 # P1
        f[:,17] = 1.0 # S1
        temp1 = (1.0 + 0.1158*cosn - 0.0029*cos2n)**2
        temp2 = (0.1554*sinn - 0.0029*sin2n)**2
        f[:,18] = np.sqrt(temp1 + temp2) # K1
        f[:,19] = 1.0 # psi1
        f[:,20] = 1.0 # phi1
        f[:,21] = 1.0 # theta1
        f[:,22] = np.sqrt((1.0+0.169*cosn)**2 + (0.227*sinn)**2) # J1
        temp1 = (1.0 + 0.640*cosn + 0.134*cos2n)**2
        temp2 = (0.640*sinn + 0.134*sin2n)**2
        f[:,23] = np.sqrt(temp1 + temp2) # OO1
        temp1 = (1.0 - 0.03731*cosn + 0.00052*cos2n)**2
        temp2 = (0.03731*sinn - 0.00052*sin2n)**2
        f[:,24] = np.sqrt(temp1 + temp2) # 2N2
        f[:,25] = f[:,24] # mu2
        f[:,26] = f[:,24] # N2
        f[:,27] = f[:,24] # nu2
        f[:,28] = 1.0 # M2a
        f[:,29] = f[:,24] # M2
        f[:,30] = 1.0 # M2b
        f[:,31] = 1.0 # lambda2
        Ltmp1 = 1.0 - 0.25*np.cos(2*p*dtr) - \
            0.11*np.cos((2.0*p - n)*dtr) - 0.04*cosn
        Ltmp2 = 0.25*np.sin(2*p*dtr) + \
            0.11*np.sin((2.0*p - n)*dtr) + 0.04*sinn
        f[:,32] = np.sqrt(Ltmp1**2 + Ltmp2**2) # L2
        f[:,33] = 1.0 # T2
        f[:,34] = 1.0 # S2
        f[:,35] = 1.0 # R2
        temp1 = (1.0 + 0.2852*cosn + 0.0324*cos2n)**2
        temp2 = (0.3108*sinn + 0.0324*sin2n)**2
        f[:,36] = np.sqrt(temp1 + temp2) # K2
        f[:,37] = np.sqrt((1.0 + 0.436*cosn)**2 + (0.436*sinn)**2) # eta2
        f[:,38] = f[:,29]**2 # MNS2
        f[:,39] = f[:,29] # 2SM2
        f[:,40] = 1.0 # M3 (wrong)
        f[:,41] = f[:,18]*f[:,29] # MK3
        f[:,42] = 1.0 # S3
        f[:,43] = f[:,29]**2 # MN4
        f[:,44] = f[:,43] # M4
        f[:,45] = f[:,43] # MS4
        f[:,46] = f[:,29]*f[:,36] # MK4
        f[:,47] = 1.0 # S4
        f[:,48] = 1.0 # S5
        f[:,49] = f[:,29]**3 # M6
        f[:,50] = 1.0 # S6
        f[:,51] = 1.0 # S7
        f[:,52] = 1.0 # S8
        # shallow water constituents
        f[:,53] = f[:,29]**4 # m8
        f[:,54] = f[:,29]*f[:,36] # mks2
        f[:,55] = f[:,4] # msqm
        f[:,56] = f[:,4] # mtm
        f[:,57] = f[:,29]**2 # n4
        f[:,58] = f[:,29] # eps2
        # mean sea level
        f[:,59] = 1.0 # Z0

        # nodal angles
        u[:,0] = 0.0 # Sa
        u[:,1] = 0.0 # Ssa
        u[:,2] = 0.0 # Mm
        u[:,3] = 0.0 # MSf
        u[:,4] = -23.7*sinn + 2.7*sin2n - 0.4*sin3n # Mf
        temp1 = -(0.203*sinn + 0.040*sin2n)
        temp2 = (1.0 + 0.203*cosn + 0.040*cos2n)
        u[:,5] = np.arctan(temp1/temp2)/dtr # Mt
        u[:,6] = 0.0 # alpha1
        u[:,7] = np.arctan(0.189*sinn/(1.0 + 0.189*cosn))/dtr # 2Q1
        u[:,8] = u[:,7] # sigma1
        u[:,9] = u[:,7] # q1
        u[:,10] = u[:,7] # rho1
        u[:,11] = 10.8*sinn - 1.3*sin2n + 0.2*sin3n # O1
        u[:,12] = 0.0 # tau1
        u[:,13] = np.arctan2(Mtmp2,Mtmp1)/dtr # M1
        u[:,14] = np.arctan(-0.221*sinn/(1.0+0.221*cosn))/dtr # chi1
        u[:,15] = 0.0 # pi1
        u[:,16] = 0.0 # P1
        u[:,17] = 0.0 # S1
        temp1 = (-0.1554*sinn + 0.0029*sin2n)
        temp2 = (1.0 + 0.1158*cosn - 0.0029*cos2n)
        u[:,18] = np.arctan(temp1/temp2)/dtr # K1
        u[:,19] = 0.0 # psi1
        u[:,20] = 0.0 # phi1
        u[:,21] = 0.0 # theta1
        u[:,22] = np.arctan(-0.227*sinn/(1.0+0.169*cosn))/dtr # J1
        temp1 = -(0.640*sinn + 0.134*sin2n)
        temp2 = (1.0 + 0.640*cosn + 0.134*cos2n)
        u[:,23] = np.arctan(temp1/temp2)/dtr # OO1
        temp1 = (-0.03731*sinn + 0.00052*sin2n)
        temp2 = (1.0 - 0.03731*cosn + 0.00052*cos2n)
        u[:,24] = np.arctan(temp1/temp2)/dtr # 2N2
        u[:,25] = u[:,24] # mu2
        u[:,26] = u[:,24] # N2
        u[:,27] = u[:,24] # nu2
        u[:,28] = 0.0 # M2a
        u[:,29] = u[:,24] # M2
        u[:,30] = 0.0 # M2b
        u[:,31] = 0.0 # lambda2
        u[:,32] = np.arctan(-Ltmp2/Ltmp1)/dtr # L2
        u[:,33] = 0.0 # T2
        u[:,34] = 0.0 # S2
        u[:,35] = 0.0 # R2
        temp1 = -(0.3108*sinn+0.0324*sin2n)
        temp2 = (1.0 + 0.2852*cosn + 0.0324*cos2n)
        u[:,36] = np.arctan(temp1/temp2)/dtr # K2
        u[:,37] = np.arctan(-0.436*sinn/(1.0 + 0.436*cosn))/dtr # eta2
        u[:,38] = u[:,29]*2.0 # MNS2
        u[:,39] = u[:,29] # 2SM2
        u[:,40] = 1.50*u[:,29] # M3
        u[:,41] = u[:,29] + u[:,18] # MK3
        u[:,42] = 0.0 # S3
        u[:,43] = 2.0*u[:,29] # MN4
        u[:,44] = u[:,43] # M4
        u[:,45] = u[:,29] # MS4
        u[:,46] = u[:,29] + u[:,36] # MK4
        u[:,47] = 0.0 # S4
        u[:,48] = 0.0 # S5
        u[:,49] = 3.0*u[:,29] # M6
        u[:,50] = 0.0 # S6
        u[:,51] = 0.0 # S7
        u[:,52] = 0.0 # S8
        # mean sea level
        u[:,59] = 0.0 # Z0

    elif kwargs['corrections'] in ('FES',):
        # additional astronomical terms for FES models
        II = np.arccos(0.913694997 - 0.035692561*np.cos(n*dtr))
        at1 = np.arctan(1.01883*np.tan(n*dtr/2.0))
        at2 = np.arctan(0.64412*np.tan(n*dtr/2.0))
        xi = -at1 - at2 + n*dtr
        xi[xi > np.pi] -= 2.0*np.pi
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

        # nodal factors
        f[:,0] = 1.0 # Sa
        f[:,1] = 1.0 # Ssa
        f[:,2] = (2.0/3.0 - np.power(np.sin(II),2.0))/0.5021 # Mm
        f[:,3] = 1.0 # MSf
        f[:,4] = np.power(np.sin(II),2.0)/0.1578  # Mf
        f[:,7] = np.sin(II)*(np.cos(II/2.0)**2)/0.38 # 2Q1
        f[:,8] = f[:,7] # sigma1
        f[:,9] = f[:,7] # q1
        f[:,10] = f[:,7] # rho1
        f[:,11] = f[:,7] # O1
        if (kwargs['M1'] == 'Doodson'):
            # A. T. Doodson's coefficients for M1 tides
            Mtmp1 = 2.0*np.cos(p*dtr) + 0.4*np.cos((p-n)*dtr)
            Mtmp2 = np.sin(p*dtr) + 0.2*np.sin((p-n)*dtr)
        elif (kwargs['M1'] == 'Ray'):
            # R. Ray's coefficients for M1 tides
            Mtmp1 = 1.36*np.cos(p*dtr) + 0.267*np.cos((p-n)*dtr)
            Mtmp2 = 0.64*np.sin(p*dtr) + 0.135*np.sin((p-n)*dtr)
        f[:,13] = np.sqrt(Mtmp1**2 + Mtmp2**2) # M1
        f[:,14] = np.sin(2.0*II) / 0.7214 # chi1
        f[:,15] = 1.0 # pi1
        f[:,16] = 1.0 # P1
        f[:,17] = 1.0 # S1
        temp1 = 0.8965*np.power(np.sin(2.0*II),2.0)
        temp2 = 0.6001*np.sin(2.0*II)*np.cos(nu)
        f[:,18] = np.sqrt(temp1 + temp2 + 0.1006) # K1
        f[:,19] = 1.0 # psi1
        f[:,20] = 1.0 # phi1
        f[:,21] = f[:,14] # theta1
        f[:,22] = f[:,14] # J1
        f[:,23] = np.sin(II)*np.power(np.sin(II/2.0),2.0)/0.01640 # OO1
        f[:,24] = np.power(np.cos(II/2.0),4.0)/0.9154 # 2N2
        f[:,25] = f[:,24] # mu2
        f[:,26] = f[:,24] # N2
        f[:,27] = f[:,24] # nu2
        f[:,28] = 1.0 # M2a
        f[:,29] = f[:,24] # M2
        f[:,30] = 1.0 # M2b
        f[:,31] = f[:,29] # lambda2
        f[:,32] = f[:,29]*Ra1 # L2
        f[:,33] = 1.0 # T2
        f[:,34] = 1.0 # S2
        f[:,35] = 1.0 # R2
        temp1 = 19.0444 * np.power(np.sin(II),4.0)
        temp2 = 2.7702 * np.power(np.sin(II),2.0) * np.cos(2.0*nu)
        f[:,36] = np.sqrt(temp1 + temp2 + 0.0981) # K2
        f[:,37] = np.power(np.sin(II),2.0)/0.1565 # eta2
        f[:,38] = f[:,29]**2 # MNS2
        f[:,39] = f[:,29] # 2SM2
        f[:,40] = np.power(np.cos(II/2.0), 6.0) / 0.8758 # M3
        f[:,41] = f[:,18]*f[:,29] # MK3
        f[:,42] = 1.0 # S3
        f[:,43] = f[:,29]**2 # MN4
        f[:,44] = f[:,43] # M4
        f[:,45] = f[:,29] # MS4
        f[:,46] = f[:,29]*f[:,36] # MK4
        f[:,47] = 1.0 # S4
        f[:,48] = 1.0 # S5
        f[:,49] = f[:,29]**3 # M6
        f[:,50] = 1.0 # S6
        f[:,51] = 1.0 # S7
        f[:,52] = 1.0 # S8
        # shallow water constituents
        f[:,53] = f[:,29]**4 # m8
        f[:,54] = f[:,29]*f[:,36] # mks2
        f[:,55] = f[:,4] # msqm
        f[:,56] = f[:,4] # mtm
        f[:,57] = f[:,29]**2 # n4
        f[:,58] = f[:,29] # eps2
        # mean sea level
        f[:,59] = 1.0 # Z0

        # nodal angles
        u[:,0] = 0.0 # Sa
        u[:,1] = 0.0 # Ssa
        u[:,2] = 0.0 # Mm
        u[:,3] = (2.0*xi - 2.0*nu)/dtr # MSf
        u[:,4] = -2.0*xi/dtr # Mf
        u[:,7] = (2.0*xi - nu)/dtr # 2Q1
        u[:,8] = u[:,7] # sigma1
        u[:,9] = u[:,7] # q1
        u[:,10] = u[:,7] # rho1
        u[:,11] = u[:,7] # O1
        u[:,13] = np.arctan2(Mtmp2,Mtmp1)/dtr # M1
        u[:,14] = -nu/dtr # chi1
        u[:,15] = 0.0 # pi1
        u[:,16] = 0.0 # P1
        u[:,17] = 0.0 # S1
        u[:,18] = -nu_prime/dtr # K1
        u[:,19] = 0.0 # psi1
        u[:,20] = 0.0 # phi1
        u[:,21] = -nu/dtr # theta1
        u[:,22] = u[:,21] # J1
        u[:,23] = (-2.0*xi - nu)/dtr # OO1
        u[:,24] = (2.0*xi - 2.0*nu)/dtr # 2N2
        u[:,25] = u[:,24] # mu2
        u[:,26] = u[:,24] # N2
        u[:,27] = u[:,24] # nu2
        u[:,29] = u[:,24] # M2
        u[:,31] = (2.0*xi - 2.0*nu)/dtr # lambda2
        u[:,32] = (2.0*xi - 2.0*nu - R)/dtr # L2
        u[:,33] = 0.0 # T2
        u[:,34] = 0.0 # S2
        u[:,35] = 0.0 # R2
        u[:,36] = -2.0*nu_sec/dtr # K2
        u[:,37] = -2.0*nu/dtr # eta2
        u[:,38] = (4.0*xi - 4.0*nu)/dtr # mns2
        u[:,39] = (2.0*xi - 2.0*nu)/dtr # 2SM2
        u[:,40] = (3.0*xi - 3.0*nu)/dtr # M3
        u[:,41] = (2.0*xi - 2.0*nu - 2.0*nu_prime)/dtr # MK3
        u[:,42] = 0.0 # S3
        u[:,43] = (4.0*xi - 4.0*nu)/dtr # MN4
        u[:,44] = (4.0*xi - 4.0*nu)/dtr # M4
        u[:,45] = (2.0*xi - 2.0*nu)/dtr  # MS4
        u[:,46] = (2.0*xi - 2.0*nu - 2.0*nu_sec)/dtr # MK4
        u[:,47] = 0.0 # S4
        u[:,48] = 0.0 # S5
        u[:,49] = (6.0*xi - 6.0*nu)/dtr # M6
        u[:,50] = 0.0 # S6
        u[:,51] = 0.0 # S7
        u[:,52] = 0.0 # S8
        # shallow water constituents
        u[:,53] = (8.0*xi - 8.0*nu)/dtr # m8
        u[:,54] = (2.0*xi - 2.0*nu - 2.0*nu_sec)/dtr # mks2
        u[:,55] = u[:,4] # msqm
        u[:,56] = u[:,4] # mtm
        u[:,57] = (4.0*xi - 4.0*nu)/dtr # n4
        u[:,58] = u[:,29] # eps2
        # mean sea level
        u[:,59] = 0.0 # Z0

    elif kwargs['corrections'] in ('GOT',):
        # nodal factors
        f[:,9] = 1.009 + 0.187*cosn - 0.015*cos2n# Q1
        f[:,11] = f[:,9]# O1
        f[:,16] = 1.0 # P1
        f[:,17] = 1.0 # S1
        f[:,18] = 1.006 + 0.115*cosn - 0.009*cos2n# K1
        f[:,26] = 1.000 - 0.037*cosn# N2
        f[:,29] = f[:,26]# M2
        f[:,34] = 1.0 # S2
        f[:,36] = 1.024 + 0.286*cosn + 0.008*cos2n# K2
        f[:,44] = f[:,29]**2# M4

        # nodal angles
        u[:,9] = 10.8*sinn - 1.3*sin2n# Q1
        u[:,11] = u[:,9]# O1
        u[:,16] = 0.0 # P1
        u[:,17] = 0.0 # S1
        u[:,18] = -8.9*sinn + 0.7*sin2n# K1
        u[:,26] = -2.1*sinn# N2
        u[:,29] = u[:,26]# M2
        u[:,34] = 0.0 # S2
        u[:,36] = -17.7*sinn + 0.7*sin2n# K2
        u[:,44] = -4.2*sinn# M4

    # number of constituents of interest
    nc = len(constituents)
    # nodal factor corrections for given constituents
    pu = np.zeros((nt,nc))
    # nodal angle corrections for given constituents
    pf = np.zeros((nt,nc))
    # equilibrium arguments for given constituents
    G = np.zeros((nt,nc))
    for i,c in enumerate(constituents):
        # map between given constituents and supported in tidal program
        assert c.lower() in cindex, f'Unsupported constituent {c.lower()}'
        j, = [j for j,val in enumerate(cindex) if (val == c.lower())]
        pu[:,i] = u[:,j]*dtr
        pf[:,i] = f[:,j]
        G[:,i] = arg[:,j]

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
        use nodal corrections from OTIS/ATLAS or GOT models

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
    ASTRO5 = True if kwargs['corrections'] in ('GOT', 'FES') else False
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
        II = np.arccos(0.913694997 - 0.035692561*np.cos(n*dtr))
        at1 = np.arctan(1.01883*np.tan(n*dtr/2.0))
        at2 = np.arctan(0.64412*np.tan(n*dtr/2.0))
        xi = -at1 - at2 + n*dtr
        xi[xi > np.pi] -= 2.0*np.pi
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
        f[:,9] = f[:,5] # J1
        f[:,10] = np.sin(II)*np.power(np.sin(II/2.0),2.0)/0.01640 # OO1
        f[:,11] = np.power(np.cos(II/2.0),4.0)/0.9154 # 2N2
        f[:,12] = f[:,11] # mu2
        f[:,13] = f[:,11] # nu2
        f[:,14] = f[:,11] # lambda2
        f[:,15] = f[:,11]*Ra1 # L2
        f[:,18] = f[:,11] # eps2
        f[:,19] = np.power(np.sin(II),2.0)/0.1565 # eta2

        # nodal angle corrections for minor constituents
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

    # return values as tuple
    return (u, f, arg)

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
        use arguments from OTIS/ATLAS or GOT models
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

    # constituents array (not all are included in tidal program)
    cindex = ['sa', 'ssa', 'mm', 'msf', 'mf', 'mt', 'alpha1', '2q1', 'sigma1',
        'q1', 'rho1', 'o1', 'tau1', 'm1', 'chi1', 'pi1', 'p1', 's1', 'k1',
        'psi1', 'phi1', 'theta1', 'j1', 'oo1', '2n2', 'mu2', 'n2', 'nu2', 'm2a',
        'm2', 'm2b', 'lambda2', 'l2', 't2', 's2', 'r2', 'k2', 'eta2', 'mns2',
        '2sm2', 'm3', 'mk3', 's3', 'mn4', 'm4', 'ms4', 'mk4', 's4', 's5', 'm6',
        's6', 's7', 's8', 'm8', 'mks2', 'msqm', 'mtm', 'n4', 'eps2', 'z0']
    # get the table of coefficients
    coefficients = _arguments_table(**kwargs)
    if isinstance(constituents, str):
        # check that given constituents is supported in tidal program
        if (constituents.lower() not in cindex) and kwargs['raise_error']:
            raise ValueError(f'Unsupported constituent {constituents}')
        elif (constituents.lower() not in cindex):
            return None
        # map between given constituents and supported in tidal program
        j, = [j for j,val in enumerate(cindex) if (val == constituents.lower())]
        # extract identifier in formalism
        if (kwargs['formalism'] == 'Cartwright'):
            # extract Cartwright number
            numbers = np.array(coefficients[:6,j])
        elif (kwargs['formalism'] == 'Doodson'):
            # convert from coefficients to Doodson number
            numbers = _to_doodson_number(coefficients[:,j], **kwargs)
    else:
        # output dictionary with Doodson numbers
        numbers = {}
        # for each input constituent
        for i,c in enumerate(constituents):
            # check that given constituents is supported in tidal program
            if (c.lower() not in cindex) and kwargs['raise_error']:
                raise ValueError(f'Unsupported constituent {c}')
            elif (c.lower() not in cindex):
                numbers[c] = None
                continue
            # map between given constituents and supported in tidal program
            j, = [j for j,val in enumerate(cindex) if (val == c.lower())]
            # convert from coefficients to Doodson number
            if (kwargs['formalism'] == 'Cartwright'):
                # extract Cartwright number
                numbers[c] = np.array(coefficients[:6,j])
            elif (kwargs['formalism'] == 'Doodson'):
                # convert from coefficients to Doodson number
                numbers[c] = _to_doodson_number(coefficients[:,j], **kwargs)
    # return the Doodson or Cartwright number
    return numbers

def _arguments_table(**kwargs):
    """
    Arguments table for tidal constituents [1]_ [2]_

    Parameters
    ----------
    corrections: str, default 'OTIS'
        use arguments from OTIS/ATLAS or GOT models

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
    coef = np.zeros((7, 60))
    coef[:,0] = [0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0] # Sa
    coef[:,1] = [0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0] # Ssa
    coef[:,2] = [0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0] # Mm
    coef[:,3] = [0.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0] # MSf
    coef[:,4] = [0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0] # Mf
    coef[:,5] = [0.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0] # Mt
    coef[:,6] = [1.0, -4.0, 2.0, 1.0, 0.0, 0.0, -1.0] # alpha1
    coef[:,7] = [1.0, -3.0, 0.0, 2.0, 0.0, 0.0, -1.0] # 2q1
    coef[:,8] = [1.0, -3.0, 2.0, 0.0, 0.0, 0.0, -1.0] # sigma1
    coef[:,9] = [1.0, -2.0, 0.0, 1.0, 0.0, 0.0, -1.0]# q1
    coef[:,10] = [1.0, -2.0, 2.0, -1.0, 0.0, 0.0, -1.0] # rho1
    coef[:,11] = [1.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0] # o1
    coef[:,12] = [1.0, -1.0, 2.0, 0.0, 0.0, 0.0, 1.0] # tau1
    coef[:,13] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0] # m1
    coef[:,14] = [1.0, 0.0, 2.0, -1.0, 0.0, 0.0, 1.0] # chi1
    coef[:,15] = [1.0, 1.0, -3.0, 0.0, 0.0, 1.0, -1.0] # pi1
    coef[:,16] = [1.0, 1.0, -2.0, 0.0, 0.0, 0.0, -1.0] # p1
    if kwargs['corrections'] in ('OTIS', 'ATLAS', 'TMD3', 'netcdf'):
        coef[:,17] = [1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 1.0] # s1
    elif kwargs['corrections'] in ('GOT', 'FES'):
        coef[:,17] = [1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 2.0] # s1 (Doodson's phase)
    coef[:,18] = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0] # k1
    coef[:,19] = [1.0, 1.0, 1.0, 0.0, 0.0, -1.0, 1.0] # psi1
    coef[:,20] = [1.0, 1.0, 2.0, 0.0, 0.0, 0.0, 1.0] # phi1
    coef[:,21] = [1.0, 2.0, -2.0, 1.0, 0.0, 0.0, 1.0] # theta1
    coef[:,22] = [1.0, 2.0, 0.0, -1.0, 0.0, 0.0, 1.0] # j1
    coef[:,23] = [1.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0] # oo1
    coef[:,24] = [2.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0] # 2n2
    coef[:,25] = [2.0, -2.0, 2.0, 0.0, 0.0, 0.0, 0.0] # mu2
    coef[:,26] = [2.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0] # n2
    coef[:,27] = [2.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0] # nu2
    coef[:,28] = [2.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0] # m2a
    coef[:,29] = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # m2
    coef[:,30] = [2.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0] # m2b
    coef[:,31] = [2.0, 1.0, -2.0, 1.0, 0.0, 0.0, 2.0] # lambda2
    coef[:,32] = [2.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0] # l2
    coef[:,33] = [2.0, 2.0, -3.0, 0.0, 0.0, 1.0, 0.0] # t2
    coef[:,34] = [2.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0] # s2
    coef[:,35] = [2.0, 2.0, -1.0, 0.0, 0.0, -1.0, 2.0] # r2
    coef[:,36] = [2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0] # k2
    coef[:,37] = [2.0, 3.0, 0.0, 0.0, 0.0, -1.0, 0.0] # eta2
    coef[:,38] = [2.0, -3.0, 2.0, 1.0, 0.0, 0.0, 0.0] # mns2
    coef[:,39] = [2.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0] # 2sm2
    coef[:,40] = [3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # m3
    coef[:,41] = coef[:,18] + coef[:,29] # mk3
    coef[:,42] = [3.0, 3.0, -3.0, 0.0, 0.0, 0.0, 0.0] # s3
    coef[:,43] = coef[:,26] + coef[:,29] # mn4
    coef[:,44] = [4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # m4
    coef[:,45] = coef[:,29] + coef[:,34] # ms4
    coef[:,46] = coef[:,29] + coef[:,36] # mk4
    coef[:,47] = [4.0, 4.0, -4.0, 0.0, 0.0, 0.0, 0.0] # s4
    coef[:,48] = [5.0, 5.0, -5.0, 0.0, 0.0, 0.0, 0.0] # s5
    coef[:,49] = [6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # m6
    coef[:,50] = [6.0, 6.0, -6.0, 0.0, 0.0, 0.0, 0.0] # s6
    coef[:,51] = [7.0, 7.0, -7.0, 0.0, 0.0, 0.0, 0.0] # s7
    coef[:,52] = [8.0, 8.0, -8.0, 0.0, 0.0, 0.0, 0.0] # s8
    # shallow water constituents
    coef[:,53] = [8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # m8
    coef[:,54] = coef[:,29] + coef[:,36] - coef[:,34] # mks2
    coef[:,55] = [0.0, 4.0, -2.0, 0.0, 0.0, 0.0, 0.0] # msqm
    coef[:,56] = [0.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0] # mtm
    coef[:,57] = [4.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0] # n4
    coef[:,58] = [2.0, -3.0, 2.0, 1.0, 0.0, 0.0, 0.0] # eps2
    # mean sea level
    coef[:,59] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # z0
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
    coef = np.zeros((7, 20))
    coef[:,0] = [1.0, -3.0, 0.0, 2.0, 0.0, 0.0, -1.0] # 2q1
    coef[:,1] = [1.0, -3.0, 2.0, 0.0, 0.0, 0.0, -1.0] # sigma1
    coef[:,2] = [1.0, -2.0, 2.0, -1.0, 0.0, 0.0, -1.0] # rho1
    coef[:,3] = [1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0] # m1
    coef[:,4] = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0] # m1
    coef[:,5] = [1.0, 0.0, 2.0, -1.0, 0.0, 0.0, 1.0] # chi1
    coef[:,6] = [1.0, 1.0, -3.0, 0.0, 0.0, 1.0, -1.0] # pi1
    coef[:,7] = [1.0, 1.0, 2.0, 0.0, 0.0, 0.0, 1.0] # phi1
    coef[:,8] = [1.0, 2.0, -2.0, 1.0, 0.0, 0.0, 1.0] # theta1
    coef[:,9] = [1.0, 2.0, 0.0, -1.0, 0.0, 0.0, 1.0] # j1
    coef[:,10] = [1.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0] # oo1
    coef[:,11] = [2.0, -2.0, 0.0, 2.0, 0.0, 0.0, 0.0] # 2n2
    coef[:,12] = [2.0, -2.0, 2.0, 0.0, 0.0, 0.0, 0.0] # mu2
    coef[:,13] = [2.0, -1.0, 2.0, -1.0, 0.0, 0.0, 0.0] # nu2
    coef[:,14] = [2.0, 1.0, -2.0, 1.0, 0.0, 0.0, 2.0]# lambda2
    coef[:,15] = [2.0, 1.0, 0.0, -1.0, 0.0, 0.0, 2.0] # l2
    coef[:,16] = [2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0] # l2
    coef[:,17] = [2.0, 2.0, -3.0, 0.0, 0.0, 1.0, 0.0] # t2
    coef[:,18] = [2.0, -3.0, 2.0, 1.0, 0.0, 0.0, 0.0] # eps2
    coef[:,19] = [2.0, 3.0, 0.0, 0.0, 0.0, -1.0, 0.0] # eta2
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
