#!/usr/bin/env python
u"""
test_arguments.py (01/2024)
Verify arguments table matches prior arguments array
"""
import pytest
import numpy as np
import pyTMD.astro
from pyTMD.arguments import _arguments_table

@pytest.mark.parametrize("MJD", np.random.randint(58000,61000,size=4))
@pytest.mark.parametrize("corrections", ['OTIS', 'GOT'])
def test_arguments(MJD, corrections):
    """
    Tests the calculation of nodal corrections for tidal constituents

    Parameters
    ----------
    MJD: np.ndarray
        modified julian day of input date
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS/ATLAS or GOT models
    """

    # set function for astronomical longitudes
    ASTRO5 = True if corrections in ('GOT', 'FES') else False
    # convert from Modified Julian Dates into Ephemeris Time
    s, h, p, omega, pp = pyTMD.astro.mean_longitudes(MJD,
        ASTRO5=ASTRO5)
    # number of temporal values
    nt = len(np.atleast_1d(MJD))
    # initial time conversions
    hour = 24.0*np.mod(MJD, 1)
    # convert from hours into degrees
    t1 = 15.0*hour
    t2 = 30.0*hour
    # variable for multiples of 90 degrees (Ray technical note 2017)
    k = 90.0 + np.zeros((nt))

    # Determine equilibrium arguments
    arg = np.zeros((nt, 60))
    arg[:,0] = h - pp # Sa
    arg[:,1] = 2.0*h # Ssa
    arg[:,2] = s - p # Mm
    arg[:,3] = 2.0*s - 2.0*h # MSf
    arg[:,4] = 2.0*s # Mf
    arg[:,5] = 3.0*s - p # Mt
    arg[:,6] = t1 - 5.0*s + 3.0*h + p - 90.0 # alpha1
    arg[:,7] = t1 - 4.0*s + h + 2.0*p - 90.0 # 2Q1
    arg[:,8] = t1 - 4.0*s + 3.0*h - 90.0 # sigma1
    arg[:,9] = t1 - 3.0*s + h + p - 90.0 # q1
    arg[:,10] = t1 - 3.0*s + 3.0*h - p - 90.0 # rho1
    arg[:,11] = t1 - 2.0*s + h - 90.0 # o1
    arg[:,12] = t1 - 2.0*s + 3.0*h + 90.0 # tau1
    arg[:,13] = t1 - s + h + 90.0 # M1
    arg[:,14] = t1 - s + 3.0*h - p + 90.0 # chi1
    arg[:,15] = t1 - 2.0*h + pp - 90.0 # pi1
    arg[:,16] = t1 - h - 90.0 # p1
    if corrections in ('OTIS', 'ATLAS', 'TMD3', 'netcdf'):
        arg[:,17] = t1 + 90.0 # s1
    elif corrections in ('GOT', 'FES'):
        arg[:,17] = t1 + 180.0 # s1 (Doodson's phase)
    arg[:,18] = t1 + h + 90.0 # k1
    arg[:,19] = t1 + 2.0*h - pp + 90.0 # psi1
    arg[:,20] = t1 + 3.0*h + 90.0 # phi1
    arg[:,21] = t1 + s - h + p + 90.0 # theta1
    arg[:,22] = t1 + s + h - p + 90.0 # J1
    arg[:,23] = t1 + 2.0*s + h + 90.0 # OO1
    arg[:,24] = t2 - 4.0*s + 2.0*h + 2.0*p # 2N2
    arg[:,25] = t2 - 4.0*s + 4.0*h # mu2
    arg[:,26] = t2 - 3.0*s + 2.0*h + p # n2
    arg[:,27] = t2 - 3.0*s + 4.0*h - p # nu2
    arg[:,28] = t2 - 2.0*s + h + pp # M2a
    arg[:,29] = t2 - 2.0*s + 2.0*h # M2
    arg[:,30] = t2 - 2.0*s + 3.0*h - pp # M2b
    arg[:,31] = t2 - s + p + 180.0 # lambda2
    arg[:,32] = t2 - s + 2.0*h - p + 180.0 # L2
    arg[:,33] = t2 - h + pp # T2
    arg[:,34] = t2 # S2
    arg[:,35] = t2 + h - pp + 180.0 # R2
    arg[:,36] = t2 + 2.0*h # K2
    arg[:,37] = t2 + s + 2.0*h - pp # eta2
    arg[:,38] = t2 - 5.0*s + 4.0*h + p # MNS2
    arg[:,39] = t2 + 2.0*s - 2.0*h # 2SM2
    arg[:,40] = 1.5*arg[:,29] # M3
    arg[:,41] = arg[:,18] + arg[:,29] # MK3
    arg[:,42] = 3.0*t1 # S3
    arg[:,43] = arg[:,26] + arg[:,29] # MN4
    arg[:,44] = 2.0*arg[:,29] # M4
    arg[:,45] = arg[:,29] + arg[:,34] # MS4
    arg[:,46] = arg[:,29] + arg[:,36] # MK4
    arg[:,47] = 4.0*t1 # S4
    arg[:,48] = 5.0*t1 # S5
    arg[:,49] = 3.0*arg[:,29] # M6
    arg[:,50] = 3.0*t2 # S6
    arg[:,51] = 7.0*t1 # S7
    arg[:,52] = 4.0*t2 # S8
    # shallow water constituents
    arg[:,53] = 4.0*arg[:,29] # m8
    arg[:,54] = arg[:,29] + arg[:,36] - arg[:,34] # mks2
    arg[:,55] = 4.0*s - 2.0*h # msqm
    arg[:,56] = 3.0*s - p # mtm
    arg[:,57] = 2.0*arg[:,26] # n4
    arg[:,58] = t2 - 5.0*s + 4.0*h + p # eps2
    # mean sea level
    arg[:,59] = 0.0 # Z0

    # determine equilibrium arguments
    fargs = np.c_[t1, s, h, p, omega, pp, k]
    test = np.dot(fargs, _arguments_table(corrections=corrections))

    # validate arguments
    assert np.all(arg == test)
