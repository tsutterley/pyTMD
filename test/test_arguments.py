#!/usr/bin/env python
u"""
test_arguments.py (01/2024)
Verify arguments table matches prior arguments array
"""
import pytest
import numpy as np
import pyTMD.astro
from pyTMD.arguments import (
    doodson_number,
    _arguments_table,
    _minor_table,
    _to_doodson_number,
    _from_doodson_number
)

@pytest.mark.parametrize("corrections", ['OTIS', 'GOT'])
def test_arguments(corrections):
    """
    Tests the calculation of nodal corrections for tidal constituents

    Parameters
    ----------
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS/ATLAS or GOT models
    """
    # use a random set of modified Julian days
    MJD = np.random.randint(58000, 61000, size=10)
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
    # convert from hours solar time into mean lunar time in degrees
    tau = 15.0*hour - s + h
    # variable for multiples of 90 degrees (Ray technical note 2017)
    k = 90.0 + np.zeros((nt))

    # determine equilibrium arguments
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

    # determine equilibrium arguments using table
    fargs = np.c_[tau, s, h, p, omega, pp, k]
    test = np.dot(fargs, _arguments_table(corrections=corrections))

    # validate arguments between methods
    assert np.all(arg == test)

def test_minor():
    """
    Tests the calculation of nodal corrections for minor tidal constituents
    """
    # use a random set of modified Julian days
    MJD = np.random.randint(58000, 61000, size=10)
    # convert from Modified Julian Dates into Ephemeris Time
    s, h, p, n, pp = pyTMD.astro.mean_longitudes(MJD)
    # number of temporal values
    nt = len(np.atleast_1d(MJD))
    # initial time conversions
    hour = 24.0*np.mod(MJD, 1)
    # convert from hours into degrees
    t1 = 15.0*hour
    t2 = 30.0*hour
    # convert from hours solar time into mean lunar time in degrees
    tau = 15.0*hour - s + h
    # variable for multiples of 90 degrees (Ray technical note 2017)
    k = 90.0 + np.zeros((nt))

    # determine equilibrium tidal arguments
    arg = np.zeros((nt, 20))
    arg[:,0] = t1 - 4.0*s + h + 2.0*p - 90.0# 2Q1
    arg[:,1] = t1 - 4.0*s + 3.0*h - 90.0# sigma1
    arg[:,2] = t1 - 3.0*s + 3.0*h - p - 90.0# rho1
    arg[:,3] = t1 - s + h - p + 90.0# M12
    arg[:,4] = t1 - s + h + p + 90.0# M11
    arg[:,5] = t1 - s + 3.0*h - p + 90.0# chi1
    arg[:,6] = t1 - 2.0*h + pp - 90.0# pi1
    arg[:,7] = t1 + 3.0*h + 90.0# phi1
    arg[:,8] = t1 + s - h + p + 90.0# theta1
    arg[:,9] = t1 + s + h - p + 90.0# J1
    arg[:,10] = t1 + 2.0*s + h + 90.0# OO1
    arg[:,11] = t2 - 4.0*s + 2.0*h + 2.0*p# 2N2
    arg[:,12] = t2 - 4.0*s + 4.0*h# mu2
    arg[:,13] = t2 - 3.0*s + 4.0*h - p# nu2
    arg[:,14] = t2 - s + p + 180.0# lambda2
    arg[:,15] = t2 - s + 2.0*h - p + 180.0# L2
    arg[:,16] = t2 - s + 2.0*h + p# L2
    arg[:,17] = t2 - h + pp# t2
    arg[:,18] = t2 - 5.0*s + 4.0*h + p # eps2
    arg[:,19] = t2 + s + 2.0*h - pp # eta2

    # determine equilibrium arguments
    fargs = np.c_[tau, s, h, p, n, pp, k]
    test = np.dot(fargs, _minor_table())

    # validate arguments between methods
    assert np.all(arg == test)

def test_doodson():
    """
    Tests the calculation of Doodson numbers
    """
    # expected values
    exp = {}
    # semi-diurnal species
    exp['m2'] = 255.555
    exp['s2'] = 273.555
    exp['n2'] = 245.655
    exp['nu2'] = 247.455
    exp['mu2'] = 237.555
    exp['2n2'] = 235.755
    exp['lambda2'] = 263.655
    exp['l2'] = 265.455
    exp['k2'] = 275.555
    # diurnal species
    exp['m1'] = 155.555
    exp['s1'] = 164.555
    exp['o1'] = 145.555
    exp['oo1'] = 185.555
    exp['k1'] = 165.555
    exp['q1'] = 135.655
    exp['2q1'] = 125.755
    exp['p1'] = 163.555
    # long-period species
    exp['mm'] = 65.455
    exp['ssa'] = 57.555
    exp['msf'] = 73.555
    exp['mf'] = 75.555
    exp['msqm'] = 93.555
    exp['mtm'] = 85.455
    # short-period species
    exp['m3'] = 355.555
    exp['m4'] = 455.555
    exp['m6'] = 655.555
    exp['m8'] = 855.555
    # get observed values for constituents
    obs = doodson_number(exp.keys())
    cartwright = doodson_number(exp.keys(), formalism='Cartwright')
    # check values
    for key,val in exp.items():
        assert val == obs[key]
        # check values when entered as string
        test = doodson_number(key)
        assert val == test
        # check conversion to and from Doodson numbers
        doodson = _to_doodson_number(cartwright[key])
        # check values when entered as Cartwright
        assert val == doodson
        # check values when entered as Doodson
        coefficients = _from_doodson_number(val)
        assert np.all(cartwright[key] == coefficients)
