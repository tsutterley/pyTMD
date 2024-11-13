#!/usr/bin/env python
u"""
test_arguments.py (11/2024)
Verify arguments table matches prior arguments array
Verify nodal corrections match prior estimates

UPDATE HISTORY:
    Updated 11/2024: moved normalize_angle to math.py
    Updated 10/2024: add comparisons for formatted Doodson numbers
        add function to parse tide potential tables
    Updated 08/2024: add comparisons for nodal corrections
    Written 01/2024
"""
import pytest
import numpy as np
import pyTMD.astro
import pyTMD.math
from pyTMD.arguments import (
    nodal,
    doodson_number,
    _arguments_table,
    _minor_table,
    _parse_tide_potential_table,
    _to_doodson_number,
    _to_extended_doodson,
    _from_doodson_number,
    _from_extended_doodson,
    _ct1971_table_5,
    _ce1973_table_1
)

@pytest.mark.parametrize("corrections", ['OTIS', 'GOT'])
def test_arguments(corrections):
    """
    Tests the calculation of nodal corrections for tidal constituents

    Parameters
    ----------
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS, FES or GOT models
    """
    # use a random set of modified Julian days
    MJD = np.random.randint(58000, 61000, size=10)
    # set function for astronomical longitudes
    ASTRO5 = corrections not in ('OTIS','ATLAS','TMD3','netcdf')
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
    else:
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
    arg[:,37] = t2 + s + 2.0*h - p # eta2
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

@pytest.mark.parametrize("corrections", ['OTIS', 'GOT'])
def test_table(corrections):
    """
    Compare Doodson coefficients tables
    """
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
    if corrections in ('OTIS', 'ATLAS', 'TMD3', 'netcdf'):
        coef[:,17] = [1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 1.0] # s1
    else:
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
    coef[:,37] = [2.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0] # eta2
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

    # get Doodson coefficients
    test = pyTMD.arguments._arguments_table(corrections=corrections)
    # validate arguments between methods
    assert np.all(coef == test)

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
    arg[:,19] = t2 + s + 2.0*h - p # eta2

    # determine equilibrium arguments
    fargs = np.c_[tau, s, h, p, n, pp, k]
    test = np.dot(fargs, _minor_table())

    # validate arguments between methods
    assert np.all(arg == test)

@pytest.mark.parametrize("corrections", ['OTIS', 'FES', 'GOT'])
@pytest.mark.parametrize("M1", ['Doodson', 'Ray'])
def test_nodal(corrections, M1):
    # degrees to radians
    dtr = np.pi/180.0
    # use a random set of modified Julian days
    MJD = np.random.randint(58000, 61000, size=10)
    # set function for astronomical longitudes
    ASTRO5 = corrections not in ('OTIS','ATLAS','TMD3','netcdf')
    # convert from Modified Julian Dates into Ephemeris Time
    s, h, p, n, pp = pyTMD.astro.mean_longitudes(MJD,
        ASTRO5=ASTRO5)
    # number of temporal values
    nt = len(np.atleast_1d(MJD))

    # constituents array (not all are included in tidal program)
    cindex = ['sa', 'ssa', 'mm', 'msf', 'mf', 'mt', 'alpha1', '2q1', 'sigma1',
        'q1', 'rho1', 'o1', 'tau1', 'm1', 'chi1', 'pi1', 'p1', 's1', 'k1',
        'psi1', 'phi1', 'theta1', 'j1', 'oo1', '2n2', 'mu2', 'n2', 'nu2', 'm2a',
        'm2', 'm2b', 'lambda2', 'l2', 't2', 's2', 'r2', 'k2', 'eta2', 'mns2',
        '2sm2', 'm3', 'mk3', 's3', 'mn4', 'm4', 'ms4', 'mk4', 's4', 's5', 'm6',
        's6', 's7', 's8', 'm8', 'mks2', 'msqm', 'mtm', 'n4', 'eps2', 'z0']

    # get nodal corrections
    pu, pf = nodal(n, p, cindex, corrections=corrections, M1=M1)

    # trigonometric factors for nodal corrections
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
    if corrections in ('OTIS', 'ATLAS', 'TMD3', 'netcdf'):
        # constituents to test
        constituents = ['q1','o1','p1','k1','n2','m2','s2','k2','m4',
            'ms4','mn4','2n2','mf','mm','s1','m1']
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
        f[:,7] = np.sqrt((1.0 + 0.189*cosn)**2 + (0.189*sinn)**2) # 2Q1
        f[:,8] = f[:,7] # sigma1
        f[:,9] = f[:,7] # q1
        f[:,10] = f[:,7] # rho1
        temp1 = (1.0 + 0.189*cosn - 0.0058*cos2n)**2
        temp2 = (0.189*sinn - 0.0058*sin2n)**2
        f[:,11] = np.sqrt(temp1 + temp2) # O1
        f[:,12] = 1.0 # tau1
        if (M1 == 'Doodson'):
            # A. T. Doodson's coefficients for M1 tides
            Mtmp1 = 2.0*np.cos(p*dtr) + 0.4*np.cos((p-n)*dtr)
            Mtmp2 = np.sin(p*dtr) + 0.2*np.sin((p-n)*dtr)
        elif (M1 == 'Ray'):
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
        temp1 = -(0.3108*sinn + 0.0324*sin2n)
        temp2 = (1.0 + 0.2852*cosn + 0.0324*cos2n)
        u[:,36] = np.arctan(temp1/temp2)/dtr # K2

        temp1 = (1.0 + 0.2852*cosn + 0.0324*cos2n)**2
        temp2 = (0.3108*sinn + 0.0324*sin2n)**2

        u[:,37] = np.arctan(-0.436*sinn/(1.0 + 0.436*cosn))/dtr # eta2
        u[:,38] = u[:,29]*2.0 # MNS2
        u[:,39] = -u[:,29] # 2SM2
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
        # shallow water constituents
        u[:,53] = 4.0*u[:,29] # m8
        u[:,54] = u[:,29] + u[:,36] # mks2
        u[:,55] = u[:,4] # msqm
        u[:,56] = u[:,4] # mtm
        u[:,57] = 2.0*u[:,29] # n4
        u[:,58] = u[:,29] # eps2
        # mean sea level
        u[:,59] = 0.0 # Z0

    elif corrections in ('FES',):
        # constituents to test
        constituents = ['nu2', 'mtm', 'msf', 'j1', 'q1', 'l2', '2n2', 'k2',
            'k1', 'msqm', 'mm', 'm8', 'sa', 'lambda2', 'mn4', 'p1', 's1',
            's4', 'ssa', 'n4', 'm2', 'eps2', 'm4', 'mf', 'm3', 's2', 'n2',
            'm6', 'mu2', 'o1', 'mks2', 't2', 'ms4', 'r2', 'm1']
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
        if (M1 == 'Doodson'):
            # A. T. Doodson's coefficients for M1 tides
            Mtmp1 = 2.0*np.cos(p*dtr) + 0.4*np.cos((p-n)*dtr)
            Mtmp2 = np.sin(p*dtr) + 0.2*np.sin((p-n)*dtr)
        elif (M1 == 'Ray'):
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

    elif corrections in ('GOT',):
        # constituents to test
        constituents = ['q1','o1','k1','n2','m2','s2','k2','s1','m4','m1']
        # nodal factors
        f[:,9] = 1.009 + 0.187*cosn - 0.015*cos2n# Q1
        f[:,11] = f[:,9]# O1
        if (M1 == 'Doodson'):
            # A. T. Doodson's coefficients for M1 tides
            Mtmp1 = 2.0*np.cos(p*dtr) + 0.4*np.cos((p-n)*dtr)
            Mtmp2 = np.sin(p*dtr) + 0.2*np.sin((p-n)*dtr)
        elif (M1 == 'Ray'):
            # R. Ray's coefficients for M1 tides
            Mtmp1 = 1.36*np.cos(p*dtr) + 0.267*np.cos((p-n)*dtr)
            Mtmp2 = 0.64*np.sin(p*dtr) + 0.135*np.sin((p-n)*dtr)
        f[:,13] = np.sqrt(Mtmp1**2 + Mtmp2**2) # M1
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
        u[:,13] = np.arctan2(Mtmp2,Mtmp1)/dtr # M1
        u[:,16] = 0.0 # P1
        u[:,17] = 0.0 # S1
        u[:,18] = -8.9*sinn + 0.7*sin2n# K1
        u[:,26] = -2.1*sinn# N2
        u[:,29] = u[:,26]# M2
        u[:,34] = 0.0 # S2
        u[:,36] = -17.7*sinn + 0.7*sin2n# K2
        u[:,44] = -4.2*sinn# M4

    # validate arguments between methods
    for i, c in enumerate(cindex):
        # only verify a subset of constituents for each model type
        if c in constituents:
            # verify nodal factors
            assert np.all(np.isclose(f[:,i], pf[:,i], rtol=1e-2, atol=1e-2))
            # verify nodal angles in radians
            urad = u[:,i]*dtr
            assert np.all(np.isclose(urad, pu[:,i], rtol=1e-2, atol=1e-2))

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
    exp['mm'] = 065.455
    exp['ssa'] = 057.555
    exp['msf'] = 073.555
    exp['mf'] = 075.555
    exp['msqm'] = 093.555
    exp['mtm'] = 085.455
    exp['node'] = 055.565
    # short-period species
    exp['m3'] = 355.555
    exp['m4'] = 455.555
    exp['m6'] = 655.555
    exp['m8'] = 855.555
    # shallow water species
    exp['2so3'] = '3X1.555'
    exp['2jp3'] = '3X3.355'
    exp['kso3'] = '3X3.555'
    exp['2jk3'] = '3X4.355'
    exp['2ko3'] = '3X5.555'
    # 3rd degree terms
    exp["2q1'"] = 125.655
    exp["q1'"] = 135.555
    exp["o1'"] = 145.655
    exp["m1'"] = 155.555
    exp["k1'"] = 165.455
    exp["j1'"] = 175.555
    exp["2n2'"] = 235.655
    exp["n2'"] = 245.555
    exp["m2'"] = 255.655
    exp["l2'"] = 265.555
    exp["m3'"] = 355.555
    exp['lambda3'] = 363.655
    exp["l3"] = 365.455
    exp["l3b"] = 365.655
    exp["f3"] = 375.555
    exp["j3"] = 375.555
    exp["s3'"] = 382.555
    # get observed values for constituents
    obs = doodson_number(exp.keys())
    cartwright = doodson_number(exp.keys(), formalism='Cartwright')
    # check values
    for key,val in exp.items():
        assert val == obs[key]
        # check values when entered as string
        test = doodson_number(key)
        assert val == test
        # check conversion to Doodson numbers
        doodson = _to_doodson_number(cartwright[key])
        # check values when entered as Cartwright
        assert val == doodson
        # check values when entered as Doodson
        coefficients = _from_doodson_number(val)
        assert np.all(cartwright[key] == coefficients)

def test_extended():
    """
    Tests the calculation of UKHO Extended Doodson numbers
    """
    # expected values
    exp = {}
    # semi-diurnal species
    exp['m2'] = 'BZZZZZZ'
    exp['s2'] = 'BBXZZZZ'
    exp['n2'] = 'BYZAZZZ'
    exp['nu2'] = 'BYBYZZZ'
    exp['mu2'] = 'BXBZZZZ'
    exp['2n2'] = 'BXZBZZZ'
    exp['lambda2'] = 'BAXAZZB'
    exp['l2'] = 'BAZYZZB'
    exp['k2'] = 'BBZZZZZ'
    # diurnal species
    exp['m1'] = 'AZZZZZA'
    exp['s1'] = 'AAYZZZA'
    exp['o1'] = 'AYZZZZY'
    exp['oo1'] = 'ACZZZZA'
    exp['k1'] = 'AAZZZZA'
    exp['q1'] = 'AXZAZZY'
    exp['2q1'] = 'AWZBZZY'
    exp['p1'] = 'AAXZZZY'
    # long-period species
    exp['mm'] = 'ZAZYZZZ'
    exp['ssa'] = 'ZZBZZZZ'
    exp['msf'] = 'ZBXZZZZ'
    exp['mf'] = 'ZBZZZZZ'
    exp['msqm'] = 'ZDXZZZZ'
    exp['mtm'] = 'ZCZYZZZ'
    exp['node'] = 'ZZZZAZB'
    # short-period species
    exp['m3'] = 'CZZZZZZ'
    exp['m4'] = 'DZZZZZZ'
    exp['n4'] = 'DXZBZZZ'
    exp['m6'] = 'FZZZZZZ'
    exp['n6'] = 'FWZCZZZ'
    exp['m8'] = 'HZZZZZZ'
    exp['m10'] = 'JZZZZZZ'
    exp['m12'] = 'LZZZZZZ'
    # shallow water species
    exp['2so3'] = 'CEVZZZA'
    exp['2jp3'] = 'CEXXZZA'
    exp['kso3'] = 'CEXZZZA'
    exp['2jk3'] = 'CEYXZZZ'
    exp['2ko3'] = 'CEZZZZA'
    # get observed values for constituents
    obs = doodson_number(exp.keys(), formalism='Extended')
    # check values
    for key,val in exp.items():
        assert val == obs[key]

def test_parse_tables():
    """
    Tests the parsing of tables for tide potential coefficients
    """
    # Cartwright and Tayler (1971) table with 3rd-degree values
    # Cartwright and Edden (1973) table with updated values
    for table in [_ct1971_table_5, _ce1973_table_1]:
        # parse table
        CTE = _parse_tide_potential_table(table)
        for i, line in enumerate(CTE):
            # convert Doodson number to Cartwright numbers
            tau, s, h, p, n, pp = _from_doodson_number(line['DO'])
            assert tau == line['tau'], line
            assert s == line['s'], line
            assert h == line['h'], line
            assert p == line['p'], line
            assert n == line['n'], line
            assert pp == line['pp'], line

def test_normalize_angle():
    """
    Tests the normalization of angles to between 0 and 360 degrees
    """
    # test angles
    angles = np.array([-180, -90, 0, 90, 180, 270, 360, 450])
    # expected values
    exp = np.array([180, 270, 0, 90, 180, 270, 0, 90])
    # test normalization of angles
    test = pyTMD.math.normalize_angle(angles)
    assert np.all(exp == test)
