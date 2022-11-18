#!/usr/bin/env python
u"""
infer_minor_corrections.py (11/2022)
Return correction for minor constituents based on Richard Ray's PERTH3 code
    PERTH: PREdict Tidal Heights

CALLING SEQUENCE:
    dh = infer_minor_corrections(t,zmajor,constituents)

INPUTS:
    t: days relative to Jan 1, 1992 (48622 MJD)
    zmajor: Complex HC for given constituents/points
    constituents: tidal constituent IDs

OUTPUT:
    dh: height from minor constituents

OPTIONS:
    deltat: time correction for converting to Ephemeris Time (days)
    corrections: use nodal corrections from OTIS/ATLAS or GOT models

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
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 08/2020: change time variable names to not overwrite functions
        update nodal corrections for FES models
    Updated 07/2020: added function docstrings
        reduce list of minor constituents if in list of major values
    Updated 11/2019: output as numpy masked arrays instead of nan-filled arrays
    Updated 08/2018: added correction option ATLAS for localized OTIS solutions
    Updated 07/2018: added option to use GSFC GOT nodal corrections
        use the number of dates if calculating a tidal time series at a point
    Updated 09/2017: Rewritten in Python
"""
import copy
import warnings
import numpy as np
from pyTMD.calc_astrol_longitudes import calc_astrol_longitudes

def infer_minor_corrections(t, zmajor, constituents, **kwargs):
    """
    Calculate the tidal corrections for minor constituents inferred using
    major constituents

    Parameters
    ----------
    t: float
        days relative to 1992-01-01T00:00:00
    zmajor: complex
        Complex HC for given constituents/points
    constituents: list
        tidal constituent IDs
    deltat: float, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS/ATLAS or GOT models

    Returns
    -------
    dh: float
        Height from minor constituents

    References
    ----------
    .. [1] A. T. Doodson and H. Warburg, "Admiralty Manual of Tides", HMSO, (1941).
    .. [2] P. Schureman, "Manual of Harmonic Analysis and Prediction of Tides"
        US Coast and Geodetic Survey, Special Publication, 98, (1958).
    .. [3] M. G. G. Foreman and R. F. Henry, "The harmonic analysis of tidal model
        time series", Advances in Water Resources, 12, (1989).
    """
    # set default keyword arguments
    kwargs.setdefault('deltat', 0.0)
    kwargs.setdefault('corrections', 'OTIS')
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(DELTAT='deltat',CORRECTIONS='corrections')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn(f"""Deprecated keyword argument {old}.
                Changed to '{new}'""", DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])

    # degrees to radians
    dtr = np.pi/180.0
    # number of constituents
    npts,nc = np.shape(zmajor)
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
    S,H,P,omega,pp = calc_astrol_longitudes(MJD + kwargs['deltat'],
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
