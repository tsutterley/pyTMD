#!/usr/bin/env python
u"""
astro.py
Written by Tyler Sutterley (04/2023)
Astronomical and nutation routines

mean_longitudes is a modification of the ASTROL fortran
subroutine by Richard Ray written in 03/1999

Computes the basic astronomical mean longitudes: s, h, p, N and PP
Note N is not N', i.e. N is decreasing with time.

Formulae for the period 1990--2010 were derived by David Cartwright
MEEUS and ASTRO5 formulae are from versions of Meeus's Astronomical Algorithms

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    PyEphem: Astronomical Ephemeris for Python
        https://rhodesmill.org/pyephem/

REFERENCES:
    Jean Meeus, Astronomical Algorithms, 2nd edition, 1998.
    Oliver Montenbruck, Practical Ephemeris Calculations, 1989.

UPDATE HISTORY:
    Updated 04/2023: added low resolution solar and lunar positions
        added function with more phase angles of the sun and moon
        functions to calculate solar and lunar positions with ephemerides
        add jplephem documentation to Spacecraft and Planet Kernel segments
        fix solar ephemeride function to include SSB to sun segment
        use a higher resolution estimate of the Greenwich hour angle
    Updated 03/2023: add basic variable typing to function inputs
    Updated 10/2022: fix MEEUS solar perigee rate
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 08/2020: change time variable names to not overwrite functions
    Updated 07/2020: added function docstrings
    Updated 07/2018: added option ASTRO5 to use coefficients from Richard Ray
        for use with the GSFC Global Ocean Tides (GOT) model
        added longitude of solar perigee (PP) as an additional output
    Updated 09/2017: added option MEEUS to use additional coefficients
        from Meeus Astronomical Algorithms to calculate mean longitudes
    Updated 09/2017: Rewritten in Python
    Rewritten in Matlab by Lana Erofeeva 2003
    Written by Richard Ray 12/1990
"""
from __future__ import annotations

import logging
import numpy as np
from pyTMD.time import timescale
from pyTMD.eop import iers_polar_motion
from pyTMD.utilities import get_data_path

# attempt imports
try:
    import jplephem.spk
except (ImportError, ModuleNotFoundError) as exc:
    logging.debug("jplephem not available")

# default JPL Spacecraft and Planet ephemerides kernel
_default_kernel = get_data_path(['data','de440s.bsp'])

# PURPOSE: calculate the sum of a polynomial function of time
def polynomial_sum(coefficients: list | np.ndarray, t: np.ndarray):
    """
    Calculates the sum of a polynomial function of time

    Parameters
    ----------
    coefficients: list or np.ndarray
        leading coefficient of polynomials of increasing order
    t: np.ndarray
        delta time in units for a given astronomical longitudes calculation
    """
    # convert time to array if importing a single value
    t = np.atleast_1d(t)
    return np.sum([c * (t ** i) for i, c in enumerate(coefficients)], axis=0)

def mxv(m, v):
    """Multiply a matrix by a vector."""
    return np.einsum('ij...,j...->i...', m, v)

def mxm(m1, m2):
    """Multiply 2 matrices together
    """
    return np.einsum('ij...,jk...->ik...', m1, m2)

def mxmxm(m1, m2, m3):
    """Multiply 3 matrices together
    """
    return np.einsum('ij...,jk...,kl...->il...', m1, m2, m3)

def rotate(theta: float | np.ndarray, axis: str = 'x'):
    """
    Rotate a 3-dimensional matrix about a given axis

    Parameters
    ----------
    theta: float or np.ndarray
        Angle of rotation in radians
    axis: str
        Axis of rotation (``'x'``, ``'y'``, or ``'z'``)
    """
    # allocate for output rotation matrix
    R = np.zeros((3, 3, len(np.atleast_1d(theta))))
    if (axis.lower() == 'x'):
        # rotate about x-axis
        R[0,0,:] = 1.0
        R[1,1,:] = np.cos(theta)
        R[1,2,:] = np.sin(theta)
        R[2,1,:] = -np.sin(theta)
        R[2,2,:] = np.cos(theta)
    elif (axis.lower() == 'y'):
        # rotate about y-axis
        R[0,0,:] = np.cos(theta)
        R[0,2,:] = -np.sin(theta)
        R[1,1,:] = 1.0
        R[2,0,:] = np.sin(theta)
        R[2,2,:] = np.cos(theta)
    elif (axis.lower() == 'z'):
        # rotate about z-axis
        R[0,0,:] = np.cos(theta)
        R[0,1,:] = np.sin(theta)
        R[1,0,:] = -np.sin(theta)
        R[1,1,:] = np.cos(theta)
        R[2,2,:] = 1.0
    else:
        raise ValueError(f'Invalid axis {axis}')
    # return the rotation matrix
    return R

def scale(factor: float | np.ndarray):
    """
    Scale a 3-dimensional matrix by a given factor

    Parameters
    ----------
    factor: float or np.ndarray
        Scaling factor
    """
    # calculate output scaling matrix
    S = np.diag(np.ones((3))*factor)
    # return the scaling matrix
    return S

# PURPOSE: compute the basic astronomical mean longitudes
def mean_longitudes(
        MJD: np.ndarray,
        MEEUS: bool = False,
        ASTRO5: bool = False
    ):
    """
    Computes the basic astronomical mean longitudes:
    `S`, `H`, `P`, `N` and `PP` [1]_

    Note `N` is not `N'`, i.e. `N` is decreasing with time.

    Parameters
    ----------
    MJD: np.ndarray
        Modified Julian Day (MJD) of input date
    MEEUS: bool, default False
        use additional coefficients from Meeus Astronomical Algorithms
    ASTRO5: bool, default False
        use Meeus Astronomical coefficients as implemented in ``ASTRO5``

    Returns
    -------
    S: np.ndarray
        mean longitude of moon (degrees)
    H: np.ndarray
        mean longitude of sun (degrees)
    P: np.ndarray
        mean longitude of lunar perigee (degrees)
    N: np.ndarray
        mean longitude of ascending lunar node (degrees)
    PP: np.ndarray
        longitude of solar perigee (degrees)

    References
    ----------
    .. [1] J. Meeus, *Astronomical Algorithms*, 2nd edition, 477 pp., (1998).
    """
    circle = 360.0
    if MEEUS:
        # convert from MJD to days relative to 2000-01-01T12:00:00
        T = MJD - 51544.5
        # mean longitude of moon
        lunar_longitude = np.array([218.3164591, 13.17639647754579,
            -9.9454632e-13, 3.8086292e-20, -8.6184958e-27])
        S = polynomial_sum(lunar_longitude, T)
        # mean longitude of sun
        solar_longitude = np.array([280.46645, 0.985647360164271,
            2.2727347e-13])
        H = polynomial_sum(solar_longitude, T)
        # mean longitude of lunar perigee
        lunar_perigee = np.array([83.3532430, 0.11140352391786447,
            -7.7385418e-12, -2.5636086e-19, 2.95738836e-26])
        P = polynomial_sum(lunar_perigee, T)
        # mean longitude of ascending lunar node
        lunar_node = np.array([125.0445550, -0.052953762762491446,
            1.55628359e-12, 4.390675353e-20, -9.26940435e-27])
        N = polynomial_sum(lunar_node, T)
        # mean longitude of solar perigee (Simon et al., 1994)
        PP = 282.94 + 1.7192 * T / 36525.0
    elif ASTRO5:
        # convert from MJD to centuries relative to 2000-01-01T12:00:00
        T = (MJD - 51544.5)/36525.0
        # mean longitude of moon (p. 338)
        lunar_longitude = np.array([218.3164477, 481267.88123421, -1.5786e-3,
             1.855835e-6, -1.53388e-8])
        S = polynomial_sum(lunar_longitude, T)
        # mean longitude of sun (p. 338)
        lunar_elongation = np.array([297.8501921, 445267.1114034, -1.8819e-3,
             1.83195e-6, -8.8445e-9])
        H = polynomial_sum(lunar_longitude-lunar_elongation, T)
        # mean longitude of lunar perigee (p. 343)
        lunar_perigee = np.array([83.3532465, 4069.0137287, -1.032e-2,
            -1.249172e-5])
        P = polynomial_sum(lunar_perigee, T)
        # mean longitude of ascending lunar node (p. 144)
        lunar_node = np.array([125.04452, -1934.136261, 2.0708e-3, 2.22222e-6])
        N = polynomial_sum(lunar_node, T)
        # mean longitude of solar perigee (Simon et al., 1994)
        PP = 282.94 + 1.7192 * T
    else:
        # Formulae for the period 1990--2010 were derived by David Cartwright
        # convert from MJD to days relative to 2000-01-01T12:00:00
        # convert from Universal Time to Dynamic Time at 2000-01-01
        T = MJD - 51544.4993
        # mean longitude of moon
        S = 218.3164 + 13.17639648 * T
        # mean longitude of sun
        H = 280.4661 + 0.98564736 * T
        # mean longitude of lunar perigee
        P = 83.3535 + 0.11140353 * T
        # mean longitude of ascending lunar node
        N = 125.0445 - 0.05295377 * T
        # solar perigee at epoch 2000
        PP = 282.8
    # take the modulus of each
    S = np.mod(S, circle)
    H = np.mod(H, circle)
    P = np.mod(P, circle)
    N = np.mod(N, circle)
    # return as tuple
    return (S, H, P, N, PP)

# PURPOSE: computes the phase angles of astronomical means
def phase_angles(MJD: np.ndarray):
    """
    Computes astronomical phase angles for the sun and moon:
    `S`, `H`, `P`, `TAU`, `ZNS` and `PS` [1]_ [2]_

    Parameters
    ----------
    MJD: np.ndarray
        Modified Julian Day (MJD) of input date

    Returns
    -------
    S: np.ndarray
        mean longitude of moon (radians)
    H: np.ndarray
        mean longitude of sun (radians)
    P: np.ndarray
        mean longitude of lunar perigee (radians)
    TAU: np.ndarray
        time angle in lunar days (radians)
    ZNS: np.ndarray
        mean longitude of ascending lunar node `N'` (radians)
    PS: np.ndarray
        mean longitude of solar perigee (radians)

    References
    ----------
    .. [1] A. T. Doodson and H. Lamb, "The harmonic development of
        the tide-generating potential", *Proceedings of the Royal Society
        of London. Series A, Containing Papers of a Mathematical and
        Physical Character*, 100(704), 305--329, (1921).
        `doi: 10.1098/rspa.1921.0088 <https://doi.org/10.1098/rspa.1921.0088>`_
    .. [2] J. Meeus, *Astronomical Algorithms*, 2nd edition, 477 pp., (1998).
    """
    # degrees to radians
    dtr = np.pi/180.0
    # convert from MJD to centuries relative to 2000-01-01T12:00:00
    T = (MJD - 51544.5)/36525.0
    # 360 degrees
    circle = 360.0
    # hour of the day
    FHR = np.mod(MJD, 1)*24.0
    # calculate phase angles
    # mean longitude of moon (degrees)
    S = polynomial_sum(np.array([218.3164477, 481267.88123421,
        -1.5786e-3, 1.855835e-6, -1.53388e-8]), T)
    # time angle in lunar days (degrees)
    TAU = ((FHR*15.0) - S + polynomial_sum(np.array([280.4606184,
        36000.7700536, 3.8793e-4, -2.58e-8]), T))
    # calculate correction for mean lunar longitude (degrees)
    PR = polynomial_sum(np.array([0.0, 1.396971278,
        3.08889e-4, 2.1e-8, 7.0e-9]), T)
    S += PR
    # mean longitude of sun (degrees)
    H = polynomial_sum(np.array([280.46645, 36000.7697489,
        3.0322222e-4, 2.0e-8, -6.54e-9]), T)
    # mean longitude of lunar perigee (degrees)
    P = polynomial_sum(np.array([83.3532465, 4069.0137287,
        -1.032172222e-2, -1.24991e-5, 5.263e-8]), T)
    # mean longitude of ascending lunar node (degrees)
    ZNS = polynomial_sum(np.array([234.95544499, 1934.13626197,
        -2.07561111e-3, -2.13944e-6, 1.65e-8]), T)
    # mean longitude of solar perigee (degrees)
    PS = polynomial_sum(np.array([282.93734098, 1.71945766667,
        4.5688889e-4, -1.778e-8, -3.34e-9]), T)
    # take the modulus of each and convert to radians
    S = dtr*np.mod(S, circle)
    H = dtr*np.mod(H, circle)
    P = dtr*np.mod(P, circle)
    TAU = dtr*np.mod(TAU, circle)
    ZNS = dtr*np.mod(ZNS, circle)
    PS = dtr*np.mod(PS, circle)
    # return as tuple
    return (S, H, P, TAU, ZNS, PS)

def delaunay_arguments(MJD: np.ndarray):
    """
    Computes astronomical phase angles for the five primary Delaunay
    Arguments of Nutation: `l`, `l'`, `F`, `D`, and `N` [1]_ [2]_ [3]_

    Parameters
    ----------
    MJD: np.ndarray
        Modified Julian Day (MJD) of input date

    Returns
    -------
    l: np.ndarray
        mean anomaly of moon (radians)
    lp: np.ndarray
        mean anomaly of the sun (radians)
    F: np.ndarray
        mean argument of the moon (radians)
    D: np.ndarray
        mean elongation of the moon from the sun (radians)
    N: np.ndarray
        mean longitude of ascending lunar node (radians)

    References
    ----------
    .. [1] J. Meeus, *Astronomical Algorithms*, 2nd edition, 477 pp., (1998).
    .. [2] G. Petit and B. Luzum (eds.), *IERS Conventions (2010)*,
        International Earth Rotation and Reference Systems Service (IERS),
        `IERS Technical Note No. 36
        <https://iers-conventions.obspm.fr/content/tn36.pdf>`_
    .. [3] N. Capitaine, P. T. Wallace, and J. Chapront,
        "Expressions for IAU 2000 precession quantities",
        *Astronomy & Astrophysics*, 412, 567--586, (2003).
        `doi: 10.1051/0004-6361:20031539
        <https://doi.org/10.1051/0004-6361:20031539>`_
    """
    # arcseconds to radians
    atr = np.pi/648000.0
    # convert from MJD to centuries relative to 2000-01-01T12:00:00
    T = (MJD - 51544.5)/36525.0
    # 360 degrees
    circle = 1296000
    # mean anomaly of the moon (arcseconds)
    l = polynomial_sum(np.array([485868.249036, 1717915923.2178,
        31.8792, 0.051635, -2.447e-04]), T)
    # mean anomaly of the sun (arcseconds)
    lp = polynomial_sum(np.array([1287104.79305,  129596581.0481,
        -0.5532, 1.36e-4, -1.149e-05]), T)
    # mean argument of the moon (arcseconds)
    # (angular distance from the ascending node)
    F = polynomial_sum(np.array([335779.526232, 1739527262.8478,
        -12.7512, -1.037e-3, 4.17e-6]), T)
    # mean elongation of the moon from the sun (arcseconds)
    D = polynomial_sum(np.array([1072260.70369, 1602961601.2090,
        -6.3706, 6.593e-3, -3.169e-05]), T)
    # mean longitude of the ascending node of the moon (arcseconds)
    N = polynomial_sum(np.array([450160.398036, -6962890.5431,
        7.4722, 7.702e-3, -5.939e-05]), T)
    # take the modulus of each and convert to radians
    l = atr*np.mod(l, circle)
    lp = atr*np.mod(lp, circle)
    F = atr*np.mod(F, circle)
    D = atr*np.mod(D, circle)
    N = atr*np.mod(N, circle)
    # return as tuple
    return (l, lp, F, D, N)

def mean_obliquity(MJD: np.ndarray):
    """Mean obliquity of the ecliptic

    Parameters
    ----------
    MJD: np.ndarray
        Modified Julian Day (MJD) of input date

    References
    ----------
    .. [1] N. Capitaine, P. T. Wallace, and J. Chapront,
        "Expressions for IAU 2000 precession quantities",
        *Astronomy & Astrophysics*, 412, 567--586, (2003).
        `doi: 10.1051/0004-6361:20031539
        <https://doi.org/10.1051/0004-6361:20031539>`_
    .. [1] N. Capitaine, J. Chapront, S. Lambert, and P. T. Wallace,
        "Expressions for the Celestial Intermediate Pole and
        Celestial Ephemeris Origin consistent with the IAU 2000A
        precession-nutation model", *Astronomy & Astrophysics*,
        400, 1145--1154, (2003). `doi: 10.1051/0004-6361:20030077
        <https://doi.org/10.1051/0004-6361:20030077>`_
    """
    # arcseconds to radians
    atr = np.pi/648000.0
    # convert from MJD to centuries relative to 2000-01-01T12:00:00
    T = (MJD - 51544.5)/36525.0
    # mean obliquity of the ecliptic (arcseconds)
    epsilon0 = np.array([84381.406, -46.836769, -1.831e-4,
        2.00340e-4, -5.76e-07, -4.34e-08])
    return atr*polynomial_sum(epsilon0, T)

# PURPOSE: compute coordinates of the sun in an ECEF frame
def solar_ecef(MJD: np.ndarray):
    """
    Computes approximate positional coordinates of the sun in an
    Earth-centric, Earth-Fixed (ECEF) frame [1]_ [2]_

    Parameters
    ----------
    MJD: np.ndarray
        Modified Julian Day (MJD) of input date

    Returns
    -------
    X, Y, Z: np.ndarray
        ECEF coordinates of the sun (meters)

    References
    ----------
    .. [1] J. Meeus, *Astronomical Algorithms*, 2nd edition, 477 pp., (1998).
    .. [2] O. Montenbruck, *Practical Ephemeris Calculations*,
        146 pp., (1989).
    """
    # degrees and arcseconds to radians
    dtr = np.pi/180.0
    atr = np.pi/648000.0
    # create timescale from Modified Julian Day (MJD)
    ts = timescale(MJD=MJD)
    # convert from MJD to centuries relative to 2000-01-01T12:00:00
    T = ts.T
    # mean longitude of solar perigee (radians)
    PP = dtr*(282.94 + 1.7192 * T)
    # mean anomaly of the sun (radians)
    solar_anomaly = np.array([357.5256, 35999.049, -1.559e-4, -4.8e-7])
    M = dtr*polynomial_sum(solar_anomaly, T)
    # series expansion for mean anomaly in solar radius (meters)
    r_sun = 1e9*(149.619 - 2.499*np.cos(M) - 0.021*np.cos(2.0*M))
    # series expansion for ecliptic longitude of the sun (radians)
    lambda_sun = PP + M + atr*(6892.0*np.sin(M) + 72.0*np.sin(2.0*M))
    # ecliptic latitude is equal to 0 within 1 arcminute
    # obliquity of the J2000 ecliptic (radians)
    epsilon_j2000 = 23.43929111*dtr
    # convert to position vectors
    x = r_sun*np.cos(lambda_sun)
    y = r_sun*np.sin(lambda_sun)*np.cos(epsilon_j2000)
    z = r_sun*np.sin(lambda_sun)*np.sin(epsilon_j2000)
    # Greenwich hour angle (radians)
    rot_z = rotate(ts.gha*dtr, 'z')
    # rotate to cartesian (ECEF) coordinates
    # ignoring polar motion and length-of-day variations
    X = rot_z[0,0,:]*x + rot_z[0,1,:]*y + rot_z[0,2,:]*z
    Y = rot_z[1,0,:]*x + rot_z[1,1,:]*y + rot_z[1,2,:]*z
    Z = rot_z[2,0,:]*x + rot_z[2,1,:]*y + rot_z[2,2,:]*z
    # return the ECEF coordinates
    return (X, Y, Z)

# PURPOSE: compute coordinates of the sun in an ECEF frame
def solar_ephemerides(MJD: np.ndarray, **kwargs):
    """
    Computes positional coordinates of the sun in an Earth-centric,
    Earth-Fixed (ECEF) frame using JPL ephemerides [1]_ [2]_

    Parameters
    ----------
    MJD: np.ndarray
        Modified Julian Day (MJD) of input date
    kernel: str or pathlib.Path
        Path to JPL ephemerides kernel file

    Returns
    -------
    X, Y, Z: np.ndarray
        ECEF coordinates of the sun (meters)

    References
    ----------
    .. [1] J. Meeus, *Astronomical Algorithms*, 2nd edition, 477 pp., (1998).
    .. [2] R. S. Park, W. M. Folkner, and J. G. Williams, and D. H. Boggs,
        "The JPL Planetary and Lunar Ephemerides DE440 and DE441",
        *The Astronomical Journal*, 161(3), 105, (2021).
        `doi: 10.3847/1538-3881/abd414
        <https://doi.org/10.3847/1538-3881/abd414>`_
    """
    # set default keyword arguments
    kwargs.setdefault('kernel', _default_kernel)
    # degrees and arcseconds to radians
    dtr = np.pi/180.0
    # create timescale from Modified Julian Day (MJD)
    ts = timescale(MJD=MJD)
    # convert from MJD to Julian days (dynamical time)
    JD = ts.tt
    # read JPL ephemerides kernel
    SPK = jplephem.spk.SPK.open(kwargs['kernel'])
    # segments for computing position of the sun
    # segment 0 SOLAR SYSTEM BARYCENTER -> segment 10 SUN
    SSB_to_Sun = SPK[0, 10]
    # segment 0 SOLAR SYSTEM BARYCENTER -> segment 3 EARTH BARYCENTER
    SSB_to_EMB = SPK[0, 3]
    # segment 3 EARTH BARYCENTER -> segment 399 EARTH
    EMB_to_Earth = SPK[3, 399]
    # compute the position of the sun relative to the Earth in meters
    # Earth_to_Sun = Earth_to_EMB + EMB_to_SSB + SSB_to_Sun
    #              = -EMB_to_Earth - SSB_to_EMB + SSB_to_Sun
    x, y, z = 1e3*(SSB_to_Sun.compute(JD) - SSB_to_EMB.compute(JD) -
        EMB_to_Earth.compute(JD))
    # Greenwich hour angle (radians)
    rot_z = rotate(ts.gha*dtr, 'z')
    # rotate to cartesian (ECEF) coordinates
    # ignoring polar motion and length-of-day variations
    X = rot_z[0,0,:]*x + rot_z[0,1,:]*y + rot_z[0,2,:]*z
    Y = rot_z[1,0,:]*x + rot_z[1,1,:]*y + rot_z[1,2,:]*z
    Z = rot_z[2,0,:]*x + rot_z[2,1,:]*y + rot_z[2,2,:]*z
    # return the ECEF coordinates
    return (X, Y, Z)

# PURPOSE: compute coordinates of the moon in an ECEF frame
def lunar_ecef(MJD: np.ndarray):
    """
    Computes approximate positional coordinates of the moon in an
    Earth-centric, Earth-Fixed (ECEF) frame [1]_ [2]_

    Parameters
    ----------
    MJD: np.ndarray
        Modified Julian Day (MJD) of input date

    Returns
    -------
    X, Y, Z: np.ndarray
        ECEF coordinates of the moon (meters)

    References
    ----------
    .. [1] J. Meeus, *Astronomical Algorithms*, 2nd edition, 477 pp., (1998).
    .. [2] O. Montenbruck, *Practical Ephemeris Calculations*,
        146 pp., (1989).
    """
    # degrees and arcseconds to radians
    dtr = np.pi/180.0
    atr = np.pi/648000.0
    # create timescale from Modified Julian Day (MJD)
    ts = timescale(MJD=MJD)
    # convert from MJD to centuries relative to 2000-01-01T12:00:00
    T = ts.T
    # mean longitude of moon (p. 338)
    lunar_longitude = np.array([218.3164477, 481267.88123421, -1.5786e-3,
            1.855835e-6, -1.53388e-8])
    s = dtr*polynomial_sum(lunar_longitude, T)
    # difference between the mean longitude of sun and moon (p. 338)
    lunar_elongation = np.array([297.8501921, 445267.1114034, -1.8819e-3,
            1.83195e-6, -8.8445e-9])
    D = dtr*polynomial_sum(lunar_elongation, T)
    # mean longitude of ascending lunar node (p. 144)
    lunar_node = np.array([125.04452, -1934.136261, 2.0708e-3, 2.22222e-6])
    N = dtr*polynomial_sum(lunar_node, T)
    F = s - N
    # mean anomaly of the sun (radians)
    M = dtr*(357.5256 + 35999.049*T)
    # mean anomaly of the moon (radians)
    l = dtr*(134.96292 + 477198.86753*T)
    # series expansion for mean anomaly in moon radius (meters)
    r_moon = 1e3*(385000.0 - 20905.0*np.cos(l) - 3699.0*np.cos(2.0*D - l) -
        2956.0*np.cos(2.0*D) - 570.0*np.cos(2.0*l) +
        246.0*np.cos(2.0*l - 2.0*D) - 205.0*np.cos(M - 2.0*D) -
        171.0*np.cos(l + 2.0*D) - 152.0*np.cos(l + M - 2.0*D))
    # series expansion for ecliptic longitude of the moon (radians)
    lambda_moon = s + atr*(
        22640.0*np.sin(l) + 769.0*np.sin(2.0*l) -
        4586.0*np.sin(l - 2.0*D) + 2370.0*np.sin(2.0*D) -
        668.0*np.sin(M) - 412.0*np.sin(2.0*F) -
        212.0*np.sin(2.0*l - 2.0*D) - 206.0*np.sin(l + M - 2.0*D) +
        192.0*np.sin(l + 2.0*D) - 165.0*np.sin(M - 2.0*D) -
        148.0*np.sin(l - M) - 125.0*np.sin(D) -
        110.0*np.sin(l + M) - 55.0*np.sin(2.0*F - 2.0*D)
    )
    # series expansion for ecliptic latitude of the moon (radians)
    q = atr*(412.0*np.sin(2.0*F) + 541.0*np.sin(M))
    beta_moon = atr*(18520.0*np.sin(F + lambda_moon - s + q) -
        526.0*np.sin(F - 2*D) + 44.0*np.sin(l + F - 2.0*D) -
        31.0*np.sin(-l + F - 2.0*D) - 25.0*np.sin(-2.0*l + F) -
        23.0*np.sin(M + F - 2.0*D) + 21.0*np.sin(-l + F) +
        11.0*np.sin(-M + F - 2.0*D)
    )
    # convert to position vectors
    x = r_moon*np.cos(lambda_moon)*np.cos(beta_moon)
    y = r_moon*np.sin(lambda_moon)*np.cos(beta_moon)
    z = r_moon*np.sin(beta_moon)
    # obliquity of the J2000 ecliptic (radians)
    epsilon_j2000 = 23.43929111*dtr
    # rotate by ecliptic
    rot_x = rotate(-epsilon_j2000, 'x')
    u = rot_x[0,0,:]*x + rot_x[0,1,:]*y + rot_x[0,2,:]*z
    v = rot_x[1,0,:]*x + rot_x[1,1,:]*y + rot_x[1,2,:]*z
    w = rot_x[2,0,:]*x + rot_x[2,1,:]*y + rot_x[2,2,:]*z
    # Greenwich hour angle (radians)
    rot_z = rotate(ts.gha*dtr, 'z')
    # rotate to cartesian (ECEF) coordinates
    # ignoring polar motion and length-of-day variations
    X = rot_z[0,0,:]*u + rot_z[0,1,:]*v + rot_z[0,2,:]*w
    Y = rot_z[1,0,:]*u + rot_z[1,1,:]*v + rot_z[1,2,:]*w
    Z = rot_z[2,0,:]*u + rot_z[2,1,:]*v + rot_z[2,2,:]*w
    # return the ECEF coordinates
    return (X, Y, Z)

# PURPOSE: compute coordinates of the moon in an ECEF frame
def lunar_ephemerides(MJD: np.ndarray, **kwargs):
    """
    Computes positional coordinates of the moon in an Earth-centric,
    Earth-Fixed (ECEF) frame using JPL ephemerides [1]_ [2]_

    Parameters
    ----------
    MJD: np.ndarray
        Modified Julian Day (MJD) of input date
    kernel: str or pathlib.Path
        Path to JPL ephemerides kernel file

    Returns
    -------
    X, Y, Z: np.ndarray
        ECEF coordinates of the moon (meters)

    References
    ----------
    .. [1] J. Meeus, *Astronomical Algorithms*, 2nd edition, 477 pp., (1998).
    .. [2] R. S. Park, W. M. Folkner, and J. G. Williams, and D. H. Boggs,
        "The JPL Planetary and Lunar Ephemerides DE440 and DE441",
        *The Astronomical Journal*, 161(3), 105, (2021).
        `doi: 10.3847/1538-3881/abd414
        <https://doi.org/10.3847/1538-3881/abd414>`_
    """
    # set default keyword arguments
    kwargs.setdefault('kernel', _default_kernel)
    # degrees to radians
    dtr = np.pi/180.0
    # create timescale from Modified Julian Day (MJD)
    ts = timescale(MJD=MJD)
    # convert from MJD to Julian days (dynamical time)
    JD = ts.tt
    # read JPL ephemerides kernel
    SPK = jplephem.spk.SPK.open(kwargs['kernel'])
    # segments for computing position of the moon
    # segment 3 EARTH BARYCENTER -> segment 399 EARTH
    EMB_to_Earth = SPK[3, 399]
    # segment 3 EARTH BARYCENTER -> segment 301 MOON
    EMB_to_Moon = SPK[3, 301]
    # compute the position of the moon relative to the Earth in meters
    # Earth_to_Moon = Earth_to_EMB + EMB_to_Moon
    #               = -EMB_to_Earth + EMB_to_Moon
    x, y, z = 1e3*(EMB_to_Moon.compute(JD) - EMB_to_Earth.compute(JD))
    # Greenwich hour angle (radians)
    rot_z = rotate(ts.gha*dtr, 'z')
    # rotate to cartesian (ECEF) coordinates
    # ignoring polar motion and length-of-day variations
    X = rot_z[0,0,:]*x + rot_z[0,1,:]*y + rot_z[0,2,:]*z
    Y = rot_z[1,0,:]*x + rot_z[1,1,:]*y + rot_z[1,2,:]*z
    Z = rot_z[2,0,:]*x + rot_z[2,1,:]*y + rot_z[2,2,:]*z
    # return the ECEF coordinates
    return (X, Y, Z)

def _precession_matrix(T: float | np.ndarray):
    # arcseconds to radians
    atr = np.pi/648000.0
    # equatorial precession angles Lieske et al. (1977)
    # Capitaine et al. (2003), eqs. (4), (37), & (39).
    # obliquity of the ecliptic
    epsilon0 = 84381.406
    EPS = epsilon0 * atr
    # lunisolar precession
    phi0 = np.array([0.0, 5038.481507, -1.0790069,
        -1.14045e-3, 1.32851e-4, -9.51e-8])
    PSI = atr*polynomial_sum(phi0, T)
    # inclination of moving equator on fixed ecliptic
    omega0 = np.array([epsilon0, -2.5754e-2, 5.12623e-2,
        -7.72503e-3, -4.67e-7, 3.337e-7])
    OMEGA = atr*polynomial_sum(omega0, T)
    # planetary precession
    chi0 = np.array([0.0, 10.556403, -2.3814292,
        -1.21197e-3, 1.70663e-4, -5.60e-8])
    CHI = atr*polynomial_sum(chi0, T)
    # compute elements of precession rotation matrix
    P = np.zeros((3,3,len(np.atleast_1d(T))))
    P[0,0,:] = np.cos(CHI)*np.cos(-PSI) - \
        np.sin(-PSI)*np.sin(CHI)*np.cos(-OMEGA)
    P[0,1,:] = np.cos(CHI)*np.sin(-PSI)*np.cos(EPS) + \
        np.sin(CHI)*np.cos(-OMEGA)*np.cos(-PSI)*np.cos(EPS) - \
        np.sin(EPS)*np.sin(CHI)*np.sin(-OMEGA)
    P[0,2,:] = np.cos(CHI)*np.sin(-PSI)*np.sin(EPS) + \
        np.sin(CHI)*np.cos(-OMEGA)*np.cos(-PSI)*np.sin(EPS) + \
        np.cos(EPS)*np.sin(CHI)*np.sin(-OMEGA)
    P[1,0,:] = -np.sin(CHI)*np.cos(-PSI) - \
        np.sin(-PSI)*np.cos(CHI)*np.cos(-OMEGA)
    P[1,1,:] = -np.sin(CHI)*np.sin(-PSI)*np.cos(EPS) + \
        np.cos(CHI)*np.cos(-OMEGA)*np.cos(-PSI)*np.cos(EPS) - \
        np.sin(EPS)*np.cos(CHI)*np.sin(-OMEGA)
    P[1,2,:] = -np.sin(CHI)*np.sin(-PSI)*np.sin(EPS) + \
        np.cos(CHI)*np.cos(-OMEGA)*np.cos(-PSI)*np.sin(EPS) + \
        np.cos(EPS)*np.cos(CHI)*np.sin(-OMEGA)
    P[2,0,:] = np.sin(-PSI)*np.sin(-OMEGA)
    P[2,1,:] = -np.sin(-OMEGA)*np.cos(-PSI)*np.cos(EPS) - \
        np.sin(EPS)*np.cos(-OMEGA)
    P[2,2,:] = -np.sin(-OMEGA)*np.cos(-PSI)*np.sin(EPS) + \
        np.cos(-OMEGA)*np.cos(EPS)
    # return the rotation matrix
    return P

def _nutation_matrix(
        mean_obliquity: float | np.ndarray,
        true_obliquity: float | np.ndarray,
        psi: float | np.ndarray
    ):
    # compute elements of nutation rotation matrix
    R = np.zeros((3,3,len(np.atleast_1d(psi))))
    R[0,0,:] = np.cos(psi)
    R[0,1,:] = -np.sin(psi)*np.cos(mean_obliquity)
    R[0,2,:] = -np.sin(psi)*np.sin(mean_obliquity)
    R[1,0,:] = np.sin(psi)*np.cos(true_obliquity)
    R[1,1,:] = np.cos(psi)*np.cos(mean_obliquity)*np.cos(true_obliquity) + \
        np.sin(mean_obliquity)*np.sin(true_obliquity)
    R[1,2,:] = np.cos(psi)*np.sin(mean_obliquity)*np.cos(true_obliquity) - \
        np.cos(mean_obliquity)*np.sin(true_obliquity)
    R[2,0,:] = np.sin(psi)*np.sin(true_obliquity)
    R[2,1,:] = np.cos(psi)*np.cos(mean_obliquity)*np.sin(true_obliquity) - \
        np.sin(mean_obliquity)*np.cos(true_obliquity)
    R[2,2,:] = np.cos(psi)*np.sin(mean_obliquity)*np.sin(true_obliquity) + \
        np.cos(mean_obliquity)*np.cos(true_obliquity)
    # return the rotation matrix
    return R

def _frame_bias_matrix():
    # arcseconds to radians
    atr = np.pi/648000.0
    xi0  = -0.0166170*atr
    eta0 = -0.0068192*atr
    da0  = -0.01460*atr
    # compute elements of the frame bias matrix
    B = np.zeros((3,3))
    B[0,1] = da0
    B[0,2] = -xi0
    B[1,0] = -da0
    B[1,2] = -eta0
    B[2,0] =  xi0
    B[2,1] =  eta0
    # second-order corrections to diagonal elements
    B[0,0] = 1.0 - 0.5 * (da0**2 + xi0**2)
    B[1,1] = 1.0 - 0.5 * (da0**2 + eta0**2)
    B[2,2] = 1.0 - 0.5 * (eta0**2 + xi0**2)
    # return the rotation matrix
    return B

def _polar_motion_matrix(T):
    # arcseconds to radians
    atr = np.pi/648000.0
    # convert to MJD from centuries relative to 2000-01-01T12:00:00
    MJD = T*36525.0 + 51544.5
    sprime = -4.7e-5*T
    px, py = iers_polar_motion(MJD)
    return mxmxm(rotate(py*atr,'x'), rotate(px*atr,'y'), rotate(-sprime*atr,'z'))

def _parse_table_5_2e():
    """Parse table with expressions for Greenwich Sidereal Time
    provided in `Chapter 5 of IERS Conventions
    <https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.2e.txt>`_
    """
    table_5_2e = get_data_path(['data','tab5.2e.txt'])
    with table_5_2e.open(mode='r', encoding='utf8') as f:
        file_contents = f.readlines()
    # names and formats
    names = ('i','Cs','Cc','l','lp','F','D','Om','L_Me','L_Ve',
        'L_E','L_Ma','L_J','L_Sa','L_U','L_Ne','p_A')
    formats = ('i','f','f','i','i','i','i','i','i',
        'i','i','i','i','i','i','i','i')
    dtype = np.dtype({'names':names, 'formats':formats})
    # j = 0 terms
    n0 = 33
    j0 = np.zeros((n0), dtype=dtype)
    for i,line in enumerate(file_contents[53:53+n0]):
        j0[i] = np.array(tuple(line.split()), dtype=dtype)
    # j = 1 terms
    n1 = 1
    j1 = np.zeros((n1), dtype=dtype)
    for i,line in enumerate(file_contents[90:90+n1]):
        j1[i] = np.array(tuple(line.split()), dtype=dtype)
    # return the table
    return (j0, j1)

def _parse_table_5_3a():
    """Parse table with IAU 2000A lunisolar and planetary components
    of nutation in longitude provided in `Chapter 5 of IERS Conventions
    <https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.3a.txt>`_
    """
    table_5_3a = get_data_path(['data','tab5.3a.txt'])
    with table_5_3a.open(mode='r', encoding='utf8') as f:
        file_contents = f.readlines()
    # names and formats
    names = ('i','As','Ac','l','lp','F','D','Om','L_Me','L_Ve',
        'L_E','L_Ma','L_J','L_Sa','L_U','L_Ne','p_A')
    formats = ('i','f','f','i','i','i','i','i','i',
        'i','i','i','i','i','i','i','i')
    dtype = np.dtype({'names':names, 'formats':formats})
    # j = 0 terms
    n0 = 1320
    j0 = np.zeros((n0), dtype=dtype)
    for i,line in enumerate(file_contents[22:22+n0]):
        j0[i] = np.array(tuple(line.split()), dtype=dtype)
    # j = 1 terms
    n1 = 38
    j1 = np.zeros((n1), dtype=dtype)
    for i,line in enumerate(file_contents[1348:1348+n1]):
        j1[i] = np.array(tuple(line.split()), dtype=dtype)
    # return the table
    return (j0, j1)

def _parse_table_5_3b():
    """Parse table with IAU 2000A lunisolar and planetary components
    of nutation in obliquity provided in `Chapter 5 of IERS Conventions
    <https://iers-conventions.obspm.fr/content/chapter5/additional_info/tab5.3b.txt>`_
    """
    table_5_3b = get_data_path(['data','tab5.3b.txt'])
    with table_5_3b.open(mode='r', encoding='utf8') as f:
        file_contents = f.readlines()
    # names and formats
    names = ('i','Bc','Bs','l','lp','F','D','Om','L_Me','L_Ve',
        'L_E','L_Ma','L_J','L_Sa','L_U','L_Ne','p_A')
    formats = ('i','f','f','i','i','i','i','i','i',
        'i','i','i','i','i','i','i','i')
    dtype = np.dtype({'names':names, 'formats':formats})
    # j = 0 terms
    n0 = 1037
    j0 = np.zeros((n0), dtype=dtype)
    for i,line in enumerate(file_contents[22:22+n0]):
        j0[i] = np.array(tuple(line.split()), dtype=dtype)
    # j = 1 terms
    n1 = 19
    j1 = np.zeros((n1), dtype=dtype)
    for i,line in enumerate(file_contents[1065:1065+n1]):
        j1[i] = np.array(tuple(line.split()), dtype=dtype)
    # return the table
    return (j0, j1)
