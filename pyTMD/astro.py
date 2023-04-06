#!/usr/bin/env python
u"""
astro.py (04/2023)
Modification of ASTROL fortran subroutine by Richard Ray 03/1999

Computes the basic astronomical mean longitudes: s, h, p, N and PP
Note N is not N', i.e. N is decreasing with time.

Formulae for the period 1990--2010 were derived by David Cartwright
MEEUS and ASTRO5 formulae are from versions of Meeus's Astronomical Algorithms

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

REFERENCES:
    Jean Meeus, Astronomical Algorithms, 2nd edition, 1998.
    Oliver Montenbruck, Practical Ephemeris Calculations, 1989.

UPDATE HISTORY:
    Updated 04/2023: added low resolution solar and lunar ephemerides
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

import numpy as np

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
    return np.sum([c * (t ** i) for i,c in enumerate(coefficients)],axis=0)

# PURPOSE: compute the basic astronomical mean longitudes
def mean_longitudes(
        MJD: np.ndarray,
        MEEUS: bool = False,
        ASTRO5: bool = False
    ):
    """
    Computes the basic astronomical mean longitudes:
    `s`, `h`, `p`, `N` and `PP` [1]_

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
    s: np.ndarray
        mean longitude of moon (degrees)
    h: np.ndarray
        mean longitude of sun (degrees)
    p: np.ndarray
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
        s = polynomial_sum(lunar_longitude,T)
        # mean longitude of sun
        solar_longitude = np.array([280.46645, 0.985647360164271,
            2.2727347e-13])
        h = polynomial_sum(solar_longitude,T)
        # mean longitude of lunar perigee
        lunar_perigee = np.array([83.3532430, 0.11140352391786447,
            -7.7385418e-12, -2.5636086e-19, 2.95738836e-26])
        p = polynomial_sum(lunar_perigee,T)
        # mean longitude of ascending lunar node
        lunar_node = np.array([125.0445550, -0.052953762762491446,
            1.55628359e-12, 4.390675353e-20, -9.26940435e-27])
        N = polynomial_sum(lunar_node,T)
        # mean longitude of solar perigee (Simon et al., 1994)
        PP = 282.94 + 1.7192 * T / 36525.0
    elif ASTRO5:
        # convert from MJD to centuries relative to 2000-01-01T12:00:00
        T = (MJD - 51544.5)/36525.0
        # mean longitude of moon (p. 338)
        lunar_longitude = np.array([218.3164477, 481267.88123421, -1.5786e-3,
             1.855835e-6, -1.53388e-8])
        s = polynomial_sum(lunar_longitude,T)
        # mean longitude of sun (p. 338)
        lunar_elongation = np.array([297.8501921, 445267.1114034, -1.8819e-3,
             1.83195e-6, -8.8445e-9])
        h = polynomial_sum(lunar_longitude-lunar_elongation,T)
        # mean longitude of lunar perigee (p. 343)
        lunar_perigee = np.array([83.3532465, 4069.0137287, -1.032e-2,
            -1.249172e-5])
        p = polynomial_sum(lunar_perigee,T)
        # mean longitude of ascending lunar node (p. 144)
        lunar_node = np.array([125.04452, -1934.136261, 2.0708e-3, 2.22222e-6])
        N = polynomial_sum(lunar_node,T)
        # mean longitude of solar perigee (Simon et al., 1994)
        PP = 282.94 + 1.7192 * T
    else:
        # convert from MJD to days relative to 2000-01-01T12:00:00
        # convert from Universal Time to Dynamic Time at 2000-01-01
        T = MJD - 51544.4993
        # mean longitude of moon
        s = 218.3164 + 13.17639648 * T
        # mean longitude of sun
        h = 280.4661 + 0.98564736 * T
        # mean longitude of lunar perigee
        p = 83.3535 + 0.11140353 * T
        # mean longitude of ascending lunar node
        N = 125.0445 - 0.05295377 * T
        # solar perigee at epoch 2000
        PP = 282.8

    # take the modulus of each
    s = np.mod(s, circle)
    h = np.mod(h, circle)
    p = np.mod(p, circle)
    N = np.mod(N, circle)

    # return as tuple
    return (s, h, p, N, PP)

# PURPOSE: compute coordinates of the sun in an ECEF frame
def solar_ecef(MJD: np.ndarray):
    """
    Computes low-resolution coordinates of the sun in an
    Earth-centric, Earth-Fixed (ECEF) frame [1]_ [2]_

    Parameters
    ----------
    MJD: np.ndarray
        Modified Julian Day (MJD) of input date

    Returns
    -------
    xyz: np.ndarray
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
    # convert from MJD to centuries relative to 2000-01-01T12:00:00
    T = (MJD - 51544.5)/36525.0
    # mean longitude of solar perigee (radians)
    PP = dtr*(282.94 + 1.7192 * T)
    # mean anomaly of the sun (radians)
    M = dtr*(357.5256 + 35999.049*T)
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
    # Greenwich Mean Sidereal Time (degrees)
    GMST = np.mod(280.46061837504 + 360.9856473662862*(T*36525.0), 360.0)
    # rotate to cartesian (ECEF) coordinates
    # ignoring polar motion and length-of-day variations
    X = np.cos(GMST*dtr)*x + np.sin(GMST*dtr)*y
    Y = np.cos(GMST*dtr)*y - np.sin(GMST*dtr)*x
    # return the ECEF coordinates
    return (X, Y, z)

# PURPOSE: compute coordinates of the moon in an ECEF frame
def lunar_ecef(MJD: np.ndarray):
    """
    Computes low-resolution coordinates of the moon in an
    Earth-centric, Earth-Fixed (ECEF) frame [1]_ [2]_

    Parameters
    ----------
    MJD: np.ndarray
        Modified Julian Day (MJD) of input date

    Returns
    -------
    xyz: np.ndarray
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
    # convert from MJD to centuries relative to 2000-01-01T12:00:00
    T = (MJD - 51544.5)/36525.0
    # mean longitude of moon (p. 338)
    lunar_longitude = np.array([218.3164477, 481267.88123421, -1.5786e-3,
            1.855835e-6, -1.53388e-8])
    s = dtr*polynomial_sum(lunar_longitude,T)
    # difference between the mean longitude of sun and moon (p. 338)
    lunar_elongation = np.array([297.8501921, 445267.1114034, -1.8819e-3,
            1.83195e-6, -8.8445e-9])
    D = dtr*polynomial_sum(lunar_elongation, T)
    # mean longitude of ascending lunar node (p. 144)
    lunar_node = np.array([125.04452, -1934.136261, 2.0708e-3, 2.22222e-6])
    N = dtr*polynomial_sum(lunar_node,T)
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
    v = y*np.cos(-epsilon_j2000) + z*np.sin(-epsilon_j2000)
    w = z*np.cos(-epsilon_j2000) - y*np.sin(-epsilon_j2000)
    # Greenwich Mean Sidereal Time (degrees)
    GMST = np.mod(280.46061837504 + 360.9856473662862*(T*36525.0), 360.0)
    # rotate to cartesian (ECEF) coordinates
    # ignoring polar motion and length-of-day variations
    X = np.cos(GMST*dtr)*x + np.sin(GMST*dtr)*v
    Y = np.cos(GMST*dtr)*v - np.sin(GMST*dtr)*x
    # return the ECEF coordinates
    return (X, Y, w)
