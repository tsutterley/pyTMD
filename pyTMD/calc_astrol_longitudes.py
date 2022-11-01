#!/usr/bin/env python
u"""
calc_astrol_longitudes.py (10/2022)
Modification of ASTROL fortran subroutine by Richard Ray 03/1999

Computes the basic astronomical mean longitudes: s, h, p, N and PP
Note N is not N', i.e. N is decreasing with time.

Formulae for the period 1990--2010 were derived by David Cartwright
MEEUS and ASTRO5 formulae are from versions of Meeus's Astronomical Algorithms

CALLING SEQUENCE:
    s,h,p,N,PP = calc_astrol_longitudes(MJD, ASTRO5=True)

INPUTS:
    MJD: Modified Julian Day of input date

OUTPUTS:
    s: mean longitude of moon (degrees)
    h: mean longitude of sun (degrees)
    p: mean longitude of lunar perigee (degrees)
    N: mean longitude of ascending lunar node (degrees)
    PP: longitude of solar perigee (degrees)

OPTIONS:
    MEEUS: use additional coefficients from Meeus Astronomical Algorithms
    ASTRO5: use Meeus Astronomical coefficients as implemented in ASTRO5

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

REFERENCES:
    Jean Meeus, Astronomical Algorithms, 2nd edition, 1998.

UPDATE HISTORY:
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
import numpy as np

# PURPOSE: calculate the sum of a polynomial function of time
def polynomial_sum(coefficients, t):
    """
    Calculates the sum of a polynomial function of time

    Parameters
    ----------
    coefficients: list
        leading coefficient of polynomials of increasing order
    t: float
        delta time in units for a given astronomical longitudes calculation
    """
    # convert time to array if importing a single value
    t = np.atleast_1d(t)
    return np.sum([c * (t ** i) for i,c in enumerate(coefficients)],axis=0)

# PURPOSE: compute the basic astronomical mean longitudes
def calc_astrol_longitudes(MJD, MEEUS=False, ASTRO5=False):
    """
    Computes the basic astronomical mean longitudes: s, h, p, N and PP
    [Meeus1998]_

    Parameters
    ----------
    MJD: float
        Modified Julian Day (MJD) of input date
    MEEUS: bool, default False
        use additional coefficients from Meeus Astronomical Algorithms
    ASTRO5: bool, default False
        use Meeus Astronomical coefficients as implemented in ASTRO5

    Returns
    -------
    s: float
        mean longitude of moon (degrees)
    h: float
        mean longitude of sun (degrees)
    p: float
        mean longitude of lunar perigee (degrees)
    N: float
        mean longitude of ascending lunar node (degrees)
    PP: float
        longitude of solar perigee (degrees)

    References
    ----------
    .. [Meeus1998] J. Meeus, *Astronomical Algorithms*, 2nd edition, 477 pp., (1998).
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
        p =  83.3535 + 0.11140353 * T
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
