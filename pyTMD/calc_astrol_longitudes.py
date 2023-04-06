#!/usr/bin/env python
u"""
calc_astrol_longitudes.py (04/2023)
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
    Updated 04/2023: deprecated in favor of pyTMD.astro.mean_longitudes
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

import warnings
import numpy as np
import pyTMD.astro

# PURPOSE: compute the basic astronomical mean longitudes
def calc_astrol_longitudes(*args, **kwargs):
    """
    Computes the basic astronomical mean longitudes:
    `s`, `h`, `p`, `N` and `PP`

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
    warnings.filterwarnings("module")
    warnings.warn("Deprecated. Please use pyTMD.astro instead",
        DeprecationWarning)
    warnings.filterwarnings("ignore")
    # call updated function to not break current workflows
    return pyTMD.astro.mean_longitudes(*args, **kwargs)
