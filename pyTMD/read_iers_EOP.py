#!/usr/bin/env python
u"""
read_iers_EOP.py
Written by Tyler Sutterley (12/2022)
Provides the daily earth orientation parameters (EOP) from IERS
    http://www.usno.navy.mil/USNO/earth-orientation/eo-products/weekly
Data format: http://maia.usno.navy.mil/ser7/readme.finals

INPUTS:
    input_file: full path to IERS EOP "finals" file

OUTPUTS:
    MJD: modified julian date of EOP measurements
    x: Angular coordinate x [arcsec]
    y: Angular coordinate y [arcsec]

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

REFERENCE:
    Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
        IERS Technical Note No. 36, BKG (2010)

UPDATE HISTORY:
    Updated 12/2022: refactor daily orientation as part of eop module
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
    Updated 07/2021: added check that IERS finals file is accessible
    Updated 03/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 07/2020: added function docstrings
    Written 09/2017
"""
import warnings
import pyTMD.eop

# calculates the daily earth orientation parameters (EOP) from IERS
def read_iers_EOP(input_file):
    """
    Calculates the daily earth orientation parameters (EOP) from IERS

    Parameters
    ----------
    input_file: str
        full path to IERS EOP "finals" file

    Returns
    -------
    MJD: float
        modified Julian date of EOP measurements
    x: float
        Angular coordinate x [arcsec]
    y: float
        Angular coordinate y [arcsec]

    References
    ----------
    .. [1] G. Petit and B. Luzum (eds.), *IERS Conventions (2010)*,
        International Earth Rotation and Reference Systems Service (IERS),
        `IERS Technical Note No. 36 <https://iers-conventions.obspm.fr/content/tn36.pdf>`_
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.eop instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.eop.iers_daily_EOP(input_file)
