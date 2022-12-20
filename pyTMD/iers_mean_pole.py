#!/usr/bin/env python
u"""
iers_mean_pole.py
Written by Tyler Sutterley (12/2022)
Provides the angular coordinates of the IERS Conventional Mean Pole (CMP)
Coordinates are based on the table of values from
    ftp://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.tab
Based on IERS_CMP_YYYY.f located at ftp://tai.bipm.org/iers/convupdt/chapter7/
    http://iers-conventions.obspm.fr/

INPUTS:
    input_file: full path to mean-pole.tab file provided by IERS
    input_epoch: dates for which the angular coordinates of the Conventional
        Mean Pole are desired in decimal years
    version: Year of the conventional model.
        Limited to values of '2003', '2010', '2015'

OUTPUTS:
    x: Angular coordinate x of conventional mean pole [arcsec]
    y: Angular coordinate y of conventional mean pole [arcsec]
    flag: epoch is valid for version and version number is valid
        data will be set to fill_value if flag == False

OPTIONS:
    fill_value: value for invalid flags

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

REFERENCE:
    Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
        IERS Technical Note No. 36, BKG (2010)

UPDATE HISTORY:
    Updated 12/2022: refactor mean pole as part of eop module
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 07/2020: added function docstrings
    Updated 12/2018: using future division for python3 compatibility
    Written 09/2017
"""
from __future__ import division

import copy
import warnings
import numpy as np
import pyTMD.eop

# read table of mean pole values, calculate angular coordinates at epoch
def iers_mean_pole(input_file, input_epoch, version, **kwargs):
    """
    Calculates the angular coordinates of the IERS Conventional Mean Pole (CMP)

    Parameters
    ----------
    input_file: str
        Full path to mean-pole.tab file provided by IERS
    input_epoch: float
        Dates for the angular coordinates of the Conventional Mean Pole
        in decimal years
    version: str
        Year of the conventional model
    fill_value: float, default np.nan
        Value for invalid flags

    Returns
    -------
    x: float
        Angular coordinate x of conventional mean pole
    y: float
        Angular coordinate y of conventional mean pole
    flag: bool
        epoch is valid for version and version number is valid

    References
    ----------
    .. [1] G. Petit and B. Luzum (eds.), *IERS Conventions (2010)*,
        International Earth Rotation and Reference Systems Service (IERS),
        `IERS Technical Note No. 36 <https://iers-conventions.obspm.fr/content/tn36.pdf>`_
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.eop instead",DeprecationWarning)
    # set default keyword arguments
    kwargs.setdefault('fill_value', np.nan)
    # raise warnings for deprecated keyword argument
    deprecated_keywords = dict(FILL_VALUE='fill_value')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn(f"""Deprecated keyword argument {old}.
                Changed to '{new}'""", DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])
    # call renamed version to not break workflows
    return pyTMD.eop.iers_mean_pole(input_file, input_epoch, version, **kwargs)
