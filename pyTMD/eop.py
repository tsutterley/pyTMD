#!/usr/bin/env python
u"""
eop.py
Written by Tyler Sutterley (04/2024)
Utilities for maintaining and calculating Earth Orientation Parameters (EOP)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml
    timescale: Python tools for time and astronomical calculations
        https://pypi.org/project/timescale/

UPDATE HISTORY:
    Updated 04/2024: deprecated in favor of timescale.eop
    Updated 04/2023: using pathlib to define and expand paths
        add wrapper function for interpolating daily EOP values
        have mean pole and finals file as attributes of EOP module
    Updated 03/2023: add secular pole model from IERS 2018 conventions
    Updated 11/2022: added encoding for writing ascii files
        use f-strings for formatting verbose or ascii output
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 03/2021: replaced numpy bool/int to prevent deprecation warnings
    Written 11/2020
"""
from __future__ import annotations

import warnings
import logging
import timescale.eop
import timescale.utilities

# IERS mean pole file for 2015 conventional mean pole
_mean_pole_file = timescale.utilities.get_data_path(['data','mean-pole.tab'])
# daily finals polar motion file from IERS
_finals_file = timescale.utilities.get_data_path(['data','finals.all'])

# PURPOSE: connects to servers and downloads mean pole files
def update_mean_pole(**kwargs):
    """
    Connects to servers to download mean-pole.tab files from HPIERS servers

    - ftp://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.readme

    Servers and Mirrors

    - ftp://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.tab
    - http://hpiers.obspm.fr/eoppc/eop/eopc01/mean-pole.tab

    Parameters
    ----------
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file
    """
    warnings.warn("Deprecated. Please use timescale.eop.update_mean_pole",
        DeprecationWarning)
    timescale.eop.update_mean_pole(**kwargs)

# PURPOSE: read table of IERS pole coordinates and calculate Gaussian average
def calculate_mean_pole(**kwargs):
    """
    Calculates the mean pole coordinates x and y are obtained by a
    Gaussian-weighted average of the IERS pole coordinates

    - ftp://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.readme

    Servers and Mirrors

    - ftp://ftp.iers.org/products/eop/long-term/c01/eopc01.iau2000.1900-now.dat
    - ftp://hpiers.obspm.fr/iers/eop/eopc01/eopc01.iau2000.1900-now.dat
    - http://hpiers.obspm.fr/eoppc/eop/eopc01/eopc01.iau2000.1900-now.dat

    Parameters
    ----------
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file
    """
    warnings.warn("Deprecated. Please use timescale.eop.calculate_mean_pole",
        DeprecationWarning)
    timescale.eop.calculate_mean_pole(**kwargs)

# PURPOSE: connects to servers and downloads IERS pole coordinates files
def pull_pole_coordinates(*args, **kwargs):
    """
    Connects to servers and downloads IERS pole coordinate files

    Servers and Mirrors

    - ftp://ftp.iers.org/products/eop/long-term/c01/eopc01.iau2000.1900-now.dat
    - ftp://hpiers.obspm.fr/iers/eop/eopc01/eopc01.iau2000.1900-now.dat
    - http://hpiers.obspm.fr/eoppc/eop/eopc01/eopc01.iau2000.1900-now.dat

    Parameters
    ----------
    FILE: str
        IERS pole coordinate file to download from remote servers

            - eopc01.1846-now.dat
            - eopc01.1900-now.dat
            - eopc01.iau2000.1900-now.dat
            - eopc01.iau2000.1846-now.dat
    verbose: bool, default False
        print file information about output file
    """
    warnings.warn("Deprecated. Please use timescale.eop.pull_pole_coordinates",
        DeprecationWarning)
    timescale.eop.pull_pole_coordinates(*args, **kwargs)

# PURPOSE: connects to servers and downloads finals files
def update_finals_file(**kwargs):
    """
    Connects to servers and downloads finals EOP files

    Servers and Mirrors

    - http://maia.usno.navy.mil/ser7/
    - https://cddis.nasa.gov/archive/products/iers/
    - ftp://cddis.nasa.gov/products/iers/
    - ftp://cddis.gsfc.nasa.gov/pub/products/iers/

    Parameters
    ----------
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    timeout: int or NoneType, default 20
        timeout in seconds for blocking operations
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file
    """
    warnings.warn("Deprecated. Please use timescale.eop.pull_pole_coordinates",
        DeprecationWarning)
    timescale.eop.update_finals_file(**kwargs)

# IERS mean or secular pole conventions
_conventions = ('2003', '2010', '2015', '2018')
# read table of mean pole values, calculate angular coordinates at epoch
def iers_mean_pole(*args, **kwargs):
    """
    Calculates the angular coordinates of the IERS Conventional Mean Pole (CMP)
    or IERS Secular Pole (2018 convention)

    Parameters
    ----------
    input_epoch: np.ndarray
        Dates for the angular coordinates of the Conventional Mean Pole
        in decimal years
    convention: str, default '2018'
        IERS Mean or Secular Pole Convention

            - ``'2003'``
            - ``'2010'``
            - ``'2015'``
            - ``'2018'``
    input_file: str or pathlib.Path
        Full path to mean-pole.tab file provided by IERS
    fill_value: float, default np.nan
        Value for invalid flags

    Returns
    -------
    x: np.ndarray
        Angular coordinate x of conventional mean pole or secular pole
    y: np.ndarray
        Angular coordinate y of conventional mean pole or secular pole
    flag: np.ndarray
        epoch is valid for version and version number is valid

    References
    ----------
    .. [1] G. Petit and B. Luzum (eds.), *IERS Conventions (2010)*,
        International Earth Rotation and Reference Systems Service (IERS),
        `IERS Technical Note No. 36
        <https://iers-conventions.obspm.fr/content/tn36.pdf>`_
    """
    warnings.warn("Deprecated. Please use timescale.eop.iers_mean_pole",
        DeprecationWarning)
    return timescale.eop.iers_mean_pole(*args, **kwargs)

# PURPOSE: read daily earth orientation parameters (EOP) file from IERS
def iers_daily_EOP(**kwargs):
    """
    Read daily earth orientation parameters (EOP) file from IERS

    Parameters
    ----------
    input_file: str or Pathlib.Path
        full path to IERS EOP "finals" file

    Returns
    -------
    MJD: np.ndarray
        Modified Julian Date of EOP measurements
    x: np.ndarray
        Angular coordinate x [arcsec]
    y: np.ndarray
        Angular coordinate y [arcsec]

    References
    ----------
    .. [1] G. Petit and B. Luzum (eds.), *IERS Conventions (2010)*,
        International Earth Rotation and Reference Systems Service (IERS),
        `IERS Technical Note No. 36
        <https://iers-conventions.obspm.fr/content/tn36.pdf>`_
    """
    warnings.warn("Deprecated. Please use timescale.eop.iers_daily_EOP",
        DeprecationWarning)
    return timescale.eop.iers_daily_EOP(**kwargs)

def iers_polar_motion(*args, **kwargs):
    """
    Interpolates daily earth orientation parameters (EOP) file from IERS

    Parameters
    ----------
    MJD: np.ndarray
        Modified Julian Date for interpolated measurements
    file: str or Pathlib.Path
        default path to IERS EOP "finals" file
    k: int
        Degree of the spline fit
    s: int or float
        Positive smoothing factor for the spline fit

    Returns
    -------
    px: np.ndarray
        Angular coordinate x [arcsec]
    py: np.ndarray
        Angular coordinate y [arcsec]

    References
    ----------
    .. [1] G. Petit and B. Luzum (eds.), *IERS Conventions (2010)*,
        International Earth Rotation and Reference Systems Service (IERS),
        `IERS Technical Note No. 36
        <https://iers-conventions.obspm.fr/content/tn36.pdf>`_
    """
    warnings.warn("Deprecated. Please use timescale.eop.iers_polar_motion",
        DeprecationWarning)
    return timescale.eop.iers_polar_motion(*args, **kwargs)
