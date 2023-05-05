#!/usr/bin/env python
u"""
eop.py
Written by Tyler Sutterley (04/2023)
Utilities for maintaining and calculating Earth Orientation Parameters (EOP)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
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

import logging
import pathlib
import traceback
import numpy as np
import scipy.interpolate
import pyTMD.utilities

# IERS mean pole file for 2015 conventional mean pole
_mean_pole_file = pyTMD.utilities.get_data_path(['data','mean-pole.tab'])
# daily polar motion file from IERS
_finals_file = pyTMD.utilities.get_data_path(['data','finals.all'])

# PURPOSE: connects to servers and downloads mean pole files
def update_mean_pole(verbose: bool = False, mode: oct = 0o775):
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
    # check hash of local version of file
    LOCAL = _mean_pole_file
    HASH = pyTMD.utilities.get_hash(_mean_pole_file)

    # try downloading from Paris Observatory IERS Centers ftp servers
    HOST = ['hpiers.obspm.fr', 'iers', 'eop', 'eopc01', 'mean-pole.tab']
    try:
        pyTMD.utilities.from_ftp(HOST, timeout=20, local=LOCAL,
            hash=HASH, verbose=verbose, mode=mode)
    except Exception as exc:
        logging.debug(traceback.format_exc(exc))
        pass
    else:
        return

    # try downloading from Paris Observatory IERS Centers https servers
    HOST = ['http://hpiers.obspm.fr', 'eoppc', 'eop', 'eopc01', 'mean-pole.tab']
    try:
        pyTMD.utilities.from_http(HOST, timeout=20, local=LOCAL,
            hash=HASH, verbose=verbose, mode=mode)
    except Exception as exc:
        logging.debug(traceback.format_exc(exc))
        pass
    else:
        return

    # raise exception
    raise RuntimeError(f'Unable to download {LOCAL}')

# PURPOSE: read table of IERS pole coordinates and calculate Gaussian average
def calculate_mean_pole(verbose: bool = False, mode: oct = 0o775):
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
    # download the IERS pole coordinates file from remote servers
    FILE = 'eopc01.1900-now.dat'
    try:
        remote_buffer = pull_pole_coordinates(FILE, verbose=verbose)
    except Exception as exc:
        return

    # read contents from input file object
    file_contents = remote_buffer.read().decode('utf8').splitlines()
    header = file_contents[0][1:].split()
    nlines = len(file_contents) - 1
    data = {h:np.zeros((nlines)) for h in header}
    # extract data for all lines
    for i, line in enumerate(file_contents[1:]):
        line_contents = line.split()
        for h, l in zip(header, line_contents):
            data[h][i] = np.float64(l)
    # output mean pole coordinates
    xm = np.zeros((nlines))
    ym = np.zeros((nlines))
    # output file with mean pole coordinates
    LOCAL = _mean_pole_file
    fid = LOCAL.open(mode='w', encoding='utf8')
    logging.info(str(LOCAL))
    for i, T in enumerate(data['an']):
        # mean pole is Gaussian Weight of all dates with a = 3.40 years.
        Wi = np.exp(-0.5*((data['an'] - T)/3.4)**2)
        xm[i] = np.sum(Wi*data['x(")'])/np.sum(Wi)
        ym[i] = np.sum(Wi*data['y(")'])/np.sum(Wi)
        print(f'{T:6.2f} {xm[i]:11.7f} {ym[i]:11.7f}', file=fid)
    # close the output file
    fid.close()
    # change the permissions mode of the output mean pole file
    LOCAL.chmod(mode)

# PURPOSE: connects to servers and downloads IERS pole coordinates files
def pull_pole_coordinates(FILE: str, verbose: bool = False):
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
    # try downloading from IERS ftp server
    HOST = ['ftp.iers.org', 'products', 'eop', 'long-term', 'c01', FILE]
    try:
        buffer = pyTMD.utilities.from_ftp(HOST, verbose=verbose, timeout=20)
    except Exception as exc:
        pass
    else:
        return buffer

    # try downloading from Paris Observatory IERS Centers ftp servers
    HOST = ['hpiers.obspm.fr', 'iers', 'eop', 'eopc01', FILE]
    try:
        buffer = pyTMD.utilities.from_ftp(HOST, verbose=verbose, timeout=20)
    except Exception as exc:
        pass
    else:
        return buffer

    # try downloading from Paris Observatory IERS Centers https servers
    HOST = ['http://hpiers.obspm.fr', 'eoppc', 'eop', 'eopc01', FILE]
    try:
        buffer = pyTMD.utilities.from_http(HOST, verbose=verbose, timeout=20)
    except Exception as exc:
        pass
    else:
        return buffer

    # raise exception
    raise RuntimeError(f'Unable to download {FILE}')

# PURPOSE: connects to servers and downloads finals files
def update_finals_file(
        username: str or None = None,
        password: str or None = None,
        timeout: int or None = 20,
        verbose: bool = False,
        mode: oct = 0o775
    ):
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
    # local version of file
    LOCAL = _finals_file
    HASH = pyTMD.utilities.get_hash(LOCAL)

    # try downloading from US Naval Oceanography Portal
    HOST = ['http://maia.usno.navy.mil', 'ser7', 'finals.all']
    try:
        pyTMD.utilities.from_http(HOST,
            timeout=timeout,
            local=LOCAL,
            hash=HASH,
            verbose=verbose,
            mode=mode
        )
    except Exception as exc:
        pass
    else:
        return

    # try downloading from NASA Crustal Dynamics Data Information System
    # note: anonymous ftp access will be discontinued on 2020-10-31
    # will require using the following https Earthdata server after that date
    server = []
    server.append(['cddis.nasa.gov', 'pub', 'products', 'iers', 'finals.all'])
    server.append(['cddis.gsfc.nasa.gov', 'products', 'iers', 'finals.all'])
    for HOST in server:
        try:
            pyTMD.utilities.from_ftp(HOST,
                timeout=timeout,
                local=LOCAL,
                hash=HASH,
                verbose=verbose,
                mode=mode)
        except Exception as exc:
            pass
        else:
            return

    # try downloading from NASA Crustal Dynamics Data Information System
    # using NASA Earthdata credentials stored in netrc file
    HOST = ['https://cddis.nasa.gov', 'archive', 'products', 'iers', 'finals.all']
    try:
        pyTMD.utilities.from_cddis(HOST,
            username=username,
            password=password,
            timeout=timeout,
            local=LOCAL,
            hash=HASH,
            verbose=verbose,
            mode=mode
        )
    except Exception as exc:
        pass
    else:
        return

# IERS mean or secular pole conventions
_conventions = ('2003', '2010', '2015', '2018')
# read table of mean pole values, calculate angular coordinates at epoch
def iers_mean_pole(
        input_epoch: np.ndarray,
        convention: str = '2018',
        **kwargs
    ):
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
    # set default keyword arguments
    kwargs.setdefault('file', _mean_pole_file)
    kwargs.setdefault('fill_value', np.nan)
    # verify IERS model version
    assert convention in _conventions, "Incorrect IERS model convention"
    # read the conventional mean pole file
    if (convention == '2015'):
        # read mean pole file
        input_file = pathlib.Path(kwargs['file']).expanduser().absolute()
        table = np.loadtxt(input_file)
        # reduce to 1971 to end date
        ii, = np.nonzero(table[:, 0] >= 1971)
        table = np.copy(table[ii,:])
        # reduce to yearly values
        jj, = np.nonzero((table[:, 0] % 1) == 0.0)
        table = np.copy(table[jj,:])
        end_time = table[-1, 0] + 0.2
        # final shape of the table
        nrows, *_ = np.shape(table)
    else:
        end_time = np.inf
    # allocate for output arrays
    x = np.full_like(input_epoch, kwargs['fill_value'])
    y = np.full_like(input_epoch, kwargs['fill_value'])
    flag = np.zeros_like(input_epoch, dtype=bool)
    for t, epoch in enumerate(input_epoch):
        # Conventional mean pole model in IERS Conventions 2003
        if (convention == '2003') and (epoch >= 1975):
            x[t] = 0.054 + 0.00083*(epoch - 2000.0)
            y[t] = 0.357 + 0.00395*(epoch - 2000.0)
            flag[t] = True
        # Conventional mean pole model in IERS Conventions 2010
        elif (convention == '2010') and (epoch >= 1975):
            dx = epoch - 2000.0
            if (dx < 10.0):
                x[t] = 0.055974 + 1.8243e-3*dx + 1.8413e-4*dx**2 + 7.024e-6*dx**3
                y[t] = 0.346346 + 1.7896e-3*dx - 1.0729e-4*dx**2 - 0.908e-6*dx**3
            else:
                x[t] = 0.023513 + 0.0076141*dx
                y[t] = 0.358891 - 0.0006287*dx
            flag[t] = True
        # Conventional mean pole model in IERS Conventions 2015
        # must be below maximum valid date within file (e.g. 2015.2 for 2015)
        elif (convention == '2015') and (epoch >= 1975) and (epoch < end_time):
            # find epoch within mean pole table
            i = 1
            j = nrows+1
            while (j > (i+1)):
                k = (i+j)//2
                if (epoch < table[k, 0]):
                    j = k
                else:
                    i = k
            # calculate differential from point in table
            dx = epoch - table[i, 0]
            if (i == (nrows-1)):
                x[t] = table[i, 1] + dx*(table[nrows-1, 1]-table[nrows-2, 1])
                y[t] = table[i, 2] + dx*(table[nrows-1, 1]-table[nrows-2, 2])
            else:
                x[t] = table[i, 1] + dx*(table[i+1, 1]-table[i, 1])
                y[t] = table[i, 2] + dx*(table[i+1, 2]-table[i, 2])
            flag[t] = True
        # Secular pole model in IERS Conventions 2018
        elif (convention == '2018'):
            # calculate secular pole
            x[t] = 0.055 + 0.001677*(epoch - 2000.0)
            y[t] = 0.3205 + 0.00346*(epoch - 2000.0)
            flag[t] = True
    # return mean/secular pole values
    return (x, y, flag)

# PURPOSE: read daily earth orientation parameters (EOP) file from IERS
def iers_daily_EOP(input_file: str | pathlib.Path = _finals_file):
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
    # tilde-expansion of input file
    input_file = pathlib.Path(input_file).expanduser().absolute()
    # check that IERS finals file is accessible
    if not input_file.exists():
        raise FileNotFoundError(input_file)
    # read data file splitting at line breaks
    with input_file.open(mode='r', encoding='utf8') as f:
        file_contents = f.read().splitlines()
    # number of data lines
    n_lines = len(file_contents)
    dinput = {}
    dinput['MJD'] = np.zeros((n_lines))
    dinput['x'] = np.zeros((n_lines))
    dinput['y'] = np.zeros((n_lines))
    # for each line in the file
    flag = 'I'
    counter = 0
    while (flag == 'I'):
        line = file_contents[counter]
        i = 2+2+2+1; j = i+8
        dinput['MJD'][counter] = np.float64(line[i:j])
        i = j+1
        flag = line[i]
        i += 2; j = i+9
        dinput['x'][counter] = np.float64(line[i:j])
        i = j+10; j = i+9
        dinput['y'][counter] = np.float64(line[i:j])
        counter += 1
    # reduce to data values
    dinput['MJD'] = dinput['MJD'][:counter]
    dinput['x'] = dinput['x'][:counter]
    dinput['y'] = dinput['y'][:counter]
    # return the date, flag and polar motion values
    return dinput

def iers_polar_motion(
        MJD: float | np.ndarray,
        file: str | pathlib.Path = _finals_file,
        **kwargs
    ):
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
    # set default parameters
    kwargs.setdefault('k', 3)
    kwargs.setdefault('s', 0)
    # read IERS daily polar motion values
    EOP = pyTMD.eop.iers_daily_EOP(file)
    # interpolate daily polar motion values to MJD using cubic splines
    xSPL = scipy.interpolate.UnivariateSpline(EOP['MJD'], EOP['x'], **kwargs)
    ySPL = scipy.interpolate.UnivariateSpline(EOP['MJD'], EOP['y'], **kwargs)
    px = xSPL(MJD)
    py = ySPL(MJD)
    return (px, py)
