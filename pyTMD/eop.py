#!/usr/bin/env python
u"""
eop.py
Written by Tyler Sutterley (11/2022)
Utilities for maintaining and calculating Earth Orientation Parameters (EOP)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 11/2022: added encoding for writing ascii files
        use f-strings for formatting verbose or ascii output
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 03/2021: replaced numpy bool/int to prevent deprecation warnings
    Written 11/2020
"""
import os
import logging
import numpy as np
import pyTMD.utilities

# PURPOSE: connects to servers and downloads mean pole files
def update_mean_pole(verbose=False, mode=0o775):
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
    # local version of file
    FILE = 'mean-pole.tab'
    LOCAL = pyTMD.utilities.get_data_path(['data',FILE])
    HASH = pyTMD.utilities.get_hash(LOCAL)

    # try downloading from Paris Observatory IERS Centers ftp servers
    HOST = ['hpiers.obspm.fr','iers','eop','eopc01',FILE]
    try:
        pyTMD.utilities.from_ftp(HOST, timeout=20, local=LOCAL,
            hash=HASH, verbose=verbose, mode=mode)
    except Exception as e:
        pass
    else:
        return

    # try downloading from Paris Observatory IERS Centers https servers
    HOST = ['http://hpiers.obspm.fr','eoppc','eop','eopc01',FILE]
    try:
        pyTMD.utilities.from_http(HOST, timeout=20, local=LOCAL,
            hash=HASH, verbose=verbose, mode=mode)
    except Exception as e:
        pass
    else:
        return

    # raise exception
    raise RuntimeError(f'Unable to download {FILE}')

# PURPOSE: read table of IERS pole coordinates and calculate Gaussian average
def calculate_mean_pole(verbose=False, mode=0o775):
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
    except Exception as e:
        return

    # read contents from input file object
    file_contents = remote_buffer.read().decode('utf8').splitlines()
    header = file_contents[0][1:].split()
    nlines = len(file_contents) - 1
    data = {h:np.zeros((nlines)) for h in header}
    # extract data for all lines
    for i,line in enumerate(file_contents[1:]):
        line_contents = line.split()
        for h,l in zip(header,line_contents):
            data[h][i] = np.float64(l)
    # output mean pole coordinates
    xm = np.zeros((nlines))
    ym = np.zeros((nlines))
    # output file with mean pole coordinates
    LOCAL = pyTMD.utilities.get_data_path(['data','mean-pole.tab'])
    fid = open(LOCAL, mode='w', encoding='utf8')
    logging.info(LOCAL)
    for i,T in enumerate(data['an']):
        # mean pole is Gaussian Weight of all dates with a = 3.40 years.
        Wi = np.exp(-0.5*((data['an']-T)/3.4)**2)
        xm[i] = np.sum(Wi*data['x(")'])/np.sum(Wi)
        ym[i] = np.sum(Wi*data['y(")'])/np.sum(Wi)
        print('{0:6.2f} {1:11.7f} {2:11.7f}'.format(T,xm[i],ym[i]),file=fid)
    # close the output file
    fid.close()
    # change the permissions mode of the output mean pole file
    os.chmod(LOCAL, mode)

# PURPOSE: connects to servers and downloads IERS pole coordinates files
def pull_pole_coordinates(FILE, verbose=False):
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
    HOST = ['ftp.iers.org','products','eop','long-term','c01',FILE]
    try:
        buffer = pyTMD.utilities.from_ftp(HOST, verbose=verbose, timeout=20)
    except Exception as e:
        pass
    else:
        return buffer

    # try downloading from Paris Observatory IERS Centers ftp servers
    HOST = ['hpiers.obspm.fr','iers','eop','eopc01',FILE]
    try:
        buffer = pyTMD.utilities.from_ftp(HOST, verbose=verbose, timeout=20)
    except Exception as e:
        pass
    else:
        return buffer

    # try downloading from Paris Observatory IERS Centers https servers
    HOST = ['http://hpiers.obspm.fr','eoppc','eop','eopc01',FILE]
    try:
        buffer = pyTMD.utilities.from_http(HOST, verbose=verbose, timeout=20)
    except Exception as e:
        pass
    else:
        return buffer

    # raise exception
    raise RuntimeError(f'Unable to download {FILE}')

# PURPOSE: connects to servers and downloads finals files
def update_finals_file(username=None, password=None, verbose=False, mode=0o775):
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
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file
    """
    # local version of file
    LOCAL = pyTMD.utilities.get_data_path(['data','finals.all'])
    HASH = pyTMD.utilities.get_hash(LOCAL)

    # try downloading from US Naval Oceanography Portal
    HOST = ['http://maia.usno.navy.mil','ser7','finals.all']
    try:
        pyTMD.utilities.from_http(HOST, timeout=5, local=LOCAL,
            hash=HASH, verbose=verbose, mode=mode)
    except Exception as e:
        pass
    else:
        return

    # try downloading from NASA Crustal Dynamics Data Information System
    # note: anonymous ftp access will be discontinued on 2020-10-31
    # will require using the following https Earthdata server after that date
    server = []
    server.append(['cddis.nasa.gov','pub','products','iers','finals.all'])
    server.append(['cddis.gsfc.nasa.gov','products','iers','finals.all'])
    for HOST in server:
        try:
            pyTMD.utilities.from_ftp(HOST, timeout=20, local=LOCAL,
                hash=HASH, verbose=verbose, mode=mode)
        except Exception as e:
            pass
        else:
            return

    # try downloading from NASA Crustal Dynamics Data Information System
    # using NASA Earthdata credentials stored in netrc file
    HOST = ['https://cddis.nasa.gov','archive','products','iers','finals.all']
    try:
        pyTMD.utilities.from_cddis(HOST, username=username, password=password,
            timeout=20, local=LOCAL, hash=HASH, verbose=verbose, mode=mode)
    except Exception as e:
        pass
    else:
        return

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
    # set default keyword arguments
    kwargs.setdefault('fill_value', np.nan)
    # verify IERS model version
    assert version in ('2003','2010','2015'), "Incorrect IERS model version"
    # read mean pole file
    table = np.loadtxt(os.path.expanduser(input_file))
    # reduce to 1971 to end date
    ii, = np.nonzero(table[:,0] >= 1971)
    table = np.copy(table[ii,:])
    # reduce to yearly values
    jj, = np.nonzero((table[:,0] % 1) == 0.0)
    table = np.copy(table[jj,:])
    end_time = table[-1,0] + 0.2
    # final shape of the table
    nrows, ncols = np.shape(table)
    # allocate for output arrays
    x = np.full_like(input_epoch, kwargs['fill_value'])
    y = np.full_like(input_epoch, kwargs['fill_value'])
    flag = np.zeros_like(input_epoch, dtype=bool)
    for t,epoch in enumerate(input_epoch):
        # Conventional mean pole model in IERS Conventions 2003
        if (version == '2003') and (epoch >= 1975) and (epoch < 2004):
            x[t] = 0.054 + 0.00083*(epoch-2000.0)
            y[t] = 0.357 + 0.00395*(epoch-2000.0)
            flag[t] = True
        # Conventional mean pole model in IERS Conventions 2010
        elif (version == '2010') and (epoch >= 1975) and (epoch < 2011):
            dx = epoch-2000.0
            if (dx < 10.0):
                x[t] = 0.055974 + 1.8243e-3*dx + 1.8413e-4*dx**2 + 7.024e-6*dx**3
                y[t] = 0.346346 + 1.7896e-3*dx + 1.0729e-4*dx**2 + 0.908e-6*dx**3
            else:
                x[t] = 0.023513 + 0.0076141*dx
                y[t] = 0.358891 - 0.0006287*dx
            flag[t] = True
        # Conventional mean pole model in IERS Conventions 2015
        # must be below maximum valid date within file (e.g. 2015.2 for 2015)
        elif (version == '2015') and (epoch >= 1975) and (epoch < end_time):
            # find epoch within mean pole table
            i = 1
            j = nrows+1
            while (j > (i+1)):
                k = (i+j)//2
                if (epoch < table[k,0]):
                    j = k
                else:
                    i = k
            # calculate differential from point in table
            dx = epoch - table[i,0]
            if (i == (nrows-1)):
                x[t] = table[i,1] + dx*(table[nrows-1,1]-table[nrows-2,1])
                y[t] = table[i,2] + dx*(table[nrows-1,1]-table[nrows-2,2])
            else:
                x[t] = table[i,1] + dx*(table[i+1,1]-table[i,1])
                y[t] = table[i,2] + dx*(table[i+1,2]-table[i,2])
            flag[t] = True
    # return mean pole values
    return (x, y, flag)

# PURPOSE: read daily earth orientation parameters (EOP) file from IERS
def iers_daily_EOP(input_file):
    """
    Read daily earth orientation parameters (EOP) file from IERS

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
    # tilde-expansion of input file
    input_file = os.path.expanduser(input_file)
    # check that IERS finals file is accessible
    if not os.access(input_file, os.F_OK):
        raise FileNotFoundError(input_file)
    # read data file splitting at line breaks
    with open(input_file, mode='r', encoding='utf8') as f:
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
