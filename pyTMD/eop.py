#!/usr/bin/env python
u"""
eop.py
Written by Tyler Sutterley (11/2022)
Utilities for maintaining Earth Orientation Parameter (EOP) files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 11/2022: added encoding for writing ascii files
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
    raise RuntimeError('Unable to download {0}'.format(FILE))

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
    raise RuntimeError('Unable to download {0}'.format(FILE))

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
