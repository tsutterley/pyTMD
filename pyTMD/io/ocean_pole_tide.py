#!/usr/bin/env python
u"""
ocean_pole_tide.py
Written by Tyler Sutterley (05/2023)

Reads ocean pole load tide coefficients provided by IERS
http://maia.usno.navy.mil/conventions/2010/2010_official/chapter7/tn36_c7.pdf
http://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/icc7.pdf

IERS 0.5x0.5 map of ocean pole tide coefficients:
ftp://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/additional_info/
    opoleloadcoefcmcor.txt.gz

OUTPUTS:
    ur: radial ocean pole tide coefficients
    un: north ocean pole tide coefficients
    ue: east ocean pole tide coefficients
    glon: ocean grid longitude
    glat: ocean grid latitude

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

REFERENCES:
    S. Desai, "Observing the pole tide with satellite altimetry", Journal of
        Geophysical Research: Oceans, 107(C11), 2002. doi: 10.1029/2001JC001224
    S. Desai, J. Wahr and B. Beckley "Revisiting the pole tide for and from
        satellite altimetry", Journal of Geodesy, 89(12), p1233-1243, 2015.
        doi: 10.1007/s00190-015-0848-7

UPDATE HISTORY:
    Updated 05/2023: add default for ocean pole tide file
    Updated 04/2023: using pathlib to define and expand paths
    Updated 03/2023: add basic variable typing to function inputs
    Updated 12/2022: refactor ocean pole tide read programs under io
    Updated 04/2022: updated docstrings to numpy documentation format
        use longcomplex data format to be windows compliant
    Updated 07/2021: added check that ocean pole tide file is accessible
    Updated 03/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 08/2020: output north load and east load deformation components
    Updated 07/2020: added function docstrings
    Updated 12/2018: Compatibility updates for Python3
    Written 09/2017
"""
from __future__ import annotations

import re
import gzip
import pathlib
import numpy as np
from pyTMD.utilities import get_data_path

# ocean pole tide file from Desai (2002) and IERS conventions
_ocean_pole_tide_file = get_data_path(['data','opoleloadcoefcmcor.txt.gz'])

# PURPOSE: read real and imaginary ocean pole tide coefficients
def ocean_pole_tide(input_file: str | pathlib.Path = _ocean_pole_tide_file):
    """
    Read real and imaginary ocean pole tide coefficients

    Parameters
    ----------
    input_file: str
        IERS map of ocean pole tide coefficients

    Returns
    -------
    ur: np.ndarray
        radial ocean pole tide coefficients
    un: np.ndarray
        north ocean pole tide coefficients
    ue: np.ndarray
        east ocean pole tide coefficients
    glon: np.ndarray
        ocean grid longitude
    glat: np.ndarray
        ocean grid latitude

    References
    ----------
    .. [1] S. Desai, "Observing the pole tide with satellite altimetry", *Journal of
        Geophysical Research: Oceans*, 107(C11), (2002).
        `doi: 10.1029/2001JC001224 <https://doi.org/10.1029/2001JC001224>`_
    .. [2] S. Desai, J. Wahr and B. Beckley "Revisiting the pole tide for and from
        satellite altimetry", *Journal of Geodesy*, 89(12), p1233-1243, (2015).
        `doi: 10.1007/s00190-015-0848-7 <https://doi.org/10.1007/s00190-015-0848-7>`_
    """
    # convert input file to tilde-expanded pathlib object
    input_file = pathlib.Path(input_file).expanduser()
    # check that ocean pole tide file is accessible
    if not input_file.exists():
        raise FileNotFoundError(str(input_file))

    # read GZIP ocean pole tide file
    with gzip.open(input_file, 'rb') as f:
        file_contents = f.read().splitlines()

    # counts the number of lines in the header
    count = 0
    # Reading over header text
    HEADER = True
    while HEADER:
        # file line at count
        line = file_contents[count]
        # find --------- within line to set HEADER flag to False when found
        HEADER = not bool(re.match(rb'---------',line))
        # add 1 to counter
        count += 1

    # grid parameters and dimensions
    dlon,dlat = (0.50,0.50)
    glon = np.arange(dlon/2.0,360+dlon/2.0,dlon)
    glat = np.arange(90.0-dlat/2.0,-90.0-dlat/2.0,-dlat)
    nlon = len(glon)
    nlat = len(glat)
    # allocate for output grid maps
    ur = np.zeros((nlon,nlat),dtype=np.longcomplex)
    un = np.zeros((nlon,nlat),dtype=np.longcomplex)
    ue = np.zeros((nlon,nlat),dtype=np.longcomplex)
    # read lines of file and add to output variables
    for i,line in enumerate(file_contents[count:]):
        ln,lt,urr,uri,unr,uni,uer,uei = np.array(line.split(), dtype='f8')
        ilon = int(ln/dlon)
        ilat = int((90.0-lt)/dlat)
        ur[ilon,ilat] = urr + 1j*uri
        un[ilon,ilat] = unr + 1j*uni
        ue[ilon,ilat] = uer + 1j*uei

    # extend matrix for bilinear interpolation
    glon = extend_array(glon,dlon)
    ur = extend_matrix(ur)
    un = extend_matrix(un)
    ue = extend_matrix(ue)
    # return values
    return (ur, un, ue, glon, glat)

# PURPOSE: Extend a longitude array
def extend_array(input_array: np.ndarray, step_size: float):
    """
    Extends a longitude array

    Parameters
    ----------
    input_array: np.ndarray
        array to extend
    step_size: float
        step size between elements of array

    Returns
    -------
    temp: np.ndarray
        extended array
    """
    n = len(input_array)
    temp = np.zeros((n+2), dtype=input_array.dtype)
    # extended array [x-1,x0,...,xN,xN+1]
    temp[0] = input_array[0] - step_size
    temp[1:-1] = input_array[:]
    temp[-1] = input_array[-1] + step_size
    return temp

# PURPOSE: Extend a global matrix
def extend_matrix(input_matrix: np.ndarray):
    """
    Extends a global matrix

    Parameters
    ----------
    input_matrix: np.ndarray
        matrix to extend

    Returns
    -------
    temp: np.ndarray
        extended matrix
    """
    nx,ny = np.shape(input_matrix)
    temp = np.zeros((nx+2,ny), dtype=input_matrix.dtype)
    temp[0,:] = input_matrix[-1,:]
    temp[1:-1,:] = input_matrix[:,:]
    temp[-1,:] = input_matrix[0,:]
    return temp
