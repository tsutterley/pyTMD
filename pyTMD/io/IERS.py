#!/usr/bin/env python
u"""
IERS.py
Written by Tyler Sutterley (08/2024)

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
    Updated 08/2024: convert outputs to be in -180:180 longitude convention
        added function to interpolate ocean pole tide values to coordinates
        renamed from ocean_pole_tide to IERS
    Updated 06/2024: use np.clongdouble instead of np.longcomplex
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
import warnings
import numpy as np
import scipy.interpolate
from pyTMD.utilities import get_data_path

# ocean pole tide file from Desai (2002) and IERS conventions
_ocean_pole_tide_file = get_data_path(['data','opoleloadcoefcmcor.txt.gz'])

# PURPOSE: extract ocean pole tide values from Desai (2002) at coordinates
def extract_coefficients(
        lon: np.ndarray, lat: np.ndarray,
        **kwargs
    ):
    """
    Reads ocean pole tide file from [1]_ [2]_ and spatially interpolates
    to input coordinates

    Parameters
    ----------
    lon: np.ndarray
        longitude to interpolate
    lat: np.ndarray
        latitude to interpolate
    model_file: str
        IERS map of ocean pole tide coefficients
    method: str, default 'spline'
        Interpolation method

            - ``'spline'``: scipy bivariate spline interpolation
            - ``'linear'``, ``'nearest'``: scipy regular grid interpolations

    Returns
    -------
    ur: np.ndarray
        radial ocean pole tide coefficients
    un: np.ndarray
        north ocean pole tide coefficients
    ue: np.ndarray
        east ocean pole tide coefficients

    References
    ----------
    .. [1] S. Desai, "Observing the pole tide with satellite altimetry", *Journal of
        Geophysical Research: Oceans*, 107(C11), (2002).
        `doi: 10.1029/2001JC001224 <https://doi.org/10.1029/2001JC001224>`_
    .. [2] S. Desai, J. Wahr and B. Beckley "Revisiting the pole tide for and from
        satellite altimetry", *Journal of Geodesy*, 89(12), p1233-1243, (2015).
        `doi: 10.1007/s00190-015-0848-7 <https://doi.org/10.1007/s00190-015-0848-7>`_
    """
    # default keyword arguments
    kwargs.setdefault('model_file', _ocean_pole_tide_file)
    kwargs.setdefault('method', 'spline')
    # number of points
    npts = len(np.atleast_1d(lat))
    # read ocean pole tide map from Desai (2002)
    Umap = {}
    Umap['R'], Umap['N'], Umap['E'], ilon, ilat = read_binary_file(**kwargs)
    # interpolate ocean pole tide map from Desai (2002)
    Uint = {}
    if (kwargs['method'] == 'spline'):
        # use scipy bivariate splines to interpolate to output points
        for key,val in Umap.items():
            Uint[key] = np.zeros((npts), dtype=np.clongdouble)
            f1 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
                val[:,::-1].real, kx=1, ky=1)
            f2 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
                val[:,::-1].imag, kx=1, ky=1)
            Uint[key].real = f1.ev(lon, lat)
            Uint[key].imag = f2.ev(lon, lat)
    else:
        # use scipy regular grid to interpolate values for a given method
        for key,val in Umap.items():
            Uint[key] = np.zeros((npts), dtype=np.clongdouble)
            r1 = scipy.interpolate.RegularGridInterpolator((ilon,ilat[::-1]),
                val[:,::-1], bounds_error=False, method=kwargs['method'])
            Uint[key][:] = r1.__call__(np.c_[lon, lat])
    # return the interpolated values
    return (Uint['R'], Uint['N'], Uint['E'])

# PURPOSE: read real and imaginary ocean pole tide coefficients
def read_binary_file(**kwargs):
    """
    Read real and imaginary ocean pole tide coefficients from [1]_ [2]_

    Parameters
    ----------
    model_file: str or pathlib.Path
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
    # default keyword arguments
    kwargs.setdefault('model_file', _ocean_pole_tide_file)
    # convert input file to tilde-expanded pathlib object
    input_file = pathlib.Path(kwargs['model_file']).expanduser()
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
    U = {}
    U['R'] = np.zeros((nlon,nlat), dtype=np.clongdouble)
    U['N'] = np.zeros((nlon,nlat), dtype=np.clongdouble)
    U['E'] = np.zeros((nlon,nlat), dtype=np.clongdouble)
    # read lines of file and add to output variables
    for i,line in enumerate(file_contents[count:]):
        ln,lt,urr,uri,unr,uni,uer,uei = np.array(line.split(), dtype='f8')
        ilon = int(ln/dlon)
        ilat = int((90.0-lt)/dlat)
        U['R'][ilon,ilat] = urr + 1j*uri
        U['N'][ilon,ilat] = unr + 1j*uni
        U['E'][ilon,ilat] = uer + 1j*uei
    # shift ocean pole tide grid to -180:180
    longitudes = np.copy(glon)
    for key, val in U.items():
        U[key], glon = _shift(val, longitudes, lon0=180.0,
            cyclic=360.0, direction='west')
    # extend matrix for bilinear interpolation
    glon = _extend_array(glon, dlon)
    # pad ends of matrix for interpolation
    for key, val in U.items():
        U[key] = _extend_matrix(val)
    # return values
    return (U['R'], U['N'], U['E'], glon, glat)

# PURPOSE: deprecated function to read ocean pole tide coefficients
def ocean_pole_tide(input_file: str | pathlib.Path = _ocean_pole_tide_file):
    warnings.warn("Deprecated. Please use pyTMD.io.IERS instead",
        DeprecationWarning)
    # pass the input file to the read_binary_file function
    return read_binary_file(model_file=input_file)

# PURPOSE: Extend a longitude array
def _extend_array(input_array: np.ndarray, step_size: float):
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
def _extend_matrix(input_matrix: np.ndarray):
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

# PURPOSE: shift a grid east or west
def _shift(
        input_matrix: np.ndarray,
        ilon: np.ndarray,
        lon0: int | float = 180,
        cyclic: int | float = 360,
        direction: str = 'west'
    ):
    """
    Shift global grid east or west to a new base longitude

    Parameters
    ----------
    input_matrix: np.ndarray
        input matrix to shift
    ilon: np.ndarray
        longitude of tidal model
    lon0: int or float, default 180
        Starting longitude for shifted grid
    cyclic: int or float, default 360
        width of periodic domain
    direction: str, default 'west'
        Direction to shift grid

            - ``'west'``
            - ``'east'``

    Returns
    -------
    temp: np.ndarray
        shifted matrix
    lon: np.ndarray
        shifted longitude
    """
    # find the starting index if cyclic
    offset = 0 if (np.fabs(ilon[-1]-ilon[0]-cyclic) > 1e-4) else 1
    i0 = np.argmin(np.fabs(ilon - lon0))
    # shift longitudinal values
    lon = np.zeros(ilon.shape, ilon.dtype)
    lon[0:-i0] = ilon[i0:]
    lon[-i0:] = ilon[offset: i0+offset]
    # add or remove the cyclic
    if (direction == 'east'):
        lon[-i0:] += cyclic
    elif (direction == 'west'):
        lon[0:-i0] -= cyclic
    # allocate for shifted data
    if np.ma.isMA(input_matrix):
        temp = np.ma.zeros(input_matrix.shape,input_matrix.dtype)
    else:
        temp = np.zeros(input_matrix.shape, input_matrix.dtype)
    # shift data values
    temp[:-i0,:] = input_matrix[i0:, :]
    temp[-i0:,:] = input_matrix[offset: i0+offset, :]
    # return the shifted values
    return (temp, lon)
