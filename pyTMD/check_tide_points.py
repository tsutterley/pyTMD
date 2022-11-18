#!/usr/bin/env python
u"""
check_tide_points.py
Written by Tyler Sutterley (11/2022)
Check if points are within a tide model domain

OTIS format tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/
Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
or Finite Element Solution (FES) models provided by AVISO

INPUTS:
    x: x-coordinates in projection EPSG
    y: y-coordinates in projection EPSG

OPTIONS:
    DIRECTORY: working data directory for tide models
    MODEL: Tide model to use
    ATLAS_FORMAT: ATLAS tide model format (OTIS, netcdf)
    GZIP: Tide model files are gzip compressed
    DEFINITION_FILE: Tide model definition file for use
    EPSG: input coordinate system
        default: 3031 Polar Stereographic South, WGS84
    METHOD: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, nearest: scipy regular grid interpolations

OUTPUTS:
    valid: array describing if input coordinate is within model domain

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    model.py: retrieves tide model parameters for named tide models
    convert_ll_xy.py: convert lat/lon points to and from projected coordinates
    read_tide_model.py: extract tidal harmonic constants from OTIS tide models
    read_netcdf_model.py: extract tidal harmonic constants from netcdf models
    read_GOT_model.py: extract tidal harmonic constants from GSFC GOT models
    read_FES_model.py: extract tidal harmonic constants from FES tide models
    bilinear_interp.py: bilinear interpolation of data to coordinates

UPDATE HISTORY:
    Updated 11/2022: place some imports within try/except statements
        use f-strings for formatting verbose or ascii output
    Updated 05/2022: added ESR netCDF4 formats to list of model types
        updated keyword arguments to read tide model programs
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 09/2021: refactor to use model class for files and attributes
    Updated 07/2021: added check that tide model directory is accessible
    Updated 06/2021: add try/except for input projection strings
    Written 05/2021
"""
from __future__ import print_function

import os
import warnings
import numpy as np
import scipy.interpolate
import pyTMD.model
import pyTMD.convert_ll_xy
import pyTMD.read_tide_model
import pyTMD.read_netcdf_model
import pyTMD.read_GOT_model
import pyTMD.read_FES_model
from pyTMD.bilinear_interp import bilinear_interp

# attempt imports
try:
    import pyproj
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("pyproj not available")
    warnings.warn("Some functions will throw an exception if called")
# ignore warnings
warnings.filterwarnings("ignore")

# PURPOSE: compute tides at points and times using tide model algorithms
def check_tide_points(x, y, DIRECTORY=None, MODEL=None,
    ATLAS_FORMAT='netcdf', GZIP=False, DEFINITION_FILE=None,
    EPSG=3031, METHOD='spline'):
    """
    Check if points are within a tide model domain

    Parameters
    ----------
    x: float
        x-coordinates in projection EPSG
    y: float
        y-coordinates in projection EPSG
    DIRECTORY: str or NoneType, default None
        working data directory for tide models
    MODEL:  str or NoneType, default None
        Tide model to use
    ATLAS_FORMAT: str, default 'netcdf'
        ATLAS tide model format

            - ``'OTIS'``
            - ``'netcdf'``
    GZIP: bool, default False
        Tide model files are gzip compressed
    DEFINITION_FILE: str or NoneType, default None
        Tide model definition file for use
    EPSG: int, default: 3031 (Polar Stereographic South, WGS84)
        Input coordinate system
    METHOD: str
        interpolation method

            - ```bilinear```: quick bilinear interpolation
            - ```spline```: scipy bivariate spline interpolation
            - ```linear```, ```nearest```: scipy regular grid interpolations

    Returns
    -------
    valid: bool
        array describing if input coordinate is within model domain
    """

    # check that tide directory is accessible
    try:
        os.access(DIRECTORY, os.F_OK)
    except:
        raise FileNotFoundError("Invalid tide directory")

    # get parameters for tide model
    if DEFINITION_FILE is not None:
        model = pyTMD.model(DIRECTORY).from_file(DEFINITION_FILE)
    else:
        model = pyTMD.model(DIRECTORY, format=ATLAS_FORMAT,
            compressed=GZIP).elevation(MODEL)

    # input shape of data
    idim = np.shape(x)
    # converting x,y from EPSG to latitude/longitude
    try:
        # EPSG projection code string or int
        crs1 = pyproj.CRS.from_epsg(int(EPSG))
    except (ValueError,pyproj.exceptions.CRSError):
        # Projection SRS string
        crs1 = pyproj.CRS.from_string(EPSG)
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon,lat = transformer.transform(np.atleast_1d(x).flatten(),
        np.atleast_1d(y).flatten())

    # read tidal constants and interpolate to grid points
    if model.format in ('OTIS','ATLAS','ESR'):
        # if reading a single OTIS solution
        xi,yi,hz,mz,iob,dt = pyTMD.read_tide_model.read_tide_grid(model.grid_file)
        # invert model mask
        mz = np.logical_not(mz)
        # adjust dimensions of input coordinates to be iterable
        # run wrapper function to convert coordinate systems of input lat/lon
        X,Y = pyTMD.convert_ll_xy(lon,lat,model.projection,'F')
    elif (model.format == 'netcdf'):
        # if reading a netCDF OTIS atlas solution
        xi,yi,hz = pyTMD.read_netcdf_model.read_netcdf_grid(model.grid_file,
            compressed=model.compressed, type=model.type)
        # copy bathymetry mask
        mz = np.copy(hz.mask)
        # copy latitude and longitude and adjust longitudes
        X,Y = np.copy([lon,lat]).astype(np.float64)
        lt0, = np.nonzero(X < 0)
        X[lt0] += 360.0
    elif (model.format == 'GOT'):
        # if reading a NASA GOT solution
        hc,xi,yi,c = pyTMD.read_GOT_model.read_GOT_grid(model.model_file[0],
            compressed=model.compressed)
        # copy tidal constituent mask
        mz = np.copy(hc.mask)
        # copy latitude and longitude and adjust longitudes
        X,Y = np.copy([lon,lat]).astype(np.float64)
        lt0, = np.nonzero(X < 0)
        X[lt0] += 360.0
    elif (model.format == 'FES'):
        # if reading a FES netCDF solution
        hc,xi,yi = pyTMD.read_FES_model.read_netcdf_file(model.model_file[0],
            compressed=model.compressed, type=model.type,
            version=model.version)
        # copy tidal constituent mask
        mz = np.copy(hc.mask)
        # copy latitude and longitude and adjust longitudes
        X,Y = np.copy([lon,lat]).astype(np.float64)
        lt0, = np.nonzero(X < 0)
        X[lt0] += 360.0

    # interpolate masks
    if (METHOD == 'bilinear'):
        # replace invalid values with nan
        mz1 = bilinear_interp(xi,yi,mz,X,Y)
        mask = np.floor(mz1).astype(mz.dtype)
    elif (METHOD == 'spline'):
        f1=scipy.interpolate.RectBivariateSpline(xi,yi,mz.T,kx=1,ky=1)
        mask = np.floor(f1.ev(X,Y)).astype(mz.dtype)
    else:
        # use scipy regular grid to interpolate values
        r1 = scipy.interpolate.RegularGridInterpolator((yi,xi),mz,
            method=METHOD,bounds_error=False,fill_value=1)
        mask = np.floor(r1.__call__(np.c_[y,x])).astype(mz.dtype)

    # reshape to original dimensions
    valid = np.logical_not(mask).reshape(idim).astype(mz.dtype)
    # replace points outside model domain with invalid
    valid &= (X >= xi.min()) & (X <= xi.max())
    valid &= (Y >= yi.min()) & (Y <= yi.max())
    # return the valid mask
    return valid