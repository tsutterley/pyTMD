#!/usr/bin/env python
u"""
check_tide_points.py
Written by Tyler Sutterley (04/2023)
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
    convert_crs.py: convert points to and from Coordinates Reference Systems
    io/model.py: retrieves tide model parameters for named tide models
    io/OTIS.py: extract tidal harmonic constants from OTIS tide models
    io/ATLAS.py: extract tidal harmonic constants from ATLAS netcdf models
    io/GOT.py: extract tidal harmonic constants from GSFC GOT models
    io/FES.py: extract tidal harmonic constants from FES tide models
    interpolate.py: interpolation routines for spatial data

UPDATE HISTORY:
    Updated 04/2023: using pathlib to define and expand paths
    Updated 03/2023: add basic variable typing to function inputs
    Updated 12/2022: refactored tide read programs under io
        refactored bilinear interpolation routine
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
from __future__ import print_function, annotations

import logging
import pathlib
import numpy as np
import scipy.interpolate

import pyTMD.io
import pyTMD.io.model
import pyTMD.convert_crs
import pyTMD.interpolate

# attempt imports
try:
    import pyproj
except (ImportError, ModuleNotFoundError) as exc:
    logging.debug("pyproj not available")

# PURPOSE: compute tides at points and times using tide model algorithms
def check_tide_points(x: np.ndarray, y: np.ndarray,
        DIRECTORY: str | pathlib.Path | None = None,
        MODEL: str | None = None,
        ATLAS_FORMAT: str = 'netcdf',
        GZIP: bool = False,
        DEFINITION_FILE: str | pathlib.Path | None = None,
        EPSG: str | int = 3031,
        METHOD: str = 'spline'
    ):
    """
    Check if points are within a tide model domain

    Parameters
    ----------
    x: np.ndarray
        x-coordinates in projection EPSG
    y: np.ndarray
        y-coordinates in projection EPSG
    DIRECTORY: str or NoneType, default None
        working data directory for tide models
    MODEL: str or NoneType, default None
        Tide model to use
    ATLAS_FORMAT: str, default 'netcdf'
        ATLAS tide model format

            - ``'OTIS'``
            - ``'netcdf'``
    GZIP: bool, default False
        Tide model files are gzip compressed
    DEFINITION_FILE: str or NoneType, default None
        Tide model definition file for use
    EPSG: str or int, default: 3031 (Polar Stereographic South, WGS84)
        Input coordinate system
    METHOD: str, default 'spline'
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
    if DIRECTORY is not None:
        DIRECTORY = pathlib.Path(DIRECTORY).expanduser()
        if not DIRECTORY.exists():
            raise FileNotFoundError("Invalid tide directory")

    # get parameters for tide model
    if DEFINITION_FILE is not None:
        model = pyTMD.io.model(DIRECTORY).from_file(
            pathlib.Path(DEFINITION_FILE).expanduser())
    else:
        model = pyTMD.io.model(DIRECTORY, format=ATLAS_FORMAT,
            compressed=GZIP).elevation(MODEL)

    # input shape of data
    idim = np.shape(x)
    # converting x,y from input coordinate reference system
    try:
        # EPSG projection code string or int
        crs1 = pyproj.CRS.from_epsg(int(EPSG))
    except (ValueError,pyproj.exceptions.CRSError):
        # Projection SRS string
        crs1 = pyproj.CRS.from_string(EPSG)
    # convert to latitude and longitude
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon, lat = transformer.transform(
        np.atleast_1d(x).flatten(), np.atleast_1d(y).flatten()
    )

    # read tidal constants and interpolate to grid points
    if model.format in ('OTIS','ATLAS','ESR'):
        # if reading a single OTIS solution
        xi, yi, hz, mz, iob, dt = pyTMD.io.OTIS.read_otis_grid(
            pathlib.Path(model.grid_file).expanduser())
        # invert model mask
        mz = np.logical_not(mz)
        # adjust dimensions of input coordinates to be iterable
        # run wrapper function to convert coordinate systems of input lat/lon
        X, Y = pyTMD.convert_crs(lon, lat, model.projection, 'F')
    elif (model.format == 'netcdf'):
        # if reading a netCDF OTIS atlas solution
        xi, yi, hz = pyTMD.io.ATLAS.read_netcdf_grid(
            pathlib.Path(model.grid_file).expanduser(),
            compressed=model.compressed, type=model.type)
        # copy bathymetry mask
        mz = np.copy(hz.mask)
        # copy latitude and longitude and adjust longitudes
        X,Y = np.copy([lon,lat]).astype(np.float64)
        lt0, = np.nonzero(X < 0)
        X[lt0] += 360.0
    elif (model.format == 'GOT'):
        # if reading a NASA GOT solution
        hc, xi, yi, c = pyTMD.io.GOT.read_ascii_file(
            pathlib.Path(model.model_file[0]).expanduser(),
            compressed=model.compressed)
        # copy tidal constituent mask
        mz = np.copy(hc.mask)
        # copy latitude and longitude and adjust longitudes
        X, Y = np.copy([lon,lat]).astype(np.float64)
        lt0, = np.nonzero(X < 0)
        X[lt0] += 360.0
    elif (model.format == 'FES'):
        # if reading a FES netCDF solution
        hc, xi, yi = pyTMD.io.FES.read_netcdf_file(
            pathlib.Path(model.model_file[0]).expanduser(),
            compressed=model.compressed, type=model.type,
            version=model.version)
        # copy tidal constituent mask
        mz = np.copy(hc.mask)
        # copy latitude and longitude and adjust longitudes
        X, Y = np.copy([lon,lat]).astype(np.float64)
        lt0, = np.nonzero(X < 0)
        X[lt0] += 360.0

    # interpolate masks
    if (METHOD == 'bilinear'):
        # replace invalid values with nan
        mz1 = pyTMD.interpolate.bilinear(xi, yi, mz, X, Y)
        mask = np.floor(mz1).astype(mz.dtype)
    elif (METHOD == 'spline'):
        f1 = scipy.interpolate.RectBivariateSpline(xi, yi, mz.T,
            kx=1, ky=1)
        mask = np.floor(f1.ev(X, Y)).astype(mz.dtype)
    else:
        # use scipy regular grid to interpolate values
        r1 = scipy.interpolate.RegularGridInterpolator((yi, xi), mz,
            method=METHOD, bounds_error=False, fill_value=1)
        mask = np.floor(r1.__call__(np.c_[y, x])).astype(mz.dtype)

    # reshape to original dimensions
    valid = np.logical_not(mask).reshape(idim).astype(mz.dtype)
    # replace points outside model domain with invalid
    valid &= (X >= xi.min()) & (X <= xi.max())
    valid &= (Y >= yi.min()) & (Y <= yi.max())
    # return the valid mask
    return valid
