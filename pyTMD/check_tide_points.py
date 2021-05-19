#!/usr/bin/env python
u"""
check_tide_points.py
Written by Tyler Sutterley (05/2021)
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
    MODEL: Tide model to use in correction
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
    convert_ll_xy.py: convert lat/lon points to and from projected coordinates
    read_tide_model.py: extract tidal harmonic constants from OTIS tide models
    read_netcdf_model.py: extract tidal harmonic constants from netcdf models
    read_GOT_model.py: extract tidal harmonic constants from GSFC GOT models
    read_FES_model.py: extract tidal harmonic constants from FES tide models
    bilinear_interp.py: bilinear interpolation of data to coordinates

UPDATE HISTORY:
    Written 05/2021
"""
from __future__ import print_function

import os
import pyproj
import numpy as np
import scipy.interpolate
import pyTMD.convert_ll_xy
import pyTMD.read_tide_model
import pyTMD.read_netcdf_model
import pyTMD.read_GOT_model
import pyTMD.read_FES_model
from pyTMD.bilinear_interp import bilinear_interp

# PURPOSE: compute tides at points and times using tide model algorithms
def check_tide_points(x,y,DIRECTORY=None,MODEL=None,EPSG=3031,METHOD='spline'):
    """
    Check if points are within a tide model domain

    Arguments
    ---------
    x: x-coordinates in projection EPSG
    y: y-coordinates in projection EPSG

    Keyword arguments
    -----------------
    DIRECTORY: working data directory for tide models
    MODEL: Tide model to use in correction
    EPSG: input coordinate system
        default: 3031 Polar Stereographic South, WGS84
    METHOD: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, nearest: scipy regular grid interpolations

    Returns
    -------
    valid: array describing if input coordinate is within model domain
    """

    # select between tide models
    if (MODEL == 'CATS0201'):
        grid_file = os.path.join(DIRECTORY,'cats0201_tmd','grid_CATS')
        model_format = 'OTIS'
        model_EPSG = '4326'
    elif (MODEL == 'CATS2008'):
        grid_file = os.path.join(DIRECTORY,'CATS2008','grid_CATS2008')
        model_format = 'OTIS'
        model_EPSG = 'CATS2008'
    elif (MODEL == 'TPXO9-atlas'):
        grid_file = os.path.join(DIRECTORY,'TPXO9_atlas','grid_tpxo9_atlas.nc.gz')
        model_format = 'netcdf'
    elif (MODEL == 'TPXO9-atlas-v2'):
        grid_file = os.path.join(DIRECTORY,'TPXO9_atlas_v2','grid_tpxo9_atlas_30_v2.nc.gz')
        model_format = 'netcdf'
    elif (MODEL == 'TPXO9-atlas-v3'):
        grid_file = os.path.join(DIRECTORY,'TPXO9_atlas_v3','grid_tpxo9_atlas_30_v3.nc.gz')
        model_format = 'netcdf'
    elif (MODEL == 'TPXO9-atlas-v4'):
        grid_file = os.path.join(DIRECTORY,'TPXO9_atlas_v4','grid_tpxo9_atlas_30_v4')
        model_format = 'OTIS'
        model_EPSG = '4326'
    elif (MODEL == 'TPXO9.1'):
        grid_file = os.path.join(DIRECTORY,'TPXO9.1','DATA','grid_tpxo9')
        model_format = 'OTIS'
        model_EPSG = '4326'
    elif (MODEL == 'TPXO8-atlas'):
        grid_file = os.path.join(DIRECTORY,'tpxo8_atlas','grid_tpxo8atlas_30_v1')
        model_format = 'ATLAS'
        model_EPSG = '4326'
    elif (MODEL == 'TPXO7.2'):
        grid_file = os.path.join(DIRECTORY,'TPXO7.2_tmd','grid_tpxo7.2')
        model_format = 'OTIS'
        model_EPSG = '4326'
    elif (MODEL == 'AODTM-5'):
        grid_file = os.path.join(DIRECTORY,'aodtm5_tmd','grid_Arc5km')
        model_format = 'OTIS'
        model_EPSG = 'PSNorth'
    elif (MODEL == 'AOTIM-5'):
        grid_file = os.path.join(DIRECTORY,'aotim5_tmd','grid_Arc5km')
        model_format = 'OTIS'
        model_EPSG = 'PSNorth'
    elif (MODEL == 'AOTIM-5-2018'):
        grid_file = os.path.join(DIRECTORY,'Arc5km2018','grid_Arc5km2018')
        model_format = 'OTIS'
        model_EPSG = 'PSNorth'
    elif (MODEL == 'GOT4.7'):
        model_directory = os.path.join(DIRECTORY,'GOT4.7','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        model_format = 'GOT'
    elif (MODEL == 'GOT4.8'):
        model_directory = os.path.join(DIRECTORY,'got4.8','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        model_format = 'GOT'
    elif (MODEL == 'GOT4.10'):
        model_directory = os.path.join(DIRECTORY,'GOT4.10c','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        model_format = 'GOT'
    elif (MODEL == 'FES2014'):
        model_directory = os.path.join(DIRECTORY,'fes2014','ocean_tide')
        model_files = ['2n2.nc.gz','eps2.nc.gz','j1.nc.gz','k1.nc.gz',
            'k2.nc.gz','l2.nc.gz','la2.nc.gz','m2.nc.gz','m3.nc.gz','m4.nc.gz',
            'm6.nc.gz','m8.nc.gz','mf.nc.gz','mks2.nc.gz','mm.nc.gz',
            'mn4.nc.gz','ms4.nc.gz','msf.nc.gz','msqm.nc.gz','mtm.nc.gz',
            'mu2.nc.gz','n2.nc.gz','n4.nc.gz','nu2.nc.gz','o1.nc.gz','p1.nc.gz',
            'q1.nc.gz','r2.nc.gz','s1.nc.gz','s2.nc.gz','s4.nc.gz','sa.nc.gz',
            'ssa.nc.gz','t2.nc.gz']
        model_format = 'FES'

    # input shape of data
    idim = np.shape(x)
    # converting x,y from EPSG to latitude/longitude
    crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(EPSG))
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon,lat = transformer.transform(np.atleast_1d(x).flatten(),
        np.atleast_1d(y).flatten())

    # read tidal constants and interpolate to grid points
    if model_format in ('OTIS','ATLAS'):
        # if reading a single OTIS solution
        xi,yi,hz,mz,iob,dt = pyTMD.read_tide_model.read_tide_grid(grid_file)
        # invert model mask
        mz = np.logical_not(mz)
        # adjust dimensions of input coordinates to be iterable
        # run wrapper function to convert coordinate systems of input lat/lon
        X,Y = pyTMD.convert_ll_xy(lon,lat,model_EPSG,'F')
    elif (model_format == 'netcdf'):
        # if reading a netCDF OTIS atlas solution
        xi,yi,hz = pyTMD.read_netcdf_model.read_netcdf_grid(grid_file,
            GZIP=True, TYPE='z')
        # copy bathymetry mask
        mz = np.copy(hz.mask)
        # copy latitude and longitude and adjust longitudes
        X,Y = np.copy([lon,lat]).astype(np.float64)
        lt0, = np.nonzero(X < 0)
        X[lt0] += 360.0
    elif (model_format == 'GOT'):
        # if reading a NASA GOT solution
        input_file = os.path.join(model_directory,model_files[0])
        hc,xi,yi,c = pyTMD.read_GOT_model.read_GOT_grid(input_file, GZIP=True)
        # copy tidal constituent mask
        mz = np.copy(hc.mask)
        # copy latitude and longitude and adjust longitudes
        X,Y = np.copy([lon,lat]).astype(np.float64)
        lt0, = np.nonzero(X < 0)
        X[lt0] += 360.0
    elif (model_format == 'FES'):
        # if reading a FES netCDF solution
        input_file = os.path.join(model_directory,model_files[0])
        hc,xi,yi = pyTMD.read_FES_model.read_netcdf_file(input_file,GZIP=True,
            TYPE='z',VERSION='FES2014')
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