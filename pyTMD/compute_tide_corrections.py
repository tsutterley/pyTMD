#!/usr/bin/env python
u"""
compute_tide_corrections.py
Written by Tyler Sutterley (08/2020)
Calculates tidal elevations for correcting elevation or imagery data

Uses OTIS format tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/
or Global Tide Model (GOT) solutions provided by Richard Ray at GSFC

INPUTS:
    x: x-coordinates in projection EPSG
    y: y-coordinates in projection EPSG
    delta_time: seconds since EPOCH

OPTIONS:
    DIRECTORY: working data directory for tide models
    MODEL: Tide model to use in correction
    EPOCH: time period for calculating delta times
        default: J2000 (seconds since 2000-01-01T00:00:00)
    TYPE: input data type
        drift: drift buoys or satellite/airborne altimetry (time per data point)
        grid: spatial grids or images (single time per image)
    TIME: time type if need to compute leap seconds to convert to UTC
        GPS: leap seconds needed
        TAI: leap seconds needed (TAI = GPS + 19 seconds)
        UTC: no leap seconds needed
    EPSG: input coordinate system
        default: 3031 Polar Stereographic South, WGS84
    METHOD: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, cubic, nearest: scipy griddata interpolations
    FILL_VALUE: output invalid value (default NaN)

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
    time.py: utilities for calculating time operations
    utilities: download and management utilities for syncing files
    calc_astrol_longitudes.py: computes the basic astronomical mean longitudes
    calc_delta_time.py: calculates difference between universal and dynamic time
    convert_ll_xy.py: convert lat/lon points to and from projected coordinates
    infer_minor_corrections.py: return corrections for minor constituents
    load_constituent.py: loads parameters for a given tidal constituent
    load_nodal_corrections.py: load the nodal corrections for tidal constituents
    predict_tide.py: predict tides at single times using harmonic constants
    predict_tide_drift.py: predict tidal elevations using harmonic constants
    read_tide_model.py: extract tidal harmonic constants from OTIS tide models
    read_netcdf_model.py: extract tidal harmonic constants from netcdf models
    read_GOT_model.py: extract tidal harmonic constants from GSFC GOT models
    read_FES_model.py: extract tidal harmonic constants from FES tide models

UPDATE HISTORY:
    Updated 08/2020: using builtin time operations
    Updated 07/2020: added function docstrings, FES2014 and TPX09-atlas-v2
        use merged delta time files combining biannual, monthly and daily files
    Updated 03/2020: added TYPE, TIME, FILL_VALUE and METHOD options
    Written 03/2020
"""
from __future__ import print_function

import os
import pyproj
import datetime
import numpy as np
import pyTMD.time
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.infer_minor_corrections import infer_minor_corrections
from pyTMD.predict_tide import predict_tide
from pyTMD.predict_tide_drift import predict_tide_drift
from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.read_netcdf_model import extract_netcdf_constants
from pyTMD.read_GOT_model import extract_GOT_constants

#-- PURPOSE: convert value to numpy arrays if single point
def point_to_array(val):
    """
    Convert a value to a numpy array if originally a single point
    """
    return np.array([val]) if (np.ndim(val) == 0) else np.copy(val)

#-- PURPOSE: compute tides at points and times using tide model algorithms
def compute_tide_corrections(x, y, delta_time, DIRECTORY=None, MODEL=None,
    EPSG=3031, EPOCH=(2000,1,1,0,0,0), TYPE='drift', TIME='UTC',
    METHOD='spline', FILL_VALUE=np.nan):
    """
    Compute tides at points and times using tidal harmonics

    Arguments
    ---------
    x: x-coordinates in projection EPSG
    y: y-coordinates in projection EPSG
    delta_time: seconds since EPOCH

    Keyword arguments
    -----------------
    DIRECTORY: working data directory for tide models
    MODEL: Tide model to use in correction
    EPOCH: time period for calculating delta times
        default: J2000 (seconds since 2000-01-01T00:00:00)
    TYPE: input data type
        drift: drift buoys or satellite/airborne altimetry (time per data point)
        grid: spatial grids or images (single time per image)
    TIME: time type if need to compute leap seconds to convert to UTC
        GPS: leap seconds needed
        TAI: leap seconds needed (TAI = GPS + 19 seconds)
        UTC: no leap seconds needed
    EPSG: input coordinate system
        default: 3031 Polar Stereographic South, WGS84
    METHOD: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, cubic, nearest: scipy griddata interpolations
    FILL_VALUE: output invalid value (default NaN)

    Returns
    -------
    tide: tidal elevation at coordinates and time in meters
    """

    #-- select between tide models
    if (MODEL == 'CATS0201'):
        grid_file = os.path.join(DIRECTORY,'cats0201_tmd','grid_CATS')
        model_file = os.path.join(DIRECTORY,'cats0201_tmd','h0_CATS02_01')
        model_format = 'OTIS'
        model_EPSG = '4326'
        model_type = 'z'
    elif (MODEL == 'CATS2008'):
        grid_file = os.path.join(DIRECTORY,'CATS2008','grid_CATS2008')
        model_file = os.path.join(DIRECTORY,'CATS2008','hf.CATS2008.out')
        model_format = 'OTIS'
        model_EPSG = 'CATS2008'
        model_type = 'z'
    elif (MODEL == 'CATS2008_load'):
        grid_file = os.path.join(DIRECTORY,'CATS2008a_SPOTL_Load','grid_CATS2008a_opt')
        model_file = os.path.join(DIRECTORY,'CATS2008a_SPOTL_Load','h_CATS2008a_SPOTL_load')
        model_format = 'OTIS'
        model_EPSG = 'CATS2008'
        model_type = 'z'
    elif (MODEL == 'TPXO9-atlas'):
        model_directory = os.path.join(DIRECTORY,'TPXO9_atlas')
        grid_file = 'grid_tpxo9_atlas.nc.gz'
        model_files = ['h_q1_tpxo9_atlas_30.nc.gz','h_o1_tpxo9_atlas_30.nc.gz',
            'h_p1_tpxo9_atlas_30.nc.gz','h_k1_tpxo9_atlas_30.nc.gz',
            'h_n2_tpxo9_atlas_30.nc.gz','h_m2_tpxo9_atlas_30.nc.gz',
            'h_s2_tpxo9_atlas_30.nc.gz','h_k2_tpxo9_atlas_30.nc.gz',
            'h_m4_tpxo9_atlas_30.nc.gz','h_ms4_tpxo9_atlas_30.nc.gz',
            'h_mn4_tpxo9_atlas_30.nc.gz','h_2n2_tpxo9_atlas_30.nc.gz']
        model_format = 'netcdf'
        model_type = 'z'
        SCALE = 1.0/1000.0
    elif (MODEL == 'TPXO9-atlas-v2'):
        model_directory = os.path.join(tide_dir,'TPXO9_atlas_v2')
        grid_file = 'grid_tpxo9_atlas_v2.nc.gz'
        model_files = ['h_q1_tpxo9_atlas_30_v2.nc.gz','h_o1_tpxo9_atlas_30_v2.nc.gz',
            'h_p1_tpxo9_atlas_30_v2.nc.gz','h_k1_tpxo9_atlas_30_v2.nc.gz',
            'h_n2_tpxo9_atlas_30_v2.nc.gz','h_m2_tpxo9_atlas_30_v2.nc.gz',
            'h_s2_tpxo9_atlas_30_v2.nc.gz','h_k2_tpxo9_atlas_30_v2.nc.gz',
            'h_m4_tpxo9_atlas_30_v2.nc.gz','h_ms4_tpxo9_atlas_30_v2.nc.gz',
            'h_mn4_tpxo9_atlas_30_v2.nc.gz','h_2n2_tpxo9_atlas_30_v2.nc.gz']
        model_format = 'netcdf'
        TYPE = 'z'
        SCALE = 1.0/1000.0
    elif (MODEL == 'TPXO9.1'):
        grid_file = os.path.join(DIRECTORY,'TPXO9.1','DATA','grid_tpxo9')
        model_file = os.path.join(DIRECTORY,'TPXO9.1','DATA','h_tpxo9.v1')
        model_format = 'OTIS'
        model_EPSG = '4326'
        model_type = 'z'
    elif (MODEL == 'TPXO8-atlas'):
        grid_file = os.path.join(DIRECTORY,'tpxo8_atlas','grid_tpxo8atlas_30_v1')
        model_file = os.path.join(DIRECTORY,'tpxo8_atlas','hf.tpxo8_atlas_30_v1')
        model_format = 'ATLAS'
        model_EPSG = '4326'
        model_type = 'z'
    elif (MODEL == 'TPXO7.2'):
        grid_file = os.path.join(DIRECTORY,'TPXO7.2_tmd','grid_tpxo7.2')
        model_file = os.path.join(DIRECTORY,'TPXO7.2_tmd','h_tpxo7.2')
        model_format = 'OTIS'
        model_EPSG = '4326'
        model_type = 'z'
    elif (MODEL == 'TPXO7.2_load'):
        grid_file = os.path.join(DIRECTORY,'TPXO7.2_load','grid_tpxo6.2')
        model_file = os.path.join(DIRECTORY,'TPXO7.2_load','h_tpxo7.2_load')
        model_format = 'OTIS'
        model_EPSG = '4326'
        model_type = 'z'
    elif (MODEL == 'AODTM-5'):
        grid_file = os.path.join(DIRECTORY,'aodtm5_tmd','grid_Arc5km')
        model_file = os.path.join(DIRECTORY,'aodtm5_tmd','h0_Arc5km.oce')
        model_format = 'OTIS'
        model_EPSG = 'PSNorth'
        model_type = 'z'
    elif (MODEL == 'AOTIM-5'):
        grid_file = os.path.join(DIRECTORY,'aotim5_tmd','grid_Arc5km')
        model_file = os.path.join(DIRECTORY,'aotim5_tmd','h_Arc5km.oce')
        model_format = 'OTIS'
        model_EPSG = 'PSNorth'
        model_type = 'z'
    elif (MODEL == 'AOTIM-5-2018'):
        grid_file = os.path.join(DIRECTORY,'Arc5km2018','grid_Arc5km2018')
        model_file = os.path.join(DIRECTORY,'Arc5km2018','h_Arc5km2018')
        model_format = 'OTIS'
        model_EPSG = 'PSNorth'
        model_type = 'z'
    elif (MODEL == 'GOT4.7'):
        model_directory = os.path.join(DIRECTORY,'GOT4.7','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        model_format = 'GOT'
        SCALE = 1.0/100.0
    elif (MODEL == 'GOT4.7_load'):
        model_directory = os.path.join(DIRECTORY,'GOT4.7','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        model_format = 'GOT'
        SCALE = 1.0/1000.0
    elif (MODEL == 'GOT4.8'):
        model_directory = os.path.join(DIRECTORY,'got4.8','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        model_format = 'GOT'
        SCALE = 1.0/100.0
    elif (MODEL == 'GOT4.8_load'):
        model_directory = os.path.join(DIRECTORY,'got4.8','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        model_format = 'GOT'
        SCALE = 1.0/1000.0
    elif (MODEL == 'GOT4.10'):
        model_directory = os.path.join(DIRECTORY,'GOT4.10c','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        model_format = 'GOT'
        SCALE = 1.0/100.0
    elif (MODEL == 'GOT4.10_load'):
        model_directory = os.path.join(DIRECTORY,'GOT4.10c','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        model_format = 'GOT'
        SCALE = 1.0/1000.0
    elif (MODEL == 'FES2014'):
        model_directory = os.path.join(DIRECTORY,'fes2014','ocean_tide')
        model_files = ['2n2.nc.gz','eps2.nc.gz','j1.nc.gz','k1.nc.gz',
            'k2.nc.gz','l2.nc.gz','la2.nc.gz','m2.nc.gz','m3.nc.gz','m4.nc.gz',
            'm6.nc.gz','m8.nc.gz','mf.nc.gz','mks2.nc.gz','mm.nc.gz',
            'mn4.nc.gz','ms4.nc.gz','msf.nc.gz','msqm.nc.gz','mtm.nc.gz',
            'mu2.nc.gz','n2.nc.gz','n4.nc.gz','nu2.nc.gz','o1.nc.gz','p1.nc.gz',
            'q1.nc.gz','r2.nc.gz','s1.nc.gz','s2.nc.gz','s4.nc.gz','sa.nc.gz',
            'ssa.nc.gz','t2.nc.gz']
        c = ['2n2','eps2','j1','k1','k2','l2','lambda2','m2','m3','m4','m6','m8',
            'mf','mks2','mm','mn4','ms4','msf','msqm','mtm','mu2','n2','n4',
            'nu2','o1','p1','q1','r2','s1','s2','s4','sa','ssa','t2']
        model_format = 'FES'
        TYPE = 'z'
        SCALE = 1.0/100.0
    elif (MODEL == 'FES2014_load'):
        model_directory = os.path.join(DIRECTORY,'fes2014','load_tide')
        model_files = ['2n2.nc.gz','eps2.nc.gz','j1.nc.gz','k1.nc.gz',
            'k2.nc.gz','l2.nc.gz','la2.nc.gz','m2.nc.gz','m3.nc.gz','m4.nc.gz',
            'm6.nc.gz','m8.nc.gz','mf.nc.gz','mks2.nc.gz','mm.nc.gz',
            'mn4.nc.gz','ms4.nc.gz','msf.nc.gz','msqm.nc.gz','mtm.nc.gz',
            'mu2.nc.gz','n2.nc.gz','n4.nc.gz','nu2.nc.gz','o1.nc.gz','p1.nc.gz',
            'q1.nc.gz','r2.nc.gz','s1.nc.gz','s2.nc.gz','s4.nc.gz','sa.nc.gz',
            'ssa.nc.gz','t2.nc.gz']
        c = ['2n2','eps2','j1','k1','k2','l2','lambda2','m2','m3','m4','m6',
            'm8','mf','mks2','mm','mn4','ms4','msf','msqm','mtm','mu2','n2',
            'n4','nu2','o1','p1','q1','r2','s1','s2','s4','sa','ssa','t2']
        model_format = 'FES'
        model_type = 'z'
        SCALE = 1.0/100.0
    else:
        raise Exception("Unlisted tide model")

    #-- converting x,y from EPSG to latitude/longitude
    proj1 = pyproj.Proj("+init=EPSG:{0:d}".format(EPSG))
    proj2 = pyproj.Proj("+init=EPSG:{0:d}".format(4326))
    lon,lat = pyproj.transform(proj1, proj2, x.flatten(), y.flatten())

    #-- convert delta time from point to array
    delta_time = point_to_array(delta_time)
    #-- calculate leap seconds if specified
    if (TIME.upper() == 'GPS'):
        GPS_Time = pyTMD.time.convert_delta_time(delta_time, epoch1=EPOCH,
            epoch2=(1980,1,6,0,0,0), scale=1.0)
        leap_seconds = pyTMD.time.count_leap_seconds(GPS_Time)
    elif (TIME.upper() == 'TAI'):
        #-- TAI time is ahead of GPS time by 19 seconds
        GPS_Time = pyTMD.time.convert_delta_time(delta_time-19.0, epoch1=EPOCH,
            epoch2=(1980,1,6,0,0,0), scale=1.0)
        leap_seconds = pyTMD.time.count_leap_seconds(GPS_Time)
    else:
        leap_seconds = 0.0

    #-- convert time to days relative to Jan 1, 1992 (48622mjd)
    t = pyTMD.time.convert_delta_time(delta_time - leap_seconds, epoch1=EPOCH,
        epoch2=(1992,1,1,0,0,0), scale=(1.0/86400.0))

    #-- read tidal constants and interpolate to grid points
    if model_format in ('OTIS','ATLAS'):
        amp,ph,D,c = extract_tidal_constants(lon, lat, grid_file, model_file,
            model_EPSG, TYPE=model_type, METHOD=METHOD, GRID=model_format)
        deltat = np.zeros_like(t)
    elif (model_format == 'netcdf'):
        amp,ph,D,c = extract_netcdf_constants(lon, lat, model_directory,
            grid_file, model_files, TYPE=model_type, METHOD=METHOD, SCALE=SCALE)
        deltat = np.zeros_like(t)
    elif (model_format == 'GOT'):
        amp,ph = extract_GOT_constants(lon, lat, model_directory, model_files,
            METHOD=METHOD, SCALE=SCALE)
        #-- interpolate delta times from calendar dates to tide time
        delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
        deltat = calc_delta_time(delta_file, t)
    elif (model_format == 'FES'):
        amp,ph = extract_FES_constants(lon, lat, model_directory, model_files,
            TYPE=model_type, VERSION=MODEL, METHOD=METHOD, SCALE=SCALE)
        deltat = np.zeros_like(t)

    #-- calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    #-- calculate constituent oscillation
    hc = amp*np.exp(cph)

    #-- predict tidal elevations at time and infer minor corrections
    if (TYPE.lower() == 'grid'):
        ny,nx = np.shape(x); nt = len(t)
        tide = np.ma.zeros((ny,nx,nt),fill_value=FILL_VALUE)
        tide.mask = np.zeros((ny,nx,nt),dtype=np.bool)
        for i in range(nt):
            TIDE = predict_tide(t[i], hc, c,
                DELTAT=deltat[i], CORRECTIONS=model_format)
            MINOR = infer_minor_corrections(t[i], hc, c,
                DELTAT=deltat[i], CORRECTIONS=model_format)
            #-- add major and minor components and reform grid
            tide[:,:,i] = np.reshape((TIDE+MINOR), (ny,nx))
            tide.mask[:,:,i] = np.reshape((TIDE.mask | MINOR.mask), (ny,nx))
    else:
        npts = len(t)
        tide = np.ma.zeros((npts), fill_value=FILL_VALUE)
        tide.mask = np.any(hc.mask,axis=1)
        tide.data[:] = predict_tide_drift(t, hc, c,
            DELTAT=deltat, CORRECTIONS=model_format)
        minor = infer_minor_corrections(t, hc, c,
            DELTAT=deltat, CORRECTIONS=model_format)
        tide.data[:] += minor.data[:]
    #-- replace invalid values with fill value
    tide.data[tide.mask] = tide.fill_value

    #-- return the tide correction
    return tide
