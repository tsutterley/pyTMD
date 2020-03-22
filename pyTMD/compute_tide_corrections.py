#!/usr/bin/env python
u"""
compute_tide_corrections.py
Written by Tyler Sutterley (03/2020)
Calculates tidal elevations for correcting elevation data

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
    LEAPS: need to compute leap seconds to convert to UTC (True/False)
    EPSG: input coordinate system
        default: 3031 Polar Stereographic South, WGS84

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        http://www.numpy.org
        http://www.scipy.org/NumPy_for_Matlab_Users
    scipy: Scientific Tools for Python
        http://www.scipy.org/
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    calc_astrol_longitudes.py: computes the basic astronomical mean longitudes
    calc_delta_time.py: calculates difference between universal and dynamic time
    convert_xy_ll.py: convert lat/lon points to and from projected coordinates
    infer_minor_corrections.py: return corrections for 16 minor constituents
    load_constituent.py: loads parameters for a given tidal constituent
    load_nodal_corrections.py: load the nodal corrections for tidal constituents
    predict_tide_drift.py: predict tidal elevations using harmonic constants
    read_tide_model.py: extract tidal harmonic constants from OTIS tide models
    read_netcdf_model.py: extract tidal harmonic constants from netcdf models
    read_GOT_model.py: extract tidal harmonic constants from GSFC GOT models

UPDATE HISTORY:
    Written 03/2020
"""
from __future__ import print_function

import os
import datetime
import numpy as np
from pyTMD.count_leap_seconds import count_leap_seconds
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.infer_minor_corrections import infer_minor_corrections
from pyTMD.predict_tide_drift import predict_tide_drift
from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.read_netcdf_model import extract_netcdf_constants
from pyTMD.read_GOT_model import extract_GOT_constants

#-- PURPOSE: convert delta_time from time since EPOCH1 to time since EPOCH2
def convert_delta_time(delta_time, EPOCH1=None, EPOCH2=None, SCALE=(1./86400.)):
    """
    Convert delta time from seconds since EPOCH to time since EPOCH2
    """
    epoch1 = datetime.datetime(*EPOCH1)
    epoch2 = datetime.datetime(*EPOCH2)
    delta_time_epochs = (epoch2 - epoch1).total_seconds()
    #-- subtract difference in time and rescale to output units
    return SCALE*(delta_time - delta_time_epochs)

#-- PURPOSE: compute tides at points and times using tide model algorithms
def compute_tide_corrections(x, y, delta_time, DIRECTORY=None, MODEL=None,
    EPSG=3031, EPOCH=(2000,1,1,0,0,0), LEAPS=True):

    #-- select between tide models
    if (MODEL == 'CATS0201'):
        grid_file = os.path.join(tide_dir,'cats0201_tmd','grid_CATS')
        model_file = os.path.join(tide_dir,'cats0201_tmd','h0_CATS02_01')
        model_format = 'OTIS'
        model_EPSG = '4326'
        type = 'z'
    elif (MODEL == 'CATS2008'):
        grid_file = os.path.join(tide_dir,'CATS2008','grid_CATS2008')
        model_file = os.path.join(tide_dir,'CATS2008','hf.CATS2008.out')
        model_format = 'OTIS'
        model_EPSG = 'CATS2008'
        type = 'z'
    elif (MODEL == 'CATS2008_load'):
        grid_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','grid_CATS2008a_opt')
        model_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','h_CATS2008a_SPOTL_load')
        model_format = 'OTIS'
        model_EPSG = 'CATS2008'
        type = 'z'
    elif (MODEL == 'TPXO9-atlas'):
        model_directory = os.path.join(tide_dir,'TPXO9_atlas')
        grid_file = 'grid_tpxo9_atlas.nc.gz'
        model_files = ['h_q1_tpxo9_atlas_30.nc.gz','h_o1_tpxo9_atlas_30.nc.gz',
            'h_p1_tpxo9_atlas_30.nc.gz','h_k1_tpxo9_atlas_30.nc.gz',
            'h_n2_tpxo9_atlas_30.nc.gz','h_m2_tpxo9_atlas_30.nc.gz',
            'h_s2_tpxo9_atlas_30.nc.gz','h_k2_tpxo9_atlas_30.nc.gz',
            'h_m4_tpxo9_atlas_30.nc.gz','h_ms4_tpxo9_atlas_30.nc.gz',
            'h_mn4_tpxo9_atlas_30.nc.gz','h_2n2_tpxo9_atlas_30.nc.gz']
        model_format = 'netcdf'
        type = 'z'
        SCALE = 1.0/1000.0
    elif (MODEL == 'TPXO9.1'):
        grid_file = os.path.join(tide_dir,'TPXO9.1','DATA','grid_tpxo9')
        model_file = os.path.join(tide_dir,'TPXO9.1','DATA','h_tpxo9.v1')
        model_format = 'OTIS'
        model_EPSG = '4326'
        type = 'z'
    elif (MODEL == 'TPXO8-atlas'):
        grid_file = os.path.join(tide_dir,'tpxo8_atlas','grid_tpxo8atlas_30_v1')
        model_file = os.path.join(tide_dir,'tpxo8_atlas','hf.tpxo8_atlas_30_v1')
        model_format = 'ATLAS'
        model_EPSG = '4326'
        type = 'z'
    elif (MODEL == 'TPXO7.2'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_tmd','grid_tpxo7.2')
        model_file = os.path.join(tide_dir,'TPXO7.2_tmd','h_tpxo7.2')
        model_format = 'OTIS'
        model_EPSG = '4326'
        type = 'z'
    elif (MODEL == 'TPXO7.2_load'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_load','grid_tpxo6.2')
        model_file = os.path.join(tide_dir,'TPXO7.2_load','h_tpxo7.2_load')
        model_format = 'OTIS'
        model_EPSG = '4326'
        type = 'z'
    elif (MODEL == 'AODTM-5'):
        grid_file = os.path.join(tide_dir,'aodtm5_tmd','grid_Arc5km')
        model_file = os.path.join(tide_dir,'aodtm5_tmd','h0_Arc5km.oce')
        model_format = 'OTIS'
        model_EPSG = 'PSNorth'
        type = 'z'
    elif (MODEL == 'AOTIM-5'):
        grid_file = os.path.join(tide_dir,'aotim5_tmd','grid_Arc5km')
        model_file = os.path.join(tide_dir,'aotim5_tmd','h_Arc5km.oce')
        model_format = 'OTIS'
        model_EPSG = 'PSNorth'
        type = 'z'
    elif (MODEL == 'AOTIM-5-2018'):
        grid_file = os.path.join(tide_dir,'Arc5km2018','grid_Arc5km2018')
        model_file = os.path.join(tide_dir,'Arc5km2018','h_Arc5km2018')
        model_format = 'OTIS'
        model_EPSG = 'PSNorth'
        type = 'z'
    elif (MODEL == 'GOT4.7'):
        model_directory = os.path.join(tide_dir,'GOT4.7','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        model_format = 'GOT'
        SCALE = 1.0/100.0
    elif (MODEL == 'GOT4.7_load'):
        model_directory = os.path.join(tide_dir,'GOT4.7','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        model_format = 'GOT'
        SCALE = 1.0/1000.0
    elif (MODEL == 'GOT4.8'):
        model_directory = os.path.join(tide_dir,'got4.8','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        model_format = 'GOT'
        SCALE = 1.0/100.0
    elif (MODEL == 'GOT4.8_load'):
        model_directory = os.path.join(tide_dir,'got4.8','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        model_format = 'GOT'
        SCALE = 1.0/1000.0
    elif (MODEL == 'GOT4.10'):
        model_directory = os.path.join(tide_dir,'GOT4.10c','grids_oceantide')
        model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
            'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        model_format = 'GOT'
        SCALE = 1.0/100.0
    elif (MODEL == 'GOT4.10_load'):
        model_directory = os.path.join(tide_dir,'GOT4.10c','grids_loadtide')
        model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
            'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
            's1load.d.gz','m4load.d.gz']
        c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
        model_format = 'GOT'
        SCALE = 1.0/1000.0

    #-- converting x,y from EPSG to latitude/longitude
    proj1 = pyproj.Proj("+init=EPSG:{0:d}".format(EPSG))
    proj2 = pyproj.Proj("+init=EPSG:{0:d}".format(4326))
    lon,lat = pyproj.transform(proj1, proj2, x, y)

    #-- calculate leap seconds if specified
    if LEAPS:
        GPS_Time = convert_delta_time(delta_time, EPOCH1=EPOCH,
            EPOCH2=(1980,1,6,0,0,0), SCALE=1.0)
        leap_seconds = count_leap_seconds(GPS_Time)
    else:
        leap_seconds = 0.0
    #-- convert time to days relative to Jan 1, 1992 (48622mjd)
    t = convert_delta_time(delta_time - leap_seconds, EPOCH1=EPOCH,
        EPOCH2=(1992,1,1,0,0,0), SCALE=(1.0/86400.0))

    #-- read tidal constants and interpolate to grid points
    if model_format in ('OTIS','ATLAS'):
        amp,ph,D,c = extract_tidal_constants(lon, lat, grid_file, model_file,
            model_EPSG, type, METHOD=METHOD, GRID=model_format)
        deltat = np.zeros_like(t)
    elif model_format in ('netcdf'):
        amp,ph,D,c = extract_netcdf_constants(lon, lat, model_directory,
            grid_file, model_files, type, METHOD=METHOD, SCALE=SCALE)
        deltat = np.zeros_like(t)
    elif (model_format == 'GOT'):
        amp,ph = extract_GOT_constants(lon, lat, model_directory, model_files,
            METHOD=METHOD, SCALE=SCALE)
        #-- convert time to Modified Julian Days for calculating deltat
        delta_file = os.path.join(tide_dir,'deltat.data')
        deltat = calc_delta_time(delta_file, t + 48622.0)

    #-- calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    #-- calculate constituent oscillation
    hc = amp*np.exp(cph)

    #-- predict tidal elevations at time and infer minor corrections
    fill_value = np.nan
    tide = np.ma.zeros_like(delta_time,fill_value=fill_value)
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
