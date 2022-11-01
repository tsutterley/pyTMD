#!/usr/bin/env python
u"""
compute_tide_corrections.py
Written by Tyler Sutterley (11/2022)
Calculates tidal elevations for correcting elevation or imagery data

Uses OTIS format tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/
Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
or Finite Element Solution (FES) models provided by AVISO

INPUTS:
    x: x-coordinates in projection EPSG
    y: y-coordinates in projection EPSG
    delta_time: seconds since EPOCH or datetime array

OPTIONS:
    DIRECTORY: working data directory for tide models
    MODEL: Tide model to use in correction
    ATLAS_FORMAT: ATLAS tide model format (OTIS, netcdf)
    GZIP: Tide model files are gzip compressed
    DEFINITION_FILE: Tide model definition file for use as correction
    EPOCH: time period for calculating delta times
        default: J2000 (seconds since 2000-01-01T00:00:00)
    TYPE: input data type
        None: determined from input variable dimensions
        drift: drift buoys or satellite/airborne altimetry (time per data point)
        grid: spatial grids or images (single time for all data points)
        time series: time series at a single point (multiple times)
    TIME: input time standard or input type
        GPS: leap seconds needed
        LORAN: leap seconds needed (LORAN = GPS + 9 seconds)
        TAI: leap seconds needed (TAI = GPS + 19 seconds)
        UTC: no leap seconds needed
        datetime: numpy datatime array in UTC
    EPSG: input coordinate system
        default: 3031 Polar Stereographic South, WGS84
    METHOD: interpolation method
        bilinear: quick bilinear interpolation
        spline: scipy bivariate spline interpolation
        linear, nearest: scipy regular grid interpolations
    EXTRAPOLATE: extrapolate with nearest-neighbors
    CUTOFF: Extrapolation cutoff in kilometers
        set to np.inf to extrapolate for all points
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
    model.py: retrieves tide model parameters for named tide models
    spatial: utilities for reading, writing and operating on spatial data
    utilities.py: download and management utilities for syncing files
    calc_astrol_longitudes.py: computes the basic astronomical mean longitudes
    calc_delta_time.py: calculates difference between universal and dynamic time
    convert_ll_xy.py: convert lat/lon points to and from projected coordinates
    infer_minor_corrections.py: return corrections for minor constituents
    load_constituent.py: loads parameters for a given tidal constituent
    load_nodal_corrections.py: load the nodal corrections for tidal constituents
    predict_tide.py: predict tides at single times using harmonic constants
    predict_tidal_ts.py: predict tidal time series using harmonic constants
    predict_tide_drift.py: predict tidal elevations using harmonic constants
    read_tide_model.py: extract tidal harmonic constants from OTIS tide models
    read_netcdf_model.py: extract tidal harmonic constants from netcdf models
    read_GOT_model.py: extract tidal harmonic constants from GSFC GOT models
    read_FES_model.py: extract tidal harmonic constants from FES tide models
    bilinear_interp.py: bilinear interpolation of data to coordinates
    nearest_extrap.py: nearest-neighbor extrapolation of data to coordinates

UPDATE HISTORY:
    Updated 11/2022: place some imports within try/except statements
    Updated 05/2022: added ESR netCDF4 formats to list of model types
        updated keyword arguments to read tide model programs
        added option to apply flexure to heights for applicable models
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 12/2021: added function to calculate a tidal time series
        verify coordinate dimensions for each input data type
        added option for converting from LORAN times to UTC
    Updated 09/2021: refactor to use model class for files and attributes
    Updated 07/2021: can use numpy datetime arrays as input time variable
        added function for determining the input spatial variable type
        added check that tide model directory is accessible
    Updated 06/2021: added new Gr1km-v2 1km Greenland model from ESR
        add try/except for input projection strings
    Updated 05/2021: added option for extrapolation cutoff in kilometers
    Updated 03/2021: added TPXO9-atlas-v4 in binary OTIS format
        simplified netcdf inputs to be similar to binary OTIS read program
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 12/2020: added valid data extrapolation with nearest_extrap
    Updated 11/2020: added model constituents from TPXO9-atlas-v3
    Updated 08/2020: using builtin time operations.
        calculate difference in leap seconds from start of epoch
        using conversion protocols following pyproj-2 updates
    Updated 07/2020: added function docstrings, FES2014 and TPXO9-atlas-v2
        use merged delta time files combining biannual, monthly and daily files
    Updated 03/2020: added TYPE, TIME, FILL_VALUE and METHOD options
    Written 03/2020
"""
from __future__ import print_function

import os
import warnings
import numpy as np
import pyTMD.time
import pyTMD.model
import pyTMD.spatial
import pyTMD.utilities
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.infer_minor_corrections import infer_minor_corrections
from pyTMD.predict_tide import predict_tide
from pyTMD.predict_tidal_ts import predict_tidal_ts
from pyTMD.predict_tide_drift import predict_tide_drift
from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.read_netcdf_model import extract_netcdf_constants
from pyTMD.read_GOT_model import extract_GOT_constants
from pyTMD.read_FES_model import extract_FES_constants

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
def compute_tide_corrections(x, y, delta_time, DIRECTORY=None, MODEL=None,
    ATLAS_FORMAT='netcdf', GZIP=False, DEFINITION_FILE=None, EPSG=3031,
    EPOCH=(2000,1,1,0,0,0), TYPE='drift', TIME='UTC', METHOD='spline',
    EXTRAPOLATE=False, CUTOFF=10.0, APPLY_FLEXURE=False, FILL_VALUE=np.nan):
    """
    Compute tides at points and times using tidal harmonics

    Parameters
    ----------
    x: float
        x-coordinates in projection EPSG
    y: float
        y-coordinates in projection EPSG
    delta_time: float
        seconds since EPOCH or datetime array
    DIRECTORY: str or NoneType, default None
        working data directory for tide models
    MODEL: str or NoneType, default None
        Tide model to use in correction
    ATLAS_FORMAT: str, default 'netcdf'
        ATLAS tide model format (OTIS, netcdf)
        ATLAS tide model format

            - ``'OTIS'``
            - ``'netcdf'``
    GZIP: bool, default False
        Tide model files are gzip compressed
    DEFINITION_FILE: str or NoneType, default None
        Tide model definition file for use
    EPSG: int, default: 3031 (Polar Stereographic South, WGS84)
        Input coordinate system
    EPOCH: tuple, default (2000,1,1,0,0,0)
        Time period for calculating delta times
    TYPE: str or NoneType, default 'drift'
        Input data type

            - ``None``: determined from input variable dimensions
            - ``'drift'``: drift buoys or satellite/airborne altimetry
            - ``'grid'``: spatial grids or images
            - ``'time series'``: time series at a single point
    TIME: str, default 'UTC'
        Time type if need to compute leap seconds to convert to UTC

            - ``'GPS'``: leap seconds needed
            - ``'LORAN'``: leap seconds needed (LORAN = GPS + 9 seconds)
            - ``'TAI'``: leap seconds needed (TAI = GPS + 19 seconds)
            - ``'UTC'``: no leap seconds needed
            - ``'datetime'``: numpy datatime array in UTC
    METHOD: str
        Interpolation method

            - ```bilinear```: quick bilinear interpolation
            - ```spline```: scipy bivariate spline interpolation
            - ```linear```, ```nearest```: scipy regular grid interpolations

    EXTRAPOLATE: bool, default False
        Extrapolate with nearest-neighbors
    CUTOFF: float, default 10.0
        Extrapolation cutoff in kilometers

        Set to np.inf to extrapolate for all points
    APPLY_FLEXURE: bool, default False
        Apply ice flexure scaling factor to height constituents

        Only valid for models containing flexure fields
    FILL_VALUE: float, default np.nan
        Output invalid value

    Returns
    -------
    tide: float
        tidal elevation at coordinates and time in meters
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

    # determine input data type based on variable dimensions
    if not TYPE:
        TYPE = pyTMD.spatial.data_type(x, y, delta_time)
    # reform coordinate dimensions for input grids
    # or verify coordinate dimension shapes
    if (TYPE.lower() == 'grid') and (np.size(x) != np.size(y)):
        x,y = np.meshgrid(np.copy(x),np.copy(y))
    elif (TYPE.lower() == 'grid'):
        x = np.atleast_2d(x)
        y = np.atleast_2d(y)
    elif TYPE.lower() in ('time series', 'drift'):
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)

    # converting x,y from EPSG to latitude/longitude
    try:
        # EPSG projection code string or int
        crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(int(EPSG)))
    except (ValueError,pyproj.exceptions.CRSError):
        # Projection SRS string
        crs1 = pyproj.CRS.from_string(EPSG)
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon,lat = transformer.transform(x.flatten(), y.flatten())

    # assert delta time is an array
    delta_time = np.atleast_1d(delta_time)
    # calculate leap seconds if specified
    if (TIME.upper() == 'GPS'):
        GPS_Epoch_Time = pyTMD.time.convert_delta_time(0, epoch1=EPOCH,
            epoch2=(1980,1,6,0,0,0), scale=1.0)
        GPS_Time = pyTMD.time.convert_delta_time(delta_time, epoch1=EPOCH,
            epoch2=(1980,1,6,0,0,0), scale=1.0)
        # calculate difference in leap seconds from start of epoch
        leap_seconds = pyTMD.time.count_leap_seconds(GPS_Time) - \
            pyTMD.time.count_leap_seconds(np.atleast_1d(GPS_Epoch_Time))
    elif (TIME.upper() == 'LORAN'):
        # LORAN time is ahead of GPS time by 9 seconds
        GPS_Epoch_Time = pyTMD.time.convert_delta_time(-9.0, epoch1=EPOCH,
            epoch2=(1980,1,6,0,0,0), scale=1.0)
        GPS_Time = pyTMD.time.convert_delta_time(delta_time-9.0, epoch1=EPOCH,
            epoch2=(1980,1,6,0,0,0), scale=1.0)
        # calculate difference in leap seconds from start of epoch
        leap_seconds = pyTMD.time.count_leap_seconds(GPS_Time) - \
            pyTMD.time.count_leap_seconds(np.atleast_1d(GPS_Epoch_Time))
    elif (TIME.upper() == 'TAI'):
        # TAI time is ahead of GPS time by 19 seconds
        GPS_Epoch_Time = pyTMD.time.convert_delta_time(-19.0, epoch1=EPOCH,
            epoch2=(1980,1,6,0,0,0), scale=1.0)
        GPS_Time = pyTMD.time.convert_delta_time(delta_time-19.0, epoch1=EPOCH,
            epoch2=(1980,1,6,0,0,0), scale=1.0)
        # calculate difference in leap seconds from start of epoch
        leap_seconds = pyTMD.time.count_leap_seconds(GPS_Time) - \
            pyTMD.time.count_leap_seconds(np.atleast_1d(GPS_Epoch_Time))
    else:
        leap_seconds = 0.0

    # convert delta times or datetimes objects
    if (TIME.lower() == 'datetime'):
        # convert delta time array from datetime object
        # to days relative to 1992-01-01T00:00:00
        t = pyTMD.time.convert_datetime(delta_time,
            epoch=(1992,1,1,0,0,0))/86400.0
    else:
        # convert time to days relative to Jan 1, 1992 (48622mjd)
        t = pyTMD.time.convert_delta_time(delta_time - leap_seconds,
            epoch1=EPOCH, epoch2=(1992,1,1,0,0,0), scale=(1.0/86400.0))
    # delta time (TT - UT1) file
    delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])

    # read tidal constants and interpolate to grid points
    if model.format in ('OTIS','ATLAS','ESR'):
        amp,ph,D,c = extract_tidal_constants(lon, lat, model.grid_file,
            model.model_file, model.projection, type=model.type,
            method=METHOD, extrapolate=EXTRAPOLATE, cutoff=CUTOFF,
            grid=model.format, apply_flexure=APPLY_FLEXURE)
        deltat = np.zeros_like(t)
    elif (model.format == 'netcdf'):
        amp,ph,D,c = extract_netcdf_constants(lon, lat, model.grid_file,
            model.model_file, type=model.type, method=METHOD,
            extrapolate=EXTRAPOLATE, cutoff=CUTOFF, scale=model.scale,
            compressed=model.compressed)
        deltat = np.zeros_like(t)
    elif (model.format == 'GOT'):
        amp,ph,c = extract_GOT_constants(lon, lat, model.model_file,
            method=METHOD, extrapolate=EXTRAPOLATE, cutoff=CUTOFF,
            scale=model.scale, compressed=model.compressed)
        # interpolate delta times from calendar dates to tide time
        deltat = calc_delta_time(delta_file, t)
    elif (model.format == 'FES'):
        amp,ph = extract_FES_constants(lon, lat, model.model_file,
            type=model.type, version=model.version, method=METHOD,
            extrapolate=EXTRAPOLATE, cutoff=CUTOFF, scale=model.scale,
            compressed=model.compressed)
        # available model constituents
        c = model.constituents
        # interpolate delta times from calendar dates to tide time
        deltat = calc_delta_time(delta_file, t)

    # calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    # calculate constituent oscillation
    hc = amp*np.exp(cph)

    # predict tidal elevations at time and infer minor corrections
    if (TYPE.lower() == 'grid'):
        ny,nx = np.shape(x); nt = len(t)
        tide = np.ma.zeros((ny,nx,nt),fill_value=FILL_VALUE)
        tide.mask = np.zeros((ny,nx,nt),dtype=bool)
        for i in range(nt):
            TIDE = predict_tide(t[i], hc, c,
                deltat=deltat[i], corrections=model.format)
            MINOR = infer_minor_corrections(t[i], hc, c,
                deltat=deltat[i], corrections=model.format)
            # add major and minor components and reform grid
            tide[:,:,i] = np.reshape((TIDE+MINOR), (ny,nx))
            tide.mask[:,:,i] = np.reshape((TIDE.mask | MINOR.mask), (ny,nx))
    elif (TYPE.lower() == 'drift'):
        npts = len(t)
        tide = np.ma.zeros((npts), fill_value=FILL_VALUE)
        tide.mask = np.any(hc.mask,axis=1)
        tide.data[:] = predict_tide_drift(t, hc, c,
            deltat=deltat, corrections=model.format)
        minor = infer_minor_corrections(t, hc, c,
            deltat=deltat, corrections=model.format)
        tide.data[:] += minor.data[:]
    elif (TYPE.lower() == 'time series'):
        npts = len(t)
        tide = np.ma.zeros((npts), fill_value=FILL_VALUE)
        tide.mask = np.any(hc.mask,axis=1)
        tide.data[:] = predict_tidal_ts(t, hc, c,
            deltat=deltat, corrections=model.format)
        minor = infer_minor_corrections(t, hc, c,
            deltat=deltat, corrections=model.format)
        tide.data[:] += minor.data[:]
    # replace invalid values with fill value
    tide.data[tide.mask] = tide.fill_value

    # return the tide correction
    return tide
