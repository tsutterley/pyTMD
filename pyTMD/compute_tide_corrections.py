#!/usr/bin/env python
u"""
compute_tide_corrections.py
Written by Tyler Sutterley (05/2023)
Calculates tidal elevations for correcting elevation or imagery data

Ocean and Load Tides
Uses OTIS format tidal solutions provided by Ohio State University and ESR
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/
Global Tide Model (GOT) solutions provided by Richard Ray at GSFC
or Finite Element Solution (FES) models provided by AVISO

Long-Period Equilibrium Tides (LPET)
Calculates long-period equilibrium tidal elevations for correcting
elevation or imagery data from the summation of fifteen spectral lines
    https://doi.org/10.1111/j.1365-246X.1973.tb03420.x

Load Pole Tides (LPT)
Calculates radial load pole tide displacements following IERS Convention
(2010) guidelines for correcting elevation or imagery data
    https://iers-conventions.obspm.fr/chapter7.php

Ocean Pole Tides (OPT)
Calculates radial ocean pole load tide displacements following IERS Convention
(2010) guidelines for correcting elevation or imagery data
    https://iers-conventions.obspm.fr/chapter7.php

Ocean Pole Tides (SET)
Calculates radial Solid Earth tide displacements following IERS Convention
(2010) guidelines for correcting elevation or imagery data
    https://iers-conventions.obspm.fr/chapter7.php

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
    spatial: utilities for reading, writing and operating on spatial data
    utilities.py: download and management utilities for syncing files
    arguments.py: load the nodal corrections for tidal constituents
    astro.py: computes the basic astronomical mean longitudes
    convert_crs.py: convert points to and from Coordinates Reference Systems
    load_constituent.py: loads parameters for a given tidal constituent
    predict.py: predict tide values using harmonic constants
    io/model.py: retrieves tide model parameters for named tide models
    io/OTIS.py: extract tidal harmonic constants from OTIS tide models
    io/ATLAS.py: extract tidal harmonic constants from netcdf models
    io/GOT.py: extract tidal harmonic constants from GSFC GOT models
    io/FES.py: extract tidal harmonic constants from FES tide models
    interpolate.py: interpolation routines for spatial data

UPDATE HISTORY:
    Updated 05/2023: use timescale class for time conversion operations
        use defaults from eop module for pole tide and EOP files
        add option for using higher resolution ephemerides from JPL
    Updated 04/2023: added function for radial solid earth tides
        using pathlib to define and expand paths
    Updated 03/2023: add basic variable typing to function inputs
        added function for long-period equilibrium tides
        added function for radial load pole tides
        added function for radial ocean pole tides
    Updated 12/2022: refactored tide read and prediction programs
    Updated 11/2022: place some imports within try/except statements
        use f-strings for formatting verbose or ascii output
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
from __future__ import print_function, annotations

import logging
import pathlib
import numpy as np
import scipy.interpolate
import pyTMD.constants
import pyTMD.eop
import pyTMD.io
import pyTMD.time
import pyTMD.io.model
import pyTMD.predict
import pyTMD.spatial
import pyTMD.utilities

# attempt imports
try:
    import pyproj
except (ImportError, ModuleNotFoundError) as exc:
    logging.critical("pyproj not available")

# PURPOSE: wrapper function for computing corrections
def compute_corrections(
        x: np.ndarray, y: np.ndarray, delta_time: np.ndarray,
        CORRECTION: str = 'ocean',
        **kwargs
    ):
    """
    Wrapper function to compute tide corrections at points and times

    Parameters
    ----------
    x: np.ndarray
        x-coordinates in projection EPSG
    y: np.ndarray
        y-coordinates in projection EPSG
    delta_time: np.ndarray
        seconds since EPOCH or datetime array
    CORRECTION: str, default 'ocean'
        Correction type to compute

            - ``'ocean'``: ocean tide from model constituents
            - ``'load'``: load tide from model constituents
            - ``'LPET'``: long-period equilibrium tide
            - ``'LPT'``: solid earth load pole tide
            - ``'OPT'``: ocean pole tide
            - ``'SET'``: solid earth tide
    **kwargs: dict
        keyword arguments for correction functions

    Returns
    -------
    correction: np.ndarray
        tidal correction at coordinates and time in meters
    """
    if CORRECTION.lower() in ('ocean', 'load'):
        return compute_tide_corrections(x, y, delta_time, **kwargs)
    elif (CORRECTION.upper() == 'LPET'):
        return compute_LPET_corrections(x, y, delta_time, **kwargs)
    elif (CORRECTION.upper() == 'LPT'):
        return compute_LPT_corrections(x, y, delta_time, **kwargs)
    elif (CORRECTION.upper() == 'OPT'):
        return compute_OPT_corrections(x, y, delta_time, **kwargs)
    elif (CORRECTION.upper() == 'SET'):
        return compute_SET_corrections(x, y, delta_time, **kwargs)
    else:
        raise ValueError(f'Unrecognized correction type: {CORRECTION}')

# PURPOSE: compute tides at points and times using tide model algorithms
def compute_tide_corrections(
        x: np.ndarray, y: np.ndarray, delta_time: np.ndarray,
        DIRECTORY: str | pathlib.Path | None = None,
        MODEL: str | None = None,
        ATLAS_FORMAT: str = 'netcdf',
        GZIP: bool = False,
        DEFINITION_FILE: str | pathlib.Path | None = None,
        EPSG: str | int = 3031,
        EPOCH: list | tuple = (2000, 1, 1, 0, 0, 0),
        TYPE: str or None = 'drift',
        TIME: str = 'UTC',
        METHOD: str = 'spline',
        EXTRAPOLATE: bool = False,
        CUTOFF: int | float=10.0,
        APPLY_FLEXURE: bool = False,
        FILL_VALUE: float = np.nan,
        **kwargs
    ):
    """
    Compute ocean or load tides at points and times from
    model constituents

    Parameters
    ----------
    x: np.ndarray
        x-coordinates in projection EPSG
    y: np.ndarray
        y-coordinates in projection EPSG
    delta_time: np.ndarray
        seconds since EPOCH or datetime array
    DIRECTORY: str or NoneType, default None
        working data directory for tide models
    MODEL: str or NoneType, default None
        Tide model to use in correction
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
    CUTOFF: int or float, default 10.0
        Extrapolation cutoff in kilometers

        Set to ``np.inf`` to extrapolate for all points
    APPLY_FLEXURE: bool, default False
        Apply ice flexure scaling factor to height constituents

        Only valid for models containing flexure fields
    FILL_VALUE: float, default np.nan
        Output invalid value

    Returns
    -------
    tide: np.ndarray
        tidal elevation at coordinates and time in meters
    """

    # check that tide directory is accessible
    if DIRECTORY is not None:
        DIRECTORY = pathlib.Path(DIRECTORY).expanduser()
        if not DIRECTORY.exists():
            raise FileNotFoundError("Invalid tide directory")

    # validate input arguments
    assert TIME in ('GPS', 'LORAN', 'TAI', 'UTC', 'datetime')
    assert METHOD in ('bilinear', 'spline', 'linear', 'nearest')

    # get parameters for tide model
    if DEFINITION_FILE is not None:
        model = pyTMD.io.model(DIRECTORY).from_file(
            pathlib.Path(DEFINITION_FILE).expanduser())
    else:
        model = pyTMD.io.model(DIRECTORY, format=ATLAS_FORMAT,
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
        crs1 = pyproj.CRS.from_epsg(int(EPSG))
    except (ValueError,pyproj.exceptions.CRSError):
        # Projection SRS string
        crs1 = pyproj.CRS.from_string(EPSG)
    # output coordinate reference system
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon, lat = transformer.transform(x.flatten(), y.flatten())

    # assert delta time is an array
    delta_time = np.atleast_1d(delta_time)
    # convert delta times or datetimes objects to timescale
    if (TIME.lower() == 'datetime'):
        timescale = pyTMD.time.timescale().from_datetime(
            delta_time.flatten())
    else:
        timescale = pyTMD.time.timescale().from_deltatime(delta_time,
            epoch=EPOCH, standard=TIME)
    # number of time points
    nt = len(timescale)

    # read tidal constants and interpolate to grid points
    if model.format in ('OTIS','ATLAS','ESR'):
        amp,ph,D,c = pyTMD.io.OTIS.extract_constants(lon, lat, model.grid_file,
            model.model_file, model.projection, type=model.type,
            method=METHOD, extrapolate=EXTRAPOLATE, cutoff=CUTOFF,
            grid=model.format, apply_flexure=APPLY_FLEXURE)
        # use delta time at 2000.0 to match TMD outputs
        deltat = np.zeros((nt), dtype=np.float64)
    elif (model.format == 'netcdf'):
        amp,ph,D,c = pyTMD.io.ATLAS.extract_constants(lon, lat, model.grid_file,
            model.model_file, type=model.type, method=METHOD,
            extrapolate=EXTRAPOLATE, cutoff=CUTOFF, scale=model.scale,
            compressed=model.compressed)
        # use delta time at 2000.0 to match TMD outputs
        deltat = np.zeros((nt), dtype=np.float64)
    elif (model.format == 'GOT'):
        amp,ph,c = pyTMD.io.GOT.extract_constants(lon, lat, model.model_file,
            method=METHOD, extrapolate=EXTRAPOLATE, cutoff=CUTOFF,
            scale=model.scale, compressed=model.compressed)
        # delta time (TT - UT1)
        deltat = timescale.tt_ut1
    elif (model.format == 'FES'):
        amp,ph = pyTMD.io.FES.extract_constants(lon, lat, model.model_file,
            type=model.type, version=model.version, method=METHOD,
            extrapolate=EXTRAPOLATE, cutoff=CUTOFF, scale=model.scale,
            compressed=model.compressed)
        # available model constituents
        c = model.constituents
        # delta time (TT - UT1)
        deltat = timescale.tt_ut1

    # calculate complex phase in radians for Euler's
    cph = -1j*ph*np.pi/180.0
    # calculate constituent oscillation
    hc = amp*np.exp(cph)

    # predict tidal elevations at time and infer minor corrections
    if (TYPE.lower() == 'grid'):
        ny,nx = np.shape(x)
        tide = np.ma.zeros((ny,nx,nt),fill_value=FILL_VALUE)
        tide.mask = np.zeros((ny,nx,nt),dtype=bool)
        for i in range(nt):
            TIDE = pyTMD.predict.map(timescale.tide[i], hc, c,
                deltat=deltat[i], corrections=model.format)
            MINOR = pyTMD.predict.infer_minor(timescale.tide[i], hc, c,
                deltat=deltat[i], corrections=model.format)
            # add major and minor components and reform grid
            tide[:,:,i] = np.reshape((TIDE+MINOR), (ny,nx))
            tide.mask[:,:,i] = np.reshape((TIDE.mask | MINOR.mask), (ny,nx))
    elif (TYPE.lower() == 'drift'):
        tide = np.ma.zeros((nt), fill_value=FILL_VALUE)
        tide.mask = np.any(hc.mask,axis=1)
        tide.data[:] = pyTMD.predict.drift(timescale.tide, hc, c,
            deltat=deltat, corrections=model.format)
        minor = pyTMD.predict.infer_minor(timescale.tide, hc, c,
            deltat=deltat, corrections=model.format)
        tide.data[:] += minor.data[:]
    elif (TYPE.lower() == 'time series'):
        nstation = len(x)
        tide = np.ma.zeros((nstation,nt), fill_value=FILL_VALUE)
        tide.mask = np.zeros((nstation,nt),dtype=bool)
        for s in range(nstation):
            TIDE = pyTMD.predict.time_series(timescale.tide, hc[s,None,:], c,
                deltat=deltat, corrections=model.format)
            MINOR = pyTMD.predict.infer_minor(timescale.tide, hc[s,None,:], c,
                deltat=deltat, corrections=model.format)
            tide.data[s,:] = TIDE.data[:] + MINOR.data[:]
            tide.mask[s,:] = (TIDE.mask | MINOR.mask)
    # replace invalid values with fill value
    tide.data[tide.mask] = tide.fill_value

    # return the ocean or load tide correction
    return tide

# PURPOSE: compute long-period equilibrium tidal elevations
def compute_LPET_corrections(
        x: np.ndarray, y: np.ndarray, delta_time: np.ndarray,
        EPSG: str | int = 3031,
        EPOCH: list | tuple = (2000, 1, 1, 0, 0, 0),
        TYPE: str or None = 'drift',
        TIME: str = 'UTC',
        **kwargs
    ):
    """
    Compute long-period equilibrium tidal elevations at points and times

    Parameters
    ----------
    x: np.ndarray
        x-coordinates in projection EPSG
    y: np.ndarray
        y-coordinates in projection EPSG
    delta_time: np.ndarray
        seconds since EPOCH or datetime array
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
    FILL_VALUE: float, default np.nan
        Output invalid value

    Returns
    -------
    tide_lpe: np.ndarray
        long-period equilibrium tide at coordinates and time in meters
    """

    # validate input arguments
    assert TIME in ('GPS', 'LORAN', 'TAI', 'UTC', 'datetime')
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
        crs1 = pyproj.CRS.from_epsg(int(EPSG))
    except (ValueError,pyproj.exceptions.CRSError):
        # Projection SRS string
        crs1 = pyproj.CRS.from_string(EPSG)
    # output coordinate reference system
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon, lat = transformer.transform(x.flatten(), y.flatten())

    # assert delta time is an array
    delta_time = np.atleast_1d(delta_time)
    # convert delta times or datetimes objects to timescale
    if (TIME.lower() == 'datetime'):
        timescale = pyTMD.time.timescale().from_datetime(
            delta_time.flatten())
    else:
        timescale = pyTMD.time.timescale().from_deltatime(delta_time,
            epoch=EPOCH, standard=TIME)
    # number of time points
    nt = len(timescale)
    # convert tide times to dynamic time
    tide_time = timescale.tide + timescale.tt_ut1

    # predict long-period equilibrium tides at time
    if (TYPE == 'grid'):
        ny,nx = np.shape(x)
        tide_lpe = np.zeros((ny,nx,nt))
        for i in range(nt):
            lpet = pyTMD.predict.equilibrium_tide(tide_time[i], lat)
            tide_lpe[:,:,i] = np.reshape(lpet,(ny,nx))
    elif (TYPE == 'drift'):
        tide_lpe = pyTMD.predict.equilibrium_tide(tide_time, lat)
    elif (TYPE == 'time series'):
        nstation = len(x)
        tide_lpe = np.zeros((nstation,nt))
        for s in range(nstation):
            lpet = pyTMD.predict.equilibrium_tide(tide_time, lat[s])
            tide_lpe[s,:] = np.copy(lpet)

    # return the long-period equilibrium tide corrections
    return tide_lpe

# PURPOSE: compute radial load pole tide displacements
# following IERS Convention (2010) guidelines
def compute_LPT_corrections(
        x: np.ndarray, y: np.ndarray, delta_time: np.ndarray,
        h: float | np.ndarray = 0.0,
        EPSG: str | int = 3031,
        EPOCH: list | tuple = (2000, 1, 1, 0, 0, 0),
        TYPE: str or None = 'drift',
        TIME: str = 'UTC',
        ELLIPSOID: str = 'WGS84',
        CONVENTION: str = '2018',
        FILL_VALUE: float = np.nan,
        **kwargs
    ):
    """
    Compute radial load pole tide displacements at points and times
    following IERS Convention (2010) guidelines

    Parameters
    ----------
    x: np.ndarray
        x-coordinates in projection EPSG
    y: np.ndarray
        y-coordinates in projection EPSG
    delta_time: np.ndarray
        seconds since EPOCH or datetime array
    h: float or np.ndarray, default 0.0
        height coordinates in meters
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
    ELLIPSOID: str, default 'WGS84'
        Ellipsoid for calculating Earth parameters
    CONVENTION: str, default '2018'
        IERS Mean or Secular Pole Convention

            - ``'2003'``
            - ``'2010'``
            - ``'2015'``
            - ``'2018'``
    FILL_VALUE: float, default np.nan
        Output invalid value

    Returns
    -------
    Srad: np.ndarray
        solid earth pole tide at coordinates and time in meters
    """

    # validate input arguments
    assert TIME in ('GPS', 'LORAN', 'TAI', 'UTC', 'datetime')
    assert ELLIPSOID in pyTMD._ellipsoids
    assert CONVENTION in pyTMD.eop._conventions
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
        crs1 = pyproj.CRS.from_epsg(int(EPSG))
    except (ValueError,pyproj.exceptions.CRSError):
        # Projection SRS string
        crs1 = pyproj.CRS.from_string(EPSG)
    # output coordinate reference system
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon,lat = transformer.transform(x.flatten(), y.flatten())

    # assert delta time is an array
    delta_time = np.atleast_1d(delta_time)
    # convert delta times or datetimes objects to timescale
    if (TIME.lower() == 'datetime'):
        timescale = pyTMD.time.timescale().from_datetime(
            delta_time.flatten())
    else:
        timescale = pyTMD.time.timescale().from_deltatime(delta_time,
            epoch=EPOCH, standard=TIME)

    # convert dynamic time to Modified Julian Days (MJD)
    MJD = timescale.tt - 2400000.5
    # convert Julian days to calendar dates
    Y,M,D,h,m,s = pyTMD.time.convert_julian(timescale.tt, format='tuple')
    # calculate time in year-decimal format
    time_decimal = pyTMD.time.convert_calendar_decimal(Y, M, day=D,
        hour=h, minute=m, second=s)
    # number of time points
    nt = len(time_decimal)

    # degrees and arcseconds to radians
    dtr = np.pi/180.0
    atr = np.pi/648000.0
    # earth and physical parameters for ellipsoid
    units = pyTMD.constants(ELLIPSOID)
    # tidal love number appropriate for the load tide
    hb2 = 0.6207

    # flatten heights
    h = np.array(h).flatten() if np.ndim(h) else h
    # convert from geodetic latitude to geocentric latitude
    # calculate X, Y and Z from geodetic latitude and longitude
    X,Y,Z = pyTMD.spatial.to_cartesian(lon.flatten(), lat.flatten(), h=h,
        a_axis=units.a_axis, flat=units.flat)
    rr = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/dtr
    # geocentric colatitude and longitude in radians
    theta = dtr*(90.0 - latitude_geocentric)
    phi = dtr*lon.flatten()

    # compute normal gravity at spatial location and elevation of points.
    # Normal gravity at height h. p. 82, Eqn.(2-215)
    gamma_h = units.gamma_h(theta, h)
    dfactor = -hb2*atr*(units.omega**2*rr**2)/(2.0*gamma_h)

    # calculate angular coordinates of mean/secular pole at time
    mpx, mpy, fl = pyTMD.eop.iers_mean_pole(time_decimal, convention=CONVENTION)
    # read and interpolate IERS daily polar motion values
    px, py = pyTMD.eop.iers_polar_motion(MJD, k=3, s=0)
    # calculate differentials from mean/secular pole positions
    mx = px - mpx
    my = -(py - mpy)

    # calculate radial displacement at time
    if (TYPE == 'grid'):
        ny,nx = np.shape(x)
        Srad = np.ma.zeros((ny,nx,nt), fill_value=FILL_VALUE)
        Srad.mask = np.zeros((ny,nx,nt),dtype=bool)
        for i in range(nt):
            SRAD = dfactor*np.sin(2.0*theta)*(mx[i]*np.cos(phi)+my[i]*np.sin(phi))
            # reform grid
            Srad.data[:,:,i] = np.reshape(SRAD, (ny,nx))
            Srad.mask[:,:,i] = np.isnan(Srad.data[:,:,i])
    elif (TYPE == 'drift'):
        Srad = np.ma.zeros((nt), fill_value=FILL_VALUE)
        Srad.data[:] = dfactor*np.sin(2.0*theta)*(mx*np.cos(phi)+my*np.sin(phi))
        Srad.mask = np.isnan(Srad.data)
    elif (TYPE == 'time series'):
        nstation = len(x)
        Srad = np.ma.zeros((nstation,nt), fill_value=FILL_VALUE)
        Srad.mask = np.zeros((nstation,nt),dtype=bool)
        for s in range(nstation):
            SRAD = dfactor[s]*np.sin(2.0*theta[s])*(mx*np.cos(phi[s])+my*np.sin(phi[s]))
            Srad.data[s,:] = np.copy(SRAD)
            Srad.mask[s,:] = np.isnan(Srad.data[s,:])
    # replace invalid data with fill values
    Srad.data[Srad.mask] = Srad.fill_value

    # return the solid earth pole tide corrections
    return Srad

# PURPOSE: compute radial load pole tide displacements
# following IERS Convention (2010) guidelines
def compute_OPT_corrections(
        x: np.ndarray, y: np.ndarray, delta_time: np.ndarray,
        h: float | np.ndarray = 0.0,
        EPSG: str | int = 3031,
        EPOCH: list | tuple = (2000, 1, 1, 0, 0, 0),
        TYPE: str or None = 'drift',
        TIME: str = 'UTC',
        ELLIPSOID: str = 'WGS84',
        CONVENTION: str = '2018',
        METHOD: str = 'spline',
        FILL_VALUE: float = np.nan,
        **kwargs
    ):
    """
    Compute radial ocean pole tide displacements at points and times
    following IERS Convention (2010) guidelines

    Parameters
    ----------
    x: np.ndarray
        x-coordinates in projection EPSG
    y: np.ndarray
        y-coordinates in projection EPSG
    delta_time: np.ndarray
        seconds since EPOCH or datetime array
    h: float or np.ndarray, default 0.0
        height coordinates in meters
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
    ELLIPSOID: str, default 'WGS84'
        Ellipsoid for calculating Earth parameters
    CONVENTION: str, default '2018'
        IERS Mean or Secular Pole Convention

            - ``'2003'``
            - ``'2010'``
            - ``'2015'``
            - ``'2018'``
    METHOD: str
        Interpolation method

            - ```bilinear```: quick bilinear interpolation
            - ```spline```: scipy bivariate spline interpolation
            - ```linear```, ```nearest```: scipy regular grid interpolations
    FILL_VALUE: float, default np.nan
        Output invalid value

    Returns
    -------
    Urad: np.ndarray
        ocean pole tide at coordinates and time in meters
    """

    # validate input arguments
    assert TIME in ('GPS', 'LORAN', 'TAI', 'UTC', 'datetime')
    assert ELLIPSOID in pyTMD._ellipsoids
    assert CONVENTION in pyTMD.eop._conventions
    assert METHOD in ('bilinear', 'spline', 'linear', 'nearest')
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
        crs1 = pyproj.CRS.from_epsg(int(EPSG))
    except (ValueError,pyproj.exceptions.CRSError):
        # Projection SRS string
        crs1 = pyproj.CRS.from_string(EPSG)
    # output coordinate reference system
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon,lat = transformer.transform(x.flatten(), y.flatten())

    # assert delta time is an array
    delta_time = np.atleast_1d(delta_time)
    # convert delta times or datetimes objects to timescale
    if (TIME.lower() == 'datetime'):
        timescale = pyTMD.time.timescale().from_datetime(
            delta_time.flatten())
    else:
        timescale = pyTMD.time.timescale().from_deltatime(delta_time,
            epoch=EPOCH, standard=TIME)

    # convert dynamic time to Modified Julian Days (MJD)
    MJD = timescale.tt - 2400000.5
    # convert Julian days to calendar dates
    Y,M,D,h,m,s = pyTMD.time.convert_julian(timescale.tt, format='tuple')
    # calculate time in year-decimal format
    time_decimal = pyTMD.time.convert_calendar_decimal(Y, M, day=D,
        hour=h, minute=m, second=s)
    # number of time points
    nt = len(time_decimal)

    # degrees and arcseconds to radians
    dtr = np.pi/180.0
    atr = np.pi/648000.0
    # earth and physical parameters for ellipsoid
    units = pyTMD.constants(ELLIPSOID)
    # mean equatorial gravitational acceleration [m/s^2]
    ge = 9.7803278
    # density of sea water [kg/m^3]
    rho_w = 1025.0
    # tidal love number differential (1 + kl - hl) for pole tide frequencies
    gamma = 0.6870 + 0.0036j

    # flatten heights
    h = np.array(h).flatten() if np.ndim(h) else h
    # convert from geodetic latitude to geocentric latitude
    # calculate X, Y and Z from geodetic latitude and longitude
    X,Y,Z = pyTMD.spatial.to_cartesian(lon.flatten(), lat.flatten(), h=h,
        a_axis=units.a_axis, flat=units.flat)
    rr = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/dtr

    # pole tide displacement scale factor
    Hp = np.sqrt(8.0*np.pi/15.0)*(units.omega**2*units.a_axis**4)/units.GM
    K = 4.0*np.pi*units.G*rho_w*Hp*units.a_axis/(3.0*ge)
    K1 = 4.0*np.pi*units.G*rho_w*Hp*units.a_axis**3/(3.0*units.GM)

    # calculate angular coordinates of mean/secular pole at time
    mpx, mpy, fl = pyTMD.eop.iers_mean_pole(time_decimal, convention=CONVENTION)
    # read and interpolate IERS daily polar motion values
    px, py = pyTMD.eop.iers_polar_motion(MJD, k=3, s=0)
    # calculate differentials from mean/secular pole positions
    mx = px - mpx
    my = -(py - mpy)

    # read ocean pole tide map from Desai (2002)
    iur, iun, iue, ilon, ilat = pyTMD.io.ocean_pole_tide()
    # interpolate ocean pole tide map from Desai (2002)
    if (METHOD == 'spline'):
        # use scipy bivariate splines to interpolate to output points
        f1 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
            iur[:,::-1].real, kx=1, ky=1)
        f2 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
            iur[:,::-1].imag, kx=1, ky=1)
        UR = np.zeros((len(latitude_geocentric)), dtype=np.longcomplex)
        UR.real = f1.ev(lon.flatten(), latitude_geocentric)
        UR.imag = f2.ev(lon.flatten(), latitude_geocentric)
    else:
        # use scipy regular grid to interpolate values for a given method
        r1 = scipy.interpolate.RegularGridInterpolator((ilon,ilat[::-1]),
            iur[:,::-1], method=METHOD)
        UR = r1.__call__(np.c_[lon.flatten(), latitude_geocentric])

    # calculate radial displacement at time
    if (TYPE == 'grid'):
        ny,nx = np.shape(x)
        Urad = np.ma.zeros((ny,nx,nt),fill_value=FILL_VALUE)
        Urad.mask = np.zeros((ny,nx,nt),dtype=bool)
        for i in range(nt):
            URAD = K*atr*np.real((mx[i]*gamma.real + my[i]*gamma.imag)*UR.real +
                (my[i]*gamma.real - mx[i]*gamma.imag)*UR.imag)
            # reform grid
            Urad.data[:,:,i] = np.reshape(URAD, (ny,nx))
            Urad.mask[:,:,i] = np.isnan(Urad.data[:,:,i])
    elif (TYPE == 'drift'):
        Urad = np.ma.zeros((nt),fill_value=FILL_VALUE)
        Urad.data[:] = K*atr*np.real((mx*gamma.real + my*gamma.imag)*UR.real +
            (my*gamma.real - mx*gamma.imag)*UR.imag)
        Urad.mask = np.isnan(Urad.data)
    elif (TYPE == 'time series'):
        nstation = len(x)
        Urad = np.ma.zeros((nstation,nt),fill_value=FILL_VALUE)
        Urad.mask = np.zeros((nstation,nt),dtype=bool)
        for s in range(nstation):
            URAD = K*atr*np.real((mx*gamma.real + my*gamma.imag)*UR.real[s] +
                (my*gamma.real - mx*gamma.imag)*UR.imag[s])
            Urad.data[s,:] = np.copy(URAD)
            Urad.mask[s,:] = np.isnan(Urad.data[s,:])
    # replace invalid data with fill values
    Urad.data[Urad.mask] = Urad.fill_value

    # return the ocean pole tide corrections
    return Urad

# PURPOSE: compute solid earth tidal elevations
def compute_SET_corrections(
        x: np.ndarray, y: np.ndarray, delta_time: np.ndarray,
        h: float | np.ndarray = 0.0,
        EPSG: str | int = 3031,
        EPOCH: list | tuple = (2000, 1, 1, 0, 0, 0),
        TYPE: str or None = 'drift',
        TIME: str = 'UTC',
        ELLIPSOID: str = 'WGS84',
        TIDE_SYSTEM='tide_free',
        EPHEMERIDES='approximate',
        **kwargs
    ):
    """
    Compute solid earth tidal elevations at points and times
    following IERS Convention (2010) guidelines

    Parameters
    ----------
    x: np.ndarray
        x-coordinates in projection EPSG
    y: np.ndarray
        y-coordinates in projection EPSG
    delta_time: np.ndarray
        seconds since EPOCH or datetime array
    h: float or np.ndarray, default 0.0
        height coordinates in meters
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
    ELLIPSOID: str, default 'WGS84'
        Ellipsoid for calculating Earth parameters
    TIDE_SYSTEM: str, default 'tide_free'
        Permanent tide system for the output solid Earth tide

            - ``'tide_free'``: no permanent direct and indirect tidal potentials
            - ``'mean_tide'``: permanent tidal potentials (direct and indirect)
    EPHEMERIDES: str, default 'approximate'
        Ephemerides for calculating Earth parameters

            - ``'approximate'``: approximate lunisolar parameters
            - ``'JPL'``: computed from JPL ephmerides kernel

    Returns
    -------
    tide_se: np.ndarray
        solid earth tide at coordinates and time in meters
    """

    # validate input arguments
    assert TIME in ('GPS', 'LORAN', 'TAI', 'UTC', 'datetime')
    assert TIDE_SYSTEM in ('mean_tide', 'tide_free')
    assert EPHEMERIDES in ('approximate', 'JPL')
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
        crs1 = pyproj.CRS.from_epsg(int(EPSG))
    except (ValueError,pyproj.exceptions.CRSError):
        # Projection SRS string
        crs1 = pyproj.CRS.from_string(EPSG)
    # output coordinate reference system
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon, lat = transformer.transform(x.flatten(), y.flatten())

    # assert delta time is an array
    delta_time = np.atleast_1d(delta_time)
    # convert delta times or datetimes objects to timescale
    if (TIME.lower() == 'datetime'):
        timescale = pyTMD.time.timescale().from_datetime(
            delta_time.flatten())
    else:
        timescale = pyTMD.time.timescale().from_deltatime(delta_time,
            epoch=EPOCH, standard=TIME)
    # convert tide times to dynamical time
    tide_time = timescale.tide + timescale.tt_ut1
    # number of time points
    nt = len(timescale)

    # earth and physical parameters for ellipsoid
    units = pyTMD.constants(ELLIPSOID)

    # convert input coordinates to cartesian
    X, Y, Z = pyTMD.spatial.to_cartesian(lon, lat, h=h,
        a_axis=units.a_axis, flat=units.flat)
    # compute ephemerides for lunisolar coordinates
    if (EPHEMERIDES.lower() == 'approximate'):
        # get low-resolution solar and lunar ephemerides
        SX, SY, SZ = pyTMD.astro.solar_ecef(timescale.MJD)
        LX, LY, LZ = pyTMD.astro.lunar_ecef(timescale.MJD)
    elif (EPHEMERIDES.upper() == 'JPL'):
        # compute solar and lunar ephemerides from JPL kernel
        SX, SY, SZ = pyTMD.astro.solar_ephemerides(timescale.MJD)
        LX, LY, LZ = pyTMD.astro.lunar_ephemerides(timescale.MJD)

    # calculate radial displacement at time
    if (TYPE == 'grid'):
        ny,nx = np.shape(x)
        tide_se = np.zeros((ny,nx,nt))
        # convert coordinates to column arrays
        XYZ = np.c_[X, Y, Z]
        for i in range(nt):
            # reshape time to match spatial
            t = tide_time[i] + np.ones((ny*nx))
            # convert coordinates to column arrays
            SXYZ = np.repeat(np.c_[SX[i], SY[i], SZ[i]], ny*nx, axis=0)
            LXYZ = np.repeat(np.c_[LX[i], LY[i], LZ[i]], ny*nx, axis=0)
            # predict solid earth tides (cartesian)
            dxi = pyTMD.predict.solid_earth_tide(t,
                XYZ, SXYZ, LXYZ, a_axis=units.a_axis,
                tide_system=TIDE_SYSTEM)
            # calculate radial component of solid earth tides
            dln,dlt,drad = pyTMD.spatial.to_geodetic(
                X + dxi[:,0], Y + dxi[:,1], Z + dxi[:,2],
                a_axis=units.a_axis, flat=units.flat)
            # remove effects of original topography
            # (if added when computing cartesian coordinates)
            tide_se[:,:,i] = np.reshape(drad - h, (ny,nx))
    elif (TYPE == 'drift'):
        # convert coordinates to column arrays
        XYZ = np.c_[X, Y, Z]
        SXYZ = np.c_[SX, SY, SZ]
        LXYZ = np.c_[LX, LY, LZ]
        # predict solid earth tides (cartesian)
        dxi = pyTMD.predict.solid_earth_tide(tide_time,
            XYZ, SXYZ, LXYZ, a_axis=units.a_axis,
            tide_system=TIDE_SYSTEM)
        # calculate radial component of solid earth tides
        dln,dlt,drad = pyTMD.spatial.to_geodetic(
            X + dxi[:,0], Y + dxi[:,1], Z + dxi[:,2],
            a_axis=units.a_axis, flat=units.flat)
        # remove effects of original topography
        # (if added when computing cartesian coordinates)
        tide_se = drad - h
    elif (TYPE == 'time series'):
        nstation = len(x)
        tide_se = np.zeros((nstation,nt))
        # convert coordinates to column arrays
        SXYZ = np.c_[SX, SY, SZ]
        LXYZ = np.c_[LX, LY, LZ]
        for s in range(nstation):
            # convert coordinates to column arrays
            XYZ = np.repeat(np.c_[X[s], Y[s], Z[s]], nt, axis=0)
            # predict solid earth tides (cartesian)
            dxi = pyTMD.predict.solid_earth_tide(tide_time,
                XYZ, SXYZ, LXYZ, a_axis=units.a_axis,
                tide_system=TIDE_SYSTEM)
            # calculate radial component of solid earth tides
            dln,dlt,drad = pyTMD.spatial.to_geodetic(
                X + dxi[:,0], Y + dxi[:,1], Z + dxi[:,2],
                a_axis=units.a_axis, flat=units.flat)
            # remove effects of original topography
            # (if added when computing cartesian coordinates)
            tide_se[s,:] = drad - h

    # return the solid earth tide corrections
    return tide_se
