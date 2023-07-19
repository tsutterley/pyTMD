#!/usr/bin/env python
u"""
compute_SET_displacements.py
Written by Tyler Sutterley (05/2023)
Calculates radial solid earth tide displacements for an input file
    following IERS Convention (2010) guidelines
    https://iers-conventions.obspm.fr/chapter7.php

INPUTS:
    csv file with columns for spatial and temporal coordinates
    HDF5 file with variables for spatial and temporal coordinates
    netCDF4 file with variables for spatial and temporal coordinates
    geotiff file with bands in spatial coordinates

COMMAND LINE OPTIONS:
    -F X, --format X: input and output data format
        csv (default)
        netCDF4
        HDF5
        geotiff
    -v X, --variables X: variable names of data in csv, HDF5 or netCDF4 file
        for csv files: the order of the columns within the file
        for HDF5 and netCDF4 files: time, y, x and data variable names
    -H X, --header X: number of header lines for csv files
    --delimiter X: Delimiter for csv or ascii files
    -t X, --type X: input data type
        drift: drift buoys or satellite/airborne altimetry (time per data point)
        grid: spatial grids or images (single time for all data points)
        time series: time series at a single point
    -e X, --epoch X: Reference epoch of input time (default Modified Julian Day)
        days since 1858-11-17T00:00:00
    -d X, --deltatime X: Input delta time for files without date information
        can be set to 0 to use exact calendar date from epoch
    -s X, --standard X: Input time standard for delta times or input time type
        UTC: Coordinate Universal Time
        GPS: GPS Time
        LORAN: Long Range Navigator Time
        TAI: International Atomic Time
        datetime: formatted datetime string in UTC
    -P X, --projection X: spatial projection as EPSG code or PROJ4 string
        4326: latitude and longitude coordinates on WGS84 reference ellipsoid
    -E X, --ellipsoid X: Ellipsoid for calculating astronomical parameters
    -p X, --tide-system X: Permanent tide system for output values
        tide_free: no permanent direct and indirect tidal potentials
        mean_tide: permanent tidal potentials (direct and indirect)
    -c X, --ephemerides X: method for calculating lunisolar ephemerides
        approximate: low-resolution ephemerides (default)
        JPL: computed solar and lunar ephemerides from JPL kernels
    -V, --verbose: Verbose output of processing run
    -M X, --mode X: Permission mode of output file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org/
    netCDF4: Python interface to the netCDF C library
         https://unidata.github.io/netcdf4-python/netCDF4/index.html
    gdal: Pythonic interface to the Geospatial Data Abstraction Library (GDAL)
        https://pypi.python.org/pypi/GDAL
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    spatial.py: utilities for reading and writing spatial data
    utilities.py: download and management utilities for syncing files
    predict.py: calculates solid Earth tides

REFERENCES:
    P. M. Mathews, B. A. Buffett, T. A. Herring and I. I Shapiro,
        "Forced nutations of the Earth: Influence of inner core dynamics:
        1. Theory", Journal of Geophysical Research: Solid Earth,
        96(B5), 8219--8242, (1991). doi: 10.1029/90JB01955
    P. M. Mathews, V. Dehant and J. M. Gipson,
        "Tidal station displacements", Journal of Geophysical
        Research: Solid Earth, 102(B9), 20469--20477, (1997).
        doi: 10.1029/97JB01515
    J. C. Ries, R. J. Eanes, C. K. Shum and M. M. Watkins,
        "Progress in the determination of the gravitational
        coefficient of the Earth", Geophysical Research Letters,
        19(6), 529--531, (1992). doi: 10.1029/92GL00259
    J. M. Wahr, "Body tides on an elliptical, rotating, elastic
        and oceanless Earth", Geophysical Journal of the Royal
        Astronomical Society, 64(3), 677--703, (1981).
        doi: 10.1111/j.1365-246X.1981.tb02690.x

UPDATE HISTORY:
    Updated 05/2023: use timescale class for time conversion operations
        add option for using higher resolution ephemerides from JPL
    Updated 04/2023: using pathlib to define and expand paths
    Written 04/2023
"""
from __future__ import print_function

import sys
import logging
import pathlib
import argparse
import numpy as np
import pyTMD

# attempt imports
try:
    import pyproj
except (ImportError, ModuleNotFoundError) as exc:
    logging.critical("pyproj not available")

# PURPOSE: try to get the projection information for the input file
def get_projection(attributes, PROJECTION):
    # coordinate reference system string from file
    try:
        crs = pyproj.CRS.from_string(attributes['projection'])
    except (ValueError,KeyError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    # EPSG projection code
    try:
        crs = pyproj.CRS.from_epsg(int(PROJECTION))
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    # coordinate reference system string
    try:
        crs = pyproj.CRS.from_string(PROJECTION)
    except (ValueError,pyproj.exceptions.CRSError):
        pass
    else:
        return crs
    # no projection can be made
    raise pyproj.exceptions.CRSError

# PURPOSE: compute the solid earth tide radial displacements following
# IERS conventions (2010)
def compute_SET_displacements(input_file, output_file,
    FORMAT='csv',
    VARIABLES=['time','lat','lon','data'],
    HEADER=0,
    DELIMITER=',',
    TYPE='drift',
    TIME_UNITS='days since 1858-11-17T00:00:00',
    TIME=None,
    TIME_STANDARD='UTC',
    PROJECTION='4326',
    ELLIPSOID='WGS84',
    TIDE_SYSTEM='tide_free',
    EPHEMERIDES='approximate',
    VERBOSE=False,
    MODE=0o775):

    # create logger for verbosity level
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    # invalid value
    fill_value = -9999.0
    # output netCDF4 and HDF5 file attributes
    # will be added to YAML header in csv files
    attrib = {}
    # latitude
    attrib['lat'] = {}
    attrib['lat']['long_name'] = 'Latitude'
    attrib['lat']['units'] = 'Degrees_North'
    # longitude
    attrib['lon'] = {}
    attrib['lon']['long_name'] = 'Longitude'
    attrib['lon']['units'] = 'Degrees_East'
    # solid earth tides
    attrib['tide_earth'] = {}
    attrib['tide_earth']['long_name'] = 'Solid_Earth_Tide'
    attrib['tide_earth']['description'] = ('Solid_earth_tides_in_the_'
        f'{TIDE_SYSTEM}_system')
    attrib['tide_earth']['reference'] = 'https://doi.org/10.1029/97JB01515'
    attrib['tide_earth']['units'] = 'meters'
    attrib['tide_earth']['_FillValue'] = fill_value
    # time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'days since 1992-01-01T00:00:00'
    attrib['time']['calendar'] = 'standard'

    # read input file to extract time, spatial coordinates and data
    if (FORMAT == 'csv'):
        parse_dates = (TIME_STANDARD.lower() == 'datetime')
        dinput = pyTMD.spatial.from_ascii(input_file, columns=VARIABLES,
            delimiter=DELIMITER, header=HEADER, parse_dates=parse_dates)
    elif (FORMAT == 'netCDF4'):
        field_mapping = pyTMD.spatial.default_field_mapping(VARIABLES)
        dinput = pyTMD.spatial.from_netCDF4(input_file,
            field_mapping=field_mapping)
    elif (FORMAT == 'HDF5'):
        field_mapping = pyTMD.spatial.default_field_mapping(VARIABLES)
        dinput = pyTMD.spatial.from_HDF5(input_file,
            field_mapping=field_mapping)
    elif (FORMAT == 'geotiff'):
        dinput = pyTMD.spatial.from_geotiff(input_file)
        # copy global geotiff attributes for projection and grid parameters
        for att_name in ['projection','wkt','spacing','extent']:
            attrib[att_name] = dinput['attributes'][att_name]
    # update time variable if entered as argument
    if TIME is not None:
        dinput['time'] = np.copy(TIME)

    # converting x,y from projection to latitude/longitude
    crs1 = get_projection(dinput['attributes'], PROJECTION)
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    if (TYPE == 'grid'):
        ny, nx = (len(dinput['y']), len(dinput['x']))
        gridx, gridy = np.meshgrid(dinput['x'], dinput['y'])
        lon, lat = transformer.transform(gridx, gridy)
    elif (TYPE == 'drift'):
        lon, lat = transformer.transform(dinput['x'], dinput['y'])
    elif (TYPE == 'time series'):
        nstation = len(dinput['y'].flatten())
        lon, lat = transformer.transform(dinput['x'], dinput['y'])

    # extract time units from netCDF4 and HDF5 attributes or from TIME_UNITS
    try:
        time_string = dinput['attributes']['time']['units']
        epoch1, to_secs = pyTMD.time.parse_date_string(time_string)
    except (TypeError, KeyError, ValueError):
        epoch1, to_secs = pyTMD.time.parse_date_string(TIME_UNITS)

    # convert delta times or datetimes objects to timescale
    if (TIME_STANDARD.lower() == 'datetime'):
        timescale = pyTMD.time.timescale().from_datetime(
            dinput['time'].flatten())
    else:
        # convert time to seconds
        delta_time = to_secs*dinput['time'].flatten()
        timescale = pyTMD.time.timescale().from_deltatime(delta_time,
            epoch=epoch1, standard=TIME_STANDARD)
    # convert tide times to dynamical time
    tide_time = timescale.tide + timescale.tt_ut1
    # number of time points
    nt = len(timescale)

    # earth and physical parameters for ellipsoid
    units = pyTMD.constants(ELLIPSOID)

    # flatten heights
    h = dinput['data'].flatten() if ('data' in dinput.keys()) else 0.0
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
            dln, dlt, drad = pyTMD.spatial.to_geodetic(
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
        dln, dlt, drad = pyTMD.spatial.to_geodetic(
            X + dxi[:,0], Y + dxi[:,1], Z + dxi[:,2],
            a_axis=units.a_axis, flat=units.flat)
        # remove effects of original topography
        # (if added when computing cartesian coordinates)
        tide_se = drad - h
    elif (TYPE == 'time series'):
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
            dln, dlt, drad = pyTMD.spatial.to_geodetic(
                X + dxi[:,0], Y + dxi[:,1], Z + dxi[:,2],
                a_axis=units.a_axis, flat=units.flat)
            # remove effects of original topography
            # (if added when computing cartesian coordinates)
            tide_se[s,:] = drad - h

    # output to file
    output = dict(time=timescale.tide, lon=lon, lat=lat, tide_se=tide_se)
    if (FORMAT == 'csv'):
        pyTMD.spatial.to_ascii(output, attrib, output_file,
            delimiter=DELIMITER, header=False,
            columns=['time','lat','lon','tide_se'])
    elif (FORMAT == 'netCDF4'):
        pyTMD.spatial.to_netCDF4(output, attrib, output_file, data_type=TYPE)
    elif (FORMAT == 'HDF5'):
        pyTMD.spatial.to_HDF5(output, attrib, output_file)
    elif (FORMAT == 'geotiff'):
        pyTMD.spatial.to_geotiff(output, attrib, output_file,
            varname='tide_se')
    # change the permissions level to MODE
    output_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates solid earth load tide displacements for
            an input file following IERS Convention (2010) guidelines
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = pyTMD.utilities.convert_arg_line_to_args
    # command line options
    # input and output file
    parser.add_argument('infile',
        type=pathlib.Path, nargs='?',
        help='Input file to run')
    parser.add_argument('outfile',
        type=pathlib.Path, nargs='?',
        help='Computed output file')
    # input and output data format
    parser.add_argument('--format','-F',
        type=str, default='csv', choices=('csv','netCDF4','HDF5','geotiff'),
        help='Input and output data format')
    # variable names (for csv names of columns)
    parser.add_argument('--variables','-v',
        type=str, nargs='+', default=['time','lat','lon','data'],
        help='Variable names of data in input file')
    # number of header lines for csv files
    parser.add_argument('--header','-H',
        type=int, default=0,
        help='Number of header lines for csv files')
    # delimiter for csv or ascii files
    parser.add_argument('--delimiter',
        type=str, default=',',
        help='Delimiter for csv or ascii files')
    # input data type
    # drift: drift buoys or satellite/airborne altimetry (time per data point)
    # grid: spatial grids or images (single time for all data points)
    # time series: station locations with multiple time values
    parser.add_argument('--type','-t',
        type=str, default='drift',
        choices=('drift','grid','time series'),
        help='Input data type')
    # time epoch (default Modified Julian Days)
    # in form "time-units since yyyy-mm-dd hh:mm:ss"
    parser.add_argument('--epoch','-e',
        type=str, default='days since 1858-11-17T00:00:00',
        help='Reference epoch of input time')
    # input delta time for files without date information
    parser.add_argument('--deltatime','-d',
        type=float, nargs='+',
        help='Input delta time for files without date variables')
    # input time standard definition
    parser.add_argument('--standard','-s',
        type=str, choices=('UTC','GPS','TAI','LORAN','datetime'), default='UTC',
        help='Input time standard for delta times')
    # spatial projection (EPSG code or PROJ4 string)
    parser.add_argument('--projection','-P',
        type=str, default='4326',
        help='Spatial projection as EPSG code or PROJ4 string')
    # ellipsoid for calculating astronomical parameters
    parser.add_argument('--ellipsoid','-E',
        type=str, choices=pyTMD._ellipsoids, default='WGS84',
        help='Ellipsoid for calculating astronomical parameters')
    # permanent tide system for output values
    parser.add_argument('--tide-system','-p',
        type=str, choices=('tide_free','mean_tide'), default='tide_free',
        help='Permanent tide system for output values')
    # method for calculating lunisolar ephemerides
    parser.add_argument('--ephemerides','-c',
        type=str, choices=('approximate','JPL'), default='approximate',
        help='Method for calculating lunisolar ephemerides')
    # verbose output of processing run
    # print information about each input and output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    # permissions mode of the local files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of output file')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # set output file from input filename if not entered
    if not args.outfile:
        vars = (args.infile.stem,'solid_earth_tide',args.infile.suffix)
        args.outfile = args.infile.with_name('{0}_{1}{2}'.format(*vars))

    # run solid earth tide program for input file
    compute_SET_displacements(args.infile, args.outfile,
        FORMAT=args.format,
        VARIABLES=args.variables,
        HEADER=args.header,
        DELIMITER=args.delimiter,
        TYPE=args.type,
        TIME_UNITS=args.epoch,
        TIME=args.deltatime,
        TIME_STANDARD=args.standard,
        PROJECTION=args.projection,
        ELLIPSOID=args.ellipsoid,
        TIDE_SYSTEM=args.tide_system,
        EPHEMERIDES=args.ephemerides,
        VERBOSE=args.verbose,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
