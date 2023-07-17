#!/usr/bin/env python
u"""
compute_LPET_elevations.py
Written by Tyler Sutterley (04/2023)
Calculates long-period equilibrium tidal elevations for an input file

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
    predict.py: calculates long-period equilibrium ocean tides

UPDATE HISTORY:
    Updated 05/2023: use timescale class for time conversion operations
    Updated 04/2023: check if datetime before converting to seconds
        using pathlib to define and expand paths
    Updated 02/2023: added functionality for time series type
    Updated 01/2023: added default field mapping for reading from netCDF4/HDF5
        added data type keyword for netCDF4 output
    Updated 12/2022: single implicit import of pyTMD tools
    Updated 11/2022: place some imports within try/except statements
        use f-strings for formatting verbose or ascii output
    Updated 10/2022: added delimiter option and datetime parsing for ascii files
    Updated 04/2022: use argparse descriptions within documentation
    Updated 01/2022: added option for changing the time standard
    Updated 11/2021: add function for attempting to extract projection
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: can use prefix files to define command line arguments
    Updated 11/2020: added options to read from and write to geotiff image files
    Updated 10/2020: using argparse to set command line parameters
    Written 09/2020
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

# PURPOSE: read csv, netCDF or HDF5 data
# compute long-period equilibrium tides at points and times
def compute_LPET_elevations(input_file, output_file,
    FORMAT='csv',
    VARIABLES=['time','lat','lon','data'],
    HEADER=0,
    DELIMITER=',',
    TYPE='drift',
    TIME_UNITS='days since 1858-11-17T00:00:00',
    TIME_STANDARD='UTC',
    TIME=None,
    PROJECTION='4326',
    VERBOSE=False,
    MODE=0o775):

    # create logger for verbosity level
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

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
    # long-period equilibrium tides
    attrib['tide_lpe'] = {}
    attrib['tide_lpe']['long_name'] = 'Equilibrium_Tide'
    attrib['tide_lpe']['description'] = ('Long-period_equilibrium_tidal_'
        'elevation_from_the_summation_of_fifteen_tidal_spectral_lines')
    attrib['tide_lpe']['reference'] = ('https://doi.org/10.1111/'
        'j.1365-246X.1973.tb03420.x')
    attrib['tide_lpe']['units'] = 'meters'
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
    # flatten latitudes
    phi = lat.flatten()

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
    # number of time points
    nt = len(timescale)
    # convert tide times to dynamical time
    tide_time = timescale.tide + timescale.tt_ut1

    # predict long-period equilibrium tides at time
    if (TYPE == 'grid'):
        tide_lpe = np.zeros((ny,nx,nt))
        for i in range(nt):
            lpet = pyTMD.predict.equilibrium_tide(tide_time[i], phi)
            tide_lpe[:,:,i] = np.reshape(lpet,(ny,nx))
    elif (TYPE == 'drift'):
        tide_lpe = pyTMD.predict.equilibrium_tide(tide_time, phi)
    elif (TYPE == 'time series'):
        tide_lpe = np.zeros((nstation,nt))
        for s in range(nstation):
            lpet = pyTMD.predict.equilibrium_tide(tide_time, phi[s])
            tide_lpe[s,:] = np.copy(lpet)

    # output to file
    output = dict(time=timescale.tide, lon=lon, lat=lat, tide_lpe=tide_lpe)
    if (FORMAT == 'csv'):
        pyTMD.spatial.to_ascii(output, attrib, output_file,
            delimiter=DELIMITER, header=False,
            columns=['time','lat','lon','tide_lpe'])
    elif (FORMAT == 'netCDF4'):
        pyTMD.spatial.to_netCDF4(output, attrib, output_file, data_type=TYPE)
    elif (FORMAT == 'HDF5'):
        pyTMD.spatial.to_HDF5(output, attrib, output_file)
    elif (FORMAT == 'geotiff'):
        pyTMD.spatial.to_geotiff(output, attrib, output_file,
            varname='tide_lpe')
    # change the permissions level to MODE
    output_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates long-period equilibrium tidal elevations for
            an input file
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
        vars = (args.infile.stem,'lpe_tide',args.infile.suffix)
        args.outfile = args.infile.with_name('{0}_{1}{2}'.format(*vars))

    # run long period equilibrium tide program for input file
    compute_LPET_elevations(args.infile, args.outfile,
        FORMAT=args.format,
        VARIABLES=args.variables,
        HEADER=args.header,
        DELIMITER=args.delimiter,
        TYPE=args.type,
        TIME_UNITS=args.epoch,
        TIME=args.deltatime,
        TIME_STANDARD=args.standard,
        PROJECTION=args.projection,
        VERBOSE=args.verbose,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
