#!/usr/bin/env python
u"""
compute_OPT_displacements.py
Written by Tyler Sutterley (05/2023)
Calculates radial ocean pole load tide displacements for an input file
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
    -E X, --ellipsoid X: Ellipsoid for calculating load pole tide parameters
    -c X, --convention X: IERS mean or secular pole convention
        2003
        2010
        2015
        2018
    -I X, --interpolate X: Interpolation method
        spline
        linear
        nearest
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
    eop.py: utilities for calculating Earth Orientation Parameters (EOP)
    io/ocean_pole_tide.py: read ocean pole load tide map from IERS

REFERENCES:
    S. Desai, "Observing the pole tide with satellite altimetry", Journal of
        Geophysical Research: Oceans, 107(C11), 2002. doi: 10.1029/2001JC001224
    S. Desai, J. Wahr and B. Beckley "Revisiting the pole tide for and from
        satellite altimetry", Journal of Geodesy, 89(12), p1233-1243, 2015.
        doi: 10.1007/s00190-015-0848-7

UPDATE HISTORY:
    Updated 05/2023: use timescale class for time conversion operations
        use defaults from eop module for pole tide and EOP files
    Updated 04/2023: check if datetime before converting to seconds
        using pathlib to define and expand paths
    Updated 03/2023: added option for changing the IERS mean pole convention
    Updated 02/2023: added functionality for time series type
    Updated 01/2023: added default field mapping for reading from netCDF4/HDF5
        added data type keyword for netCDF4 output
    Updated 12/2022: single implicit import of pyTMD tools
        use constants class for ellipsoidal parameters
    Updated 11/2022: place some imports within try/except statements
        use f-strings for formatting verbose or ascii output
    Updated 10/2022: added delimiter option and datetime parsing for ascii files
    Updated 04/2022: use longcomplex data format to be windows compliant
        use argparse descriptions within sphinx documentation
    Updated 01/2022: added option for changing the time standard
    Updated 11/2021: add function for attempting to extract projection
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: can use prefix files to define command line arguments
    Updated 04/2021: fix arguments to program (using internal EOP files)
        Add missing delta time argument for files without time variables
        Thanks to Karen Alley for pointing these out!
    Updated 03/2021: use cartesian coordinate conversion routine in spatial
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 12/2020: merged time conversion routines into module
    Updated 11/2020: use internal mean pole and finals EOP files
        added options to read from and write to geotiff image files
    Updated 10/2020: using argparse to set command line parameters
    Updated 09/2020: can use HDF5 and netCDF4 as inputs and outputs
    Updated 08/2020: replaced griddata with scipy regular grid interpolators
    Updated 10/2017: use mean pole coordinates from calc_mean_iers_pole.py
    Written 10/2017 for public release
"""
from __future__ import print_function

import sys
import logging
import pathlib
import argparse
import numpy as np
import scipy.interpolate
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

# PURPOSE: compute the ocean pole load tide radial displacements following
# IERS conventions (2010) and using data from Desai (2002)
def compute_OPT_displacements(input_file, output_file,
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
    CONVENTION='2018',
    METHOD='spline',
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
    # ocean pole tides
    attrib['tide_oc_pole'] = {}
    attrib['tide_oc_pole']['long_name'] = 'Ocean_Pole_Tide'
    attrib['tide_oc_pole']['description'] = ('Ocean_pole_tide_radial_'
        'displacements_time_due_to_polar_motion')
    attrib['tide_oc_pole']['reference'] = ('ftp://tai.bipm.org/iers/conv2010/'
        'chapter7/opoleloadcoefcmcor.txt.gz')
    attrib['tide_oc_pole']['units'] = 'meters'
    attrib['tide_oc_pole']['_FillValue'] = fill_value
    # Modified Julian Days
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'days since 1858-11-17T00:00:00'
    attrib['time']['description'] = 'Modified Julian Days'
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

    # convert dynamic time to Modified Julian Days (MJD)
    MJD = timescale.tt - 2400000.5
    # convert Julian days to calendar dates
    Y,M,D,h,m,s = pyTMD.time.convert_julian(timescale.tt, format='tuple')
    # calculate time in year-decimal format
    time_decimal = pyTMD.time.convert_calendar_decimal(Y,M,day=D,
        hour=h,minute=m,second=s)
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
    h = dinput['data'].flatten() if ('data' in dinput.keys()) else 0.0
    # convert from geodetic latitude to geocentric latitude
    # calculate X, Y and Z from geodetic latitude and longitude
    X,Y,Z = pyTMD.spatial.to_cartesian(lon.flatten(), lat.flatten(), h=h,
        a_axis=units.a_axis, flat=units.flat)
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
        UR = np.zeros((len(latitude_geocentric)),dtype=np.longcomplex)
        UR.real = f1.ev(lon.flatten(), latitude_geocentric)
        UR.imag = f2.ev(lon.flatten(), latitude_geocentric)
    else:
        # use scipy regular grid to interpolate values for a given method
        r1 = scipy.interpolate.RegularGridInterpolator((ilon, ilat[::-1]),
            iur[:,::-1], method=METHOD)
        UR = r1.__call__(np.c_[lon.flatten(), latitude_geocentric])

    # calculate radial displacement at time
    if (TYPE == 'grid'):
        Urad = np.ma.zeros((ny,nx,nt),fill_value=fill_value)
        Urad.mask = np.zeros((ny,nx,nt),dtype=bool)
        for i in range(nt):
            URAD = K*atr*np.real((mx[i]*gamma.real + my[i]*gamma.imag)*UR.real +
                (my[i]*gamma.real - mx[i]*gamma.imag)*UR.imag)
            # reform grid
            Urad.data[:,:,i] = np.reshape(URAD, (ny,nx))
            Urad.mask[:,:,i] = np.isnan(Urad.data[:,:,i])
    elif (TYPE == 'drift'):
        Urad = np.ma.zeros((nt),fill_value=fill_value)
        Urad.data[:] = K*atr*np.real((mx*gamma.real + my*gamma.imag)*UR.real +
            (my*gamma.real - mx*gamma.imag)*UR.imag)
        Urad.mask = np.isnan(Urad.data)
    elif (TYPE == 'time series'):
        Urad = np.ma.zeros((nstation,nt),fill_value=fill_value)
        Urad.mask = np.zeros((nstation,nt),dtype=bool)
        for s in range(nstation):
            URAD = K*atr*np.real((mx*gamma.real + my*gamma.imag)*UR.real[s] +
                (my*gamma.real - mx*gamma.imag)*UR.imag[s])
            Urad.data[s,:] = np.copy(URAD)
            Urad.mask[s,:] = np.isnan(Urad.data[s,:])
    # replace invalid data with fill values
    Urad.data[Urad.mask] = Urad.fill_value

    # output to file
    output = dict(time=timescale.MJD, lon=lon, lat=lat, tide_oc_pole=Urad)
    if (FORMAT == 'csv'):
        pyTMD.spatial.to_ascii(output, attrib, output_file,
            delimiter=DELIMITER, header=False,
            columns=['time','lat','lon','tide_oc_pole'])
    elif (FORMAT == 'netCDF4'):
        pyTMD.spatial.to_netCDF4(output, attrib, output_file, data_type=TYPE)
    elif (FORMAT == 'HDF5'):
        pyTMD.spatial.to_HDF5(output, attrib, output_file)
    elif (FORMAT == 'geotiff'):
        pyTMD.spatial.to_geotiff(output, attrib, output_file,
            varname='tide_oc_pole')
    # change the permissions level to MODE
    output_file.chmod(mode=MODE)

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Calculates radial ocean pole load tide displacements for
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
    # ellipsoid for calculating load pole tide parameters
    parser.add_argument('--ellipsoid','-E',
        type=str, choices=pyTMD._ellipsoids, default='WGS84',
        help='Ellipsoid for calculating load pole tide parameters')
    # Earth orientation parameters
    parser.add_argument('--convention','-c',
        type=str, choices=pyTMD.eop._conventions, default='2018',
        help='IERS mean or secular pole convention')
    # interpolation method
    parser.add_argument('--interpolate','-I',
        metavar='METHOD', type=str, default='spline',
        choices=('spline','linear','nearest'),
        help='Spatial interpolation method')
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
        vars = (args.infile.stem,'ocean_pole_tide',args.infile.suffix)
        args.outfile = args.infile.with_name('{0}_{1}{2}'.format(*vars))

    # run ocean pole tide program for input file
    compute_OPT_displacements(args.infile, args.outfile,
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
        CONVENTION=args.convention,
        METHOD=args.interpolate,
        VERBOSE=args.verbose,
        MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
