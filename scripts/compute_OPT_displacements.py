#!/usr/bin/env python
u"""
compute_OPT_displacements.py
Written by Tyler Sutterley (09/2020)
Calculates radial ocean pole load tide displacements for an input file
    following IERS Convention (2010) guidelines
    http://maia.usno.navy.mil/conventions/2010officialinfo.php
    http://maia.usno.navy.mil/conventions/chapter7.php

INPUTS:
    csv file with columns for spatial and temporal coordinates
    HDF5 file with variables for spatial and temporal coordinates
    netCDF4 file with variables for spatial and temporal coordinates

COMMAND LINE OPTIONS:
    -D X, --directory=X: Working data directory
    --format=X: input and output data format
        csv (default)
        netCDF4
        HDF5
    --variables=X: variable names of data in csv, HDF5 or netCDF4 file
        for csv files: the order of the columns within the file
        for HDF5 and netCDF4 files: time, y, x and data variable names
    --epoch=X: Reference epoch of input time (default Modified Julian Day)
        days since 1858-11-17T00:00:00
    --projection=X: spatial projection as EPSG code or PROJ4 string
        4326: latitude and longitude coordinates on WGS84 reference ellipsoid
    -V, --verbose: Verbose output of processing run
    -M X, --mode=X: Permission mode of output file

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
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    time.py: utilities for calculating time operations
    spatial.py: utilities for reading and writing spatial data
    utilities: download and management utilities for syncing files
    convert_julian.py: returns the calendar date and time given a Julian date
    convert_calendar_decimal.py: converts from calendar dates into decimal years
    iers_mean_pole.py: provides the angular coordinates of IERS Mean Pole
    read_iers_EOP.py: read daily earth orientation parameters from IERS
    read_ocean_pole_tide.py: read ocean pole load tide map from IERS

UPDATE HISTORY:
    Updated 09/2020: can use HDF5 and netCDF4 as inputs and outputs
    Updated 08/2020: replaced griddata with scipy regular grid interpolators
    Updated 10/2017: use mean pole coordinates from calc_mean_iers_pole.py
    Written 10/2017 for public release
"""
from __future__ import print_function

import sys
import os
import getopt
import pyproj
import numpy as np
import pyTMD.time
import pyTMD.spatial
import scipy.interpolate
from pyTMD.utilities import get_data_path
from pyTMD.convert_julian import convert_julian
from pyTMD.convert_calendar_decimal import convert_calendar_decimal
from pyTMD.iers_mean_pole import iers_mean_pole
from pyTMD.read_iers_EOP import read_iers_EOP
from pyTMD.read_ocean_pole_tide import read_ocean_pole_tide

#-- PURPOSE: compute the ocean pole load tide radial displacements following
#-- IERS conventions (2010) and using data from Desai (2002)
def compute_OPT_displacements(tide_dir, input_file, output_file,
    FORMAT='csv', VARIABLES=['time','lat','lon','data'],
    TIME_UNITS='days since 1858-11-17T00:00:00', PROJECTION='4326',
    METHOD=None, VERBOSE=False, MODE=0o775):

    #-- invalid value
    fill_value = -9999.0
    #-- output netCDF4 and HDF5 file attributes
    #-- will be added to YAML header in csv files
    attrib = {}
    #-- latitude
    attrib['lat'] = {}
    attrib['lat']['long_name'] = 'Latitude'
    attrib['lat']['units'] = 'Degrees_North'
    #-- longitude
    attrib['lon'] = {}
    attrib['lon']['long_name'] = 'Longitude'
    attrib['lon']['units'] = 'Degrees_East'
    #-- ocean pole tides
    attrib['tide_oc_pole'] = {}
    attrib['tide_oc_pole']['long_name'] = 'Ocean_Pole_Tide'
    attrib['tide_oc_pole']['description'] = ('Ocean_pole_tide_radial_'
        'displacements_time_due_to_polar_motion')
    attrib['tide_oc_pole']['reference'] = ('ftp://tai.bipm.org/iers/conv2010/'
        'chapter7/opoleloadcoefcmcor.txt.gz')
    attrib['tide_oc_pole']['units'] = 'meters'
    attrib['tide_oc_pole']['_FillValue'] = fill_value
    #-- Modified Julian Days
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'days since 1858-11-17T00:00:00'
    attrib['time']['description'] = 'Modified Julian Days'
    attrib['time']['calendar'] = 'standard'

    #-- read input file to extract time, spatial coordinates and data
    if (FORMAT == 'csv'):
        dinput = pyTMD.spatial.from_ascii(input_file, columns=VARIABLES,
            header=0, verbose=VERBOSE)
    elif (FORMAT == 'netCDF4'):
        dinput = pyTMD.spatial.from_netCDF4(input_file, timename=VARIABLES[0],
            xname=VARIABLES[2], yname=VARIABLES[1], varname=VARIABLES[3],
            verbose=VERBOSE)
    elif (FORMAT == 'HDF5'):
        dinput = pyTMD.spatial.from_HDF5(input_file, timename=VARIABLES[0],
            xname=VARIABLES[2], yname=VARIABLES[1], varname=VARIABLES[3],
            verbose=VERBOSE)

    #-- converting x,y from projection to latitude/longitude
    #-- could try to extract projection attributes from netCDF4 and HDF5 files
    try:
        crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(int(PROJECTION)))
    except (ValueError,pyproj.exceptions.CRSError):
        crs1 = pyproj.CRS.from_string(PROJECTION)
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon,lat = transformer.transform(dinput['x'].flatten(),dinput['y'].flatten())

    #-- extract time units from netCDF4 and HDF5 attributes or from TIME_UNITS
    try:
        time_string = dinput['attributes']['time']['units']
    except (TypeError, KeyError):
        epoch1,to_secs = pyTMD.time.parse_date_string(TIME_UNITS)
    else:
        epoch1,to_secs = pyTMD.time.parse_date_string(time_string)
    #-- convert dates to Modified Julian days (days since 1858-11-17T00:00:00)
    MJD = pyTMD.time.convert_delta_time(to_secs*dinput['time'].flatten(),
        epoch1=epoch1, epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
    #-- add offset to convert to Julian days and then convert to calendar dates
    Y,M,D,h,m,s = convert_julian(2400000.5 + MJD, FORMAT='tuple')
    #-- calculate time in year-decimal format
    time_decimal = convert_calendar_decimal(Y,M,DAY=D,HOUR=h,MINUTE=m,SECOND=s)
    #-- number of data points
    n_time = len(time_decimal)

    #-- degrees to radians and arcseconds to radians
    dtr = np.pi/180.0
    atr = np.pi/648000.0
    #-- earth and physical parameters (IERS and WGS84)
    G = 6.67428e-11#-- universal constant of gravitation [m^3/(kg*s^2)]
    GM = 3.986004418e14#-- geocentric gravitational constant [m^3/s^2]
    a_axis = 6378136.6#-- WGS84 equatorial radius of the Earth [m]
    flat = 1.0/298.257223563#-- flattening of the WGS84 ellipsoid
    omega = 7.292115e-5#-- mean rotation rate of the Earth [radians/s]
    rho_w = 1025.0#-- density of sea water [kg/m^3]
    ge = 9.7803278#-- mean equatorial gravitational acceleration [m/s^2]
    #-- Linear eccentricity and first numerical eccentricity
    lin_ecc = np.sqrt((2.0*flat - flat**2)*a_axis**2)
    ecc1 = lin_ecc/a_axis
    #-- tidal love number differential (1 + kl - hl) for pole tide frequencies
    gamma = 0.6870 + 0.0036j

    #-- convert from geodetic latitude to geocentric latitude
    #-- geodetic latitude in radians
    latitude_geodetic_rad = lat*dtr
    #-- prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.0*np.sin(latitude_geodetic_rad)**2.0)
    #-- calculate X, Y and Z from geodetic latitude and longitude
    X = (N+dinput['data'])*np.cos(latitude_geodetic_rad)*np.cos(lon*dtr)
    Y = (N+dinput['data'])*np.cos(latitude_geodetic_rad)*np.sin(lon*dtr)
    Z = (N * (1.0 - ecc1**2.0) + dinput['data']) * np.sin(latitude_geodetic_rad)
    #-- calculate geocentric latitude and convert to degrees
    latitude_geocentric = np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/dtr

    #-- pole tide displacement scale factor
    Hp = np.sqrt(8.0*np.pi/15.0)*(omega**2*a_axis**4)/GM
    K = 4.0*np.pi*G*rho_w*Hp*a_axis/(3.0*ge)
    K1 = 4.0*np.pi*G*rho_w*Hp*a_axis**3/(3.0*GM)

    #-- pole tide files (mean and daily)
    mean_pole_file = os.path.join(tide_dir,'mean_pole_2017-10-23.tab')
    pole_tide_file = os.path.join(tide_dir,'finals_all_2017-09-01.tab')
    #-- calculate angular coordinates of mean pole at time
    mpx,mpy,fl = iers_mean_pole(mean_pole_file,time_decimal,'2015')
    #-- read IERS daily polar motion values
    EOP = read_iers_EOP(pole_tide_file)
    #-- interpolate daily polar motion values to t1 using cubic splines
    xSPL = scipy.interpolate.UnivariateSpline(EOP['MJD'],EOP['x'],k=3,s=0)
    ySPL = scipy.interpolate.UnivariateSpline(EOP['MJD'],EOP['y'],k=3,s=0)
    px = xSPL(MJD)
    py = ySPL(MJD)
    #-- calculate differentials from mean pole positions
    mx = px - mpx
    my = -(py - mpy)

    #-- read ocean pole tide map from Desai (2002)
    ocean_pole_tide_file = get_data_path(['data','opoleloadcoefcmcor.txt.gz'])
    iur,iun,iue,ilon,ilat = read_ocean_pole_tide(ocean_pole_tide_file)
    #-- interpolate ocean pole tide map from Desai (2002)
    if (METHOD == 'spline'):
        #-- use scipy bivariate splines to interpolate to output points
        f1 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
            iur[:,::-1].real, kx=1, ky=1)
        f2 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
            iur[:,::-1].imag, kx=1, ky=1)
        UR = np.zeros((n_time),dtype=np.complex128)
        UR.real = f1.ev(lon,latitude_geocentric)
        UR.imag = f2.ev(lon,latitude_geocentric)
    else:
        #-- use scipy regular grid to interpolate values for a given method
        r1 = scipy.interpolate.RegularGridInterpolator((ilon,ilat[::-1]),
            iur[:,::-1], method=METHOD)
        UR = r1.__call__(np.c_[lon,latitude_geocentric])

    #-- calculate radial displacement at time
    Urad = np.ma.zeros((n_time),fill_value=fill_value)
    Urad.data[:] = K*atr*np.real((mx*gamma.real + my*gamma.imag)*UR.real +
        (my*gamma.real - mx*gamma.imag)*UR.imag)
    #-- replace fill values
    Urad.mask = np.isnan(Urad.data)
    Urad.data[Urad.mask] = Urad.fill_value

    #-- output to file
    output = dict(time=MJD,lon=lon,lat=lat,tide_oc_pole=Urad)
    if (FORMAT == 'csv'):
        pyTMD.spatial.to_ascii(output, attrib, output_file, delimiter=',',
            columns=['time','lat','lon','tide_oc_pole'], verbose=VERBOSE)
    elif (FORMAT == 'netCDF4'):
        pyTMD.spatial.to_netCDF4(output, attrib, output_file, verbose=VERBOSE)
    elif (FORMAT == 'HDF5'):
        pyTMD.spatial.to_HDF5(output, attrib, output_file, verbose=VERBOSE)
    #-- change the permissions level to MODE
    os.chmod(output_file, MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' --format=X\t\tInput and output data format')
    print('\tcsv\n\tnetCDF4\n\tHDF5')
    print(' --variables=X\t\tVariable names of data in input file')
    print(' --epoch=X\t\tReference epoch of input time')
    print('\tin form "time-units since yyyy-mm-dd hh:mm:ss"')
    print(' --projection=X\t\tSpatial projection as EPSG code or PROJ4 string')
    print(' -V, --verbose\t\tVerbose output of processing run')
    print(' -M X, --mode=X\t\tPermission mode of output file\n')

#-- Main program that calls compute_OPT_displacements()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','variables=','epoch=','projection=',
        'format=','verbose','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:], 'hD:VM:', long_options)

    #-- command line options
    #-- set data directory containing the pole tide files
    tide_dir = None
    #-- input and output data format
    FORMAT = 'csv'
    #-- variable names (for csv names of columns)
    VARIABLES = ['time','lat','lon','data']
    #-- time epoch (default Modified Julian Days)
    TIME_UNITS = 'days since 1858-11-17T00:00:00'
    #-- spatial projection (EPSG code or PROJ4 string)
    PROJECTION = '4326'
    #-- verbose output of processing run
    VERBOSE = False
    #-- permissions mode of the local files (number in octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-D","--directory"):
            tide_dir = os.path.expanduser(arg)
        elif opt in ("--format",):
            FORMAT = arg
        elif opt in ("--variables",):
            VARIABLES = arg.split(',')
        elif opt in ("--epoch",):
            TIME_UNITS = arg
        elif opt in ("--projection",):
            PROJECTION = arg
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- enter input and output files as system argument
    if not arglist:
        raise Exception('No System Arguments Listed')

    #-- tilde-expand input and output files
    #-- set output file from input filename if not entered
    input_file=os.path.expanduser(arglist[0])
    try:
        output_file=os.path.expanduser(arglist[1])
    except IndexError:
        fileBasename,fileExtension = os.path.splitext(input_file)
        args = (fileBasename,'ocean_pole_tide',fileExtension)
        output_file = '{0}_{1}{2}'.format(*args)
    #-- set base directory from the input file if not currently set
    tide_dir = os.path.dirname(input_file) if (tide_dir is None) else tide_dir

    #-- run load pole tide program for input file
    compute_OPT_displacements(tide_dir, input_file, output_file, FORMAT=FORMAT,
        VARIABLES=VARIABLES, TIME_UNITS=TIME_UNITS, PROJECTION=PROJECTION,
        METHOD='spline', VERBOSE=VERBOSE, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
