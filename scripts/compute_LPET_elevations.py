#!/usr/bin/env python
u"""
compute_LPET_elevations.py
Written by Tyler Sutterley (08/2020)
Calculates long-period equilibrium tidal elevations for an input file

INPUTS:
    csv file with columns for spatial and temporal coordinates
    HDF5 file with variables for spatial and temporal coordinates
    netCDF4 file with variables for spatial and temporal coordinates

COMMAND LINE OPTIONS:
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
    calc_delta_time.py: calculates difference between universal and dynamic time
    compute_equilibrium_tide.py: calculates long-period equilibrium ocean tides

UPDATE HISTORY:
    Written 09/2020
"""
from __future__ import print_function

import sys
import os
import getopt
import pyproj
import numpy as np
import pyTMD.time
import pyTMD.spatial
from pyTMD.utilities import get_data_path
from pyTMD.convert_julian import convert_julian
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.compute_equilibrium_tide import compute_equilibrium_tide

#-- PURPOSE: read csv, netCDF or HDF5 data
#-- compute long-period equilibrium tides at points and times
def compute_LPET_elevations(input_file, output_file,
    FORMAT='csv', VARIABLES=['time','lat','lon','data'],
    TIME_UNITS='days since 1858-11-17T00:00:00', PROJECTION='4326',
    VERBOSE=False, MODE=0o775):

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
    #-- long-period equilibrium tides
    attrib['tide_lpe'] = {}
    attrib['tide_lpe']['long_name'] = 'Equilibrium_Tide'
    attrib['tide_lpe']['description'] = ('Long-period_equilibrium_tidal_'
        'elevation_from_the_summation_of_fifteen_tidal_spectral_lines')
    attrib['tide_lpe']['reference'] = ('https://doi.org/10.1111/'
        'j.1365-246X.1973.tb03420.x')
    attrib['tide_lpe']['units'] = 'meters'
    #-- time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'days since 1992-01-01T00:00:00'
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
    #-- convert time from units to days since 1992-01-01T00:00:00
    tide_time = pyTMD.time.convert_delta_time(to_secs*dinput['time'].flatten(),
        epoch1=epoch1, epoch2=(1992,1,1,0,0,0), scale=1.0/86400.0)
    #-- interpolate delta times from calendar dates to tide time
    delta_file = get_data_path(['data','merged_deltat.data'])
    deltat = calc_delta_time(delta_file, tide_time)

    #-- predict long-period equilibrium tides at time
    tide_lpe = compute_equilibrium_tide(tide_time + deltat, lat)

    #-- output to file
    output = dict(time=tide_time,lon=lon,lat=lat,tide_lpe=tide_lpe)
    if (FORMAT == 'csv'):
        pyTMD.spatial.to_ascii(output, attrib, output_file, delimiter=',',
            columns=['time','lat','lon','tide_lpe'], verbose=VERBOSE)
    elif (FORMAT == 'netCDF4'):
        pyTMD.spatial.to_netCDF4(output, attrib, output_file, verbose=VERBOSE)
    elif (FORMAT == 'HDF5'):
        pyTMD.spatial.to_HDF5(output, attrib, output_file, verbose=VERBOSE)
    #-- change the permissions level to MODE
    os.chmod(output_file, MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' --format=X\t\tInput and output data format')
    print('\tcsv\n\tnetCDF4\n\tHDF5')
    print(' --variables=X\t\tVariable names of data in input file')
    print(' --epoch=X\t\tReference epoch of input time')
    print('\tin form "time-units since yyyy-mm-dd hh:mm:ss"')
    print(' --projection=X\t\tSpatial projection as EPSG code or PROJ4 string')
    print(' -V, --verbose\t\tVerbose output of processing run')
    print(' -M X, --mode=X\t\tPermission mode of output file\n')

#-- Main program that calls compute_LPET_elevations()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','variables=','epoch=','projection=','format=',
        'verbose','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:], 'hVM:', long_options)

    #-- command line options
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
        args = (fileBasename,'lpe_tide',fileExtension)
        output_file = '{0}_{1}{2}'.format(*args)

    #-- run long period equilibrium tide program for input file
    compute_LPET_elevations(input_file, output_file, FORMAT=FORMAT,
        VARIABLES=VARIABLES, TIME_UNITS=TIME_UNITS, PROJECTION=PROJECTION,
        VERBOSE=VERBOSE, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
