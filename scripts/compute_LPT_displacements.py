#!/usr/bin/env python
u"""
compute_LPT_displacements.py
Written by Tyler Sutterley (10/2017)
Calculates radial pole load tide displacements for an input csv file
    following IERS Convention (2010) guidelines
    http://maia.usno.navy.mil/conventions/2010officialinfo.php
    http://maia.usno.navy.mil/conventions/chapter7.php

INPUTS:
    csv file with columns:
        Modified Julian Day (days since 1858-11-17 at 00:00:00)
        latitude: degrees
        longitude: degrees
        elevation (height above or below WGS84 ellipsoid)

COMMAND LINE OPTIONS:
    -D X, --directory=X: Working data directory
    -M X, --mode=X: Permission mode of output file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/

PROGRAM DEPENDENCIES:
    convert_julian.py: returns the calendar date and time given a Julian date
    convert_calendar_decimal.py: converts from calendar dates into decimal years
    iers_mean_pole.py: provides the angular coordinates of IERS Mean Pole
    read_iers_EOP.py: read daily earth orientation parameters from IERS

UPDATE HISTORY:
    Updated 10/2017: use mean pole coordinates from calc_mean_iers_pole.py
    Written 10/2017 for public release
"""
from __future__ import print_function

import sys
import os
import getopt
import numpy as np
import scipy.interpolate
from pyTMD.convert_julian import convert_julian
from pyTMD.convert_calendar_decimal import convert_calendar_decimal
from pyTMD.iers_mean_pole import iers_mean_pole
from pyTMD.read_iers_EOP import read_iers_EOP
from pyTMD.read_ocean_pole_tide import read_ocean_pole_tide

#-- PURPOSE: compute the pole load tide radial displacements following
#-- IERS conventions (2010)
def compute_LPT_displacements(tide_dir, input_file, output_file, MODE=0o775):

    #-- read input *.csv file to extract MJD, latitude, longitude and elevation
    dtype = dict(names=('MJD','lat','lon','h'),formats=('f','f','f','f'))
    dinput = np.loadtxt(input_file, delimiter=',', dtype=dtype)
    #-- convert from MJD to calendar dates, then to year-decimal
    YY,MM,DD,HH,MN,SS = convert_julian(dinput['MJD'] + 2400000.5,FORMAT='tuple')
    tdec = convert_calendar_decimal(YY,MM,DAY=DD,HOUR=HH,MINUTE=MN,SECOND=SS)

    #-- degrees to radians
    dtr = np.pi/180.0
    atr = np.pi/648000.0
    #-- earth and physical parameters (IERS and WGS84)
    GM = 3.98004418e14#-- geocentric gravitational constant [m^3/s^2]
    a_axis = 6378136.6#-- semimajor axis of the WGS84 ellipsoid [m]
    flat = 1.0/298.257223563#-- flattening of the WGS84 ellipsoid
    b_axis = (1.0 -flat)*a_axis#-- semiminor axis of the WGS84 ellipsoid [m]
    omega = 7.292115e-5#-- mean rotation rate of the Earth [radians/s]
    #-- tidal love number appropriate for the load tide
    hb2 = 0.6207
    #-- Linear eccentricity, first and second numerical eccentricity
    lin_ecc = np.sqrt((2.0*flat - flat**2)*a_axis**2)
    ecc1 = lin_ecc/a_axis
    ecc2 = lin_ecc/b_axis
    #-- m parameter [omega^2*a^2*b/(GM)]. p. 70, Eqn.(2-137)
    m = omega**2*((1 -flat)*a_axis**3)/GM
    #-- flattening components
    f_2 = -flat + (5.0/2.0)*m + (1.0/2.0)*flat**2.0 - (26.0/7.0)*flat*m + \
        (15.0/4.0)*m**2.0
    f_4 = -(1.0/2.0)*flat**2.0 + (5.0/2.0)*flat*m

    #-- convert from geodetic latitude to geocentric latitude
    #-- geodetic latitude in radians
    latitude_geodetic_rad = dinput['lat']*dtr
    #-- prime vertical radius of curvature
    N = a_axis/np.sqrt(1.0 - ecc1**2.0*np.sin(latitude_geodetic_rad)**2.0)
    #-- calculate X, Y and Z from geodetic latitude and longitude
    X = (N+dinput['h'])*np.cos(latitude_geodetic_rad)*np.cos(dinput['lon']*dtr)
    Y = (N+dinput['h'])*np.cos(latitude_geodetic_rad)*np.sin(dinput['lon']*dtr)
    Z = (N * (1.0 - ecc1**2.0) + dinput['h']) * np.sin(latitude_geodetic_rad)
    rr = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
    #-- calculate geocentric latitude and convert to degrees
    latitude_geocentric = np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/dtr

    #-- colatitude and longitude in radians
    theta = dtr*(90.0 - latitude_geocentric)
    phi = dinput['lon']*dtr

    #-- compute normal gravity at spatial location and elevation of points.
    #-- normal gravity at the equator. p. 79, Eqn.(2-186)
    gamma_a = (GM/(a_axis*b_axis)) * (1.0-(3.0/2.0)*m - (3.0/14.0)*ecc2**2.0*m)
    #-- Normal gravity. p. 80, Eqn.(2-199)
    gamma_0 = gamma_a*(1.0 + f_2*np.cos(theta)**2.0 +
        f_4*np.sin(np.pi*latitude_geocentric/180.0)**4.0)
    #-- Normal gravity at height h. p. 82, Eqn.(2-215)
    gamma_h = gamma_0*(1.0-(2.0/a_axis)*(1.0+flat+m-2.0*flat*np.cos(theta)**2.0)
        *dinput['h'] + (3.0/a_axis**2.0)*dinput['h']**2.0)

    #-- pole tide files (mean and daily)
    mean_pole_file = os.path.join(tide_dir,'mean_pole_2017-10-23.tab')
    pole_tide_file = os.path.join(tide_dir,'finals_all_2017-09-01.tab')
    #-- calculate angular coordinates of mean pole at time tdec
    mpx,mpy,fl = iers_mean_pole(mean_pole_file,tdec,'2015')
    #-- read IERS daily polar motion values
    EOP = read_iers_EOP(pole_tide_file)
    #-- interpolate daily polar motion values to t1 using cubic splines
    xSPL = scipy.interpolate.UnivariateSpline(EOP['MJD'],EOP['x'],k=3,s=0)
    ySPL = scipy.interpolate.UnivariateSpline(EOP['MJD'],EOP['y'],k=3,s=0)
    px = xSPL(dinput['MJD'])
    py = ySPL(dinput['MJD'])
    #-- calculate differentials from mean pole positions
    mx = px - mpx
    my = -(py - mpy)

    #-- calculate radial displacement at time
    dfactor = -hb2*atr*(omega**2*rr**2)/(2.0*gamma_h)
    Srad = dfactor*np.sin(2.0*theta)*(mx*np.cos(1.0*phi) + my*np.sin(1.0*phi))

    #-- output to file
    with open(output_file) as f:
        for d,lt,ln,s in zip(dinput['MJD'],dinput['lat'],dinput['lon'],Srad):
            print('{0:g},{1:g},{2:g},{3:f}'.format(d,lt,ln,s), file=f)
    #-- change the permissions level to MODE
    os.chmod(output_file, MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' -M X, --mode=X\t\tPermission mode of output file\n')

#-- Main program that calls compute_LPT_displacements()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:], 'hD:M:', long_options)

    #-- set data directory containing the pole tide files
    tide_dir = None
    #-- permissions mode of the local files (number in octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-D","--directory"):
            tide_dir = os.path.expanduser(arg)
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- enter input and output files as system argument
    if not arglist:
        raise Exception('No System Arguments Listed')

    #-- tilde-expand input and output csv files
    input_file = os.path.expanduser(arglist[0])
    output_file = os.path.expanduser(arglist[1])
    #-- set base directory from the input file if not current set
    tide_dir = os.path.dirname(input_file) if (tide_dir is None) else tide_dir

    #-- run load pole tide program for input *.csv file
    compute_LPT_displacements(tide_dir, input_file, output_file, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
