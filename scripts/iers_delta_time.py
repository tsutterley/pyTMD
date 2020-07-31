#!/usr/bin/env python
u"""
iers_delta_time.py
Written by Tyler Sutterley (07/2020)

Connects to the IERS ftp server to download Bulletin-A files
    https://datacenter.iers.org/productMetadata.php?id=6
    ftp://ftp.iers.org/products/eop/rapid/bulletina/
Reads the IERS Bulletin-A files and calculates the daily delta times
Delta times are the difference between universal time and dynamical time

COMMAND LINE OPTIONS:
    -D X, --directory=X: working data directory
    -V, --verbose: Output information about each read file
    -M X, --mode=X: Permission mode of output file

OUTPUTS:
    iers_deltat.data: daily deltat file from IERS Bulletins

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    convert_julian.py: returns the calendar date and time given a Julian date
    count_leap_seconds.py: determines the number of leap seconds for a GPS time

UPDATE HISTORY:
    Written 07/2020
"""
from __future__ import print_function

import sys
import os
import re
import io
import getopt
import ftplib
import posixpath
import numpy as np
from pyTMD.convert_julian import convert_julian
from pyTMD.count_leap_seconds import count_leap_seconds

#-- PURPOSE: check internet connection
def check_connection():
    #-- attempt to connect to ftp host for IERS bulletins
    try:
        f = ftplib.FTP('ftp.iers.org')
        f.login()
        f.voidcmd("NOOP")
    except IOError:
        raise RuntimeError('Check internet connection')
    else:
        return True

#-- PURPOSE: connects to IERS servers and finds Bulletin-A files
def iers_delta_time(tide_dir, VERBOSE=False, MODE=0o775):
    #-- connect to ftp host for IERS bulletins
    ftp = ftplib.FTP('ftp.iers.org')
    ftp.login()
    #-- remote path
    remote_dir = posixpath.join('products','eop','rapid','bulletina')
    #-- regular expression pattern for finding files
    rx = re.compile('bulletina-(.*?)-(\d+).txt$',re.VERBOSE)
    #-- open output file
    fid = open(os.path.join(tide_dir,'iers_deltat.data'),'w')
    file_format = ' {0:4.0f} {1:2.0f} {2:2.0f} {3:7.4f}'
    #-- for each subdirectory
    for sub in sorted(ftp.nlst(remote_dir)):
        #-- find Bulletin-A files in ftp subdirectory
        files = [rx.search(f).group(0) for f in ftp.nlst(sub) if rx.search(f)]
        #-- for each Bulletin-A file
        for f in sorted(files):
            print(f) if VERBOSE else None
            #-- copy remote file contents to BytesIO object
            i = io.BytesIO()
            ftp.retrbinary('RETR {0}'.format(posixpath.join(sub,f)),i.write)
            i.seek(0)
            #-- read Bulletin-A file from BytesIO object
            YY,MM,DD,DELTAT = read_iers_bulletin_a(i)
            #-- print delta time for week to output file
            for Y,M,D,T in zip(YY,MM,DD,DELTAT):
                print(file_format.format(Y,M,D,T),file=fid)
            #-- close the bytesIO object
            i.close()
    #-- close the output file
    fid.close()
    #-- change the permissions mode
    os.chmod(os.path.join(tide_dir,'iers_deltat.data'),MODE)

#-- PURPOSE: reads IERS Bulletin-A and calculates the delta times
def read_iers_bulletin_a(fileID):
    """
    Read IERS Bulletin-A and calculate delta times (TT - UT1)

    Arguments
    ---------
    fileID: file object for Bulletin-A file

    Returns
    -------
    YEAR: calendar year
    MONTH: calendar month
    DAY: day of the month
    DELTAT: difference between universal time and dynamical time
    """
    #-- read contents from input file object
    file_contents = fileID.read().decode('utf-8').splitlines()

    #-- parse header text to find time offsets
    #-- TT-TAI
    TT_TAI = 0
    #-- TAI-UTC
    TAI_UTC = 0
    #-- counts the number of lines in the header
    count = 0
    HEADER = False
    #-- Reading over header text
    while not HEADER:
        #-- file line at count
        l = file_contents[count]
        #-- check if line contains time offsets
        if re.search('TT\s\=\sTAI',l):
            TT_TAI = np.float(re.findall('(\d+\.\d+)',l).pop())
        if re.search('TAI-UTC',l):
            TAI_UTC = np.float(re.findall('=\s(\d+\.\d+)',l).pop())
        #-- find line to set HEADER flag to True
        HEADER = bool(re.search('COMBINED\sEARTH\sORIENTATION\sPARAMETERS:',l))
        #-- add 1 to counter
        count += 1

    #-- convert variables to numpy arrays
    MJD = np.zeros((7))
    UT1_UTC = np.zeros((7))
    valid = 0
    #-- for each day in the week
    for i in range(7):
        try:
            #-- split numerical instances from data line
            line_contents = file_contents[count+i+4].split()
            #-- years are not always complete in the bulletin file
            #-- Modified Julian Day (days since 1858-11-17T00:00:00)
            MJD[i] = np.float(line_contents[3])
            #-- difference between UT1 and UTC times
            UT1_UTC[i] = np.float(line_contents[8])
        except IndexError:
            pass
        else:
            valid += 1

    #-- calculate components for delta time
    #-- TAI time is ahead of GPS by 19 seconds
    TAI_GPS = 19.0
    #-- calculate calendar dates from Modified Julian days
    Y,M,D,h,m,s = convert_julian(MJD[:valid]+2400000.5,FORMAT='tuple')
    #-- number of leap seconds between GPS and UTC
    #-- this finds the daily correction for weeks with leap seconds
    GPS_UTC = calc_GPS_to_UTC(Y,M,D,h,m,s+TAI_UTC-TAI_GPS)
    #-- calculate delta time (TT - UT1) -->
    #-- (TT-TAI) + (TAI-GPS) + (GPS-UTC) - (UT1-UTC)
    DELTAT = TT_TAI + TAI_GPS + GPS_UTC - UT1_UTC[:valid]

    #-- return dates and delta times
    return (Y,M,D,DELTAT)

#-- PURPOSE: calculate the number of leap seconds between GPS time
#-- (seconds since Jan 6, 1980 00:00:00) and UTC
def calc_GPS_to_UTC(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND):
    GPS = 367.*YEAR - np.floor(7.*(YEAR + np.floor((MONTH+9.)/12.))/4.) - \
        np.floor(3.*(np.floor((YEAR + (MONTH - 9.)/7.)/100.) + 1.)/4.) + \
        np.floor(275.*MONTH/9.) + DAY + 1721028.5 - 2444244.5
    GPS_Time = GPS*86400.0 + HOUR*3600.0 + MINUTE*60.0 + SECOND
    return count_leap_seconds(GPS_Time)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' -V, --verbose\t\tOutput information about each read file')
    print(' -M X, --mode=X\t\tPermission mode of output file\n')

#-- Main program that calls iers_delta_time()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','verbose','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:], 'hD:VM:', long_options)

    #-- directory with tide data
    tide_dir = os.getcwd()
    #-- verbosity settings
    VERBOSE = False
    #-- permissions mode of the local files (number in octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-D","--directory"):
            tide_dir = os.path.expanduser(arg)
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- check internet connection before attempting to run program
    if check_connection():
        iers_delta_time(tide_dir, VERBOSE=VERBOSE, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
