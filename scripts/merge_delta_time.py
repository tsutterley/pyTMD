#!/usr/bin/env python
u"""
merge_delta_time.py
Written by Tyler Sutterley (07/2020)

Connects to servers to download historic_deltat.data and deltat.data files
Reads IERS Bulletin-A produced iers_deltat.data files from iers_delta_time.py
Creates a merged delta time file combining the historic, monthly and daily files

Long-term Delta T:
    https://www.usno.navy.mil/USNO/earth-orientation/eo-products/long-term

Servers and Mirrors:
    http://maia.usno.navy.mil/ser7/
    ftp://cddis.nasa.gov/products/iers/
    ftp://cddis.gsfc.nasa.gov/pub/products/iers/

COMMAND LINE OPTIONS:
    -D X, --directory=X: working data directory
    -V, --verbose: Output information about each read file
    -M X, --mode=X: Permission mode of output file

OUTPUTS:
    historical_deltat.data: biannual historical deltat file
    deltat.data: monthly deltat file
    merged_deltat.data: merged deltat file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    convert_calendar_decimal.py: converts from calendar dates into decimal years

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
import shutil
import hashlib
import posixpath
import numpy as np
if sys.version_info[0] == 2:
    import urllib2
else:
    import urllib.request as urllib2
from pyTMD.convert_calendar_decimal import convert_calendar_decimal

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

#-- PURPOSE: Download delta time files and merge into a single
def merge_delta_time(tide_dir, VERBOSE=False, MODE=0o775):
    #-- retrieve history and monthly delta time files
    for FILE in ['deltat.data','historic_deltat.data']:
        pull_deltat_file(tide_dir,FILE,VERBOSE=VERBOSE,MODE=MODE)
    #-- read historic delta time file
    historic = np.loadtxt(os.path.join(tide_dir,'historic_deltat.data'),
        skiprows=2)
    HY = np.floor(historic[:,0])
    HM = 12.0*np.mod(historic[:,0],1.0) + 1.0
    HD = np.ones_like(historic[:,0])
    #-- read modern monthly delta time file
    monthly = np.loadtxt(os.path.join(tide_dir,'deltat.data'))
    monthly_time = convert_calendar_decimal(monthly[:,0],monthly[:,1],
        DAY=monthly[:,2])
    #-- read modern daily delta time file from IERS Bulletin A files
    daily = np.loadtxt(os.path.join(tide_dir,'iers_deltat.data'))
    daily_time = convert_calendar_decimal(daily[:,0],daily[:,1],
        DAY=daily[:,2])
    #-- write to new merged file
    fid = open(os.path.join(tide_dir,'merged_deltat.data'),'w')
    file_format = ' {0:4.0f} {1:2.0f} {2:2.0f} {3:7.4f}'
    #-- use historical values for times prior to monthly
    ind1, = np.nonzero(historic[:,0] < monthly_time[0])
    for i in ind1:
        args = (HY[i],HM[i],HD[i],historic[i,1])
        print(file_format.format(*args),file=fid)
    #-- use monthly values for times prior to daily
    ind2, = np.nonzero(monthly_time < daily_time[0])
    for i in ind2:
        args = (monthly[i,0],monthly[i,1],monthly[i,2],monthly[i,3])
        print(file_format.format(*args),file=fid)
    #-- use daily values for all times available
    for i,dt in enumerate(daily_time):
        args = (daily[i,0],daily[i,1],daily[i,2],daily[i,3])
        print(file_format.format(*args),file=fid)
    #-- close the merged file and change the permissions mode
    fid.close()
    os.chmod(os.path.join(tide_dir,'merged_deltat.data'),MODE)

#-- PURPOSE: connects to servers and downloads delta time files
def pull_deltat_file(tide_dir, FILE, VERBOSE=False, MODE=0o775):
    #-- local version of file
    LOCAL = os.path.join(tide_dir,FILE)
    #-- check if local file exists
    local_hash = ''
    if os.access(LOCAL,os.F_OK):
        #-- generate checksum hash for local file
        #-- open the local_file in binary read mode
        with open(LOCAL, 'rb') as local_buffer:
            local_hash = hashlib.md5(local_buffer.read()).hexdigest()
    #-- chunked transfer encoding size
    CHUNK = 16 * 1024

    #-- try downloading from US Naval Oceanography Portal
    REMOTE = posixpath.join('http://maia.usno.navy.mil','ser7',FILE)
    try:
        #-- Create and submit request.
        request = urllib2.Request(REMOTE)
        response = urllib2.urlopen(request)
    except:
        pass
    else:
        #-- copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO(response.read())
        remote_buffer.seek(0)
        #-- generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        #-- compare checksums
        if (local_hash != remote_hash):
            #-- print file information
            if VERBOSE:
                print('{0} -->\n\t{1}'.format(REMOTE,LOCAL))
            #-- store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with open(LOCAL, 'wb') as f:
                shutil.copyfileobj(remote_buffer, f, CHUNK)
            #-- change the permissions mode
            os.chmod(LOCAL,MODE)
        #-- leave function
        return

    #-- try downloading from NASA Crustal Dynamics Data Information System
    #-- note: anonymous ftp access will be discontinued on 2020-10-31
    #-- will have to change to an Earthdata login over https or ftp
    server = []
    server.append(['cddis.nasa.gov','pub','products','iers',FILE])
    server.append(['cddis.gsfc.nasa.gov','products','iers',FILE])
    for HOST in server:
        try:
            #-- try to connect to ftp host
            ftp = ftplib.FTP(HOST[0])
        except socket.gaierror:
            pass
        else:
            ftp.login()
            #-- remote path
            REMOTE = posixpath.join(*HOST[1:])
            #-- copy remote file contents to bytesIO object
            remote_buffer = io.BytesIO()
            ftp.retrbinary('RETR {0}'.format(REMOTE), remote_buffer.write)
            remote_buffer.seek(0)
            #-- generate checksum hash for remote file
            remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
            #-- compare checksums
            if (local_hash != remote_hash):
                #-- print file information
                if VERBOSE:
                    print('{0} -->\n\t{1}'.format(posixpath.join(*HOST),LOCAL))
                #-- store bytes to file using chunked transfer encoding
                remote_buffer.seek(0)
                with open(LOCAL, 'wb') as f:
                    shutil.copyfileobj(remote_buffer, f, CHUNK)
                #-- change the permissions mode
                os.chmod(LOCAL,MODE)
            #-- leave function
            return

    #-- return exception that no server could be connected
    raise Exception('All Server Connection Error')

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' -V, --verbose\t\tOutput information about each read file')
    print(' -M X, --mode=X\t\tPermission mode of output file\n')

#-- Main program that calls merge_delta_time()
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
        merge_delta_time(tide_dir, VERBOSE=VERBOSE, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
