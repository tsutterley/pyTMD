#!/usr/bin/env python
u"""
usap_cats_tides.py
Written by Tyler Sutterley (07/2020)
Download Circum-Antarctic Tidal Simulations from the US Antarctic Program
CATS2008: https://www.usap-dc.org/view/dataset/601235

CALLING SEQUENCE:
    python usap_cats_tides.py --tide=CATS2008

COMMAND LINE OPTIONS:
    --help: list the command line options
    --directory=X: working data directory
    --tide=X: tide model to download
        CATS2008
    -M X, --mode=X: Local permissions mode of the files downloaded

PYTHON DEPENDENCIES:
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    utilities: download and management utilities for syncing files

UPDATE HISTORY:
    Written 08/2020
"""
from __future__ import print_function

import sys
import os
import re
import time
import getopt
import zipfile
import posixpath
import pyTMD.utilities

#-- PURPOSE: Download Circum-Antarctic Tidal Simulations from USAP
def usap_cats_tides(MODEL,DIRECTORY=None,MODE=0o775):
    #-- remote subdirectories for each model
    REMOTE = {}
    REMOTE['CATS2008'] = ['601235','2019-12-19T23:26:43.6Z',
        'CATS2008.zip?dataset_id=601235']
    #-- local subdirectory for each model
    LOCAL = {}
    LOCAL['CATS2008'] = 'CATS2008'
    #-- recursively create directories if non-existent
    if not os.access(os.path.join(DIRECTORY,LOCAL[MODEL]), os.F_OK):
        os.makedirs(os.path.join(DIRECTORY,LOCAL[MODEL]), MODE)

    #-- download CATS2008 zip file and read as virtual file object
    HOST = ['https://www.usap-dc.org','dataset','usap-dc',*REMOTE[MODEL]]
    #-- download zipfile from host
    zfile = zipfile.ZipFile(pyTMD.utilities.from_http(HOST))
    print('{0} -->\n'.format(posixpath.join(*HOST)))
    #-- extract each member
    for m in zfile.filelist:
        #-- strip directories from member filename
        m.filename = posixpath.basename(m.filename)
        print('\t{0}\n'.format(os.path.join(DIRECTORY,LOCAL[MODEL],m.filename)))
        #-- extract file
        zfile.extract(m, path=os.path.join(DIRECTORY,LOCAL[MODEL]))
        #-- change permissions mode
        os.chmod(os.path.join(DIRECTORY,LOCAL[MODEL],m.filename), MODE)
    #-- close the zipfile object
    zfile.close()

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' --tide=X\t\tCircum-Antarctic Tidal Simulation to download')
    print('\tCATS2008')
    print(' -M X, --mode=X\t\tPermission mode of files downloaded\n')

#-- Main program that calls usap_cats_tides()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','tide=','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hD:M:',long_options)

    #-- command line parameters
    DIRECTORY = os.getcwd()
    MODELS = ['CATS2008']
    #-- permissions mode of the local directories and files (number in octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-D","--directory"):
            DIRECTORY = os.path.expanduser(arg)
        elif opt in ("--tide",):
            MODELS = arg.upper().split(',')
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- check internet connection before attempting to run program
    if pyTMD.utilities.check_connection('https://www.usap-dc.org'):
        for m in MODELS:
            usap_cats_tides(m,DIRECTORY=DIRECTORY,MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
