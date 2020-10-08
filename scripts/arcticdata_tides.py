#!/usr/bin/env python
u"""
arcticdata_tides.py
Written by Tyler Sutterley (10/2020)
Download Arctic Ocean Tide Models from the NSF ArcticData archive
AODTM-5: https://arcticdata.io/catalog/view/doi:10.18739/A2901ZG3N
AOTIM-5: https://arcticdata.io/catalog/view/doi:10.18739/A2S17SS80
AOTIM-5-2018: https://arcticdata.io/catalog/view/doi:10.18739/A21R6N14K

CALLING SEQUENCE:
    python arcticdata_tides.py --tide=AOTIM-5-2018

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: working data directory
    -T X, --tide X: Arctic tide model to download
        AODTM-5
        AOTIM-5
        AOTIM-5-2018
    -M X, --mode X: Local permissions mode of the files downloaded

PYTHON DEPENDENCIES:
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    utilities: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 10/2020: using argparse to set command line parameters
    Written 08/2020
"""
from __future__ import print_function

import sys
import os
import re
import time
import zipfile
import argparse
import posixpath
import pyTMD.utilities

#-- PURPOSE: Download Arctic Ocean Tide Models from the NSF ArcticData archive
def arcticdata_tides(MODEL,DIRECTORY=None,MODE=0o775):
    #-- doi for each model
    DOI = {}
    DOI['AODTM-5'] = '10.18739/A2901ZG3N'
    DOI['AOTIM-5'] = '10.18739/A2S17SS80'
    DOI['AOTIM-5-2018'] = '10.18739/A21R6N14K'
    #-- local subdirectory for each model
    LOCAL = {}
    LOCAL['AODTM-5'] = 'aodtm5_tmd'
    LOCAL['AOTIM-5'] = 'aotim5_tmd'
    LOCAL['AOTIM-5-2018'] = 'Arc5km2018'
    #-- recursively create directories if non-existent
    if not os.access(os.path.join(DIRECTORY,LOCAL[MODEL]), os.F_OK):
        os.makedirs(os.path.join(DIRECTORY,LOCAL[MODEL]), MODE)

    #-- build host url for model
    resource_map_doi = 'resource_map_doi:{0}'.format(DOI[MODEL])
    HOST = ['https://arcticdata.io','metacat','d1','mn','v2','packages',
        pyTMD.utilities.quote_plus(posixpath.join('application','bagit-097')),
        pyTMD.utilities.quote_plus(resource_map_doi)]
    #-- download zipfile from host
    zfile = zipfile.ZipFile(pyTMD.utilities.from_http(HOST))
    print('{0} -->\n'.format(posixpath.join(*HOST)))
    #-- find model files within zip file
    rx = re.compile('(grid|h[0]?|UV[0]?|Model|xy)_(.*?)',re.VERBOSE)
    members = [m for m in zfile.filelist if rx.search(m.filename)]
    #-- extract each member
    for m in members:
        #-- strip directories from member filename
        m.filename = posixpath.basename(m.filename)
        print('\t{0}\n'.format(os.path.join(DIRECTORY,LOCAL[MODEL],m.filename)))
        #-- extract file
        zfile.extract(m, path=os.path.join(DIRECTORY,LOCAL[MODEL]))
        #-- change permissions mode
        os.chmod(os.path.join(DIRECTORY,LOCAL[MODEL],m.filename), MODE)
    #-- close the zipfile object
    zfile.close()

#-- Main program that calls arcticdata_tides()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Download Arctic Ocean Tide Models from the NSF ArcticData
            archive
            """
    )
    #-- command line parameters
    #-- working data directory for location of tide models
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- Arctic Ocean tide model to download
    parser.add_argument('--tide','-T',
        metavar='TIDE', type=str, nargs='+', default=['AOTIM-5-2018'],
        choices=('AODTM-5','AOTIM-5','AOTIM-5-2018'),
        help='Arctic Ocean tide model to download')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of the files downloaded')
    args = parser.parse_args()

    #-- check internet connection before attempting to run program
    if pyTMD.utilities.check_connection('https://arcticdata.io'):
        for m in args.tide:
            arcticdata_tides(m,DIRECTORY=args.directory,MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
