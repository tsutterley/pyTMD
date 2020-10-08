#!/usr/bin/env python
u"""
usap_cats_tides.py
Written by Tyler Sutterley (10/2020)
Download Circum-Antarctic Tidal Simulations from the US Antarctic Program
CATS2008: https://www.usap-dc.org/view/dataset/601235

CALLING SEQUENCE:
    python usap_cats_tides.py --tide=CATS2008

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: working data directory
    -T X, --tide X: Circum-Antarctic tide model to download
        CATS2008
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

#-- Main program that calls usap_cats_tides()
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Download Circum-Antarctic Tidal Simulations from the
            US Antarctic Program
            """
    )
    #-- command line parameters
    #-- working data directory for location of tide models
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- Antarctic Ocean tide model to download
    parser.add_argument('--tide','-T',
        metavar='TIDE', type=str, nargs='+', default=['CATS2008'],
        choices=('CATS2008'),
        help='Circum-Antarctic tide model to download')
    #-- permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of the files downloaded')
    args = parser.parse_args()

    #-- check internet connection before attempting to run program
    if pyTMD.utilities.check_connection('https://www.usap-dc.org'):
        for m in args.tide:
            usap_cats_tides(m,DIRECTORY=args.directory,MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
