#!/usr/bin/env python
u"""
usap_cats_tides.py
Written by Tyler Sutterley (04/2023)
Download Circum-Antarctic Tidal Simulations from the US Antarctic Program
CATS2008: https://www.usap-dc.org/view/dataset/601235

NOTE: USAP now requires a captcha to download datasets

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
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 04/2023: using pathlib to define and expand paths
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 04/2022: use argparse descriptions within documentation
    Updated 10/2021: using python logging for handling verbose output
    Updated 08/2021: USAP now requires captchas for dataset downloads
    Updated 07/2021: can use prefix files to define command line arguments
    Updated 10/2020: using argparse to set command line parameters
    Written 08/2020
"""
from __future__ import print_function

import sys
import re
import time
import logging
import pathlib
import zipfile
import warnings
import argparse
import posixpath
import webbrowser
import pyTMD.utilities

# PURPOSE: Download Circum-Antarctic Tidal Simulations from USAP
def usap_cats_tides(MODEL,DIRECTORY=None,MODE=0o775):

    # create logger for verbosity level
    logger = pyTMD.utilities.build_logger(__name__,level=logging.INFO)

    # remote subdirectories for each model
    REMOTE = {}
    REMOTE['CATS2008'] = ['601235','2019-12-19T23:26:43.6Z',
        'CATS2008.zip?dataset_id=601235']
    # local subdirectory for each model
    LOCAL = {}
    LOCAL['CATS2008'] = 'CATS2008'
    # recursively create directories if non-existent
    DIRECTORY = pathlib.Path(DIRECTORY).expanduser().absolute()
    local_dir = DIRECTORY.joinpath(LOCAL[MODEL])
    local_dir.mkdir(MODE, parents=True, exist_ok=True)

    # USAP now requires a captcha to download datasets
    # use a manual download until USAP allows some sort of verification
    DATASET = {}
    DATASET['CATS2008'] = ['https://www.usap-dc.org','view','dataset','601235']
    # open USAP url in a new browser window
    webbrowser.open_new_tab(posixpath.join(*DATASET[MODEL]))
    pyTMD.utilities.file_opener(local_dir)
    return

    # download CATS2008 zip file and read as virtual file object
    HOST = ['https://www.usap-dc.org','dataset','usap-dc',*REMOTE[MODEL]]
    # download zipfile from host
    zfile = zipfile.ZipFile(pyTMD.utilities.from_http(HOST))
    logger.info('{0} -->\n'.format(posixpath.join(*HOST)))
    # extract each member
    for m in zfile.filelist:
        # strip directories from member filename
        m.filename = posixpath.basename(m.filename)
        local_file = local_dir.joinpath(m.filename)
        logger.info(f'\t{str(local_file)}\n')
        # extract file
        zfile.extract(m, path=local_dir)
        # change permissions mode
        local_file.chmod(mode=MODE)
    # close the zipfile object
    zfile.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Download Circum-Antarctic Tidal Simulations from the
            US Antarctic Program
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = pyTMD.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory for location of tide models
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # Antarctic Ocean tide model to download
    parser.add_argument('--tide','-T',
        metavar='TIDE', type=str, nargs='+', default=['CATS2008'],
        choices=('CATS2008',),
        help='Circum-Antarctic tide model to download')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permissions mode of the files downloaded')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # warn user that USAP requires a reCAPTCHA check
    warnings.filterwarnings("module")
    warnings.warn("Deprecated. USAP now requires captcha", DeprecationWarning)
    warnings.filterwarnings("ignore")
    # check internet connection before attempting to run program
    if pyTMD.utilities.check_connection('https://www.usap-dc.org'):
        for m in args.tide:
            usap_cats_tides(m, DIRECTORY=args.directory, MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
