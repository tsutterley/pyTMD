#!/usr/bin/env python
u"""
arcticdata_tides.py
Written by Tyler Sutterley (05/2023)
Download Arctic Ocean Tide Models from the NSF ArcticData archive

AODTM-5: https://arcticdata.io/catalog/view/doi:10.18739/A2901ZG3N
AOTIM-5: https://arcticdata.io/catalog/view/doi:10.18739/A2S17SS80
AOTIM-5-2018: https://arcticdata.io/catalog/view/doi:10.18739/A21R6N14K
Arc2kmTM: https://arcticdata.io/catalog/view/doi:10.18739/A2D21RK6K
Gr1kmTM: https://arcticdata.io/catalog/view/doi:10.18739/A2B853K18

CALLING SEQUENCE:
    python arcticdata_tides.py --tide=Gr1kmTM

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: working data directory
    -T X, --tide X: Arctic tide model to download
        AODTM-5
        AOTIM-5
        AOTIM-5-2018
        Arc2kmTM
        Gr1kmTM
    -t X, --timeout X: timeout in seconds for blocking operations
    -M X, --mode X: Local permissions mode of the files downloaded

PYTHON DEPENDENCIES:
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 05/2023: added option to change connection timeout
    Updated 04/2023: using pathlib to define and expand paths
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 06/2022: added Greenland 1km model (Gr1kmTM) to list of models
    Updated 04/2022: use argparse descriptions within documentation
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: can use prefix files to define command line arguments
    Updated 10/2020: using argparse to set command line parameters
    Written 08/2020
"""
from __future__ import print_function, annotations

import re
import logging
import pathlib
import zipfile
import argparse
import posixpath
import pyTMD.utilities

# PURPOSE: Download Arctic Ocean Tide Models from the NSF ArcticData archive
def arcticdata_tides(MODEL: str,
    DIRECTORY: str | pathlib.Path | None = None,
    TIMEOUT: int | None = None,
    MODE: oct = 0o775):

    # create logger for verbosity level
    logger = pyTMD.utilities.build_logger(__name__,level=logging.INFO)

    # digital object identifier (doi) for each Arctic tide model
    DOI = {}
    DOI['AODTM-5'] = '10.18739/A2901ZG3N'
    DOI['AOTIM-5'] = '10.18739/A2S17SS80'
    DOI['AOTIM-5-2018'] = '10.18739/A21R6N14K'
    DOI['Arc2kmTM'] = '10.18739/A2D21RK6K'
    DOI['Gr1kmTM'] = '10.18739/A2B853K18'
    # local subdirectory for each Arctic tide model
    LOCAL = {}
    LOCAL['AODTM-5'] = 'aodtm5_tmd'
    LOCAL['AOTIM-5'] = 'aotim5_tmd'
    LOCAL['AOTIM-5-2018'] = 'Arc5km2018'
    LOCAL['Arc2kmTM'] = 'Arc2kmTM'
    LOCAL['Gr1kmTM'] = 'Gr1kmTM'

    # recursively create directories if non-existent
    DIRECTORY = pathlib.Path(DIRECTORY).expanduser().absolute()
    local_dir = DIRECTORY.joinpath(LOCAL[MODEL])
    local_dir.mkdir(MODE, parents=True, exist_ok=True)

    # build host url for model
    resource_map_doi = f'resource_map_doi:{DOI[MODEL]}'
    HOST = ['https://arcticdata.io','metacat','d1','mn','v2','packages',
        pyTMD.utilities.quote_plus(posixpath.join('application','bagit-097')),
        pyTMD.utilities.quote_plus(resource_map_doi)]
    # download zipfile from host
    logger.info(f'{posixpath.join(*HOST)} -->\n')
    zfile = zipfile.ZipFile(pyTMD.utilities.from_http(HOST, timeout=TIMEOUT))
    # find model files within zip file
    rx = re.compile(r'(grid|h[0]?|UV[0]?|Model|xy)_(.*?)',re.VERBOSE)
    members = [m for m in zfile.filelist if rx.search(m.filename)]
    # extract each member
    for m in members:
        # strip directories from member filename
        m.filename = posixpath.basename(m.filename)
        local_file = local_dir.joinpath(m.filename)
        logger.info(str(local_file))
        # extract file
        zfile.extract(m, path=local_dir)
        # change permissions mode
        local_file.chmod(mode=MODE)
    # close the zipfile object
    zfile.close()

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Download Arctic Ocean Tide Models from the NSF ArcticData
            archive
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = pyTMD.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory for location of tide models
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # Arctic Ocean tide model to download
    parser.add_argument('--tide','-T',
        metavar='TIDE', type=str, nargs='+', default=['Gr1kmTM'],
        choices=('AODTM-5','AOTIM-5','AOTIM-5-2018','Arc2kmTM','Gr1kmTM'),
        help='Arctic Ocean tide model to download')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
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

    # check internet connection before attempting to run program
    if pyTMD.utilities.check_connection('https://arcticdata.io'):
        for m in args.tide:
            arcticdata_tides(m,
                DIRECTORY=args.directory,
                TIMEOUT=args.timeout,
                MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
