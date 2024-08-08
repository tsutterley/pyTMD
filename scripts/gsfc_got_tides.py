#!/usr/bin/env python
u"""
gsfc_got_tides.py
Written by Tyler Sutterley (08/2024)
Download GSFC Global Ocean Tide (GOT) models

CALLING SEQUENCE:
    python gsfc_got_tides.py --tide=GOT5.6

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: working data directory
    -T X, --tide X: GOT tide model to download
        GOT4.8
        GOT4.10
        GOT5.5
        GOT5.5D
        GOT5.6
    --format: GOT tide model format to download
        ascii
        netCDF
    -G, --gzip: compress output ascii and netCDF4 tide files
    -t X, --timeout X: timeout in seconds for blocking operations
    -M X, --mode X: Local permissions mode of the files downloaded

PYTHON DEPENDENCIES:
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 08/2024: keep prime nomenclature for 3rd degree tides
    Written 07/2024
"""
from __future__ import print_function, annotations

import re
import gzip
import shutil
import logging
import pathlib
import tarfile
import argparse
import posixpath
import pyTMD.utilities

# PURPOSE: Download Arctic Ocean Tide Models from the NSF ArcticData archive
def gsfc_got_tides(MODEL: str,
    DIRECTORY: str | pathlib.Path | None = None,
    FORMAT: str = 'netcdf',
    GZIP: bool = False,
    TIMEOUT: int | None = None,
    MODE: oct = 0o775):

    # create logger for verbosity level
    logger = pyTMD.utilities.build_logger(__name__,level=logging.INFO)
    # if compressing the output files
    opener = gzip.open if GZIP else open

    # url for each tide model tarfile
    PATH = {}
    PATH['GOT4.8'] = ['2022-07','got4.8.tar.gz']
    PATH['GOT4.10'] = ['2023-12','got4.10c.tar.gz']
    PATH['GOT5.5'] = ['2024-07','GOT5.5.tar%201.gz']
    PATH['GOT5.5D'] = ['2024-07','GOT5.5D.tar%201.gz']
    PATH['GOT5.6'] = ['2024-07','GOT5.6.tar%201.gz']

    # recursively create directories if non-existent
    DIRECTORY = pathlib.Path(DIRECTORY).expanduser().absolute()
    DIRECTORY.mkdir(MODE, parents=True, exist_ok=True)

    # build host url for model
    URL = ['https://earth.gsfc.nasa.gov','sites','default','files',*PATH[MODEL]]
    # download tarfile from host
    logger.info(f'{posixpath.join(*URL)} -->\n')
    fileobj = pyTMD.utilities.from_http(URL, timeout=TIMEOUT)
    # open the tar file
    tar = tarfile.open(name=PATH[MODEL][-1], fileobj=fileobj, mode='r:gz')
    # read tar file and extract all files
    member_files = [m for m in tar.getmembers() if tarfile.TarInfo.isfile(m)]
    for m in member_files:
        # extract file contents to new file
        base, sfx = posixpath.splitext(m.name)
        # skip files that are not in the desired format
        if (sfx == '.nc') and (FORMAT == 'ascii'):
            continue
        elif (sfx == '.d') and (FORMAT == 'netcdf'):
            continue
        elif re.match(r'^._', posixpath.basename(m.name)):
            continue
        # output file name
        output = f'{m.name}.gz' if sfx in ('.d','.nc') and GZIP else m.name
        local_file = DIRECTORY.joinpath(*posixpath.split(output))
        # check if the local file exists
        if local_file.exists() and newer(m.mtime, local_file.stat().st_mtime):
            # check the modification time of the local file
            # if remote file is newer: overwrite the local file
            continue
        # print the file being transferred
        logger.info(f'\t{str(local_file)}')
        # recursively create output directory if non-existent
        local_file.parent.mkdir(mode=MODE, parents=True, exist_ok=True)
        # extract file to local directory
        with tar.extractfile(m) as fi,opener(local_file, 'wb') as fo:
            shutil.copyfileobj(fi, fo)
        # get last modified date of remote file within tar file
        # keep remote modification time of file and local access time
        pathlib.os.utime(local_file, (local_file.stat().st_atime, m.mtime))
        local_file.chmod(mode=MODE)

# PURPOSE: compare the modification time of two files
def newer(t1: int, t2: int) -> bool:
    return (pyTMD.utilities.even(t1) <= pyTMD.utilities.even(t2))

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Download Global Ocean Tide models from NASA
            Goddard Space Flight Center (GSFC)
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = pyTMD.utilities.convert_arg_line_to_args
    # command line parameters
    # working data directory for location of tide models
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # Global Ocean Tide model to download
    parser.add_argument('--tide','-T',
        metavar='TIDE', type=str, nargs='+', default=['GOT5.5'],
        choices=('GOT4.8','GOT4.10','GOT5.5','GOT5.5D','GOT5.6'),
        help='Global Ocean Tide model to download')
    # Global Ocean Tide model format to download
    parser.add_argument('--format',
        type=str, default='netcdf',
        choices=('ascii','netcdf'),
        help='Global Ocean Tide model format to download')
    # compress output ascii and netCDF4 tide files with gzip
    parser.add_argument('--gzip','-G',
        default=False, action='store_true',
        help='Compress output ascii and netCDF tide files')
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
    if pyTMD.utilities.check_connection('https://earth.gsfc.nasa.gov'):
        for m in args.tide:
            gsfc_got_tides(m,
                DIRECTORY=args.directory,
                FORMAT=args.format,
                GZIP=args.gzip,
                TIMEOUT=args.timeout,
                MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
