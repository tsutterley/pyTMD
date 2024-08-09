#!/usr/bin/env python
u"""
aviso_fes_tides.py
Written by Tyler Sutterley (07/2024)
Downloads the FES (Finite Element Solution) global tide model from AVISO
Decompresses the model tar files into the constituent files and auxiliary files
    https://www.aviso.altimetry.fr/data/products/auxiliary-products/
        global-tide-fes.html
    https://www.aviso.altimetry.fr/en/data/data-access.html

CALLING SEQUENCE:
    python aviso_fes_tides.py --user <username> --tide FES2014
    where <username> is your AVISO data dissemination server username

COMMAND LINE OPTIONS:
    --help: list the command line options
    --directory X: working data directory
    -U X, --user: username for AVISO FTP servers (email)
    -P X, --password: password for AVISO FTP servers
    -N X, --netrc X: path to .netrc file for authentication
    --tide X: FES tide model to download
        FES1999
        FES2004
        FES2012
        FES2014
        FES2022
    --load: download load tide model outputs (FES2014)
    --currents: download tide model current outputs (FES2012 and FES2014)
    --extrapolated: Download extrapolated tide model outputs (FES2022)
    -G, --gzip: compress output ascii and netCDF4 tide files
    -t X, --timeout X: timeout in seconds for blocking operations
    --log: output log of files downloaded
    -M X, --mode X: Local permissions mode of the files downloaded

PYTHON DEPENDENCIES:
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 07/2024: added list and download for FES2022 tide model
        compare modification times with remote to not overwrite files
    Updated 05/2023: added option to change connection timeout
    Updated 04/2023: using pathlib to define and expand paths
        added option to include AVISO FTP password as argument
    Updated 11/2022: added encoding for writing ascii files
        use f-strings for formatting verbose or ascii output
    Updated 04/2022: use argparse descriptions within documentation
    Updated 10/2021: using python logging for handling verbose output
    Updated 07/2021: can use prefix files to define command line arguments
    Updated 05/2021: use try/except for retrieving netrc credentials
    Updated 04/2021: set a default netrc file and check access
    Updated 10/2020: using argparse to set command line parameters
    Updated 07/2020: add gzip option to compress output ascii and netCDF4 files
    Updated 06/2020: added netrc option for alternative authentication
    Updated 05/2019: new authenticated ftp host (changed 2018-05-31)
    Written 09/2017
"""
from __future__ import print_function, annotations

import sys
import os
import io
import re
import gzip
import lzma
import netrc
import shutil
import logging
import tarfile
import getpass
import pathlib
import argparse
import builtins
import posixpath
import calendar, time
import ftplib
import pyTMD.utilities

# PURPOSE: download local AVISO FES files with ftp server
def aviso_fes_tides(MODEL: str,
        DIRECTORY: str | pathlib.Path | None = None,
        USER: str = '',
        PASSWORD: str = '',
        LOAD: bool = False,
        CURRENTS: bool = False,
        EXTRAPOLATED: bool = False,
        GZIP: bool = False,
        TIMEOUT: int | None = None,
        LOG: bool = False,
        MODE: oct = 0o775
    ):

    # connect and login to AVISO ftp server
    f = ftplib.FTP('ftp-access.aviso.altimetry.fr', timeout=TIMEOUT)
    f.login(USER, PASSWORD)

    # create log file with list of downloaded files (or print to terminal)
    if LOG:
        # format: AVISO_FES_tides_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = DIRECTORY.joinpath(f'AVISO_FES_tides_{today}.log')
        fid = LOGFILE.open(mode='w', encoding='utf8')
        logger = pyTMD.utilities.build_logger(__name__,stream=fid,
            level=logging.INFO)
        logger.info(f'AVISO FES Sync Log ({today})')
        logger.info(f'\tMODEL: {MODEL}')
    else:
        # standard output (terminal output)
        logger = pyTMD.utilities.build_logger(__name__,level=logging.INFO)

    # download the FES tide model files
    if MODEL in ('FES1999','FES2004','FES2012','FES2014'):
        aviso_fes_tar(MODEL, f, logger,
            DIRECTORY=DIRECTORY,
            LOAD=LOAD,
            CURRENTS=CURRENTS,
            GZIP=GZIP,
            MODE=MODE)
    elif MODEL in ('FES2022',):
        aviso_fes_list(MODEL, f, logger,
            DIRECTORY=DIRECTORY,
            LOAD=LOAD,
            CURRENTS=CURRENTS,
            EXTRAPOLATED=EXTRAPOLATED,
            GZIP=GZIP,
            MODE=MODE)

    # close the ftp connection
    f.quit()
    # close log file and set permissions level to MODE
    if LOG:
        LOGFILE.chmod(mode=MODE)

# PURPOSE: download local AVISO FES files with ftp server
# by downloading tar files and extracting contents
def aviso_fes_tar(MODEL, f, logger,
        DIRECTORY: str | pathlib.Path | None = None,
        LOAD: bool = False,
        CURRENTS: bool = False,
        GZIP: bool = False,
        MODE: oct = 0o775
    ):

    # check if local directory exists and recursively create if not
    localpath = pathlib.Path(DIRECTORY).joinpath(MODEL.lower()).expanduser()
    localpath.mkdir(MODE, parents=True, exist_ok=True)

    # path to remote directory for FES
    FES = {}
    # mode for reading tar files
    TAR = {}
    # flatten file structure
    FLATTEN = {}

    # 1999 model
    FES['FES1999']=[]
    FES['FES1999'].append(['fes1999_fes2004','readme_fes1999.html'])
    FES['FES1999'].append(['fes1999_fes2004','fes1999.tar.gz'])
    TAR['FES1999'] = [None,'r:gz']
    FLATTEN['FES1999'] = [None,True]
    # 2004 model
    FES['FES2004']=[]
    FES['FES2004'].append(['fes1999_fes2004','readme_fes2004.html'])
    FES['FES2004'].append(['fes1999_fes2004','fes2004.tar.gz'])
    TAR['FES2004'] = [None,'r:gz']
    FLATTEN['FES2004'] = [None,True]
    # 2012 model
    FES['FES2012']=[]
    FES['FES2012'].append(['fes2012_heights','readme_fes2012_heights_v1.1'])
    FES['FES2012'].append(['fes2012_heights','fes2012_heights_v1.1.tar.lzma'])
    TAR['FES2012'] = []
    TAR['FES2012'].extend([None,'r:xz'])
    FLATTEN['FES2012'] = []
    FLATTEN['FES2012'].extend([None,True])
    if CURRENTS:
        subdir = 'fes2012_currents'
        FES['FES2012'].append([subdir,'readme_fes2012_currents_v1.1'])
        FES['FES2012'].append([subdir,'fes2012_currents_v1.1_block1.tar.lzma'])
        FES['FES2012'].append([subdir,'fes2012_currents_v1.1_block2.tar.lzma'])
        FES['FES2012'].append([subdir,'fes2012_currents_v1.1_block3.tar.lzma'])
        FES['FES2012'].append([subdir,'fes2012_currents_v1.1_block4.tar.lzma'])
        TAR['FES2012'].extend([None,'r:xz','r:xz','r:xz','r:xz'])
        FLATTEN['FES2012'].extend([None,False,False,False,False])
    # 2014 model
    FES['FES2014']=[]
    FES['FES2014'].append(['fes2014_elevations_and_load',
        'readme_fes2014_elevation_and_load_v1.2.txt'])
    FES['FES2014'].append(['fes2014_elevations_and_load',
        'fes2014b_elevations','ocean_tide.tar.xz'])
    TAR['FES2014'] = []
    TAR['FES2014'].extend([None,'r'])
    FLATTEN['FES2014'] = []
    FLATTEN['FES2014'].extend([None,False])
    if LOAD:
        FES['FES2014'].append(['fes2014_elevations_and_load',
            'fes2014a_loadtide','load_tide.tar.xz'])
        TAR['FES2014'].extend(['r'])
        FLATTEN['FES2014'].extend([False])
    if CURRENTS:
        subdir = 'fes2014a_currents'
        FES['FES2014'].append([subdir,'readme_fes2014_currents_v1.2.txt'])
        FES['FES2014'].append([subdir,'eastward_velocity.tar.xz'])
        FES['FES2014'].append([subdir,'northward_velocity.tar.xz'])
        TAR['FES2014'].extend(['r'])
        FLATTEN['FES2014'].extend([False])

    # for each file for a model
    for remotepath,tarmode,flatten in zip(FES[MODEL],TAR[MODEL],FLATTEN[MODEL]):
        # download file from ftp and decompress tar files
        ftp_download(logger, f, remotepath, localpath,
            TARMODE=tarmode,
            FLATTEN=flatten,
            GZIP=GZIP,
            MODE=MODE
        )

# PURPOSE: download local AVISO FES files with ftp server
# by downloading individual files
def aviso_fes_list(MODEL, f, logger,
        DIRECTORY: str | pathlib.Path | None = None,
        LOAD: bool = False,
        CURRENTS: bool = False,
        EXTRAPOLATED: bool = False,
        GZIP: bool = False,
        MODE: oct = 0o775
    ):
    # validate local directory
    DIRECTORY = pathlib.Path(DIRECTORY).expanduser().absolute()

    # path to remote directory for FES
    FES = {}
    # 2022 model
    FES['FES2022'] = []
    FES['FES2022'].append(['fes2022b','ocean_tide'])
    if LOAD:
        FES['FES2022'].append(['fes2022b','load_tide'])
    if EXTRAPOLATED:
        FES['FES2022'].append(['fes2022b','ocean_tide_extrapolated'])

    # for each model file type
    for subdir in FES[MODEL]:
        local_dir = DIRECTORY.joinpath(*subdir)
        file_list = ftp_list(f, subdir, basename=True, sort=True)
        for fi in file_list:
            remote_path = [*subdir, fi]
            LZMA = fi.endswith('.xz')
            ftp_download(logger, f, remote_path, local_dir,
                LZMA=LZMA,
                GZIP=GZIP,
                CHUNK=32768,
                MODE=MODE
            )

# PURPOSE: List a directory on a ftp host
def ftp_list(ftp, remote_path, basename=False, pattern=None, sort=False):
    # list remote path
    output = ftp.nlst(posixpath.join('auxiliary','tide_model',*remote_path))
    # reduce to basenames
    if basename:
        output = [posixpath.basename(i) for i in output]
    # reduce using regular expression pattern
    if pattern:
        i = [i for i,f in enumerate(output) if re.search(pattern,f)]
        # reduce list of listed items
        output = [output[indice] for indice in i]
    # sort the list
    if sort:
        i = [i for i,j in sorted(enumerate(output), key=lambda i: i[1])]
        # sort list of listed items
        output = [output[indice] for indice in i]
    # return the list of items
    return output

# PURPOSE: pull file from a remote ftp server and decompress if tar file
def ftp_download(logger, ftp, remote_path, local_dir,
        LZMA=None,
        TARMODE=None,
        FLATTEN=None,
        GZIP=False,
        CHUNK=8192,
        MODE=0o775
    ):
    # remote and local directory for data product
    remote_file = posixpath.join('auxiliary','tide_model',*remote_path)
    # if compressing the output file
    opener = gzip.open if GZIP else open

    # Printing files transferred
    remote_ftp_url = posixpath.join('ftp://', ftp.host, remote_file)
    logger.info(f'{remote_ftp_url} -->')
    if TARMODE:
        # copy remote file contents to bytesIO object
        fileobj = io.BytesIO()
        ftp.retrbinary(f'RETR {remote_file}', fileobj.write, blocksize=CHUNK)
        fileobj.seek(0)
        # open the tar file
        tar = tarfile.open(name=remote_path[-1], fileobj=fileobj, mode=TARMODE)
        # read tar file and extract all files
        member_files = [m for m in tar.getmembers() if tarfile.TarInfo.isfile(m)]
        for m in member_files:
            member = posixpath.basename(m.name) if FLATTEN else m.name
            base, sfx = posixpath.splitext(m.name)
            # extract file contents to new file
            output = f'{member}.gz' if sfx in ('.asc','.nc') and GZIP else member
            local_file = local_dir.joinpath(*posixpath.split(output))
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
    elif LZMA:
        # get last modified date of remote file and convert into unix time
        mdtm = ftp.sendcmd(f'MDTM {remote_file}')
        mtime = calendar.timegm(time.strptime(mdtm[4:],"%Y%m%d%H%M%S"))
        # output file name for compressed and uncompressed cases
        stem = posixpath.basename(posixpath.splitext(remote_file)[0])
        base, sfx = posixpath.splitext(stem)
        # extract file contents to new file
        output = f'{stem}.gz' if sfx in ('.asc','.nc') and GZIP else stem
        local_file = local_dir.joinpath(output)
        # check if the local file exists
        if local_file.exists() and newer(mtime,local_file.stat().st_mtime):
            # check the modification time of the local file
            # if remote file is newer: overwrite the local file
            return
        # print the file being transferred
        logger.info(f'\t{str(local_file)}')
        # recursively create output directory if non-existent
        local_file.parent.mkdir(mode=MODE, parents=True, exist_ok=True)
        # copy remote file contents to bytesIO object
        fileobj = io.BytesIO()
        ftp.retrbinary(f'RETR {remote_file}', fileobj.write, blocksize=CHUNK)
        fileobj.seek(0)
        # decompress lzma file and extract contents to local directory
        with lzma.open(fileobj) as fi,opener(local_file, 'wb') as fo:
            shutil.copyfileobj(fi, fo)
        # get last modified date of remote file within tar file
        # keep remote modification time of file and local access time
        pathlib.os.utime(local_file, (local_file.stat().st_atime, mtime))
        local_file.chmod(mode=MODE)
    else:
        # copy readme and uncompressed files directly
        stem = posixpath.basename(remote_file)
        base, sfx = posixpath.splitext(stem)
        # output file name for compressed and uncompressed cases
        output = f'{stem}.gz' if sfx in ('.asc','.nc') and GZIP else stem
        local_file = local_dir.joinpath(output)
        # get last modified date of remote file and convert into unix time
        mdtm = ftp.sendcmd(f'MDTM {remote_file}')
        mtime = calendar.timegm(time.strptime(mdtm[4:],"%Y%m%d%H%M%S"))
        # check if the local file exists
        if local_file.exists() and newer(mtime, local_file.stat().st_mtime):
            # check the modification time of the local file
            # if remote file is newer: overwrite the local file
            return
        # print the file being transferred
        logger.info(f'\t{str(local_file)}\n')
        # recursively create output directory if non-existent
        local_file.parent.mkdir(mode=MODE, parents=True, exist_ok=True)
        # copy remote file contents to local file
        with opener(local_file, 'wb') as f:
            ftp.retrbinary(f'RETR {remote_file}', f.write, blocksize=CHUNK)
        # keep remote modification time of file and local access time
        pathlib.os.utime(local_file, (local_file.stat().st_atime, mtime))
        local_file.chmod(mode=MODE)

# PURPOSE: compare the modification time of two files
def newer(t1: int, t2: int) -> bool:
    return (pyTMD.utilities.even(t1) <= pyTMD.utilities.even(t2))

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Downloads the FES (Finite Element Solution) global tide
            model from AVISO.  Decompresses the model tar files into the
            constituent files and auxiliary files.
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = pyTMD.utilities.convert_arg_line_to_args
    # command line parameters
    # AVISO FTP credentials
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('AVISO_USERNAME'),
        help='Username for AVISO Login')
    parser.add_argument('--password','-W',
        type=str, default=os.environ.get('AVISO_PASSWORD'),
        help='Password for AVISO Login')
    parser.add_argument('--netrc','-N',
        type=pathlib.Path, default=pathlib.Path().home().joinpath('.netrc'),
        help='Path to .netrc file for authentication')
    # working data directory
    parser.add_argument('--directory','-D',
        type=pathlib.Path, default=pathlib.Path.cwd(),
        help='Working data directory')
    # FES tide models
    choices = ['FES1999','FES2004','FES2012','FES2014','FES2022']
    parser.add_argument('--tide','-T',
        metavar='TIDE', type=str, nargs='+',
        default=['FES2022'], choices=choices,
        help='FES tide model to download')
    # download FES load tides
    parser.add_argument('--load',
        default=False, action='store_true',
        help='Download load tide model outputs')
    # download FES tidal currents
    parser.add_argument('--currents',
        default=False, action='store_true',
        help='Download tide model current outputs')
    # download extrapolate FES tidal data
    parser.add_argument('--extrapolated',
        default=False, action='store_true',
        help='Download extrapolated tide model outputs')
    # compress output ascii and netCDF4 tide files with gzip
    parser.add_argument('--gzip','-G',
        default=False, action='store_true',
        help='Compress output ascii and netCDF4 tide files')
    # connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    # Output log file in form
    # AVISO_FES_tides_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    # permissions mode of the local directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files downloaded')
    # return the parser
    return parser

# This is the main part of the program that calls the individual functions
def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # AVISO FTP Server hostname
    HOST = 'ftp-access.aviso.altimetry.fr'
    # get authentication
    if not args.user and not args.netrc.exists():
        # check that AVISO credentials were entered
        args.user = builtins.input(f'Username for {HOST}: ')
        # enter password securely from command-line
        args.password = getpass.getpass(f'Password for {args.user}@{HOST}: ')
    elif args.netrc.exists():
        args.user,_,args.password = netrc.netrc(args.netrc).authenticators(HOST)
    elif args.user and not args.password:
        # enter password securely from command-line
        args.password = getpass.getpass(f'Password for {args.user}@{HOST}: ')

    # check internet connection before attempting to run program
    if pyTMD.utilities.check_ftp_connection(HOST,args.user,args.password):
        for m in args.tide:
            aviso_fes_tides(m,
                DIRECTORY=args.directory,
                USER=args.user,
                PASSWORD=args.password,
                LOAD=args.load,
                CURRENTS=args.currents,
                EXTRAPOLATED=args.extrapolated,
                GZIP=args.gzip,
                TIMEOUT=args.timeout,
                LOG=args.log,
                MODE=args.mode)

# run main program
if __name__ == '__main__':
    main()
