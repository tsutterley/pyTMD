#!/usr/bin/env python
u"""
aviso_fes_tides.py
Written by Tyler Sutterley (07/2020)
Downloads the FES (Finite Element Solution) global tide model from AVISO
Decompresses the model tar files into the constituent files and auxiliary files
    https://www.aviso.altimetry.fr/data/products/auxiliary-products/
        global-tide-fes.html
    https://www.aviso.altimetry.fr/en/data/data-access.html

CALLING SEQUENCE:
    python aviso_fes_tides.py --user=<username> --tide=fes2014
    where <username> is your AVISO data dissemination server username

COMMAND LINE OPTIONS:
    --help: list the command line options
    --directory=X: working data directory
    --user=X: username for AVISO FTP servers (email)
    -N X, --netrc=X: path to .netrc file for authentication
    --tide=X: FES tide model to download
        FES1999
        FES2004
        FES2012
        FES2014
    --load: download load tide model outputs (fes2014)
    --currents: download tide model current outputs (fes2012 and fes2014)
    --log: output log of files downloaded
    -M X, --mode=X: Local permissions mode of the files downloaded

PYTHON DEPENDENCIES:
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

UPDATE HISTORY:
    Updated 07/2020: add gzip option to compress output ascii and netCDF4 files
    Updated 06/2020: added netrc option for alternative authentication
    Updated 05/2019: new authenticated ftp host (changed 2018-05-31)
    Written 09/2017
"""
from __future__ import print_function

import sys
import os
import io
import gzip
import netrc
import getopt
import shutil
import tarfile
import getpass
import builtins
import posixpath
import calendar, time
import ftplib

#-- PURPOSE: check internet connection
def check_connection(USER, PASSWORD):
    #-- attempt to connect to ftp host for AVISO products
    try:
        f = ftplib.FTP('ftp-access.aviso.altimetry.fr')
        f.login(USER, PASSWORD)
        f.voidcmd("NOOP")
    except IOError:
        raise RuntimeError('Check internet connection')
    except ftplib.error_perm:
        raise RuntimeError('Check login credentials')
    else:
        return True

#-- PURPOSE: download local AVISO FES files with ftp server
def aviso_fes_tides(MODEL, DIRECTORY=None, USER='', PASSWORD='', LOAD=False,
    CURRENTS=False, GZIP=False, LOG=False, MODE=None):
    #-- connect and login to AVISO ftp server
    f = ftplib.FTP('ftp-access.aviso.altimetry.fr',timeout=1000)
    f.login(USER, PASSWORD)
    #-- check if local directory exists and recursively create if not
    localpath = os.path.join(DIRECTORY,MODEL.lower())
    os.makedirs(localpath,MODE) if not os.path.exists(localpath) else None

    #-- create log file with list of downloaded files (or print to terminal)
    if LOG:
        #-- format: AVISO_FES_tides_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'AVISO_FES_tides_{0}.log'.format(today)
        fid = open(os.path.join(DIRECTORY,LOGFILE),'w')
        print('AVISO FES Sync Log ({0})'.format(today), file=fid)
        print('\tMODEL: {0}'.format(MODEL))
    else:
        #-- standard output (terminal output)
        fid = sys.stdout

    #-- path to remote directory for FES
    FES = {}
    #-- mode for reading tar files
    TAR = {}
    #-- flatten file structure
    FLATTEN = {}

    #-- 1999 model
    FES['FES1999']=[]
    FES['FES1999'].append(['fes1999_fes2004','readme_fes1999.html'])
    FES['FES1999'].append(['fes1999_fes2004','fes1999.tar.gz'])
    TAR['FES1999'] = [None,'r:gz']
    FLATTEN['FES1999'] = [None,True]
    #-- 2004 model
    FES['FES2004']=[]
    FES['FES2004'].append(['fes1999_fes2004','readme_fes2004.html'])
    FES['FES2004'].append(['fes1999_fes2004','fes2004.tar.gz'])
    TAR['FES2004'] = [None,'r:gz']
    FLATTEN['FES2004'] = [None,True]
    #-- 2012 model
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
    #-- 2014 model
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

    #-- for each file for a model
    for remotepath,tarmode,flatten in zip(FES[MODEL],TAR[MODEL],FLATTEN[MODEL]):
        #-- download file from ftp and decompress tar files
        ftp_download_file(fid,f,remotepath,localpath,tarmode,flatten,GZIP,MODE)
    #-- close the ftp connection
    f.quit()
    #-- close log file and set permissions level to MODE
    if LOG:
        fid.close()
        os.chmod(os.path.join(LOGDIR,LOGFILE), MODE)

#-- PURPOSE: pull file from a remote ftp server and decompress if tar file
def ftp_download_file(fid,ftp,remote_path,local_dir,tarmode,flatten,GZIP,MODE):
    #-- remote and local directory for data product
    remote_file = posixpath.join('auxiliary','tide_model',*remote_path)

    #-- Printing files transferred
    print('{0}{1}/{2} --> '.format('ftp://',ftp.host,remote_file),file=fid)
    if tarmode:
        #-- copy remote file contents to bytesIO object
        fileobj = io.BytesIO()
        ftp.retrbinary('RETR {0}'.format(remote_file), fileobj.write)
        fileobj.seek(0)
        #-- open the AOD1B monthly tar file
        tar = tarfile.open(name=remote_path[-1], fileobj=fileobj, mode=tarmode)
        #-- read tar file and extract all files
        member_files = [m for m in tar.getmembers() if tarfile.TarInfo.isfile(m)]
        for m in member_files:
            member = posixpath.basename(m.name) if flatten else m.name
            fileBasename,fileExtension = posixpath.splitext(m.name)
            #-- extract file contents to new file
            if fileExtension in ('.asc','.nc') and GZIP:
                local_file = os.path.join(local_dir,
                    *posixpath.split('{0}.gz'.format(member)))
                print('\t{0}'.format(local_file),file=fid)
                #-- recursively create output directory if non-existent
                if not os.access(os.path.dirname(local_file),os.F_OK):
                    os.makedirs(os.path.dirname(local_file),MODE)
                #-- extract file to compressed gzip format in local directory
                with tar.extractfile(m) as fi,gzip.open(local_file, 'wb') as fo:
                    shutil.copyfileobj(fi, fo)
            else:
                local_file = os.path.join(local_dir,*posixpath.split(member))
                print('\t{0}'.format(local_file),file=fid)
                #-- recursively create output directory if non-existent
                if not os.access(os.path.dirname(local_file),os.F_OK):
                    os.makedirs(os.path.dirname(local_file),MODE)
                #-- extract file to local directory
                with tar.extractfile(m) as fi,open(local_file, 'wb') as fo:
                    shutil.copyfileobj(fi, fo)
            #-- get last modified date of remote file within tar file
            #-- keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, m.mtime))
            os.chmod(local_file, MODE)
        print()
    else:
        #-- copy readme and uncompressed files directly
        local_file = os.path.join(local_dir,remote_path[-1])
        print('\t{0}\n'.format(local_file),file=fid)
        #-- copy remote file contents to local file
        with open(local_file, 'wb') as f:
            ftp.retrbinary('RETR {0}'.format(remote_file), f.write)
        #-- get last modified date of remote file and convert into unix time
        mdtm = ftp.sendcmd('MDTM {0}'.format(remote_file))
        remote_mtime = calendar.timegm(time.strptime(mdtm[4:],"%Y%m%d%H%M%S"))
        #-- keep remote modification time of file and local access time
        os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
        os.chmod(local_file, MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' -U X, --user=X\t\tUsername for AVISO FTP servers')
    print(' -N X, --netrc=X\tpath to .netrc file for authentication')
    print(' --tide=X\t\tFES tide model to download')
    print('\tFES1999\n\tFES2004\n\tFES2012\n\tFES2014')
    print(' --load\t\t\tDownload load tide model outputs')
    print(' --currents\t\tDownload tide model current outputs')
    print(' -G, --gzip\t\tCompress output ascii and netCDF4 tide files')
    print(' -M X, --mode=X\t\tPermission mode of files downloaded')
    print(' -l, --log\t\tOutput log file')
    today = time.strftime('%Y-%m-%d',time.localtime())
    LOGFILE = 'AVISO_FES_tides_{0}.log'.format(today)
    print('    Log file format: {}\n'.format(LOGFILE))

#-- Main program that calls aviso_fes_tides()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','user=','netrc=','tide=','load',
        'currents','gzip','log','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hD:U:N:M:l',long_options)

    #-- command line parameters
    DIRECTORY = os.getcwd()
    USER = ''
    NETRC = None
    MODELS = ['FES2014']
    LOAD = False
    CURRENTS = False
    GZIP = False
    LOG = False
    #-- permissions mode of the local directories and files (number in octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-D","--directory"):
            DIRECTORY = os.path.expanduser(arg)
        elif opt in ("-U","--user"):
            USER = arg
        elif opt in ("-N","--netrc"):
            NETRC = os.path.expanduser(arg)
        elif opt in ("--tide",):
            MODELS = arg.upper().split(',')
        elif opt in ("--load",):
            LOAD = True
        elif opt in ("--currents",):
            CURRENTS = True
        elif opt in ("-G","--gzip",):
            GZIP = True
        elif opt in ("-l","--log"):
            LOG = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- AVISO FTP Server hostname
    HOST = 'ftp.aviso.altimetry.fr'
    #-- get AVISO FTP Server credentials
    if not USER and not NETRC:
        #-- check that AVISO FTP Server credentials were entered
        USER = builtins.input('Username for {0}: '.format(HOST))
        #-- enter password securely from command-line
        PASSWORD = getpass.getpass('Password for {0}@{1}: '.format(USER,URS))
    elif NETRC:
        USER,LOGIN,PASSWORD = netrc.netrc(NETRC).authenticators(HOST)
    else:
        #-- enter password securely from command-line
        PASSWORD = getpass.getpass('Password for {0}@{1}: '.format(USER,HOST))

    #-- check internet connection before attempting to run program
    if check_connection(USER,PASSWORD):
        for m in MODELS:
            aviso_fes_tides(m,DIRECTORY=DIRECTORY,USER=USER,PASSWORD=PASSWORD,
                LOAD=LOAD,CURRENTS=CURRENTS,GZIP=GZIP,LOG=LOG,MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
