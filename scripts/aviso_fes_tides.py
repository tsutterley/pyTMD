#!/usr/bin/env python
u"""
aviso_fes_tides.py
Written by Tyler Sutterley (06/2020)
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
        fes1999
        fes2004
        fes2012
        fes2014
    --load: download load tide model outputs (fes2014)
    --currents: download tide model current outputs (fes2012 and fes2014)
    --log: output log of files downloaded
    -M X, --mode=X: Local permissions mode of the files downloaded

UPDATE HISTORY:
    Updated 06/2020: added netrc option for alternative authentication
    Updated 05/2019: new authenticated ftp host (changed 2018-05-31)
    Written 09/2017
"""
from __future__ import print_function

import sys
import os
import io
import netrc
import getopt
import shutil
import tarfile
import getpass
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
    CURRENTS=False, LOG=False, MODE=None):
    #-- connect and login to AVISO ftp server
    f = ftplib.FTP('ftp-access.aviso.altimetry.fr')
    f.login(USER, PASSWORD)
    #-- check if local directory exists and recursively create if not
    local_dir = os.path.join(DIRECTORY,MODEL)
    os.makedirs(local_dir,MODE) if not os.path.exists(local_dir) else None

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
    TAR = {}
    PATH = {}
    #-- 1999 model
    FES['fes1999']=[]
    FES['fes1999'].append(['fes1999_fes2004','readme_fes1999.html'])
    FES['fes1999'].append(['fes1999_fes2004','fes1999.tar.gz'])
    TAR['fes1999'] = [None,'r:gz']
    PATH['fes1999'] = [local_dir,local_dir]
    #-- 2004 model
    FES['fes2004']=[]
    FES['fes2004'].append(['fes1999_fes2004','readme_fes2004.html'])
    FES['fes2004'].append(['fes1999_fes2004','fes2004.tar.gz'])
    TAR['fes2004'] = [None,'r:gz']
    PATH['fes2004'] = [local_dir,DIRECTORY]
    #-- 2012 model
    FES['fes2012']=[]
    FES['fes2012'].append(['fes2012_heights','readme_fes2012_heights_v1.1'])
    FES['fes2012'].append(['fes2012_heights','fes2012_heights_v1.1.tar.lzma'])
    TAR['fes2012'] = [None,'r:xz']
    PATH['fes2012'] = [local_dir,DIRECTORY]
    if CURRENTS:
        FES['fes2012'].append(['fes2012_currents','readme_fes2012_currents_v1.1'])
        FES['fes2012'].append(['fes2012_currents','fes2012_currents_v1.1_block1.tar.lzma'])
        FES['fes2012'].append(['fes2012_currents','fes2012_currents_v1.1_block2.tar.lzma'])
        FES['fes2012'].append(['fes2012_currents','fes2012_currents_v1.1_block3.tar.lzma'])
        FES['fes2012'].append(['fes2012_currents','fes2012_currents_v1.1_block4.tar.lzma'])
        TAR['fes2012'].extend([None,'r:xz','r:xz','r:xz','r:xz'])
        PATH['fes2012'].extend([local_dir,DIRECTORY,DIRECTORY,DIRECTORY,DIRECTORY])
    #-- 2014 model
    FES['fes2014']=[]
    FES['fes2014'].append(['fes2014_elevations_and_load','readme_fes2014_elevation_and_load_v1.2.txt'])
    FES['fes2014'].append(['fes2014_elevations_and_load','fes2014b_elevations','ocean_tide.tar.xz'])
    TAR['fes2014'] = [None,'r']
    PATH['fes2004'] = [local_dir,DIRECTORY]
    if LOAD:
        FES['fes2014'].append(['fes2014_elevations_and_load','fes2014b_elevations','load_tide.tar.xz'])
        TAR['fes2014'].extend(['r'])
        PATH['fes2012'].extend([DIRECTORY])
    if CURRENTS:
        FES['fes2014'].append(['fes2014a_currents','readme_fes2014_currents_v1.2.txt'])
        FES['fes2014'].append(['fes2014a_currents','eastward_velocity.tar.xz'])
        FES['fes2014'].append(['fes2014a_currents','northward_velocity.tar.xz'])
        TAR['fes2014'].extend([None,'r','r'])
        PATH['fes2014'].extend([local_dir,DIRECTORY,DIRECTORY])

    #-- for each file for a model
    for remote_path,tarmode in zip(FES[MODEL],TAR[MODEL]):
        #-- download file from ftp and decompress tar files
        ftp_download_file(fid,f,remote_path,local_dir,tarmode,MODE)
    #-- close the ftp connection
    f.quit()
    #-- close log file and set permissions level to MODE
    if LOG:
        fid.close()
        os.chmod(os.path.join(LOGDIR,LOGFILE), MODE)

#-- PURPOSE: pull file from a remote ftp server and decompress if tar file
def ftp_download_file(fid, ftp, remote_path, local_dir, tarmode, MODE):
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
        for member in member_files:
            local_file = os.path.join(local_dir,posixpath.basename(member.name))
            print('\t{0}'.format(local_file),file=fid)
            #-- extract file contents to new file
            with tar.extractfile(member) as f_in,open(local_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            #-- get last modified date of remote file within tar file
            remote_mtime = tarfile.TarInfo.isfile(member)
            #-- keep remote modification time of file and local access time
            os.utime(local_file, (os.stat(local_file).st_atime, remote_mtime))
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
    print('\tfes1999\n\tfes2004\n\tfes2012\n\tfes2014')
    print(' --load\t\t\tDownload load tide model outputs')
    print(' --currents\t\tDownload tide model current outputs')
    print(' -M X, --mode=X\t\tPermission mode of files downloaded')
    print(' -l, --log\t\tOutput log file')
    today = time.strftime('%Y-%m-%d',time.localtime())
    LOGFILE = 'AVISO_FES_tides_{0}.log'.format(today)
    print('    Log file format: {}\n'.format(LOGFILE))

#-- Main program that calls aviso_fes_tides()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','user=','netrc=','tide=','load',
        'currents','log','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hD:U:N:M:l',long_options)

    #-- command line parameters
    DIRECTORY = os.getcwd()
    USER = ''
    NETRC = None
    MODELS = ['fes2014']
    LOAD = False
    CURRENTS = False
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
            MODELS = arg.lower().split(',')
        elif opt in ("--load",):
            LOAD = True
        elif opt in ("--currents",):
            CURRENTS = True
        elif opt in ("-l","--log"):
            LOG = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- AVISO FTP Server hostname
    HOST = 'ftp.aviso.altimetry.fr'
    #-- get AVISO FTP Server credentials
    if not USER:
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
            LOAD=LOAD,CURRENTS=CURRENTS,LOG=LOG,MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
