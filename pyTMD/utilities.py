#!/usr/bin/env python
u"""
utilities.py
Written by Tyler Sutterley (06/2023)
Download and management utilities for syncing time and auxiliary files

PYTHON DEPENDENCIES:
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

UPDATE HISTORY:
    Updated 06/2023: add functions to retrieve and revoke Earthdata tokens
    Updated 05/2023: add reify decorator for evaluation of properties
        make urs a keyword argument in CCDIS list and download functions
        add case for JPL kernel file download where local path is defined
    Updated 04/2023: using pathlib to define and expand paths
        added function to download ephemeride files from JPL SSD server
    Updated 03/2023: add basic variable typing to function inputs
    Updated 01/2023: updated SSL context to fix some deprecation warnings
    Updated 11/2022: added list program for IERS Bulletin-A https server
        use f-strings for formatting verbose or ascii output
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 10/2021: build python logging instance for handling verbose output
    Updated 09/2021: added generic list from Apache http server
    Updated 08/2021: added function to open a file path
    Updated 07/2021: add parser for converting file files to arguments
    Updated 03/2021: added sha1 option for retrieving file hashes
    Updated 01/2021: added username and password to ftp functions
        added ftp connection check
    Updated 12/2020: added file object keyword for downloads if verbose
        add url split function for creating url location lists
    Updated 11/2020: normalize source and destination paths in copy
        make context an optional keyword argument in from_http
    Updated 09/2020: copy from http and https to bytesIO object in chunks
        use netrc credentials if not entered from CDDIS functions
        generalize build opener function for different Earthdata instances
    Updated 08/2020: add GSFC CDDIS opener, login and download functions
    Written 08/2020
"""
from __future__ import print_function, division, annotations

import sys
import os
import re
import io
import ssl
import json
import netrc
import ftplib
import shutil
import base64
import socket
import getpass
import inspect
import hashlib
import logging
import pathlib
import builtins
import warnings
import posixpath
import subprocess
import lxml.etree
import calendar, time
import dateutil.parser
if sys.version_info[0] == 2:
    from urllib import quote_plus
    from cookielib import CookieJar
    import urllib2
else:
    from urllib.parse import quote_plus
    from http.cookiejar import CookieJar
    import urllib.request as urllib2

# PURPOSE: get absolute path within a package from a relative path
def get_data_path(relpath: list | str | pathlib.Path):
    """
    Get the absolute path within a package from a relative path

    Parameters
    ----------
    relpath: list, str or pathlib.Path
        relative path
    """
    # current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    filepath = pathlib.Path(filename).absolute().parent
    if isinstance(relpath, list):
        # use *splat operator to extract from list
        return filepath.joinpath(*relpath)
    elif isinstance(relpath, (str, pathlib.Path)):
        return filepath.joinpath(relpath)

class reify(object):
    """Class decorator that puts the result of the method it
    decorates into the instance"""
    def __init__(self, wrapped):
        self.wrapped = wrapped
        self.__name__ = wrapped.__name__
        self.__doc__ = wrapped.__doc__

    def __get__(self, inst, objtype=None):
        if inst is None:
            return self
        val = self.wrapped(inst)
        setattr(inst, self.wrapped.__name__, val)
        return val

# PURPOSE: platform independent file opener
def file_opener(filename: str | pathlib.Path):
    """
    Platform independent file opener

    Parameters
    ----------
    filename: str or pathlib.Path
        path to file
    """
    filename = pathlib.Path(filename).expanduser()
    if (sys.platform == "win32"):
        os.startfile(filename, "explore")
    elif (sys.platform == "darwin"):
        subprocess.call(["open", filename])
    else:
        subprocess.call(["xdg-open", filename])

# PURPOSE: get the hash value of a file
def get_hash(
        local: str | io.IOBase | pathlib.Path,
        algorithm: str = 'MD5'
    ):
    """
    Get the hash value from a local file or ``BytesIO`` object

    Parameters
    ----------
    local: obj, str or pathlib.Path
        BytesIO object or path to file
    algorithm: str, default 'MD5'
        hashing algorithm for checksum validation

            - ``'MD5'``: Message Digest
            - ``'sha1'``: Secure Hash Algorithm
    """
    # check if open file object or if local file exists
    if isinstance(local, io.IOBase):
        if (algorithm == 'MD5'):
            return hashlib.md5(local.getvalue()).hexdigest()
        elif (algorithm == 'sha1'):
            return hashlib.sha1(local.getvalue()).hexdigest()
    elif isinstance(local, (str, pathlib.Path)):
        # generate checksum hash for local file
        local = pathlib.Path(local).expanduser()
        # if file currently doesn't exist, return empty string
        if not local.exists():
            return ''
        # open the local_file in binary read mode
        with local.open(mode='rb') as local_buffer:
            # generate checksum hash for a given type
            if (algorithm == 'MD5'):
                return hashlib.md5(local_buffer.read()).hexdigest()
            elif (algorithm == 'sha1'):
                return hashlib.sha1(local_buffer.read()).hexdigest()
    else:
        return ''

# PURPOSE: get the git hash value
def get_git_revision_hash(
        refname: str = 'HEAD',
        short: bool = False
    ):
    """
    Get the ``git`` hash value for a particular reference

    Parameters
    ----------
    refname: str, default HEAD
        Symbolic reference name
    short: bool, default False
        Return the shorted hash value
    """
    # get path to .git directory from current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    basepath = pathlib.Path(filename).absolute().parent.parent
    gitpath = basepath.joinpath('.git')
    # build command
    cmd = ['git', f'--git-dir={gitpath}', 'rev-parse']
    cmd.append('--short') if short else None
    cmd.append(refname)
    # get output
    with warnings.catch_warnings():
        return str(subprocess.check_output(cmd), encoding='utf8').strip()

# PURPOSE: get the current git status
def get_git_status():
    """Get the status of a ``git`` repository as a boolean value
    """
    # get path to .git directory from current file path
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    basepath = pathlib.Path(filename).absolute().parent.parent
    gitpath = basepath.joinpath('.git')
    # build command
    cmd = ['git', f'--git-dir={gitpath}', 'status', '--porcelain']
    with warnings.catch_warnings():
        return bool(subprocess.check_output(cmd))

# PURPOSE: recursively split a url path
def url_split(s: str):
    """
    Recursively split a url path into a list

    Parameters
    ----------
    s: str
        url string
    """
    head, tail = posixpath.split(s)
    if head in ('http:','https:','ftp:','s3:'):
        return s,
    elif head in ('', posixpath.sep):
        return tail,
    return url_split(head) + (tail,)

# PURPOSE: convert file lines to arguments
def convert_arg_line_to_args(arg_line):
    """
    Convert file lines to arguments

    Parameters
    ----------
    arg_line: str
        line string containing a single argument and/or comments
    """
    # remove commented lines and after argument comments
    for arg in re.sub(r'\#(.*?)$',r'',arg_line).split():
        if not arg.strip():
            continue
        yield arg

# PURPOSE: build a logging instance with a specified name
def build_logger(name: str, **kwargs):
    """
    Builds a logging instance with the specified name

    Parameters
    ----------
    name: str
        name of the logger
    format: str
        event description message format
    level: int, default logging.CRITICAL
        lowest-severity log message logger will handle
    propagate: bool, default False
        events logged will be passed to higher level handlers
    stream: obj or NoneType, default None
        specified stream to initialize StreamHandler
    """
    # set default arguments
    kwargs.setdefault('format', '%(levelname)s:%(name)s:%(message)s')
    kwargs.setdefault('level', logging.CRITICAL)
    kwargs.setdefault('propagate',False)
    kwargs.setdefault('stream',None)
    # build logger
    logger = logging.getLogger(name)
    logger.setLevel(kwargs['level'])
    logger.propagate = kwargs['propagate']
    # create and add handlers to logger
    if not logger.handlers:
        # create handler for logger
        handler = logging.StreamHandler(stream=kwargs['stream'])
        formatter = logging.Formatter(kwargs['format'])
        handler.setFormatter(formatter)
        # add handler to logger
        logger.addHandler(handler)
    return logger

# PURPOSE: convert Roman numerals to (Arabic) integers
def roman_to_int(roman: str):
    """
    Converts a string from Roman numerals into an integer (Arabic)

    Parameters
    ----------
    roman: str
        Roman numeral string
    """
    # mapping between Roman and Arabic numerals
    roman_map = {'i':1, 'v':5, 'x':10, 'l':50, 'c':100, 'd':500, 'm':1000}
    # verify case
    roman = roman.lower()
    output = 0
    # iterate through roman numerals in string and calculate total
    for i,s in enumerate(roman):
        if (i > 0) and (roman_map[s] > roman_map[roman[i-1]]):
            output += roman_map[s] - 2*roman_map[roman[i-1]]
        else:
            output += roman_map[s]
    # return the integer value
    return output

# PURPOSE: returns the Unix timestamp value for a formatted date string
def get_unix_time(
        time_string: str,
        format: str = '%Y-%m-%d %H:%M:%S'
    ):
    """
    Get the Unix timestamp value for a formatted date string

    Parameters
    ----------
    time_string: str
        formatted time string to parse
    format: str, default '%Y-%m-%d %H:%M:%S'
        format for input time string
    """
    try:
        parsed_time = time.strptime(time_string.rstrip(), format)
    except (TypeError, ValueError):
        pass
    else:
        return calendar.timegm(parsed_time)
    # try parsing with dateutil
    try:
        parsed_time = dateutil.parser.parse(time_string.rstrip())
    except (TypeError, ValueError):
        return None
    else:
        return parsed_time.timestamp()

# PURPOSE: output a time string in isoformat
def isoformat(time_string: str):
    """
    Reformat a date string to ISO formatting

    Parameters
    ----------
    time_string: str
        formatted time string to parse
    """
    # try parsing with dateutil
    try:
        parsed_time = dateutil.parser.parse(time_string.rstrip())
    except (TypeError, ValueError):
        return None
    else:
        return parsed_time.isoformat()

# PURPOSE: rounds a number to an even number less than or equal to original
def even(value: float):
    """
    Rounds a number to an even number less than or equal to original

    Parameters
    ----------
    value: float
        number to be rounded
    """
    return 2*int(value//2)

# PURPOSE: rounds a number upward to its nearest integer
def ceil(value: float):
    """
    Rounds a number upward to its nearest integer

    Parameters
    ----------
    value: float
        number to be rounded upward
    """
    return -int(-value//1)

# PURPOSE: make a copy of a file with all system information
def copy(
        source: str | pathlib.Path,
        destination: str | pathlib.Path,
        move: bool = False,
        **kwargs
    ):
    """
    Copy or move a file with all system information

    Parameters
    ----------
    source: str
        source file
    destination: str
        copied destination file
    move: bool, default False
        remove the source file
    """
    source = pathlib.Path(source).expanduser().absolute()
    destination = pathlib.Path(destination).expanduser().absolute()
    # log source and destination
    logging.info(f'{str(source)} -->\n\t{str(destination)}')
    shutil.copyfile(source, destination)
    shutil.copystat(source, destination)
    # remove the original file if moving
    if move:
        source.unlink()

# PURPOSE: check ftp connection
def check_ftp_connection(
        HOST: str,
        username: str | None = None,
        password: str | None = None
    ):
    """
    Check internet connection with ftp host

    Parameters
    ----------
    HOST: str
        remote ftp host
    username: str or NoneType
        ftp username
    password: str or NoneType
        ftp password
    """
    # attempt to connect to ftp host
    try:
        f = ftplib.FTP(HOST)
        f.login(username, password)
        f.voidcmd("NOOP")
    except IOError:
        raise RuntimeError('Check internet connection')
    except ftplib.error_perm:
        raise RuntimeError('Check login credentials')
    else:
        return True

# PURPOSE: list a directory on a ftp host
def ftp_list(
        HOST: str | list,
        username: str | None = None,
        password: str | None = None,
        timeout: int | None = None,
        basename: bool = False,
        pattern: str | None = None,
        sort: bool = False
    ):
    """
    List a directory on a ftp host

    Parameters
    ----------
    HOST: str or list
        remote ftp host path split as list
    username: str or NoneType
        ftp username
    password: str or NoneType
        ftp password
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    basename: bool, default False
        return the file or directory basename instead of the full path
    pattern: str or NoneType, default None
        regular expression pattern for reducing list
    sort: bool, default False
        sort output list

    Returns
    -------
    output: list
        items in a directory
    mtimes: list
        last modification times for items in the directory
    """
    # verify inputs for remote ftp host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try to connect to ftp host
    try:
        ftp = ftplib.FTP(HOST[0],timeout=timeout)
    except (socket.gaierror,IOError):
        raise RuntimeError(f'Unable to connect to {HOST[0]}')
    else:
        ftp.login(username,password)
        # list remote path
        output = ftp.nlst(posixpath.join(*HOST[1:]))
        # get last modified date of ftp files and convert into unix time
        mtimes = [None]*len(output)
        # iterate over each file in the list and get the modification time
        for i,f in enumerate(output):
            try:
                # try sending modification time command
                mdtm = ftp.sendcmd(f'MDTM {f}')
            except ftplib.error_perm:
                # directories will return with an error
                pass
            else:
                # convert the modification time into unix time
                mtimes[i] = get_unix_time(mdtm[4:], format="%Y%m%d%H%M%S")
        # reduce to basenames
        if basename:
            output = [posixpath.basename(i) for i in output]
        # reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(output) if re.search(pattern,f)]
            # reduce list of listed items and last modified times
            output = [output[indice] for indice in i]
            mtimes = [mtimes[indice] for indice in i]
        # sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(output), key=lambda i: i[1])]
            # sort list of listed items and last modified times
            output = [output[indice] for indice in i]
            mtimes = [mtimes[indice] for indice in i]
        # close the ftp connection
        ftp.close()
        # return the list of items and last modified times
        return (output, mtimes)

# PURPOSE: download a file from a ftp host
def from_ftp(
        HOST: str | list,
        username: str | None = None,
        password: str | None = None,
        timeout: int | None = None,
        local: str | pathlib.Path | None = None,
        hash: str = '',
        chunk: int = 8192,
        verbose: bool = False,
        fid=sys.stdout,
        mode: oct = 0o775
    ):
    """
    Download a file from a ftp host

    Parameters
    ----------
    HOST: str or list
        remote ftp host path
    username: str or NoneType
        ftp username
    password: str or NoneType
        ftp password
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    local: str, pathlib.Path or NoneType, default None
        path to local file
    hash: str, default ''
        MD5 hash of local file
    chunk: int, default 8192
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default sys.stdout
        open file object to print if verbose
    mode: oct, default 0o775
        permissions mode of output local file

    Returns
    -------
    remote_buffer: obj
        BytesIO representation of file
    """
    # create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    # verify inputs for remote ftp host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try downloading from ftp
    try:
        # try to connect to ftp host
        ftp = ftplib.FTP(HOST[0], timeout=timeout)
    except (socket.gaierror,IOError):
        raise RuntimeError(f'Unable to connect to {HOST[0]}')
    else:
        ftp.login(username,password)
        # remote path
        ftp_remote_path = posixpath.join(*HOST[1:])
        # copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        ftp.retrbinary(f'RETR {ftp_remote_path}',
            remote_buffer.write, blocksize=chunk)
        remote_buffer.seek(0)
        # save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        # generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        # get last modified date of remote file and convert into unix time
        mdtm = ftp.sendcmd(f'MDTM {ftp_remote_path}')
        remote_mtime = get_unix_time(mdtm[4:], format="%Y%m%d%H%M%S")
        # compare checksums
        if local and (hash != remote_hash):
            # convert to absolute path
            local = pathlib.Path(local).expanduser().absolute()
            # create directory if non-existent
            local.parent.mkdir(mode=mode, parents=True, exist_ok=True)
            # print file information
            args = (posixpath.join(*HOST), str(local))
            logging.info('{0} -->\n\t{1}'.format(*args))
            # store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with local.open(mode='wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            # change the permissions mode
            local.chmod(mode)
            # keep remote modification time of file and local access time
            os.utime(local, (local.stat().st_atime, remote_mtime))
        # close the ftp connection
        ftp.close()
        # return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

# default ssl context
_default_ssl_context = ssl.SSLContext(ssl.PROTOCOL_TLS)

# PURPOSE: check internet connection
def check_connection(HOST: str, context=_default_ssl_context):
    """
    Check internet connection with http host

    Parameters
    ----------
    HOST: str
        remote http host
    context: obj, default ssl.SSLContext(ssl.PROTOCOL_TLS)
        SSL context for ``urllib`` opener object
    """
    # attempt to connect to http host
    try:
        urllib2.urlopen(HOST, timeout=20, context=context)
    except urllib2.HTTPError as exc:
        logging.debug(exc.code)
        raise RuntimeError(exc.reason) from exc
    except urllib2.URLError as exc:
        logging.debug(exc.reason)
        raise RuntimeError('Check internet connection') from exc
    else:
        return True

# PURPOSE: list a directory on an Apache http Server
def http_list(
        HOST: str | list,
        timeout: int | None = None,
        context = _default_ssl_context,
        parser = lxml.etree.HTMLParser(),
        format: str = '%Y-%m-%d %H:%M',
        pattern: str = '',
        sort: bool = False
    ):
    """
    List a directory on an Apache http Server

    Parameters
    ----------
    HOST: str or list
        remote http host path
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default ssl.SSLContext(ssl.PROTOCOL_TLS)
        SSL context for ``urllib`` opener object
    parser: obj, default lxml.etree.HTMLParser()
        HTML parser for ``lxml``
    format: str, default '%Y-%m-%d %H:%M'
        format for input time string
    pattern: str, default ''
        regular expression pattern for reducing list
    sort: bool, default False
        sort output list

    Returns
    -------
    colnames: list
        column names in a directory
    collastmod: list
        last modification times for items in the directory
    """
    # verify inputs for remote http host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try listing from http
    try:
        # Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request, timeout=timeout, context=context)
    except urllib2.HTTPError as exc:
        logging.debug(exc.code)
        raise RuntimeError(exc.reason) from exc
    except urllib2.URLError as exc:
        logging.debug(exc.reason)
        msg = 'List error from {0}'.format(posixpath.join(*HOST))
        raise Exception(msg) from exc
    else:
        # read and parse request for files (column names and modified times)
        tree = lxml.etree.parse(response, parser)
        colnames = tree.xpath('//tr/td[not(@*)]//a/@href')
        # get the Unix timestamp value for a modification time
        collastmod = [get_unix_time(i,format=format)
            for i in tree.xpath('//tr/td[@align="right"][1]/text()')]
        # reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(colnames) if re.search(pattern, f)]
            # reduce list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        # sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(colnames), key=lambda i: i[1])]
            # sort list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        # return the list of column names and last modified times
        return (colnames, collastmod)

# PURPOSE: download a file from a http host
def from_http(
        HOST: str | list,
        timeout: int | None = None,
        context = _default_ssl_context,
        local: str | pathlib.Path | None = None,
        hash: str = '',
        chunk: int = 16384,
        verbose: bool = False,
        fid = sys.stdout,
        mode: oct = 0o775
    ):
    """
    Download a file from a http host

    Parameters
    ----------
    HOST: str or list
        remote http host path split as list
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default ssl.SSLContext(ssl.PROTOCOL_TLS)
        SSL context for ``urllib`` opener object
    local: str, pathlib.Path or NoneType, default None
        path to local file
    hash: str, default ''
        MD5 hash of local file
    chunk: int, default 16384
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default sys.stdout
        open file object to print if verbose
    mode: oct, default 0o775
        permissions mode of output local file

    Returns
    -------
    remote_buffer: obj
        BytesIO representation of file
    """
    # create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    # verify inputs for remote http host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try downloading from http
    try:
        # Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request, timeout=timeout, context=context)
    except:
        raise Exception('Download error from {0}'.format(posixpath.join(*HOST)))
    else:
        # copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        shutil.copyfileobj(response, remote_buffer, chunk)
        remote_buffer.seek(0)
        # save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        # generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        # compare checksums
        if local and (hash != remote_hash):
            # convert to absolute path
            local = pathlib.Path(local).expanduser().absolute()
            # create directory if non-existent
            local.parent.mkdir(mode=mode, parents=True, exist_ok=True)
            # print file information
            args = (posixpath.join(*HOST), str(local))
            logging.info('{0} -->\n\t{1}'.format(*args))
            # store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with local.open(mode='wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            # change the permissions mode
            local.chmod(mode)
        # return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

# PURPOSE: attempt to build an opener with netrc
def attempt_login(
        urs: str,
        context=_default_ssl_context,
        password_manager: bool = True,
        get_ca_certs: bool = True,
        redirect: bool = True,
        authorization_header: bool = False,
        **kwargs
    ):
    """
    attempt to build a urllib opener for NASA Earthdata

    Parameters
    ----------
    urs: str
        Earthdata login URS 3 host
    context: obj, default ssl.SSLContext(ssl.PROTOCOL_TLS)
        SSL context for ``urllib`` opener object
    password_manager: bool, default True
        Create password manager context using default realm
    get_ca_certs: bool, default True
        Get list of loaded “certification authority” certificates
    redirect: bool, default True
        Create redirect handler object
    authorization_header: bool, default False
        Add base64 encoded authorization header to opener
    username: str, default from environmental variable
        NASA Earthdata username
    password: str, default from environmental variable
        NASA Earthdata password
    retries: int, default 5
        number of retry attempts
    netrc: str, default ~/.netrc
        path to .netrc file for authentication

    Returns
    -------
    opener: obj
        OpenerDirector instance
    """
    # set default keyword arguments
    kwargs.setdefault('username', os.environ.get('EARTHDATA_USERNAME'))
    kwargs.setdefault('password', os.environ.get('EARTHDATA_PASSWORD'))
    kwargs.setdefault('retries', 5)
    kwargs.setdefault('netrc', os.path.expanduser('~/.netrc'))
    try:
        # only necessary on jupyterhub
        os.chmod(kwargs['netrc'], 0o600)
        # try retrieving credentials from netrc
        username, _, password = netrc.netrc(kwargs['netrc']).authenticators(urs)
    except Exception as exc:
        # try retrieving credentials from environmental variables
        username, password = (kwargs['username'], kwargs['password'])
        pass
    # if username or password are not available
    if not username:
        username = builtins.input(f'Username for {urs}: ')
    if not password:
        prompt = f'Password for {username}@{urs}: '
        password = getpass.getpass(prompt=prompt)
    # for each retry
    for retry in range(kwargs['retries']):
        # build an opener for urs with credentials
        opener = build_opener(username, password,
            context=context,
            password_manager=password_manager,
            get_ca_certs=get_ca_certs,
            redirect=redirect,
            authorization_header=authorization_header,
            urs=urs)
        # try logging in by check credentials
        try:
            check_credentials()
        except Exception as exc:
            pass
        else:
            return opener
        # reattempt login
        username = builtins.input(f'Username for {urs}: ')
        password = getpass.getpass(prompt=prompt)
    # reached end of available retries
    raise RuntimeError('End of Retries: Check NASA Earthdata credentials')

# PURPOSE: "login" to NASA Earthdata with supplied credentials
def build_opener(
        username: str,
        password: str,
        context=_default_ssl_context,
        password_manager: bool = True,
        get_ca_certs: bool = True,
        redirect: bool = True,
        authorization_header: bool = False,
        urs: str = 'https://urs.earthdata.nasa.gov'
    ):
    """
    Build ``urllib`` opener for NASA Earthdata with supplied credentials

    Parameters
    ----------
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    context: obj, default ssl.SSLContext(ssl.PROTOCOL_TLS)
        SSL context for ``urllib`` opener object
    password_manager: bool, default True
        Create password manager context using default realm
    get_ca_certs: bool, default True
        Get list of loaded “certification authority” certificates
    redirect: bool, default True
        Create redirect handler object
    authorization_header: bool, default False
        Add base64 encoded authorization header to opener
    urs: str, default 'https://urs.earthdata.nasa.gov'
        Earthdata login URS 3 host

    Returns
    -------
    opener: obj
        ``OpenerDirector`` instance
    """
    # https://docs.python.org/3/howto/urllib2.html#id5
    handler = []
    # create a password manager
    if password_manager:
        password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
        # Add the username and password for NASA Earthdata Login system
        password_mgr.add_password(None, urs, username, password)
        handler.append(urllib2.HTTPBasicAuthHandler(password_mgr))
    # Create cookie jar for storing cookies. This is used to store and return
    # the session cookie given to use by the data server (otherwise will just
    # keep sending us back to Earthdata Login to authenticate).
    cookie_jar = CookieJar()
    handler.append(urllib2.HTTPCookieProcessor(cookie_jar))
    # SSL context handler
    if get_ca_certs:
        context.get_ca_certs()
    handler.append(urllib2.HTTPSHandler(context=context))
    # redirect handler
    if redirect:
        handler.append(urllib2.HTTPRedirectHandler())
    # create "opener" (OpenerDirector instance)
    opener = urllib2.build_opener(*handler)
    # Encode username/password for request authorization headers
    # add Authorization header to opener
    if authorization_header:
        b64 = base64.b64encode(f'{username}:{password}'.encode())
        opener.addheaders = [("Authorization", f"Basic {b64.decode()}")]
    # Now all calls to urllib2.urlopen use our opener.
    urllib2.install_opener(opener)
    # All calls to urllib2.urlopen will now use handler
    # Make sure not to include the protocol in with the URL, or
    # HTTPPasswordMgrWithDefaultRealm will be confused.
    return opener

# PURPOSE: generate a NASA Earthdata user token
def get_token(
        HOST: str = 'https://urs.earthdata.nasa.gov/api/users/token',
        username: str | None = None,
        password: str | None = None,
        build: bool = True,
        urs: str = 'urs.earthdata.nasa.gov',
    ):
    """
    Generate a NASA Earthdata User Token

    Parameters
    ----------
    HOST: str or list
        NASA Earthdata token API host
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    build: bool, default True
        Build opener and check WebDAV credentials
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    urs: str, default 'urs.earthdata.nasa.gov'
        NASA Earthdata URS 3 host

    Returns
    -------
    token: dict
        JSON response with NASA Earthdata User Token
    """
    # attempt to build urllib2 opener and check credentials
    if build:
        attempt_login(urs,
            username=username,
            password=password,
            password_manager=False,
            get_ca_certs=False,
            redirect=False,
            authorization_header=True)
    # create post response with Earthdata token API
    try:
        request = urllib2.Request(HOST, method='POST')
        response = urllib2.urlopen(request)
    except urllib2.HTTPError as exc:
        logging.debug(exc.code)
        raise RuntimeError(exc.reason) from exc
    except urllib2.URLError as exc:
        logging.debug(exc.reason)
        raise RuntimeError('Check internet connection') from exc
    # read and return JSON response
    return json.loads(response.read())

# PURPOSE: generate a NASA Earthdata user token
def list_tokens(
        HOST: str = 'https://urs.earthdata.nasa.gov/api/users/tokens',
        username: str | None = None,
        password: str | None = None,
        build: bool = True,
        urs: str = 'urs.earthdata.nasa.gov',
    ):
    """
    List the current associated NASA Earthdata User Tokens

    Parameters
    ----------
    HOST: str
        NASA Earthdata list token API host
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    build: bool, default True
        Build opener and check WebDAV credentials
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    urs: str, default 'urs.earthdata.nasa.gov'
        NASA Earthdata URS 3 host

    Returns
    -------
    tokens: list
        JSON response with NASA Earthdata User Tokens
    """
    # attempt to build urllib2 opener and check credentials
    if build:
        attempt_login(urs,
            username=username,
            password=password,
            password_manager=False,
            get_ca_certs=False,
            redirect=False,
            authorization_header=True)
    # create get response with Earthdata list tokens API
    try:
        request = urllib2.Request(HOST)
        response = urllib2.urlopen(request)
    except urllib2.HTTPError as exc:
        logging.debug(exc.code)
        raise RuntimeError(exc.reason) from exc
    except urllib2.URLError as exc:
        logging.debug(exc.reason)
        raise RuntimeError('Check internet connection') from exc
    # read and return JSON response
    return json.loads(response.read())

# PURPOSE: revoke a NASA Earthdata user token
def revoke_token(
        token: str,
        HOST: str = f'https://urs.earthdata.nasa.gov/api/users/revoke_token',
        username: str | None = None,
        password: str | None = None,
        build: bool = True,
        urs: str = 'urs.earthdata.nasa.gov',
    ):
    """
    Generate a NASA Earthdata User Token

    Parameters
    ----------
    token: str
        NASA Earthdata token to be revoked
    HOST: str
        NASA Earthdata revoke token API host
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    build: bool, default True
        Build opener and check WebDAV credentials
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    urs: str, default 'urs.earthdata.nasa.gov'
        NASA Earthdata URS 3 host
    """
    # attempt to build urllib2 opener and check credentials
    if build:
        attempt_login(urs,
            username=username,
            password=password,
            password_manager=False,
            get_ca_certs=False,
            redirect=False,
            authorization_header=True)
    # full path for NASA Earthdata revoke token API
    url = f'{HOST}?token={token}'
    # create post response with Earthdata revoke tokens API
    try:
        request = urllib2.Request(url, method='POST')
        response = urllib2.urlopen(request)
    except urllib2.HTTPError as exc:
        logging.debug(exc.code)
        raise RuntimeError(exc.reason) from exc
    except urllib2.URLError as exc:
        logging.debug(exc.reason)
        raise RuntimeError('Check internet connection') from exc
    # verbose response
    logging.debug(f'Token Revoked: {token}')

# PURPOSE: check that entered NASA Earthdata credentials are valid
def check_credentials():
    """
    Check that entered NASA Earthdata credentials are valid
    """
    try:
        remote_path = posixpath.join('https://cddis.nasa.gov','archive')
        request = urllib2.Request(url=remote_path)
        response = urllib2.urlopen(request, timeout=20)
    except urllib2.HTTPError as exc:
        logging.debug(exc.code)
        raise RuntimeError(exc.reason) from exc
    except urllib2.URLError as exc:
        logging.debug(exc.reason)
        raise RuntimeError('Check internet connection') from exc
    else:
        return True

# PURPOSE: list a directory on GSFC CDDIS https server
def cddis_list(
        HOST: str | list,
        username: str | None = None,
        password: str | None = None,
        build: bool = True,
        timeout: int | None = None,
        urs: str = 'urs.earthdata.nasa.gov',
        parser=lxml.etree.HTMLParser(),
        pattern: str = '',
        sort: bool = False
    ):
    """
    List a directory on GSFC CDDIS archive server

    Parameters
    ----------
    HOST: str or list
        remote https host
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    build: bool, default True
        Build opener and check Earthdata credentials
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    urs: str, default 'urs.earthdata.nasa.gov'
        NASA Earthdata URS 3 host
    parser: obj, default lxml.etree.HTMLParser()
        HTML parser for ``lxml``
    pattern: str, default ''
        regular expression pattern for reducing list
    sort: bool, default False
        sort output list

    Returns
    -------
    colnames: list
        column names in a directory
    collastmod: list
        last modification times for items in the directory
    """
    # use netrc credentials
    if build and not (username or password):
        username,_,password = netrc.netrc().authenticators(urs)
    # build urllib2 opener and check credentials
    if build:
        # build urllib2 opener with credentials
        build_opener(username, password)
        # check credentials
        check_credentials()
    # verify inputs for remote https host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # Encode username/password for request authorization headers
    b64 = base64.b64encode(f'{username}:{password}'.encode())
    authorization_header = f"Basic {b64.decode()}"
    # try listing from https
    try:
        # Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        request.add_header("Authorization", authorization_header)
        tree = lxml.etree.parse(urllib2.urlopen(request, timeout=timeout), parser)
    except:
        raise Exception('List error from {0}'.format(posixpath.join(*HOST)))
    else:
        # read and parse request for files (column names and modified times)
        # find directories
        colnames = tree.xpath('//div[@class="archiveDir"]/div/a/text()')
        collastmod = [None]*(len(colnames))
        # find files
        colnames.extend(tree.xpath('//div[@class="archiveItem"]/div/a/text()'))
        # get the Unix timestamp value for a modification time
        collastmod.extend([get_unix_time(i[:19], format='%Y:%m:%d %H:%M:%S')
            for i in tree.xpath('//div[@class="archiveItem"]/div/span/text()')])
        # reduce using regular expression pattern
        if pattern:
            i = [i for i,f in enumerate(colnames) if re.search(pattern, f)]
            # reduce list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        # sort the list
        if sort:
            i = [i for i,j in sorted(enumerate(colnames), key=lambda i: i[1])]
            # sort list of column names and last modified times
            colnames = [colnames[indice] for indice in i]
            collastmod = [collastmod[indice] for indice in i]
        # return the list of column names and last modified times
        return (colnames, collastmod)

# PURPOSE: download a file from a GSFC CDDIS https server
def from_cddis(
        HOST: str | list,
        username: str | None = None,
        password: str | None = None,
        build: bool = True,
        timeout: int | None = None,
        urs: str = 'urs.earthdata.nasa.gov',
        local: str | pathlib.Path | None = None,
        hash: str = '',
        chunk: int = 16384,
        verbose: bool = False,
        fid=sys.stdout,
        mode: oct = 0o775
    ):
    """
    Download a file from GSFC CDDIS archive server

    Parameters
    ----------
    HOST: str or list
        remote https host
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    build: bool, default True
        Build opener and check Earthdata credentials
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    urs: str, default 'urs.earthdata.nasa.gov'
        NASA Earthdata URS 3 host
    local: str, pathlib.Path or NoneType, default None
        path to local file
    hash: str, default ''
        MD5 hash of local file
    chunk: int, default 16384
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    fid: obj, default sys.stdout
        open file object to print if verbose
    mode: oct, default 0o775
        permissions mode of output local file

    Returns
    -------
    remote_buffer: obj
        BytesIO representation of file
    """
    # create logger
    loglevel = logging.INFO if verbose else logging.CRITICAL
    logging.basicConfig(stream=fid, level=loglevel)
    # use netrc credentials
    if build and not (username or password):
        username,_,password = netrc.netrc().authenticators(urs)
    # build urllib2 opener and check credentials
    if build:
        # build urllib2 opener with credentials
        build_opener(username, password)
        # check credentials
        check_credentials()
    # verify inputs for remote https host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # Encode username/password for request authorization headers
    b64 = base64.b64encode(f'{username}:{password}'.encode())
    authorization_header = f"Basic {b64.decode()}"
    # try downloading from https
    try:
        # Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        request.add_header("Authorization", authorization_header)
        response = urllib2.urlopen(request, timeout=timeout)
    except:
        raise Exception('Download error from {0}'.format(posixpath.join(*HOST)))
    else:
        # copy remote file contents to bytesIO object
        remote_buffer = io.BytesIO()
        shutil.copyfileobj(response, remote_buffer, chunk)
        remote_buffer.seek(0)
        # save file basename with bytesIO object
        remote_buffer.filename = HOST[-1]
        # generate checksum hash for remote file
        remote_hash = hashlib.md5(remote_buffer.getvalue()).hexdigest()
        # compare checksums
        if local and (hash != remote_hash):
            # convert to absolute path
            local = pathlib.Path(local).expanduser().absolute()
            # create directory if non-existent
            local.parent.mkdir(mode=mode, parents=True, exist_ok=True)
            # print file information
            args = (posixpath.join(*HOST), str(local))
            logging.info('{0} -->\n\t{1}'.format(*args))
            # store bytes to file using chunked transfer encoding
            remote_buffer.seek(0)
            with local.open(mode='wb') as f:
                shutil.copyfileobj(remote_buffer, f, chunk)
            # change the permissions mode
            local.chmod(mode)
        # return the bytesIO object
        remote_buffer.seek(0)
        return remote_buffer

# PURPOSE: list a directory on IERS https Server
def iers_list(
        HOST: str | list,
        timeout: int | None = None,
        context = _default_ssl_context,
        parser = lxml.etree.HTMLParser()
    ):
    """
    List a directory on IERS Bulletin-A https server

    Parameters
    ----------
    HOST: str or list
        remote http host path
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default ssl.SSLContext(ssl.PROTOCOL_TLS)
        SSL context for ``urllib`` opener object
    parser: obj, default lxml.etree.HTMLParser()
        HTML parser for ``lxml``

    Returns
    -------
    colnames: list
        column names in a directory
    collastmod: list
        last modification times for items in the directory
    """
    # verify inputs for remote http host
    if isinstance(HOST, str):
        HOST = url_split(HOST)
    # try listing from http
    try:
        # Create and submit request.
        request = urllib2.Request(posixpath.join(*HOST))
        response = urllib2.urlopen(request, timeout=timeout, context=context)
    except urllib2.HTTPError as exc:
        logging.debug(exc.code)
        raise RuntimeError(exc.reason) from exc
    except urllib2.URLError as exc:
        logging.debug(exc.reason)
        msg = 'List error from {0}'.format(posixpath.join(*HOST))
        raise Exception(msg) from exc
    else:
        # read and parse request for files (column names and modified times)
        tree = lxml.etree.parse(response, parser)
        colnames = tree.xpath('//tr/td[@class="$tdclass"][4]//a/@href')
        # get the Unix timestamp value for a modification time
        collastmod = [get_unix_time(i,format='%Y-%m-%d')
            for i in tree.xpath('//tr/td[@class="$tdclass"][2]/span/text()')]
        # sort list of column names and last modified times in reverse order
        # return the list of column names and last modified times
        return (colnames[::-1], collastmod[::-1])

def from_jpl_ssd(
        kernel='de440s.bsp',
        timeout: int | None = None,
        context = _default_ssl_context,
        local: str | pathlib.Path | None = None,
        hash: str = '',
        chunk: int = 16384,
        verbose: bool = False,
        mode: oct = 0o775
    ):
    """
    Download `planetary ephemeride kernels`__ from the JPL Solar
    System Dynamics server

    .. __: https://ssd.jpl.nasa.gov/planets/eph_export.html

    Parameters
    ----------
    kernel: str
        JPL kernel file to download
    timeout: int or NoneType, default None
        timeout in seconds for blocking operations
    context: obj, default ssl.SSLContext(ssl.PROTOCOL_TLS)
        SSL context for ``urllib`` opener object
    hash: str, default ''
        MD5 hash of local file
    chunk: int, default 16384
        chunk size for transfer encoding
    verbose: bool, default False
        print file transfer information
    mode: oct, default 0o775
        permissions mode of output local file
    """
    # determine which kernel file to download
    if (local is None):
        # local path to kernel file
        local = get_data_path(['data',kernel])
    elif (kernel is None) and (local is not None):
        # verify inputs for remote http host
        local = pathlib.Path(local).expanduser().absolute()
        kernel = local.name
    # remote host path to kernel file
    HOST = ['https://ssd.jpl.nasa.gov','ftp','eph','planets','bsp',kernel]
    # get kernel file from remote host
    logging.info('Downloading JPL Planetary Ephemeride Kernel File')
    from_http(HOST, timeout=timeout, context=context, local=local,
        hash=hash, chunk=chunk, verbose=verbose, mode=mode)
