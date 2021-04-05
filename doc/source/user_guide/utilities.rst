============
utilities.py
============

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can download a file from a ftp or http host
 - Can download a file from CDDIS via https when NASA Earthdata credentials are supplied
 - Checks ``MD5`` or ``sha1`` hashes between local and remote files

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/utilities.py


General Methods
===============


.. method:: pyTMD.utilities.get_data_path(relpath)

    Get the absolute path within a package from a relative path

    Arguments:

        ``relpath``: local relative path as list or string


.. method:: pyTMD.utilities.get_hash(local, algorithm='MD5')

    Get the hash value from a local file or BytesIO object

    Arguments:

        ``local``: BytesIO object or path to file

    Keyword arguments:

        ``algorithm``: hashing algorithm for checksum validation

            ``'MD5'``: Message Digest

            ``'sha1'``: Secure Hash Algorithm


.. method:: pyTMD.utilities.url_split(s)

    Recursively split a url path into a list

    Arguments:

        ``s``: url string


.. method:: pyTMD.utilities.roman_to_int(local)

    Converts a string from Roman numerals into an integer (Arabic)

    Arguments:

        ``roman``: Roman numeral string


.. method:: pyTMD.utilities.get_unix_time(time_string, format='%Y-%m-%d %H:%M:%S')

    Get the Unix timestamp value for a formatted date string

    Arguments:

        ``time_string``: formatted time string to parse

    Keyword arguments:

        ``format``: format for input time string


.. method:: gravity_toolkit.utilities.even(value)

    Rounds a number to an even number less than or equal to original

    Arguments:

        ``value``: number to be rounded


.. method:: pyTMD.utilities.copy(source, destination, verbose=False, move=False)

    Copy or move a file with all system information

    Arguments:

        ``source``: source file

        ``destination``: copied destination file

    Keyword arguments:

        ``verbose``: print file transfer information

        ``move``: remove the source file


.. method:: pyTMD.utilities.check_ftp_connection(HOST,username=None,password=None)

    Check internet connection with ftp host

    Arguments:

        ``HOST``: remote ftp host

    Keyword arguments:

        ``username``: ftp username

        ``password``: ftp password


.. method:: pyTMD.utilities.ftp_list(HOST,username=None,password=None,timeout=None,basename=False,pattern=None,sort=False)

    List a directory on a ftp host

    Arguments:

        ``HOST``: remote ftp host path split as list

    Keyword arguments:

        ``username``: ftp username

        ``password``: ftp password

        ``timeout``: timeout in seconds for blocking operations

        ``basename``: return the file or directory basename instead of the full path

        ``pattern``: regular expression pattern for reducing list

        ``sort``: sort output list

    Returns:

        ``output``: list of items in a directory

        ``mtimes``: list of last modification times for items in the directory


.. method:: pyTMD.utilities.from_ftp(HOST,username=None,password=None,timeout=None,local=None,hash='',chunk=8192,verbose=False,fid=sys.stdout,mode=0o775)

    Download a file from a ftp host

    Arguments:

        ``HOST``: remote ftp host path split as list

    Keyword arguments:

        ``username``: ftp username

        ``password``: ftp password

        ``timeout``: timeout in seconds for blocking operations

        ``local``: path to local file

        ``hash``: MD5 hash of local file

        ``chunk``: chunk size for transfer encoding

        ``verbose``: print file transfer information

        ``fid``: open file object to print if verbose

        ``mode``: permissions mode of output local file


.. method:: pyTMD.utilities.check_connection(HOST)

    Check internet connection

    Arguments:

        ``HOST``: remote http host


.. method:: pyTMD.utilities.from_http(HOST,timeout=None,context=ssl.SSLContext(),local=None,hash='',chunk=16384,verbose=False,fid=sys.stdout,mode=0o775)

    Download a file from a http host

    Arguments:

        ``HOST``: remote http host path split as list

    Keyword arguments:

        ``timeout``: timeout in seconds for blocking operations

        ``context``: SSL context for url opener object

        ``local``: path to local file

        ``hash``: MD5 hash of local file

        ``chunk``: chunk size for transfer encoding

        ``verbose``: print file transfer information

        ``fid``: open file object to print if verbose

        ``mode``: permissions mode of output local file


.. method:: pyTMD.utilities.build_opener(username, password, context=ssl.SSLContext(ssl.PROTOCOL_TLS), password_manager=True, get_ca_certs=True, redirect=True, authorization_header=False, urs='https://urs.earthdata.nasa.gov')

    build urllib opener for NASA Earthdata with supplied credentials

    Arguments:

        ``username``: NASA Earthdata username

        ``password``: NASA Earthdata password

    Keyword arguments:

        ``context``: SSL context for opener object

        ``password_manager``: create password manager context using default realm

        ``get_ca_certs``: get list of loaded “certification authority” certificates

        ``redirect``: create redirect handler object

        ``authorization_header``: add base64 encoded authorization header to opener

        ``urs``: Earthdata login URS 3 host


.. method:: pyTMD.utilities.check_credentials()

    Check that entered NASA Earthdata credentials are valid


.. method:: pyTMD.utilities.cddis_list(HOST,username=None,password=None,build=True,timeout=None,parser=None,pattern='',sort=False)

    Download a file from a NASA GSFC CDDIS https server

    Arguments:

        ``HOST``: remote http host path split as list

    Keyword arguments:

        ``username``: NASA Earthdata username

        ``password``: NASA Earthdata password

        ``build``: Build opener and check NASA Earthdata password

        ``timeout``: timeout in seconds for blocking operations

        ``parser``: HTML parser for lxml

        ``pattern``: regular expression pattern for reducing list

        ``sort``: sort output list

    Returns:

        ``colnames``: list of column names in a directory

        ``collastmod``: list of last modification times for items in the directory


.. method:: pyTMD.utilities.from_cddis(HOST,username=None,password=None,build=True,timeout=None,local=None,hash='',chunk=16384,verbose=False,fid=sys.stdout,mode=0o775)

    Download a file from a NASA GSFC CDDIS https server

    Arguments:

        ``HOST``: remote http host path split as list

    Keyword arguments:

        ``username``: NASA Earthdata username

        ``password``: NASA Earthdata password

        ``build``: Build opener and check NASA Earthdata password

        ``timeout``: timeout in seconds for blocking operations

        ``local``: path to local file

        ``hash``: MD5 hash of local file

        ``chunk``: chunk size for transfer encoding

        ``verbose``: print file transfer information

        ``fid``: open file object to print if verbose

        ``mode``: permissions mode of output local file
