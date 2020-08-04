============
utilities.py
============

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can download a file from a ftp or http host
 - Checks MD5 hashes between local and remote files

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/master/pyTMD/utilities.py


General Methods
===============


.. method:: pyTMD.utilities.get_data_path(relpath)

    Get the absolute path within a package from a relative path

    Arguments:

        `relpath`: local relative path as list or string


.. method:: pyTMD.utilities.get_hash(local)

    Get the MD5 hash value from a local file

    Arguments:

        `local`: path to file


.. method:: pyTMD.utilities.ftp_list(HOST,timeout=None,basename=False,pattern=None,sort=False)

    List a directory on a ftp host

    Arguments:

        `HOST`: remote ftp host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `basename`: return the file or directory basename instead of the full path

        `pattern`: regular expression pattern for reducing list

        `sort`: sort output list


.. method:: pyTMD.utilities.from_ftp(HOST,timeout=None,local=None,hash='',chunk=16384,verbose=False,mode=0o775)

    Download a file from a ftp host

    Arguments:

        `HOST`: remote ftp host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `mode`: permissions mode of output local file


.. method:: pyTMD.utilities.from_http(HOST,timeout=None,local=None,hash='',chunk=16384,verbose=False,mode=0o775)

    Download a file from a http host

    Arguments:

        `HOST`: remote http host path split as list

    Keyword arguments:

        `timeout`: timeout in seconds for blocking operations

        `local`: path to local file

        `hash`: MD5 hash of local file

        `chunk`: chunk size for transfer encoding

        `verbose`: print file transfer information

        `mode`: permissions mode of output local file
