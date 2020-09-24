=======
time.py
=======

Utilities for calculating time operations

 - Can convert delta time from seconds since an epoch to time since a different epoch
 - Can calculate the time in days since epoch from calendar dates
 - Can count the number of leap seconds between a given GPS time and UTC
 - Syncs leap second files with NIST servers
 - Updates differences between universal time (UT) and dynamic time (TT)

Calling Sequence
================

Count the number of leap seconds between a GPS time and UTC

.. code-block:: python

    import pyTMD.time
    leap_seconds = pyTMD.time.count_leap_seconds(gps_seconds)

Convert a time from seconds since 1980-01-06T00:00:00 to Modified Julian Days (MJD)

.. code-block:: python

    import pyTMD.time
    MJD = pyTMD.time.convert_delta_time(delta_time, epoch1=(1980,1,6,0,0,0),
        epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)

Convert a calendar date into Modified Julian Days

.. code-block:: python

    import pyTMD.time
    MJD = pyTMD.time.convert_calendar_dates(YEAR,MONTH,DAY,hour=HOUR,
        minute=MINUTE,second=SECOND,epoch=(1858,11,17,0,0,0))

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/time.py


General Methods
===============


.. method:: pyTMD.time.convert_delta_time(delta_time, epoch1=None, epoch2=None, scale=1.0)

    Convert delta time from seconds since epoch1 to time since epoch2

    Arguments:

        `delta_time`: seconds since epoch1

    Keyword arguments:

        `epoch1`: epoch for input delta_time

        `epoch2`: epoch for output delta_time

        `scale`: scaling factor for converting time to output units


.. method:: pyTMD.time.convert_calendar_dates(year, month, day, hour=0.0, minute=0.0, second=0.0, epoch=None, scale=1.0)

    Calculate the time in time units since epoch from calendar dates

    Arguments:

        `year`: calendar month

        `month`: month of the year

        `day`: day of the month

    Keyword arguments:

        `hour`: hour of the day

        `minute`: minute of the hour

        `second`: second of the minute

        `epoch`: epoch for output delta_time

        `scale`: scaling factor for converting time to output units


.. method:: pyTMD.time.count_leap_seconds(GPS_Time)

    Counts the number of leap seconds between a given GPS time and UTC

    Arguments:

        `GPS_Time`: seconds since January 6, 1980 at 00:00:00


.. method:: pyTMD.time.get_leap_seconds()

    Gets a list of GPS times for when leap seconds occurred


.. method:: pyTMD.time.update_leap_seconds(verbose=False, mode=0o775)

    Connects to servers to download leap-seconds.list files from `NIST servers`__

.. __: ftp://ftp.nist.gov/pub/time/leap-seconds.list

    Keyword arguments:

        `verbose`: print file information about output file

        `mode`: permissions mode of output file


.. method:: pyTMD.time.merge_delta_time(username=None, password=None, verbose=False, mode=0o775)

    Connects to servers to download `differences between dynamic and universal time`__

.. __: https://www.usno.navy.mil/USNO/earth-orientation/eo-products/long-term

    Reads IERS Bulletin-A produced iers_deltat.data files

    Creates a merged file combining the historic, monthly and daily files

    Keyword arguments:

        `username`: NASA Earthdata username

        `password`: NASA Earthdata password

        `verbose`: print file information about output file

        `mode`: permissions mode of output file


.. method:: pyTMD.time.merge_bulletin_a_files(username=None, password=None, verbose=False, mode=0o775)

    Attempt to connects to the IERS server and the CDDIS Earthdata server to download and merge Bulletin-A files

    Reads the IERS Bulletin-A files and calculates the daily delta times

    Keyword arguments:

        `username`: NASA Earthdata username

        `password`: NASA Earthdata password

        `verbose`: print file information about output file

        `mode`: permissions mode of output file


.. method:: pyTMD.time.iers_delta_time(daily_file, verbose=False, mode=0o775)

    Connects to the IERS server to download `Bulletin-A files`__

.. __: ftp://ftp.iers.org/products/eop/rapid/bulletina

    Reads the IERS Bulletin-A files and calculates the daily delta times

    Arguments:

        `daily_file`: output daily delta time file from merged Bulletin-A files

    Keyword arguments:

        `verbose`: print file information about output file

        `mode`: permissions mode of output file


.. method:: pyTMD.time.cddis_delta_time(daily_file, username=None, password=None, verbose=False, mode=0o775)

    Connects to the NASA GSFC CDDIS https server to download `Bulletin-A files`__

.. __: https://cddis.nasa.gov/archive/products/iers/iers_bulletins/bulletin_a

    Reads the IERS Bulletin-A files and calculates the daily delta times

    Arguments:

        `daily_file`: output daily delta time file from merged Bulletin-A files

    Keyword arguments:

        `username`: NASA Earthdata username

        `password`: NASA Earthdata password

        `verbose`: print file information about output file

        `mode`: permissions mode of output file


.. method:: pyTMD.time.read_iers_bulletin_a(fileID)

    Read a weekly `IERS Bulletin-A file`__ and calculate the daily delta times (TT - UT1)

.. __: https://datacenter.iers.org/productMetadata.php?id=6

    Arguments:

        fileID: open file object for Bulletin-A file


.. method:: pyTMD.time.pull_deltat_file(FILE, username=None, password=None, verbose=False, mode=0o775)

    Connects to `servers`__ and downloads monthly or historic delta time files

.. __: ftp://cddis.nasa.gov/products/iers/

    Arguments:

        `FILE`: delta time file to download from remote servers

    Keyword arguments:

        `username`: NASA Earthdata username

        `password`: NASA Earthdata password

        `verbose`: print file information about output file

        `mode`: permissions mode of output file
