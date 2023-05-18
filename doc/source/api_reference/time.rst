====
time
====

Utilities for calculating time operations

 - Can convert delta time from seconds since an epoch to time since a different epoch
 - Can calculate the time in days since epoch from calendar dates
 - Calculates the difference between dynamic time and universal time (`TT` - `UT1`) following Richard Ray's ``PERTH3`` algorithms
 - Can count the number of leap seconds between a given GPS time and UTC
 - Syncs leap second files with NIST servers
 - Updates differences between universal time (`UT`) and dynamic time (`TT`)

Calling Sequence
----------------

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

.. autofunction:: pyTMD.time.parse

.. autofunction:: pyTMD.time.parse_date_string

.. autofunction:: pyTMD.time.split_date_string

.. autofunction:: pyTMD.time.datetime_to_list

.. autofunction:: pyTMD.time.calendar_days

.. autofunction:: pyTMD.time.convert_datetime

.. autofunction:: pyTMD.time.convert_delta_time

.. autofunction:: pyTMD.time.convert_calendar_dates

.. autofunction:: pyTMD.time.convert_calendar_decimal

.. autofunction:: pyTMD.time.convert_julian

.. autoclass:: pyTMD.time.timescale
   :members:

.. autofunction:: pyTMD.time.interpolate_delta_time

.. autofunction:: pyTMD.time.count_leap_seconds

.. autofunction:: pyTMD.time.get_leap_seconds

.. autofunction:: pyTMD.time.update_leap_seconds

.. autofunction:: pyTMD.time.merge_delta_time

.. autofunction:: pyTMD.time.append_delta_time

.. autofunction:: pyTMD.time.merge_bulletin_a_files

.. autofunction:: pyTMD.time.iers_ftp_delta_time

.. autofunction:: pyTMD.time.iers_delta_time

.. autofunction:: pyTMD.time.cddis_delta_time

.. autofunction:: pyTMD.time.read_iers_bulletin_a

.. autofunction:: pyTMD.time.update_bulletin_a

.. autofunction:: pyTMD.time.pull_deltat_file
