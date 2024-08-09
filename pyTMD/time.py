#!/usr/bin/env python
u"""
time.py
Written by Tyler Sutterley (04/2024)
Utilities for calculating time operations

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml
    timescale: Python tools for time and astronomical calculations
        https://pypi.org/project/timescale/

UPDATE HISTORY:
    Updated 04/2024: deprecated in favor of timescale.time
    Updated 02/2024: move the immutable parameters in timescale class
    Updated 10/2023: add function to output timescale to string arrays
    Updated 06/2023: improve conversion of timescale to datetime arrays
    Updated 05/2023: add timescale class for converting between time scales
        added timescale to_datetime function to create datetime arrays
        allow epoch arguments to be numpy datetime64 variables or strings
        function to convert a string with time zone information to datetime
    Updated 04/2023: using pathlib to define and expand paths
    Updated 03/2023: add basic variable typing to function inputs
    Updated 12/2022: added interpolation for delta time (TT - UT1)
        output variables for some standard epochs used within tide programs
    Updated 11/2022: use IERS https server as default for Bulletin-A files
        added download function for latest Bulletin-A file from IERS
        added function to append from existing merged delta time file
        use f-strings for formatting verbose or ascii output
    Updated 10/2022: added encoding for reading leap seconds ascii files
    Updated 08/2022: output variables to unit conversion to seconds
        and the number of days per month for both leap and standard years
    Updated 05/2022: changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 04/2021: updated NIST ftp server url for leap-seconds.list
    Updated 03/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 02/2021: NASA CDDIS anonymous ftp access discontinued
    Updated 01/2021: added ftp connection checks
        add date parser for cases when only a calendar date with no units
    Updated 12/2020: merged with convert_julian and convert_calendar_decimal
        added calendar_days routine to get number of days per month
    Updated 09/2020: added wrapper function for merging Bulletin-A files
        can parse date strings in form "time-units since yyyy-mm-dd hh:mm:ss"
    Updated 08/2020: added NASA Earthdata routines for downloading from CDDIS
    Written 07/2020
"""
from __future__ import annotations

import logging
import warnings
import numpy as np
import timescale.time

# conversion factors between time units and seconds
_to_sec = {'microseconds': 1e-6, 'microsecond': 1e-6,
           'microsec': 1e-6, 'microsecs': 1e-6,
           'milliseconds': 1e-3, 'millisecond': 1e-3,
           'millisec': 1e-3, 'millisecs': 1e-3,
           'msec': 1e-3, 'msecs': 1e-3, 'ms': 1e-3,
           'seconds': 1.0, 'second': 1.0, 'sec': 1.0,
           'secs': 1.0, 's': 1.0,
           'minutes': 60.0, 'minute': 60.0,
           'min': 60.0, 'mins': 60.0,
           'hours': 3600.0, 'hour': 3600.0,
           'hr': 3600.0, 'hrs': 3600.0, 'h': 3600.0,
           'day': 86400.0, 'days': 86400.0, 'd': 86400.0}

# standard epochs
_mjd_epoch = (1858, 11, 17, 0, 0, 0)
_ntp_epoch = (1900, 1, 1, 0, 0, 0)
_unix_epoch = (1970, 1, 1, 0, 0, 0)
_gps_epoch = (1980, 1, 6, 0, 0, 0)
_tide_epoch = (1992, 1, 1, 0, 0, 0)
_j2000_epoch = (2000, 1, 1, 12, 0, 0)
_atlas_sdp_epoch = (2018, 1, 1, 0, 0, 0)

# PURPOSE: parse a date string and convert to a datetime object in UTC
def parse(date_string: str):
    """
    Parse a date string and convert to a naive ``datetime`` object in UTC

    Parameters
    ----------
    date_string: str
        formatted time string

    Returns
    -------
    date: obj
        output ``datetime`` object
    """
    warnings.warn("Deprecated. Please use timescale.time.parse",
        DeprecationWarning)
    return timescale.time.parse(date_string)

# PURPOSE: parse a date string into epoch and units scale
def parse_date_string(date_string: str):
    """
    Parse a date string of the form

    - time-units since ``yyyy-mm-dd hh:mm:ss``
    - ``yyyy-mm-dd hh:mm:ss`` for exact calendar dates

    Parameters
    ----------
    date_string: str
        time-units since yyyy-mm-dd hh:mm:ss

    Returns
    -------
    epoch: list
        epoch of ``delta_time``
    conversion_factor: float
        multiplication factor to convert to seconds
    """
    warnings.warn("Deprecated. Please use timescale.time.parse_date_string",
        DeprecationWarning)
    return timescale.time.parse_date_string(date_string)

# PURPOSE: split a date string into units and epoch
def split_date_string(date_string: str):
    """
    Split a date string into units and epoch

    Parameters
    ----------
    date_string: str
        time-units since yyyy-mm-dd hh:mm:ss
    """
    warnings.warn("Deprecated. Please use timescale.time.split_date_string",
        DeprecationWarning)
    return timescale.time.split_date_string(date_string)

# PURPOSE: convert a datetime object into a list
def datetime_to_list(date):
    """
    Convert a ``datetime`` object into a list

    Parameters
    ----------
    date: obj
        Input ``datetime`` object to convert

    Returns
    -------
    date: list
        [year,month,day,hour,minute,second]
    """
    warnings.warn("Deprecated. Please use timescale.time.datetime_to_list",
        DeprecationWarning)
    return timescale.time.datetime_to_list(date)

# days per month in a leap and a standard year
# only difference is February (29 vs. 28)
_dpm_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
_dpm_stnd = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

# PURPOSE: gets the number of days per month for a given year
def calendar_days(year: int | float | np.ndarray) -> np.ndarray:
    """
    Calculates the number of days per month for a given year

    Parameters
    ----------
    year: int, float or np.ndarray
        calendar year

    Returns
    -------
    dpm: np.ndarray
        number of days for each month
    """
    warnings.warn("Deprecated. Please use timescale.time.calendar_days",
        DeprecationWarning)
    return timescale.time.calendar_days(year)

# PURPOSE: convert a numpy datetime array to delta times since an epoch
def convert_datetime(*args, **kwargs):
    """
    Convert a ``numpy`` ``datetime`` array to seconds since ``epoch``

    Parameters
    ----------
    date: np.ndarray
        ``numpy`` ``datetime`` array
    epoch: str, tuple, list, np.ndarray, default (1970,1,1,0,0,0)
        epoch for output ``delta_time``

    Returns
    -------
    delta_time: float
        seconds since epoch
    """
    warnings.warn("Deprecated. Please use timescale.time.convert_datetime",
        DeprecationWarning)
    return timescale.time.convert_datetime(*args, **kwargs)

# PURPOSE: convert times from seconds since epoch1 to time since epoch2
def convert_delta_time(*args, **kwargs):
    """
    Convert delta time from seconds since ``epoch1`` to time since ``epoch2``

    Parameters
    ----------
    delta_time: np.ndarray
        seconds since epoch1
    epoch1: str, tuple, list or NoneType, default None
        epoch for input ``delta_time``
    epoch2: str, tuple, list or NoneType, default None
        epoch for output ``delta_time``
    scale: float, default 1.0
        scaling factor for converting time to output units
    """
    warnings.warn("Deprecated. Please use timescale.time.convert_delta_time",
        DeprecationWarning)
    return timescale.time.convert_delta_time(*args, **kwargs)

# PURPOSE: calculate the delta time from calendar date
# http://scienceworld.wolfram.com/astronomy/JulianDate.html
def convert_calendar_dates(*args, **kwargs) -> np.ndarray:
    """
    Calculate the time in units since ``epoch`` from calendar dates

    Parameters
    ----------
    year: np.ndarray
        calendar year
    month: np.ndarray
        month of the year
    day: np.ndarray
        day of the month
    hour: np.ndarray or float, default 0.0
        hour of the day
    minute: np.ndarray or float, default 0.0
        minute of the hour
    second: np.ndarray or float, default 0.0
        second of the minute
    epoch: tuple or list, default pyTMD.time._tide_epoch
        epoch for output ``delta_time``
    scale: float, default 1.0
        scaling factor for converting time to output units

    Returns
    -------
    delta_time: np.ndarray
        time since epoch
    """
    warnings.warn("Deprecated. Please use timescale.time.convert_calendar_dates",
        DeprecationWarning)
    return timescale.time.convert_calendar_dates(*args, **kwargs)

# PURPOSE: Converts from calendar dates into decimal years
def convert_calendar_decimal(*args, **kwargs) -> np.ndarray:
    """
    Converts from calendar date into decimal years taking into
    account leap years

    Parameters
    ----------
    year: np.ndarray
        calendar year
    month: np.ndarray
        calendar month
    day: np.ndarray, float or NoneType, default None
        day of the month
    hour: np.ndarray, float or NoneType, default None
        hour of the day
    minute: np.ndarray, float or NoneType, default None
        minute of the hour
    second: np.ndarray, float or NoneType, default None
        second of the minute
    DofY: np.ndarray, float or NoneType, default None
        day of the year

    Returns
    -------
    t_date: np.ndarray
        date in decimal-year format

    References
    ----------
    .. [1] N. Dershowitz, and E. M. Reingold.
        *Calendrical Calculations*,
        Cambridge: Cambridge University Press, (2008).
    """
    warnings.warn("Deprecated. Please use timescale.time.convert_calendar_decimal",
        DeprecationWarning)
    return timescale.time.convert_calendar_decimal(*args, **kwargs)

# PURPOSE: Converts from Julian day to calendar date and time
def convert_julian(*args, **kwargs):
    """
    Converts from Julian day to calendar date and time

    Parameters
    ----------
    JD: np.ndarray
        Julian Day (days since 01-01-4713 BCE at 12:00:00)
    astype: str or NoneType, default None
        convert output to variable type
    format: str, default 'dict'
        format of output variables

            - ``'dict'``: dictionary with variable keys
            - ``'tuple'``: tuple in most-to-least-significant order
            - ``'zip'``: aggregated variable sets

    Returns
    -------
    year: np.ndarray
        calendar year
    month: np.ndarray
        calendar month
    day: np.ndarray
        day of the month
    hour: np.ndarray
        hour of the day
    minute: np.ndarray
        minute of the hour
    second: np.ndarray
        second of the minute

    References
    ----------
    .. [1] W. H. Press, *Numerical Recipes in C*,
        Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
        Cambridge University Press, (1988).
    .. [2] D. A. Hatcher, "Simple Formulae for Julian Day Numbers and
        Calendar Dates", *Quarterly Journal of the Royal Astronomical
        Society*, 25(1), 1984.
    """
    warnings.warn("Deprecated. Please use timescale.time.convert_calendar_decimal",
        DeprecationWarning)
    return timescale.time.convert_julian(*args, **kwargs)

# delta time (TT - UT1) file
_delta_file = timescale.utilities.get_data_path(['data','merged_deltat.data'])

class timescale(timescale.time.Timescale):
    """
    Class for converting between time scales

    Attributes
    ----------
    leaps: np.ndarray
        Number of leap seconds
    MJD: np.ndarray
        Modified Julian Days
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
