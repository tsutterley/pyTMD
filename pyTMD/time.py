#!/usr/bin/env python
u"""
time.py
Written by Tyler Sutterley (11/2022)
Utilities for calculating time operations

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    lxml: processing XML and HTML in Python
        https://pypi.python.org/pypi/lxml

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
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
import os
import re
import copy
import logging
import warnings
import datetime
import traceback
import numpy as np
import dateutil.parser
import pyTMD.utilities

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

# PURPOSE: parse a date string into epoch and units scale
def parse_date_string(date_string):
    """
    parse a date string of the form

    - time-units since ``yyyy-mm-dd hh:mm:ss``
    - ``yyyy-mm-dd hh:mm:ss`` for exact calendar dates

    Parameters
    ----------
    date_string: str
        time-units since yyyy-mm-dd hh:mm:ss

    Returns
    -------
    epoch: list
        epoch of delta time
    conversion_factor: float
        multiplication factor to convert to seconds
    """
    # try parsing the original date string as a date
    try:
        epoch = dateutil.parser.parse(date_string)
    except ValueError:
        pass
    else:
        # return the epoch (as list)
        return (datetime_to_list(epoch),0.0)
    # split the date string into units and epoch
    units,epoch = split_date_string(date_string)
    if units not in _to_sec.keys():
        raise ValueError(f'Invalid units: {units}')
    # return the epoch (as list) and the time unit conversion factors
    return (datetime_to_list(epoch), _to_sec[units])

# PURPOSE: split a date string into units and epoch
def split_date_string(date_string):
    """
    split a date string into units and epoch

    Parameters
    ----------
    date_string: str
        time-units since yyyy-mm-dd hh:mm:ss
    """
    try:
        units,_,epoch = date_string.split(None, 2)
    except ValueError:
        raise ValueError(f'Invalid format: {date_string}')
    else:
        return (units.lower(), dateutil.parser.parse(epoch))

# PURPOSE: convert a datetime object into a list
def datetime_to_list(date):
    """
    convert a datetime object into a list

    Parameters
    ----------
    date: datetime object

    Returns
    -------
    date: list
        [year,month,day,hour,minute,second]
    """
    return [date.year, date.month, date.day,
            date.hour, date.minute, date.second]

# days per month in a leap and a standard year
# only difference is February (29 vs. 28)
_dpm_leap = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
_dpm_stnd = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

# PURPOSE: gets the number of days per month for a given year
def calendar_days(year):
    """
    Calculates the number of days per month for a given year

    Parameters
    ----------
    year: int or float
        calendar year

    Returns
    -------
    dpm: list
        number of days for each month
    """
    # Rules in the Gregorian calendar for a year to be a leap year:
    # divisible by 4, but not by 100 unless divisible by 400
    # True length of the year is about 365.2422 days
    # Adding a leap day every four years ==> average 365.25
    # Subtracting a leap year every 100 years ==> average 365.24
    # Adding a leap year back every 400 years ==> average 365.2425
    # Subtracting a leap year every 4000 years ==> average 365.24225
    m4 = (year % 4)
    m100 = (year % 100)
    m400 = (year % 400)
    m4000 = (year % 4000)
    # find indices for standard years and leap years using criteria
    if ((m4 == 0) & (m100 != 0) | (m400 == 0) & (m4000 != 0)):
        return np.array(_dpm_leap, dtype=np.float64)
    elif ((m4 != 0) | (m100 == 0) & (m400 != 0) | (m4000 == 0)):
        return np.array(_dpm_stnd, dtype=np.float64)

# PURPOSE: convert a numpy datetime array to delta times from the UNIX epoch
def convert_datetime(date, epoch=(1970, 1, 1, 0, 0, 0)):
    """
    Convert a numpy datetime array to seconds since ``epoch``

    Parameters
    ----------
    date: obj
        numpy datetime array
    epoch: tuple, default (1970,1,1,0,0,0)
        epoch for output delta_time

    Returns
    -------
    delta_time: float
        seconds since epoch
    """
    epoch = datetime.datetime(*epoch)
    return (date - np.datetime64(epoch)) / np.timedelta64(1, 's')

# PURPOSE: convert times from seconds since epoch1 to time since epoch2
def convert_delta_time(delta_time, epoch1=None, epoch2=None, scale=1.0):
    """
    Convert delta time from seconds since ``epoch1`` to time since ``epoch2``

    Parameters
    ----------
    delta_time: float
        seconds since epoch1
    epoch1: tuple or NoneType, default None
        epoch for input delta_time
    epoch2: tuple or NoneType, default None
        epoch for output delta_time
    scale: float, default 1.0
        scaling factor for converting time to output units
    """
    epoch1 = datetime.datetime(*epoch1)
    epoch2 = datetime.datetime(*epoch2)
    delta_time_epochs = (epoch2 - epoch1).total_seconds()
    # subtract difference in time and rescale to output units
    return scale*(delta_time - delta_time_epochs)

# PURPOSE: calculate the delta time from calendar date
# http://scienceworld.wolfram.com/astronomy/JulianDate.html
def convert_calendar_dates(year, month, day, hour=0.0, minute=0.0, second=0.0,
    epoch=(1992, 1, 1, 0, 0, 0), scale=1.0):
    """
    Calculate the time in time units since ``epoch`` from calendar dates

    Parameters
    ----------
    year: float
        calendar year
    month: float
        month of the year
    day: float
        day of the month
    hour: float, default 0.0
        hour of the day
    minute: float, default 0.0
        minute of the hour
    second: float, default 0.0
        second of the minute
    epoch: tuple, default (1992,1,1,0,0,0)
        epoch for output delta_time
    scale: float, default 1.0
        scaling factor for converting time to output units

    Returns
    -------
    delta_time: float
        days since epoch
    """
    # calculate date in Modified Julian Days (MJD) from calendar date
    # MJD: days since November 17, 1858 (1858-11-17T00:00:00)
    MJD = 367.0*year - np.floor(7.0*(year + np.floor((month+9.0)/12.0))/4.0) - \
        np.floor(3.0*(np.floor((year + (month - 9.0)/7.0)/100.0) + 1.0)/4.0) + \
        np.floor(275.0*month/9.0) + day + hour/24.0 + minute/1440.0 + \
        second/86400.0 + 1721028.5 - 2400000.5
    epoch1 = datetime.datetime(1858, 11, 17, 0, 0, 0)
    epoch2 = datetime.datetime(*epoch)
    delta_time_epochs = (epoch2 - epoch1).total_seconds()
    # return the date in days since epoch
    return scale*np.array(MJD - delta_time_epochs/86400.0, dtype=np.float64)

# PURPOSE: Converts from calendar dates into decimal years
def convert_calendar_decimal(year, month, day=None, hour=None, minute=None,
    second=None, DofY=None):
    """
    Converts from calendar date into decimal years taking into
    account leap years

    Parameters
    ----------
    year: float
        calendar year
    month: float
        calendar month
    day: float or NoneType, default None
        day of the month
    hour: float or NoneType, default None
        hour of the day
    minute: float or NoneType, default None
        minute of the hour
    second: float or NoneType, default None
        second of the minute
    DofY: float or NoneType, default None
        day of the year (January 1 = 1)

    Returns
    -------
    t_date: float
        date in decimal-year format

    References
    ----------
    .. [1] Dershowitz, N. and E.M. Reingold. 2008.
        Calendrical Calculations.
        Cambridge: Cambridge University Press.
    """

    # number of dates
    n_dates = len(np.atleast_1d(year))

    # create arrays for calendar date variables
    cal_date = {}
    cal_date['year'] = np.zeros((n_dates))
    cal_date['month'] = np.zeros((n_dates))
    cal_date['day'] = np.zeros((n_dates))
    cal_date['hour'] = np.zeros((n_dates))
    cal_date['minute'] = np.zeros((n_dates))
    cal_date['second'] = np.zeros((n_dates))
    # day of the year
    cal_date['DofY'] = np.zeros((n_dates))

    # remove singleton dimensions and use year and month
    cal_date['year'][:] = np.squeeze(year)
    cal_date['month'][:] = np.squeeze(month)

    # create output date variable
    t_date = np.zeros((n_dates))

    # days per month in a leap and a standard year
    # only difference is February (29 vs. 28)
    dpm_leap = np.array(_dpm_leap, dtype=np.float64)
    dpm_stnd = np.array(_dpm_stnd, dtype=np.float64)

    # Rules in the Gregorian calendar for a year to be a leap year:
    # divisible by 4, but not by 100 unless divisible by 400
    # True length of the year is about 365.2422 days
    # Adding a leap day every four years ==> average 365.25
    # Subtracting a leap year every 100 years ==> average 365.24
    # Adding a leap year back every 400 years ==> average 365.2425
    # Subtracting a leap year every 4000 years ==> average 365.24225
    m4 = (cal_date['year'] % 4)
    m100 = (cal_date['year'] % 100)
    m400 = (cal_date['year'] % 400)
    m4000 = (cal_date['year'] % 4000)
    # find indices for standard years and leap years using criteria
    leap, = np.nonzero((m4 == 0) & (m100 != 0) | (m400 == 0) & (m4000 != 0))
    stnd, = np.nonzero((m4 != 0) | (m100 == 0) & (m400 != 0) | (m4000 == 0))

    # calculate the day of the year
    if DofY is not None:
        # if entered directly as an input
        # remove 1 so day 1 (Jan 1st) = 0.0 in decimal format
        cal_date['DofY'][:] = np.squeeze(DofY)-1
    else:
        # use calendar month and day of the month to calculate day of the year
        # month minus 1: January = 0, February = 1, etc (indice of month)
        # in decimal form: January = 0.0
        month_m1 = np.array(cal_date['month'],dtype=int) - 1

        # day of month
        if day is not None:
            # remove 1 so 1st day of month = 0.0 in decimal format
            cal_date['day'][:] = np.squeeze(day)-1.0
        else:
            # if not entering days as an input
            # will use the mid-month value
            cal_date['day'][leap] = dpm_leap[month_m1[leap]]/2.0
            cal_date['day'][stnd] = dpm_stnd[month_m1[stnd]]/2.0

        # create matrix with the lower half = 1
        # this matrix will be used in a matrix multiplication
        # to calculate the total number of days for prior months
        # the -1 will make the diagonal == 0
        # i.e. first row == all zeros and the
        # last row == ones for all but the last element
        mon_mat=np.tri(12,12,-1)
        # using a dot product to calculate total number of days
        # for the months before the input date
        # basically is sum(i*dpm)
        # where i is 1 for all months < the month of interest
        # and i is 0 for all months >= the month of interest
        # month of interest is zero as the exact days will be
        # used to calculate the date

        # calculate the day of the year for leap and standard
        # use total days of all months before date
        # and add number of days before date in month
        cal_date['DofY'][stnd] = cal_date['day'][stnd] + \
            np.dot(mon_mat[month_m1[stnd],:],dpm_stnd)
        cal_date['DofY'][leap] = cal_date['day'][leap] + \
            np.dot(mon_mat[month_m1[leap],:],dpm_leap)

    # hour of day (else is zero)
    if hour is not None:
        cal_date['hour'][:] = np.squeeze(hour)

    # minute of hour (else is zero)
    if minute is not None:
        cal_date['minute'][:] = np.squeeze(minute)

    # second in minute (else is zero)
    if second is not None:
        cal_date['second'][:] = np.squeeze(second)

    # calculate decimal date
    # convert hours, minutes and seconds into days
    # convert calculated fractional days into decimal fractions of the year
    # Leap years
    t_date[leap] = cal_date['year'][leap] + \
        (cal_date['DofY'][leap] + cal_date['hour'][leap]/24. + \
        cal_date['minute'][leap]/1440. + \
        cal_date['second'][leap]/86400.)/np.sum(dpm_leap)
    # Standard years
    t_date[stnd] = cal_date['year'][stnd] + \
        (cal_date['DofY'][stnd] + cal_date['hour'][stnd]/24. + \
        cal_date['minute'][stnd]/1440. + \
        cal_date['second'][stnd]/86400.)/np.sum(dpm_stnd)

    return t_date

# PURPOSE: Converts from Julian day to calendar date and time
def convert_julian(JD, **kwargs):
    """
    Converts from Julian day to calendar date and time

    Parameters
    ----------
    JD: float
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
    year: float
        calendar year
    month: float
        calendar month
    day: float
        day of the month
    hour: float
        hour of the day
    minute: float
        minute of the hour
    second: float
        second of the minute

    References
    ----------
    .. [1] "Numerical Recipes in C", by William H. Press,
        Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
        Cambridge University Press, 1988 (second printing).
    .. [2] Hatcher, D. A., "Simple Formulae for Julian Day Numbers and
        Calendar Dates", Quarterly Journal of the Royal Astronomical
        Society, 25(1), 1984.
    """
    # set default keyword arguments
    kwargs.setdefault('astype', None)
    kwargs.setdefault('format', 'dict')
    # raise warnings for deprecated keyword arguments
    deprecated_keywords = dict(ASTYPE='astype', FORMAT='format')
    for old,new in deprecated_keywords.items():
        if old in kwargs.keys():
            warnings.warn(f"""Deprecated keyword argument {old}.
                Changed to '{new}'""", DeprecationWarning)
            # set renamed argument to not break workflows
            kwargs[new] = copy.copy(kwargs[old])

    # convert to array if only a single value was imported
    if (np.ndim(JD) == 0):
        JD = np.atleast_1d(JD)
        single_value = True
    else:
        single_value = False

    # verify julian day
    JDO = np.floor(JD + 0.5)
    C = np.zeros_like(JD)
    # calculate C for dates before and after the switch to Gregorian
    IGREG = 2299161.0
    ind1, = np.nonzero(JDO < IGREG)
    C[ind1] = JDO[ind1] + 1524.0
    ind2, = np.nonzero(JDO >= IGREG)
    B = np.floor((JDO[ind2] - 1867216.25)/36524.25)
    C[ind2] = JDO[ind2] + B - np.floor(B/4.0) + 1525.0
    # calculate coefficients for date conversion
    D = np.floor((C - 122.1)/365.25)
    E = np.floor((365.0 * D) + np.floor(D/4.0))
    F = np.floor((C - E)/30.6001)
    # calculate day, month, year and hour
    day = np.floor(C - E + 0.5) - np.floor(30.6001*F)
    month = F - 1.0 - 12.0*np.floor(F/14.0)
    year = D - 4715.0 - np.floor((7.0 + month)/10.0)
    hour = np.floor(24.0*(JD + 0.5 - JDO))
    # calculate minute and second
    G = (JD + 0.5 - JDO) - hour/24.0
    minute = np.floor(G*1440.0)
    second = (G - minute/1440.0) * 86400.0

    # convert all variables to output type (from float)
    if kwargs['astype'] is not None:
        year = year.astype(kwargs['astype'])
        month = month.astype(kwargs['astype'])
        day = day.astype(kwargs['astype'])
        hour = hour.astype(kwargs['astype'])
        minute = minute.astype(kwargs['astype'])
        second = second.astype(kwargs['astype'])

    # if only a single value was imported initially: remove singleton dims
    if single_value:
        year = year.item(0)
        month = month.item(0)
        day = day.item(0)
        hour = hour.item(0)
        minute = minute.item(0)
        second = second.item(0)

    # return date variables in output format
    if (kwargs['format'] == 'dict'):
        return dict(year=year, month=month, day=day,
            hour=hour, minute=minute, second=second)
    elif (kwargs['format'] == 'tuple'):
        return (year, month, day, hour, minute, second)
    elif (kwargs['format'] == 'zip'):
        return zip(year, month, day, hour, minute, second)

# PURPOSE: Count number of leap seconds that have passed for each GPS time
def count_leap_seconds(GPS_Time, truncate=True):
    """
    Counts the number of leap seconds between a given GPS time and UTC

    Parameters
    ----------
    GPS_Time: float
        seconds since January 6, 1980 at 00:00:00
    truncate: bool, default True
        Reduce list of leap seconds to positive GPS times

    Returns
    -------
    n_leaps: float
        number of elapsed leap seconds
    """
    # get the valid leap seconds
    leaps = get_leap_seconds(truncate=truncate)
    # number of leap seconds prior to GPS_Time
    n_leaps = np.zeros_like(GPS_Time,dtype=np.float64)
    for i,leap in enumerate(leaps):
        count = np.count_nonzero(GPS_Time >= leap)
        if (count > 0):
            indices = np.nonzero(GPS_Time >= leap)
            n_leaps[indices] += 1.0
    # return the number of leap seconds for converting to UTC
    return n_leaps

# PURPOSE: Define GPS leap seconds
def get_leap_seconds(truncate=True):
    """
    Gets a list of GPS times for when leap seconds occurred

    Parameters
    ----------
    truncate: bool, default True
        Reduce list of leap seconds to positive GPS times

    Returns
    -------
    GPS time: float
        GPS seconds when leap seconds occurred
    """
    leap_secs = pyTMD.utilities.get_data_path(['data','leap-seconds.list'])
    # find line with file expiration as delta time
    with open(leap_secs, mode='r', encoding='utf8') as fid:
        secs, = [re.findall(r'\d+',i).pop() for i in fid.read().splitlines()
            if re.match(r'^(?=#@)',i)]
    # check that leap seconds file is still valid
    expiry = datetime.datetime(1900,1,1) + datetime.timedelta(seconds=int(secs))
    today = datetime.datetime.now()
    update_leap_seconds() if (expiry < today) else None
    # get leap seconds
    leap_UTC,TAI_UTC = np.loadtxt(pyTMD.utilities.get_data_path(leap_secs)).T
    # TAI time is ahead of GPS by 19 seconds
    TAI_GPS = 19.0
    # convert leap second epochs from NTP to GPS
    # convert from time of 2nd leap second to time of 1st leap second
    leap_GPS = convert_delta_time(leap_UTC+TAI_UTC-TAI_GPS-1,
        epoch1=(1900,1,1,0,0,0), epoch2=(1980,1,6,0,0,0))
    # return the GPS times of leap second occurance
    if truncate:
        return leap_GPS[leap_GPS >= 0].astype(np.float64)
    else:
        return leap_GPS.astype(np.float64)

# PURPOSE: connects to servers and downloads leap second files
def update_leap_seconds(timeout=20, verbose=False, mode=0o775):
    """
    Connects to servers to download leap-seconds.list files from NIST servers

    - https://www.nist.gov/pml/time-and-frequency-division/leap-seconds-faqs

    Servers and Mirrors

    - ftp://ftp.nist.gov/pub/time/leap-seconds.list
    - https://www.ietf.org/timezones/data/leap-seconds.list

    Parameters
    ----------
    timeout: int, default 20
        timeout in seconds for blocking operations
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file
    """
    # local version of file
    FILE = 'leap-seconds.list'
    LOCAL = pyTMD.utilities.get_data_path(['data',FILE])
    HASH = pyTMD.utilities.get_hash(LOCAL)

    # try downloading from NIST ftp servers
    HOST = ['ftp.nist.gov','pub','time',FILE]
    try:
        pyTMD.utilities.check_ftp_connection(HOST[0])
        pyTMD.utilities.from_ftp(HOST, timeout=timeout, local=LOCAL,
            hash=HASH, verbose=verbose, mode=mode)
    except Exception as e:
        logging.debug(traceback.format_exc())
        pass
    else:
        return

    # try downloading from Internet Engineering Task Force (IETF) mirror
    REMOTE = ['https://www.ietf.org','timezones','data',FILE]
    try:
        pyTMD.utilities.from_http(REMOTE, timeout=timeout, local=LOCAL,
            hash=HASH, verbose=verbose, mode=mode)
    except Exception as e:
        logging.debug(traceback.format_exc())
        pass
    else:
        return

# PURPOSE: Download delta time files and merge into a single
def merge_delta_time(username=None, password=None, verbose=False, mode=0o775):
    """
    Connects to servers to download historic_deltat.data and deltat.data files

    Reads IERS Bulletin-A produced iers_deltat.data files

    Creates a merged file combining the historic, monthly and daily files

    Long-term Delta T

    - https://www.usno.navy.mil/USNO/earth-orientation/eo-products/long-term

    Parameters
    ----------
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file

    Notes
    -----
    Delta times are the difference between universal time and dynamical time
    """
    # retrieve history delta time files
    pull_deltat_file('historic_deltat.data',username=username,password=password,
        verbose=verbose,mode=mode)
    # read historic delta time file
    historic_file=pyTMD.utilities.get_data_path(['data','historic_deltat.data'])
    historic = np.loadtxt(historic_file, skiprows=2)
    HY = np.floor(historic[:,0])
    HM = 12.0*np.mod(historic[:,0],1.0) + 1.0
    HD = np.ones_like(historic[:,0])
    # retrieve monthly delta time files
    pull_deltat_file('deltat.data',username=username,password=password,
        verbose=verbose,mode=mode)
    # read modern monthly delta time file
    monthly_file = pyTMD.utilities.get_data_path(['data','deltat.data'])
    monthly = np.loadtxt(monthly_file)
    monthly_time = convert_calendar_decimal(monthly[:,0],monthly[:,1],
        day=monthly[:,2])
    # retrieve daily delta time files
    merge_bulletin_a_files(username=username, password=password,
        verbose=verbose, mode=mode)
    # read modern daily delta time file from IERS Bulletin A files
    daily_file = pyTMD.utilities.get_data_path(['data','iers_deltat.data'])
    daily = np.loadtxt(daily_file)
    daily_time = convert_calendar_decimal(daily[:,0], daily[:,1],
        day=daily[:,2])
    # write to new merged file
    merged_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
    fid = open(merged_file, mode='w', encoding='utf8')
    logging.info(merged_file)
    file_format = ' {0:4.0f} {1:2.0f} {2:2.0f} {3:7.4f}'
    # use historical values for times prior to monthly
    ind1, = np.nonzero(historic[:,0] < monthly_time[0])
    for i in ind1:
        args = (HY[i],HM[i],HD[i],historic[i,1])
        print(file_format.format(*args),file=fid)
    # use monthly values for times prior to daily
    ind2, = np.nonzero(monthly_time < np.min(daily_time))
    for i in ind2:
        args = (monthly[i,0],monthly[i,1],monthly[i,2],monthly[i,3])
        print(file_format.format(*args),file=fid)
    # use daily values for all times available
    for i in np.argsort(daily_time):
        args = (daily[i,0],daily[i,1],daily[i,2],daily[i,3])
        print(file_format.format(*args),file=fid)
    # close the merged file and change the permissions mode
    fid.close()
    os.chmod(merged_file,mode)

# PURPOSE: Append Bulletin-A file to merged delta time file
def append_delta_time(verbose=False, mode=0o775):
    """
    Appends merged delta time file with values from latest Bulletin-A file

    Parameters
    ----------
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file

    Notes
    -----
    Delta times are the difference between universal time and dynamical time
    """

    # append to merged file
    merged_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])
    fid = open(merged_file, mode='a', encoding='utf8')
    logging.info(merged_file)
    file_format = ' {0:4.0f} {1:2.0f} {2:2.0f} {3:7.4f}'
    # read latest Bulletin-A file from IERS
    bulletin_file = pyTMD.utilities.get_data_path(['data','ser7.dat'])
    logging.info(bulletin_file)
    with open(bulletin_file, mode='rb') as fileID:
        YY,MM,DD,DELTAT = read_iers_bulletin_a(fileID)
    # append latest delta time values to merged file
    for Y,M,D,T in zip(YY,MM,DD,DELTAT):
        print(file_format.format(Y,M,D,T), file=fid)
    # close the merged file and change the permissions mode
    fid.close()
    os.chmod(merged_file,mode)

# PURPOSE: connect to IERS or CDDIS server and merge Bulletin-A files
def merge_bulletin_a_files(username=None,password=None,
    verbose=False,mode=0o775):
    """
    Attempt to connects to the IERS server and the CDDIS Earthdata server
    to download and merge Bulletin-A files

    Reads the IERS Bulletin-A files and calculates the daily delta times

    Servers and Mirrors

    - https://datacenter.iers.org/availableVersions.php?id=6
    - ftp://ftp.iers.org/products/eop/rapid/bulletina/
    - https://cddis.nasa.gov/archive/products/iers/iers_bulletins/bulletin_a/

    Parameters
    ----------
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file

    Notes
    -----
    Delta times are the difference between universal time and dynamical time
    """
    # if complete: replace previous version of file
    LOCAL = pyTMD.utilities.get_data_path(['data','iers_deltat.data'])
    COPY = pyTMD.utilities.get_data_path(['data','iers_deltat.temp'])
    # try connecting to IERS http servers and merge Bulletin-A files
    try:
        iers_delta_time(COPY, verbose=verbose, mode=mode)
    except Exception as e:
        logging.debug(traceback.format_exc())
        os.remove(COPY) if os.access(COPY, os.F_OK) else None
        pass
    else:
        pyTMD.utilities.copy(COPY, LOCAL, move=True)
        return

    # try connecting to IERS ftp servers and merge Bulletin-A files
    try:
        iers_ftp_delta_time(COPY, verbose=verbose, mode=mode)
    except Exception as e:
        logging.debug(traceback.format_exc())
        os.remove(COPY) if os.access(COPY, os.F_OK) else None
        pass
    else:
        pyTMD.utilities.copy(COPY, LOCAL, move=True)
        return

    # try connecting to CDDIS https servers and merge Bulletin-A files
    try:
        cddis_delta_time(COPY, username=username, password=password,
            verbose=verbose, mode=mode)
    except Exception as e:
        logging.debug(traceback.format_exc())
        os.remove(COPY) if os.access(COPY, os.F_OK) else None
        pass
    else:
        pyTMD.utilities.copy(COPY, LOCAL, move=True)
        return

# PURPOSE: connects to IERS ftp servers and finds Bulletin-A files
def iers_ftp_delta_time(daily_file, timeout=120, verbose=False, mode=0o775):
    """
    Connects to the IERS ftp server to download Bulletin-A files

    - https://datacenter.iers.org/productMetadata.php?id=6

    Reads the IERS Bulletin-A files and calculates the daily delta times

    Servers and Mirrors

    - ftp://ftp.iers.org/products/eop/rapid/bulletina/

    Parameters
    ----------
    daily_file: str
        output daily delta time file from merged Bulletin-A files
    timeout: int, default 120
        timeout in seconds for blocking operations
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file

    Notes
    -----
    Delta times are the difference between universal time and dynamical time
    """
    # connect to ftp host for IERS bulletins
    HOST = ['ftp.iers.org','products','eop','rapid','bulletina']
    pyTMD.utilities.check_ftp_connection(HOST[0])
    # regular expression pattern for finding files
    rx = re.compile(r'bulletina-(.*?)-(\d+).txt$',re.VERBOSE)
    # open output daily delta time file
    fid = open(daily_file, mode='w', encoding='utf8')
    # output file format
    file_format = ' {0:4.0f} {1:2.0f} {2:2.0f} {3:7.4f}'
    # find subdirectories
    subdirectory,_ = pyTMD.utilities.ftp_list(HOST,timeout=timeout,
        basename=True,sort=True)
    # for each subdirectory
    for SUB in subdirectory:
        # find Bulletin-A files in ftp subdirectory
        HOST.append(SUB)
        logging.info(SUB)
        bulletin_files,_ = pyTMD.utilities.ftp_list(HOST,
            timeout=timeout,basename=True,sort=True,pattern=rx)
        # for each Bulletin-A file
        for f in sorted(bulletin_files):
            logging.info(f)
            # copy remote file contents to BytesIO object
            HOST.append(f)
            remote_buffer = pyTMD.utilities.from_ftp(HOST,timeout=timeout)
            # read Bulletin-A file from BytesIO object
            YY,MM,DD,DELTAT = read_iers_bulletin_a(remote_buffer)
            # print delta time for week to output file
            for Y,M,D,T in zip(YY,MM,DD,DELTAT):
                print(file_format.format(Y,M,D,T),file=fid)
            # close the bytesIO object
            remote_buffer.close()
            # remove the file from the list
            HOST.remove(f)
        # remove the subdirectory from the list
        HOST.remove(SUB)
    # close the output file
    fid.close()
    # change the permissions mode
    os.chmod(daily_file,mode)

# PURPOSE: connects to IERS http servers and finds Bulletin-A files
def iers_delta_time(daily_file, timeout=120, verbose=False, mode=0o775):
    """
    Connects to the IERS server to download Bulletin-A files

    - https://datacenter.iers.org/productMetadata.php?id=6

    Reads the IERS Bulletin-A files and calculates the daily delta times

    Servers and Mirrors

    - https://datacenter.iers.org/availableVersions.php?id=6

    Parameters
    ----------
    daily_file: str
        output daily delta time file from merged Bulletin-A files
    timeout: int, default 120
        timeout in seconds for blocking operations
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file

    Notes
    -----
    Delta times are the difference between universal time and dynamical time
    """
    # open output daily delta time file
    fid = open(daily_file, mode='w', encoding='utf8')
    # output file format
    file_format = ' {0:4.0f} {1:2.0f} {2:2.0f} {3:7.4f}'
    # connect to http host for IERS Bulletin-A files
    HOST = 'https://datacenter.iers.org/availableVersions.php?id=6'
    bulletin_files,_ = pyTMD.utilities.iers_list(HOST, timeout=timeout)
    # for each Bulletin-A file
    for f in bulletin_files[-100:]:
        logging.info(f)
        remote_buffer = pyTMD.utilities.from_http(f, timeout=timeout)
        # read Bulletin-A file from BytesIO object
        YY,MM,DD,DELTAT = read_iers_bulletin_a(remote_buffer)
        # print delta time for week to output file
        for Y,M,D,T in zip(YY,MM,DD,DELTAT):
            print(file_format.format(Y,M,D,T), file=fid)
        # close the bytesIO object
        remote_buffer.close()
    # close the output file
    fid.close()
    # change the permissions mode
    os.chmod(daily_file, mode)

# PURPOSE: connects to CDDIS Earthdata https server and finds Bulletin-A files
def cddis_delta_time(daily_file, username=None, password=None,
    verbose=False, mode=0o775):
    """
    Connects to the CDDIS Earthdata server to download Bulletin-A files

    Reads the IERS Bulletin-A files and calculates the daily delta times

    Servers and Mirrors

    - https://cddis.nasa.gov/archive/products/iers/iers_bulletins/bulletin_a/

    Parameters
    ----------
    daily_file: str
        output daily delta time file from merged Bulletin-A files
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file

    Notes
    -----
    Delta times are the difference between universal time and dynamical time
    """
    # connect to CDDIS Earthdata host for IERS bulletins
    HOST = ['https://cddis.nasa.gov','archive','products','iers',
        'iers_bulletins','bulletin_a']
    # build NASA Earthdata opener for CDDIS and check credentials
    pyTMD.utilities.build_opener(username, password)
    pyTMD.utilities.check_credentials()
    # regular expression pattern for finding directories
    R1 = re.compile(r'volume_(.*?)$',re.VERBOSE)
    # regular expression pattern for finding files
    R2 = re.compile(r'iers_bulletina\.(.*?)_(\d+)$',re.VERBOSE)
    # open output daily delta time file
    fid = open(daily_file, mode='w', encoding='utf8')
    file_format = ' {0:4.0f} {1:2.0f} {2:2.0f} {3:7.4f}'
    # for each subdirectory
    subdirectory,mtimes=pyTMD.utilities.cddis_list(HOST,build=False,pattern=R1)
    # extract roman numerals from subdirectories
    roman = [R1.findall(s).pop() for s in subdirectory]
    # sort the list of Roman numerals
    subdirectory = [subdirectory[i] for i,j in sorted(enumerate(roman),
        key=lambda i: pyTMD.utilities.roman_to_int(i[1]))]
    # output file format
    for SUB in subdirectory:
        # find Bulletin-A files in https subdirectory
        HOST.append(SUB)
        bulletin_files,mtimes = pyTMD.utilities.cddis_list(HOST,build=False,
            sort=True,pattern=R2)
        # for each Bulletin-A file
        for f in sorted(bulletin_files):
            logging.info(f)
            # copy remote file contents to BytesIO object
            HOST.append(f)
            remote_buffer = pyTMD.utilities.from_cddis(HOST,
                build=False,timeout=20)
            # read Bulletin-A file from BytesIO object
            YY,MM,DD,DELTAT = read_iers_bulletin_a(remote_buffer)
            # print delta time for week to output file
            for Y,M,D,T in zip(YY,MM,DD,DELTAT):
                print(file_format.format(Y,M,D,T),file=fid)
            # close the bytesIO object
            remote_buffer.close()
            # remove the file from the list
            HOST.remove(f)
        # remove the subdirectory from the list
        HOST.remove(SUB)
    # close the output file
    fid.close()
    # change the permissions mode
    os.chmod(daily_file,mode)

# PURPOSE: reads IERS Bulletin-A and calculates the delta times
def read_iers_bulletin_a(fileID):
    """
    Read a weekly IERS Bulletin-A file and calculate the
    delta times (TT - UT1)

    Parameters
    ----------
    fileID: obj
        open file object for Bulletin-A file

    Returns
    -------
    Y: float,
        calendar year
    M: float
        calendar month
    D: float
        day of the month
    DELTAT: float
        difference between universal time and dynamical time

    Notes
    -----
    Delta times are the difference between universal time and dynamical time
    """
    # read contents from input file object
    file_contents = fileID.read().decode('utf8').splitlines()

    # parse header text to find time offsets
    # TT-TAI
    TT_TAI = 0
    # TAI-UTC
    TAI_UTC = 0
    # counts the number of lines in the header
    count = 0
    HEADER = False
    # Reading over header text
    while not HEADER:
        # file line at count
        l = file_contents[count]
        # check if line contains time offsets
        if re.search(r'TT\s\=\sTAI',l):
            TT_TAI = np.float64(re.findall(r'(\d+\.\d+)',l).pop())
        if re.search(r'TAI-UTC',l):
            TAI_UTC = np.float64(re.findall(r'=\s(\d+\.\d+)',l).pop())
        # find line to set HEADER flag to True
        HEADER = bool(re.search(r'COMBINED\sEARTH\sORIENTATION\sPARAMETERS:',l))
        # add 1 to counter
        count += 1

    # convert variables to numpy arrays
    MJD = np.zeros((7))
    UT1_UTC = np.zeros((7))
    valid = 0
    # for each day in the week
    for i in range(7):
        try:
            # split numerical instances from data line
            line_contents = file_contents[count+i+4].split()
            # years are not always complete in the bulletin file
            # Modified Julian Day (days since 1858-11-17T00:00:00)
            MJD[i] = np.float64(line_contents[3])
            # difference between UT1 and UTC times
            UT1_UTC[i] = np.float64(line_contents[8])
        except (IndexError,ValueError):
            pass
        else:
            valid += 1

    # calculate components for delta time
    # TAI time is ahead of GPS by 19 seconds
    TAI_GPS = 19.0
    # calculate calendar dates from Modified Julian days
    Y,M,D,h,m,s = convert_julian(MJD[:valid]+2400000.5, format='tuple')
    # calculate GPS Time (seconds since 1980-01-06T00:00:00)
    # by converting the Modified Julian days (days since 1858-11-17T00:00:00)
    GPS_Time = convert_delta_time(MJD[:valid]*8.64e4, epoch1=(1858,11,17,0,0,0),
        epoch2=(1980,1,6,0,0,0), scale=1.0) + TAI_UTC - TAI_GPS
    # number of leap seconds between GPS and UTC
    # this finds the daily correction for weeks with leap seconds
    GPS_UTC = count_leap_seconds(GPS_Time)
    # calculate delta time (TT - UT1) -->
    # (TT-TAI) + (TAI-GPS) + (GPS-UTC) - (UT1-UTC)
    DELTAT = TT_TAI + TAI_GPS + GPS_UTC - UT1_UTC[:valid]

    # return dates and delta times
    return (Y,M,D,DELTAT)

# PURPOSE: connects to servers and downloads latest Bulletin-A file
def update_bulletin_a(timeout=20, verbose=False, mode=0o775):
    """
    Connects to IERS Rapid Service/Prediction Center (RS/PC) and
    downloads latest Bulletin-A file

    - https://maia.usno.navy.mil/ser7/readme.bulla

    Servers and Mirrors

    - https://maia.usno.navy.mil/ser7/ser7.dat

    Parameters
    ----------
    timeout: int, default 20
        timeout in seconds for blocking operations
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file
    """
    # local version of file
    LOCAL = pyTMD.utilities.get_data_path(['data','ser7.dat'])
    HASH = pyTMD.utilities.get_hash(LOCAL)

    # try downloading from IERS Rapid Service/Prediction Center (RS/PC)
    REMOTE = ['https://maia.usno.navy.mil','ser7','ser7.dat']
    try:
        pyTMD.utilities.from_http(REMOTE, timeout=timeout, local=LOCAL,
            hash=HASH, verbose=verbose, mode=mode)
    except Exception as e:
        logging.debug(traceback.format_exc())
        pass
    else:
        return

# PURPOSE: connects to servers and downloads delta time files
def pull_deltat_file(FILE,username=None,password=None,verbose=False,mode=0o775):
    """
    Connects to servers and downloads delta time files

    Servers and Mirrors

    - http://maia.usno.navy.mil/ser7/
    - https://cddis.nasa.gov/archive/products/iers/
    - ftp://cddis.nasa.gov/products/iers/
    - ftp://cddis.gsfc.nasa.gov/pub/products/iers/

    Parameters
    ----------
    FILE: str
        delta time file to download from remote servers

            - deltat.data: monthly deltat file
            - historic_deltat.data: historic deltat file
    username: str or NoneType, default None
        NASA Earthdata username
    password: str or NoneType, default None
        NASA Earthdata password
    verbose: bool, default False
        print file information about output file
    mode: oct, default 0o775
        permissions mode of output file

    Notes
    -----
    Delta times are the difference between universal time and dynamical time
    """
    # local version of file
    LOCAL = pyTMD.utilities.get_data_path(['data',FILE])
    HASH = pyTMD.utilities.get_hash(LOCAL)

    # try downloading from US Naval Oceanography Portal
    HOST = ['http://maia.usno.navy.mil','ser7',FILE]
    try:
        pyTMD.utilities.from_http(HOST,timeout=5,local=LOCAL,hash=HASH,
            verbose=verbose,mode=mode)
    except Exception as e:
        logging.debug(traceback.format_exc())
        pass
    else:
        return

    # try downloading from NASA Crustal Dynamics Data Information System
    # NOTE: anonymous ftp access was discontinued on 2020-10-31
    # requires using the following https Earthdata server
    server = []
    # server.append(['cddis.nasa.gov','pub','products','iers',FILE])
    # server.append(['cddis.gsfc.nasa.gov','products','iers',FILE])
    for HOST in server:
        try:
            pyTMD.utilities.check_ftp_connection(HOST[0])
            pyTMD.utilities.from_ftp(HOST,timeout=20,local=LOCAL,hash=HASH,
                verbose=verbose,mode=mode)
        except Exception as e:
            logging.debug(traceback.format_exc())
            pass
        else:
            return

    # try downloading from NASA Crustal Dynamics Data Information System
    # using NASA Earthdata credentials stored in netrc file
    HOST = ['https://cddis.nasa.gov','archive','products','iers',FILE]
    try:
        pyTMD.utilities.from_cddis(HOST,username=username,password=password,
            timeout=20,local=LOCAL,hash=HASH,verbose=verbose,mode=mode)
    except Exception as e:
        logging.debug(traceback.format_exc())
        pass
    else:
        return
