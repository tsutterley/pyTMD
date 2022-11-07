#!/usr/bin/env python
u"""
test_time.py (09/2020)
Verify time conversion functions
"""
import os
import pytest
import warnings
import numpy as np
import pyTMD.time
import pyTMD.utilities

# parameterize calendar dates
@pytest.mark.parametrize("YEAR", np.random.randint(1992,2020,size=2))
@pytest.mark.parametrize("MONTH", np.random.randint(1,13,size=2))
# PURPOSE: verify forward and backwards time conversions
def test_julian(YEAR,MONTH):
    # days per month in a leap and a standard year
    # only difference is February (29 vs. 28)
    dpm_leap = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    dpm_stnd = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    DPM = dpm_stnd if np.mod(YEAR,4) else dpm_leap
    # calculate Modified Julian Day (MJD) from calendar date
    DAY = np.random.randint(1,DPM[MONTH-1]+1)
    HOUR = np.random.randint(0,23+1)
    MINUTE = np.random.randint(0,59+1)
    SECOND = 60.0*np.random.random_sample(1)
    MJD = pyTMD.time.convert_calendar_dates(YEAR, MONTH, DAY,
        hour=HOUR, minute=MINUTE, second=SECOND,
        epoch=(1858,11,17,0,0,0))
    # convert MJD to calendar date
    JD = np.squeeze(MJD) + 2400000.5
    YY,MM,DD,HH,MN,SS = pyTMD.time.convert_julian(JD,
        format='tuple', astype=np.float64)
    # assert dates
    assert (YY == YEAR)
    assert (MM == MONTH)
    assert (DD == DAY)
    assert (HH == HOUR)
    assert (MN == MINUTE)
    assert np.isclose(SS, SECOND, atol=1e-2)

# parameterize calendar dates
@pytest.mark.parametrize("YEAR", np.random.randint(1992,2020,size=2))
@pytest.mark.parametrize("MONTH", np.random.randint(1,13,size=2))
# PURPOSE: verify forward and backwards time conversions
def test_decimal_dates(YEAR,MONTH):
    # days per month in a leap and a standard year
    # only difference is February (29 vs. 28)
    dpm_leap = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    dpm_stnd = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    DPM = dpm_stnd if np.mod(YEAR,4) else dpm_leap
    assert (np.sum(DPM) == pyTMD.time.calendar_days(YEAR).sum())
    # calculate Modified Julian Day (MJD) from calendar date
    DAY = np.random.randint(1,DPM[MONTH-1]+1)
    HOUR = np.random.randint(0,23+1)
    MINUTE = np.random.randint(0,59+1)
    SECOND = 60.0*np.random.random_sample(1)
    # calculate year-decimal time
    tdec = pyTMD.time.convert_calendar_decimal(YEAR, MONTH, day=DAY,
        hour=HOUR, minute=MINUTE, second=SECOND)
    # day of the year 1 = Jan 1, 365 = Dec 31 (std)
    day_temp = np.mod(tdec, 1)*np.sum(DPM)
    DofY = np.floor(day_temp) + 1
    # cumulative sum of the calendar dates
    day_cumulative = np.cumsum(np.concatenate(([0],DPM))) + 1
    # finding which month date is in
    i = np.nonzero((DofY >= day_cumulative[0:-1]) & (DofY < day_cumulative[1:]))
    month_range = np.arange(1,13)
    month = month_range[i]
    # finding day of the month
    day = (DofY - day_cumulative[i]) + 1
    # convert residuals into time (hour, minute and second)
    hour_temp = np.mod(day_temp,1)*24.0
    minute_temp = np.mod(hour_temp,1)*60.0
    second = np.mod(minute_temp,1)*60.0
    # assert dates
    assert (np.floor(tdec) == YEAR)
    assert (month == MONTH)
    assert (day == DAY)
    assert (np.floor(hour_temp) == HOUR)
    assert (np.floor(minute_temp) == MINUTE)
    assert np.isclose(second, SECOND, atol=1e-5)

# PURPOSE: test UNIX time
def test_unix_time():
    # ATLAS Standard Data Epoch
    UNIX = pyTMD.utilities.get_unix_time('2018-01-01 00:00:00')
    assert (UNIX == 1514764800)

# PURPOSE: test parsing time strings
def test_parse_date_string():
    # time string for Modified Julian Days
    time_string = 'days since 1858-11-17T00:00:00'
    epoch,to_secs = pyTMD.time.parse_date_string(time_string)
    # check the epoch and the time unit conversion factors
    assert np.all(epoch == [1858,11,17,0,0,0])
    assert (to_secs == 86400.0)
    # time string for ATLAS Standard Data Epoch
    time_string = 'seconds since 2018-01-01T00:00:00'
    epoch,to_secs = pyTMD.time.parse_date_string(time_string)
    # check the epoch and the time unit conversion factors
    assert np.all(epoch == [2018,1,1,0,0,0])
    assert (to_secs == 1.0)
    # time string for unitless case
    time_string = '2000-01-01T12:00:00'
    epoch,to_secs = pyTMD.time.parse_date_string(time_string)
    # check the epoch and the time unit conversion factors
    assert np.all(epoch == [2000,1,1,12,0,0])
    assert (to_secs == 0.0)

# PURPOSE: verify forward and backwards delta time conversions
@pytest.mark.parametrize("delta_time", np.random.randint(1,31536000,size=4))
def test_delta_time(delta_time, gps_epoch=1198800018.0):
    # convert to array if single value
    delta_time = np.atleast_1d(delta_time)
    # calculate gps time from delta_time
    gps_seconds = gps_epoch + delta_time
    time_leaps = pyTMD.time.count_leap_seconds(gps_seconds)
    # compare output delta times with original values
    output_time = pyTMD.time.convert_delta_time(gps_seconds - time_leaps,
        epoch1=(1980,1,6,0,0,0), epoch2=(2018,1,1,0,0,0), scale=1.0)
    assert (delta_time == output_time)

# PURPOSE: update delta time files and values
def test_update_delta_time(username, password):
    pyTMD.time.merge_delta_time(username=username,password=password)
    # confirm delta time files
    delta_time_files = []
    delta_time_files.append('historic_deltat.data')
    delta_time_files.append('deltat.data')
    delta_time_files.append('iers_deltat.data')
    delta_time_files.append('merged_deltat.data')
    for FILE in delta_time_files:
        assert os.access(pyTMD.utilities.get_data_path(['data',FILE]),os.F_OK)
