#!/usr/bin/env python
u"""
test_time.py (08/2020)
Verify time conversion functions
"""
import warnings
import pytest
import numpy as np
import pyTMD.time
from pyTMD.convert_julian import convert_julian
from pyTMD.convert_calendar_decimal import convert_calendar_decimal

#-- parameterize calendar dates
@pytest.mark.parametrize("YEAR", np.random.randint(1992,2020,size=2))
@pytest.mark.parametrize("MONTH", np.random.randint(1,13,size=2))
#-- PURPOSE: verify forward and backwards time conversions
def test_julian(YEAR,MONTH):
    #-- days per month in a leap and a standard year
    #-- only difference is February (29 vs. 28)
    dpm_leap = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    dpm_stnd = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    DPM = dpm_stnd if np.mod(YEAR,4) else dpm_leap
    #-- calculate Modified Julian Day (MJD) from calendar date
    DAY = np.random.randint(1,DPM[MONTH-1]+1)
    HOUR = np.random.randint(0,23+1)
    MINUTE = np.random.randint(0,59+1)
    SECOND = 60.0*np.random.random_sample(1)
    MJD = pyTMD.time.convert_calendar_dates(YEAR, MONTH, DAY,
        hour=HOUR, minute=MINUTE, second=SECOND,
        epoch=(1858,11,17,0,0,0))
    #-- convert MJD to calendar date
    YY,MM,DD,HH,MN,SS = convert_julian(MJD+2400000.5, FORMAT='tuple')
    #-- assert dates
    eps = np.finfo(np.float16).eps
    assert (YY == YEAR)
    assert (MM == MONTH)
    assert (DD == DAY)
    assert (HH == HOUR)
    assert (MN == MINUTE)
    assert (np.abs(SS - SECOND) < eps)

#-- parameterize calendar dates
@pytest.mark.parametrize("YEAR", np.random.randint(1992,2020,size=2))
@pytest.mark.parametrize("MONTH", np.random.randint(1,13,size=2))
#-- PURPOSE: verify forward and backwards time conversions
def test_decimal_dates(YEAR,MONTH):
    #-- days per month in a leap and a standard year
    #-- only difference is February (29 vs. 28)
    dpm_leap = np.array([31,29,31,30,31,30,31,31,30,31,30,31])
    dpm_stnd = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    DPM = dpm_stnd if np.mod(YEAR,4) else dpm_leap
    #-- calculate Modified Julian Day (MJD) from calendar date
    DAY = np.random.randint(1,DPM[MONTH-1]+1)
    HOUR = np.random.randint(0,23+1)
    MINUTE = np.random.randint(0,59+1)
    SECOND = 60.0*np.random.random_sample(1)
    #-- calculate year-decimal time
    tdec = convert_calendar_decimal(YEAR, MONTH, DAY=DAY,
        HOUR=HOUR, MINUTE=MINUTE, SECOND=SECOND)
    #-- day of the year 1 = Jan 1, 365 = Dec 31 (std)
    day_temp = np.mod(tdec, 1)*np.sum(DPM)
    DofY = np.floor(day_temp) + 1
    #-- cumulative sum of the calendar dates
    day_cumulative = np.cumsum(np.concatenate(([0],DPM))) + 1
    #-- finding which month date is in
    i = np.nonzero((DofY >= day_cumulative[0:-1]) & (DofY < day_cumulative[1:]))
    month_range = np.arange(1,13)
    month = month_range[i]
    #-- finding day of the month
    day = (DofY - day_cumulative[i]) + 1
    #-- convert residuals into time (hour, minute and second)
    hour_temp = np.mod(day_temp,1)*24.0
    minute_temp = np.mod(hour_temp,1)*60.0
    second = np.mod(minute_temp,1)*60.0
    #-- assert dates
    eps = np.finfo(np.float16).eps
    assert (np.floor(tdec) == YEAR)
    assert (month == MONTH)
    assert (day == DAY)
    assert (np.floor(hour_temp) == HOUR)
    assert (np.floor(minute_temp) == MINUTE)
    assert (np.abs(second - SECOND) < eps)

#-- PURPOSE: verify forward and backwards delta time conversions
@pytest.mark.parametrize("delta_time", np.random.randint(1,31536000,size=4))
def test_delta_time(delta_time, gps_epoch=1198800018.0):
    #-- convert to array if single value
    if (np.ndim(delta_time) == 0):
        delta_time = np.array([delta_time])
    #-- calculate gps time from delta_time
    gps_seconds = gps_epoch + delta_time
    time_leaps = pyTMD.time.count_leap_seconds(gps_seconds)
    #-- compare output delta times with original values
    output_time = pyTMD.time.convert_delta_time(gps_seconds - time_leaps,
        epoch1=(1980,1,6,0,0,0), epoch2=(2018,1,1,0,0,0), scale=1.0)
    assert (delta_time == output_time)
