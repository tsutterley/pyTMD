#!/usr/bin/env python
u"""
test_time.py (09/2020)
Verify time conversion functions
"""
import pytest
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
        epoch=pyTMD.time._mjd_epoch)
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
    atlas_sdp_epoch = '2018-01-01 00:00:00'
    UNIX = pyTMD.utilities.get_unix_time(atlas_sdp_epoch)
    assert (UNIX == 1514764800)
    # check UNIX time conversion with delta times
    output_time = pyTMD.time.convert_delta_time(UNIX,
        epoch1=pyTMD.time._unix_epoch, epoch2=atlas_sdp_epoch,
        scale=1.0)
    assert (output_time == 0)

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
    # time string for unitless case with a time zone
    time_string = '2000-01-01T12:00:00.000-06:00'
    epoch,to_secs = pyTMD.time.parse_date_string(time_string)
    # check the epoch and the time unit conversion factors
    assert np.all(epoch == [2000,1,1,18,0,0])
    assert (to_secs == 0.0)

# PURPOSE: test isoformat
def test_isoformat():
    time_string = '2000-01-01'
    output = pyTMD.utilities.isoformat(time_string)
    validation = '2000-01-01T00:00:00'
    assert output == validation

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
        epoch1=pyTMD.time._gps_epoch, epoch2=pyTMD.time._atlas_sdp_epoch,
        scale=1.0)
    assert (delta_time == output_time)
    # compare with original values using string epochs
    output_time = pyTMD.time.convert_delta_time(gps_seconds - time_leaps,
        epoch1='1980-01-06T00:00:00', epoch2='2018-01-01T00:00:00',
        scale=1.0)
    assert (delta_time == output_time)
    # calculate and compare with timescale values
    GPS = pyTMD.time.timescale().from_deltatime(gps_seconds,
        epoch=pyTMD.time._gps_epoch, standard='GPS')
    output_time = GPS.to_deltatime(pyTMD.time._atlas_sdp_epoch,
        scale=GPS.day)
    assert np.isclose(delta_time, output_time)

# PURPOSE: test timescale class conversions and constants
def test_timescale():
    # ATLAS Standard Data Epoch
    atlas_sdp_epoch = np.datetime64('2018-01-01T00:00:00')
    # from datetime
    ATLAS = pyTMD.time.timescale().from_datetime(atlas_sdp_epoch)
    assert np.all(ATLAS.MJD == 58119)
    assert np.all(ATLAS.tide == 9497)
    delta_time_epochs = (ATLAS.to_datetime() - atlas_sdp_epoch)
    assert np.all(delta_time_epochs/np.timedelta64(1, 'ns') == 0)
    # from deltatime
    ATLAS = pyTMD.time.timescale().from_deltatime(0, epoch=(2018,1,1))
    assert np.all(ATLAS.MJD == 58119)
    assert np.all(ATLAS.tide == 9497)
    delta_time_epochs = (ATLAS.to_datetime() - atlas_sdp_epoch)
    assert np.all(delta_time_epochs/np.timedelta64(1, 'ns') == 0)
    # from MJD
    ATLAS = pyTMD.time.timescale(MJD=58119)
    assert np.all(ATLAS.ut1 == 2458119.5)
    assert np.all(ATLAS.tide == 9497)
    assert np.all((ATLAS.MJD - 51544.5) == (ATLAS.ut1 - 2451545.0))
    # from MJD hourly array
    delta_time = np.arange(0, 365, 1.0/24.0)
    ATLAS = pyTMD.time.timescale(MJD=58119 + delta_time)
    delta_time_epochs = (ATLAS.to_datetime() - atlas_sdp_epoch)
    assert np.allclose(delta_time_epochs/np.timedelta64(1, 'D'), delta_time)
    assert np.allclose(ATLAS.to_deltatime(epoch=atlas_sdp_epoch), delta_time)
    assert np.allclose(ATLAS.to_deltatime(epoch='2018-01-01'), delta_time)
    # check constants
    assert (ATLAS.century == 36525.0)
    assert (ATLAS.day == 86400.0)
    assert (ATLAS.turn == 1.0)
    assert (ATLAS.turndeg == 360.0)
    assert (ATLAS.turnasec == 1296000.0)
    assert (ATLAS.deg2asec == 3600.0)
    assert (ATLAS.deg2rad == np.pi/180.0)
    assert (ATLAS.asec2rad == np.pi/648000.0)
    assert (ATLAS.masec2rad == np.pi/0.648e12)

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
        assert pyTMD.utilities.get_data_path(['data',FILE]).exists()
