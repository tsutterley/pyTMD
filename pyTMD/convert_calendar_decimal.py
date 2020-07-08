#!/usr/bin/env python
u"""
convert_calendar_decimal.py
Written by Tyler Sutterley (07/2020)

Converts from calendar date into decimal years
Converts year, month (day, hour, minute, second)
    into decimal years taking into account leap years

CALLING SEQUENCE:
    t_date = convert_calendar_decimal(year, month)
    t_date = convert_calendar_decimal(year, month, DAY=day, \
        HOUR=hour, MINUTE=minute, SECOND=second)

INPUTS:
    year: can be a single value or an array of dates
    month: can be a single value or an array of dates

OPTION:
    DAY: can be a single value or an array of dates
    HOUR: can be a single value or an array of dates
    MINUTE: can be a single value or an array of dates
    SECOND: can be a single value or an array of dates
    DofY: day of the year (January 1 = 1)
        can be a single value or an array of dates

OUTPUTS:
    t_date: date in decimal format (years)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

NOTES:
    Dershowitz, N. and E.M. Reingold. 2008.  Calendrical Calculations.
        Cambridge: Cambridge University Press.

UPDATE HISTORY:
    Updated 07/2020: added function docstrings
    Updated 05/2015: updated comments and minor update to nonzero statement
    Updated 05/2014: added option for day of year
    Updated 04/2014: new code from convert_J2000.py
    Updated 04/2014: updated comments and improved rules
        for leap years to include mod 100 and mod 400
    Written 04/2014
"""
import numpy as np

def convert_calendar_decimal(year, month, DAY=None, HOUR=None, MINUTE=None,
    SECOND=None, DofY=None):
    """
    Converts from calendar date into decimal years taking into
    account leap years

    Arguments
    ---------
    year: calendar year
    month: calendar month

    Keyword arguments
    -----------------
    DAY: day of the month
    HOUR: hour of the day
    MINUTE: minute of the hour
    SECOND: second of the minute
    DofY: day of the year (January 1 = 1)

    Returns
    -------
    t_date: date in decimal format
    """

    #-- number of dates
    if (np.ndim(np.squeeze(year)) == 0):
        #-- single date entered
        n_dates = 1
    else:
        #-- array of dates entered
        n_dates = len(np.squeeze(year))

    #-- create arrays for calendar date variables
    cal_date = {}
    cal_date['year'] = np.zeros((n_dates))
    cal_date['month'] = np.zeros((n_dates))
    cal_date['day'] = np.zeros((n_dates))
    cal_date['hour'] = np.zeros((n_dates))
    cal_date['minute'] = np.zeros((n_dates))
    cal_date['second'] = np.zeros((n_dates))
    #-- day of the year
    cal_date['DofY'] = np.zeros((n_dates))

    #-- remove singleton dimensions and use year and month
    cal_date['year'][:] = np.squeeze(year)
    cal_date['month'][:] = np.squeeze(month)

    #-- create output date variable
    t_date = np.zeros((n_dates))

    #-- days per month in a leap and a standard year
    #-- only difference is February (29 vs. 28)
    dpm_leap=np.array([31,29,31,30,31,30,31,31,30,31,30,31], dtype=np.float)
    dpm_stnd=np.array([31,28,31,30,31,30,31,31,30,31,30,31], dtype=np.float)

    #-- Rules in the Gregorian calendar for a year to be a leap year:
    #-- divisible by 4, but not by 100 unless divisible by 400
    #-- True length of the year is about 365.2422 days
    #-- Adding a leap day every four years ==> average 365.25
    #-- Subtracting a leap year every 100 years ==> average 365.24
    #-- Adding a leap year back every 400 years ==> average 365.2425
    #-- Subtracting a leap year every 4000 years ==> average 365.24225
    m4 = (cal_date['year'] % 4)
    m100 = (cal_date['year'] % 100)
    m400 = (cal_date['year'] % 400)
    m4000 = (cal_date['year'] % 4000)
    #-- find indices for standard years and leap years using criteria
    leap, = np.nonzero((m4 == 0) & (m100 != 0) | (m400 == 0) & (m4000 != 0))
    stnd, = np.nonzero((m4 != 0) | (m100 == 0) & (m400 != 0) | (m4000 == 0))

    #-- calculate the day of the year
    if DofY is not None:
        #-- if entered directly as an input
        #-- remove 1 so day 1 (Jan 1st) = 0.0 in decimal format
        cal_date['DofY'][:] = np.squeeze(DofY)-1
    else:
        #-- use calendar month and day of the month to calculate day of the year
        #-- month minus 1: January = 0, February = 1, etc (indice of month)
        #-- in decimal form: January = 0.0
        month_m1 = np.array(cal_date['month'],dtype=np.int) - 1

        #-- day of month
        if DAY is not None:
            #-- remove 1 so 1st day of month = 0.0 in decimal format
            cal_date['day'][:] = np.squeeze(DAY)-1.0
        else:
            #-- if not entering days as an input
            #-- will use the mid-month value
            cal_date['day'][leap] = dpm_leap[month_m1[leap]]/2.0
            cal_date['day'][stnd] = dpm_stnd[month_m1[stnd]]/2.0

        #-- create matrix with the lower half = 1
        #-- this matrix will be used in a matrix multiplication
        #-- to calculate the total number of days for prior months
        #-- the -1 will make the diagonal == 0
        #-- i.e. first row == all zeros and the
        #-- last row == ones for all but the last element
        mon_mat=np.tri(12,12,-1)
        #-- using a dot product to calculate total number of days
        #-- for the months before the input date
        #-- basically is sum(i*dpm)
        #-- where i is 1 for all months < the month of interest
        #-- and i is 0 for all months >= the month of interest
        #-- month of interest is zero as the exact days will be
        #-- used to calculate the date

        #-- calculate the day of the year for leap and standard
        #-- use total days of all months before date
        #-- and add number of days before date in month
        cal_date['DofY'][stnd] = cal_date['day'][stnd] + \
            np.dot(mon_mat[month_m1[stnd],:],dpm_stnd)
        cal_date['DofY'][leap] = cal_date['day'][leap] + \
            np.dot(mon_mat[month_m1[leap],:],dpm_leap)

    #-- hour of day (else is zero)
    if HOUR is not None:
        cal_date['hour'][:] = np.squeeze(HOUR)

    #-- minute of hour (else is zero)
    if MINUTE is not None:
        cal_date['minute'][:] = np.squeeze(MINUTE)

    #-- second in minute (else is zero)
    if SECOND is not None:
        cal_date['second'][:] = np.squeeze(SECOND)

    #-- calculate decimal date
    #-- convert hours, minutes and seconds into days
    #-- convert calculated fractional days into decimal fractions of the year
    #-- Leap years
    t_date[leap] = cal_date['year'][leap] + \
        (cal_date['DofY'][leap] + cal_date['hour'][leap]/24. + \
        cal_date['minute'][leap]/1440. + \
        cal_date['second'][leap]/86400.)/np.sum(dpm_leap)
    #-- Standard years
    t_date[stnd] = cal_date['year'][stnd] + \
        (cal_date['DofY'][stnd] + cal_date['hour'][stnd]/24. + \
        cal_date['minute'][stnd]/1440. + \
        cal_date['second'][stnd]/86400.)/np.sum(dpm_stnd)

    return t_date
