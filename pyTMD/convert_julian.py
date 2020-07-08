#!/usr/bin/env python
u"""
convert_julian.py
Written by Tyler Sutterley (07/2020)

Return the calendar date and time given Julian date.

CALLING SEQUENCE:
    YEAR,MONTH,DAY,HOUR,MINUTE,SECOND = convert_julian(JD, FORMAT='tuple')

INPUTS:
    JD: Julian Day of the specified calendar date.

OUTPUTS:
    year: Number of the desired year
    month: Number of the desired month (1 = January, ..., 12 = December)
    day: Number of day of the month
    hour: hour of the day
    minute: minute of the hour
    second: second (and fractions of a second) of the minute

OPTIONS:
    ASTYPE: convert output to variable type (e.g. int).  Default is float
    FORMAT: format of output variables
        'dict': dictionary with variable keys as listed above
        'tuple': tuple with variable order YEAR,MONTH,DAY,HOUR,MINUTE,SECOND
        'zip': aggregated variable sets

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

NOTES:
    Translated from caldat in "Numerical Recipes in C", by William H. Press,
        Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
        Cambridge University Press, 1988 (second printing).
    Hatcher, D. A., "Simple Formulae for Julian Day Numbers and Calendar Dates",
        Quarterly Journal of the Royal Astronomical Society, 25(1), 1984.

UPDATE HISTORY:
    Updated 07/2020: added function docstrings
    Updated 10/2017: updated comments and formatting of algorithm
    Updated 06/2017: added option FORMAT to change the output variables format
    Updated 06/2016: added option to convert output to variable type (e.g. int)
    Updated 11/2015: extracting the values from singleton dimension cases
    Updated 03/2015: remove singleton dimensions if initially importing value
    Updated 03/2014: updated to be able to convert arrays
    Written 05/2013
"""
import numpy as np

def convert_julian(JD, ASTYPE=None, FORMAT='dict'):
    """
    Converts from Julian day to calendar date and time

    Arguments
    ---------
    JD: Julian Day (days since 01-01-4713 BCE at 12:00:00)

    Keyword arguments
    -----------------
    ASTYPE: convert output to variable type
    FORMAT: format of output variables
        'dict': dictionary with variable keys
        'tuple': tuple with variable order YEAR,MONTH,DAY,HOUR,MINUTE,SECOND
        'zip': aggregated variable sets

    Returns
    -------
    year: calendar year
    month: calendar month
    day: day of the month
    hour: hour of the day
    minute: minute of the hour
    second: second of the minute
    """

    #-- convert to array if only a single value was imported
    if (np.ndim(JD) == 0):
        JD = np.array([JD])
        SINGLE_VALUE = True
    else:
        SINGLE_VALUE = False

    JDO = np.floor(JD + 0.5)
    C = np.zeros_like(JD)
    #-- calculate C for dates before and after the switch to Gregorian
    IGREG = 2299161.0
    ind1, = np.nonzero(JDO < IGREG)
    C[ind1] = JDO[ind1] + 1524.0
    ind2, = np.nonzero(JDO >= IGREG)
    B = np.floor((JDO[ind2] - 1867216.25)/36524.25)
    C[ind2] = JDO[ind2] + B - np.floor(B/4.0) + 1525.0
    #-- calculate coefficients for date conversion
    D = np.floor((C - 122.1)/365.25)
    E = np.floor((365.0 * D) + np.floor(D/4.0))
    F = np.floor((C - E)/30.6001)
    #-- calculate day, month, year and hour
    DAY = np.floor(C - E + 0.5) - np.floor(30.6001*F)
    MONTH = F - 1.0 - 12.0*np.floor(F/14.0)
    YEAR = D - 4715.0 - np.floor((7.0+MONTH)/10.0)
    HOUR = np.floor(24.0*(JD + 0.5 - JDO))
    #-- calculate minute and second
    G = (JD + 0.5 - JDO) - HOUR/24.0
    MINUTE = np.floor(G*1440.0)
    SECOND = (G - MINUTE/1440.0) * 86400.0

    #-- convert all variables to output type (from float)
    if ASTYPE is not None:
        YEAR = YEAR.astype(ASTYPE)
        MONTH = MONTH.astype(ASTYPE)
        DAY = DAY.astype(ASTYPE)
        HOUR = HOUR.astype(ASTYPE)
        MINUTE = MINUTE.astype(ASTYPE)
        SECOND = SECOND.astype(ASTYPE)

    #-- if only a single value was imported initially: remove singleton dims
    if SINGLE_VALUE:
        YEAR = YEAR.item(0)
        MONTH = MONTH.item(0)
        DAY = DAY.item(0)
        HOUR = HOUR.item(0)
        MINUTE = MINUTE.item(0)
        SECOND = SECOND.item(0)

    #-- return date variables in output format (default python dictionary)
    if (FORMAT == 'dict'):
        return dict(year=YEAR, month=MONTH, day=DAY,
            hour=HOUR, minute=MINUTE, second=SECOND)
    elif (FORMAT == 'tuple'):
        return (YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
    elif (FORMAT == 'zip'):
        return zip(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
