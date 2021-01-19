#!/usr/bin/env python
u"""
convert_calendar_decimal.py
Written by Tyler Sutterley (12/2020)

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
    Updated 12/2020: function deprecated.  merged with pyTMD.time
    Updated 07/2020: added function docstrings
    Updated 05/2015: updated comments and minor update to nonzero statement
    Updated 05/2014: added option for day of year
    Updated 04/2014: new code from convert_J2000.py
    Updated 04/2014: updated comments and improved rules
        for leap years to include mod 100 and mod 400
    Written 04/2014
"""
import warnings
import pyTMD.time

def convert_calendar_decimal(*args,**kwargs):
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.time instead",DeprecationWarning)
    # call renamed version to not break workflows
    kwds={key.lower():val for key,val in kwargs.items()}
    return pyTMD.time.convert_calendar_decimal(*args,**kwds)