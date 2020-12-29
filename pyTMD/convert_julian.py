#!/usr/bin/env python
u"""
convert_julian.py
Written by Tyler Sutterley (12/2020)

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
    Updated 12/2020: function deprecated.  merged with pyTMD.time
    Updated 07/2020: added function docstrings
    Updated 10/2017: updated comments and formatting of algorithm
    Updated 06/2017: added option FORMAT to change the output variables format
    Updated 06/2016: added option to convert output to variable type (e.g. int)
    Updated 11/2015: extracting the values from singleton dimension cases
    Updated 03/2015: remove singleton dimensions if initially importing value
    Updated 03/2014: updated to be able to convert arrays
    Written 05/2013
"""
import warnings
import pyTMD.time

def convert_julian(*args,**kwargs):
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.time instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.time.convert_julian(*args,**kwargs)