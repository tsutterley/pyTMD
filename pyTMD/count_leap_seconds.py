#!/usr/bin/env python
u"""
count_leap_seconds.py (07/2020)
Count number of leap seconds that have passed for each GPS time
Based partially on Tiffany Summerscales's PHP conversion algorithm
    https://www.andrews.edu/~tzs/timeconv/timealgorithm.html

INPUTS:
    GPS_Time: GPS time (standard = seconds since January 6, 1980 at 00:00)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (https://numpy.org)

UPDATE HISTORY:
    Updated 07/2020: added function docstrings
    Written 10/2017
"""
import numpy as np

#-- PURPOSE: Define GPS leap seconds
def get_leaps():
    """
    Gets a list of GPS times for when leap seconds occurred
    """
    leaps = [46828800, 78364801, 109900802, 173059203, 252028804, 315187205,
        346723206, 393984007, 425520008, 457056009, 504489610, 551750411,
        599184012, 820108813, 914803214, 1025136015, 1119744016, 1167264017]
    return leaps

#-- PURPOSE: Count number of leap seconds that have passed for each GPS time
def count_leap_seconds(GPS_Time):
    """
    Counts the number of leap seconds between a given GPS time and UTC

    Arguments
    ---------
    GPS_Time: seconds since January 6, 1980 at 00:00:00

    Returns
    -------
    n_leaps: number of elapsed leap seconds
    """
    leaps = get_leaps()
    #-- number of leap seconds prior to GPS_Time
    n_leaps = np.zeros_like(GPS_Time)
    for i,leap in enumerate(leaps):
        count = np.count_nonzero(GPS_Time >= leap)
        if (count > 0):
            indices, = np.nonzero(GPS_Time >= leap)
            n_leaps[indices] += 1.0
    return n_leaps
