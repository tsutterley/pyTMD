#!/usr/bin/env python
u"""
test_leap_seconds.py (08/2020)
"""
import pytest
import pyTMD.time
import pyTMD.utilities

# PURPOSE: Define GPS leap seconds
def get_leaps():
    """
    Gets a list of GPS times for when leap seconds occurred
    """
    leaps = [46828800, 78364801, 109900802, 173059203, 252028804, 315187205,
        346723206, 393984007, 425520008, 457056009, 504489610, 551750411,
        599184012, 820108813, 914803214, 1025136015, 1119744016, 1167264017]
    return leaps

# PURPOSE: Count number of leap seconds that have passed for each GPS time
def test_leap_seconds():
    valid_gps_leaps = get_leaps()
    test_gps_leaps = pyTMD.time.get_leap_seconds()
    assert all((v==t) for v,t in zip(valid_gps_leaps,test_gps_leaps))

# PURPOSE: update leap second file
def test_update_leap_seconds():
    pyTMD.time.update_leap_seconds(verbose=False, mode=0o775)
    FILE = pyTMD.utilities.get_data_path(['data','leap-seconds.list'])
    assert FILE.exists()
