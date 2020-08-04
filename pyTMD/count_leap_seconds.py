#!/usr/bin/env python
u"""
count_leap_seconds.py (08/2020)
Count number of leap seconds that have passed for each GPS time
"""
import warnings
import pyTMD.time

def count_leap_seconds(*args,**kwargs):
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.time instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.time.count_leap_seconds(*args)
