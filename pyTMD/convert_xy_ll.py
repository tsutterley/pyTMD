#!/usr/bin/env python
u"""
convert_xy_ll.py (07/2020)
Wrapper function to convert lat/lon points to and from projected coordinates
"""
import warnings
from pyTMD.convert_ll_xy import convert_ll_xy

def convert_xy_ll(*args,**kwargs):
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use convert_ll_xy() instead",DeprecationWarning)
    # call renamed version to not break workflows
    return convert_ll_xy(*args)
