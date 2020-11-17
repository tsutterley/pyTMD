#!/usr/bin/env python
u"""
test_eop.py (09/2020)
Verify Earth Orientation Parameter (EOP) functions
"""
import os
import pytest
import warnings
import numpy as np
import pyTMD.eop
from pyTMD.iers_mean_pole import iers_mean_pole
from pyTMD.read_iers_EOP import read_iers_EOP

#-- PURPOSE: update mean pole values
def test_update_mean_pole():
    pyTMD.eop.update_mean_pole(verbose=True, mode=0o775)
    LOCAL = pyTMD.utilities.get_data_path(['data','mean-pole.tab'])
    assert os.access(LOCAL, os.F_OK)

#-- PURPOSE: calculate updated mean pole values
def test_calculate_mean_pole():
    pyTMD.eop.calculate_mean_pole(verbose=True, mode=0o775)
    LOCAL = pyTMD.utilities.get_data_path(['data','mean-pole.tab'])
    assert os.access(LOCAL, os.F_OK)

#-- PURPOSE: calculate updated mean pole values
def test_update_finals(username, password):
    pyTMD.eop.update_finals_file(username=username, password=password,
        verbose=True, mode=0o775)
    LOCAL = pyTMD.utilities.get_data_path(['data','finals.all'])
    assert os.access(LOCAL, os.F_OK)
