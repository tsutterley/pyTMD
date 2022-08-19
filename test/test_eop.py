#!/usr/bin/env python
u"""
test_eop.py (09/2020)
Verify Earth Orientation Parameter (EOP) functions
"""
import os
import pytest
import warnings
import numpy as np
import scipy.interpolate
import pyTMD.eop
import pyTMD.time
from pyTMD.iers_mean_pole import iers_mean_pole
from pyTMD.read_iers_EOP import read_iers_EOP

# PURPOSE: update mean pole values
def test_update_mean_pole():
    pyTMD.eop.update_mean_pole(verbose=True, mode=0o775)
    LOCAL = pyTMD.utilities.get_data_path(['data','mean-pole.tab'])
    assert os.access(LOCAL, os.F_OK)

# PURPOSE: calculate updated mean pole values
def test_calculate_mean_pole():
    pyTMD.eop.calculate_mean_pole(verbose=True, mode=0o775)
    LOCAL = pyTMD.utilities.get_data_path(['data','mean-pole.tab'])
    assert os.access(LOCAL, os.F_OK)

# PURPOSE: update EOP finals values
def test_update_finals(username, password):
    pyTMD.eop.update_finals_file(username=username, password=password,
        verbose=True, mode=0o775)
    LOCAL = pyTMD.utilities.get_data_path(['data','finals.all'])
    assert os.access(LOCAL, os.F_OK)

# PURPOSE: read mean pole values
@pytest.mark.parametrize("EPOCH", ['2003','2010','2015'])
def test_read_EOP(EPOCH):
    # convert dates to Modified Julian days (days since 1858-11-17T00:00:00)
    delta_time = 86400.0*np.arange(0,365)
    MJD = pyTMD.time.convert_delta_time(delta_time, epoch1=(2000,1,1,0,0,0),
        epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
    # add offset to convert to Julian days and then convert to calendar dates
    Y,M,D,h,m,s = pyTMD.time.convert_julian(2400000.5 + MJD, format='tuple')
    # calculate time in year-decimal format
    time_decimal = pyTMD.time.convert_calendar_decimal(Y,M,day=D,
        hour=h,minute=m,second=s)
    # mean and daily EOP files
    mean_pole_file = pyTMD.utilities.get_data_path(['data','mean-pole.tab'])
    pole_tide_file = pyTMD.utilities.get_data_path(['data','finals.all'])
    # calculate angular coordinates of mean pole at time
    # iterate over different IERS conventional mean pole (CMP) formulations
    mpx,mpy,fl = iers_mean_pole(mean_pole_file,time_decimal,EPOCH)
    # check flags
    assert np.all(fl)
    # read IERS daily polar motion values
    EOP = read_iers_EOP(pole_tide_file)
    # check validity
    assert np.all(np.isfinite(EOP['x'])) & np.all(np.isfinite(EOP['y']))
    # interpolate daily polar motion values to time using cubic splines
    xSPL = scipy.interpolate.UnivariateSpline(EOP['MJD'],EOP['x'],k=3,s=0)
    ySPL = scipy.interpolate.UnivariateSpline(EOP['MJD'],EOP['y'],k=3,s=0)
    px = xSPL(MJD)
    py = ySPL(MJD)
    # calculate differentials from mean pole positions
    mx = px - mpx
    my = -(py - mpy)
    # check validity of differentials
    assert np.all(np.isfinite(mx)) & np.all(np.isfinite(my))
