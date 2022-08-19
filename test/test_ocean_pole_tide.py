#!/usr/bin/env python
u"""
test_ocean_pole_tide.py (08/2020)
"""
import os
import re
import inspect
import warnings
import pytest
import numpy as np
import scipy.interpolate
from pyTMD.utilities import get_data_path
from pyTMD.iers_mean_pole import iers_mean_pole
from pyTMD.read_iers_EOP import read_iers_EOP
from pyTMD.read_ocean_pole_tide import read_ocean_pole_tide

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = os.path.dirname(os.path.abspath(filename))

# parameterize interpolation method
@pytest.mark.parametrize("METHOD", ['spline','nearest','linear'])
# test the interpolation of ocean pole tide values
def test_ocean_pole_tide(METHOD):
    # degrees to radians and arcseconds to radians
    dtr = np.pi/180.0
    atr = np.pi/648000.0
    # earth and physical parameters (IERS)
    # G = 6.67428e-11# universal constant of gravitation [m^3/(kg*s^2)]
    G = 6.673e-11# test universal constant of gravitation [m^3/(kg*s^2)]
    GM = 3.986004418e14# geocentric gravitational constant [m^3/s^2]
    a_axis = 6378136.6# equatorial radius of the Earth [m]
    flat = 1.0/298.257223563# flattening of the ellipsoid
    omega = 7.292115e-5# mean rotation rate of the Earth [arcseconds/s]
    rho_w = 1025.0# density of sea water [kg/m^3]
    ge = 9.7803278# mean equatorial gravitational acceleration [m/s^2]
    # Linear eccentricity and first numerical eccentricity
    lin_ecc = np.sqrt((2.0*flat - flat**2)*a_axis**2)
    ecc1 = lin_ecc/a_axis
    # tidal love number differential (1 + kl - hl) for pole tide frequencies
    gamma = 0.6870 + 0.0036j

    # pole tide displacement scale factor
    Hp = np.sqrt(8.0*np.pi/15.0)*(omega**2*a_axis**4)/GM
    K = 4.0*np.pi*G*rho_w*Hp*a_axis/(3.0*ge)
    K1 = 4.0*np.pi*G*rho_w*Hp*a_axis**3/(3.0*GM)
    # determine differences with values from test data
    eps = np.finfo(np.float16).eps
    assert (np.abs(Hp - 2.8577142980e+04) < eps)
    assert (np.abs(K - 5.3394043696e+03) < eps)

    # read test file for values
    ocean_pole_test_file = os.path.join(filepath,'opoleloadcmcor.test')
    names = ('MJD','xbar_p','ybar_p','x_p','y_p','m1','m2','u_radial','u_north','u_east')
    formats = ('i4','f4','f4','f4','f4','f4','f4','f4','f4','f4')
    validation = np.loadtxt(ocean_pole_test_file, skiprows=26,
        dtype=dict(names=names, formats=formats))
    file_lines = len(validation)
    # mean pole coordinates for test
    xmean = np.array([5.4e-2, 8.30e-4])
    ymean = np.array([3.57e-1, 3.95e-3])
    t0 = 5.1544e4
    # coordinates for test
    lon,lat = (232.25,-43.75)

    # read ocean pole tide map from Desai (2002)
    ocean_pole_tide_file = get_data_path(['data','opoleloadcoefcmcor.txt.gz'])
    iur,iun,iue,ilon,ilat = read_ocean_pole_tide(ocean_pole_tide_file)

    # interpolate ocean pole tide map from Desai (2002)
    if (METHOD == 'spline'):
        # use scipy bivariate splines to interpolate to output points
        f1 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
            iur[:,::-1].real, kx=1, ky=1)
        f2 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
            iur[:,::-1].imag, kx=1, ky=1)
        UR = np.zeros((file_lines),dtype=np.longcomplex)
        UR.real = f1.ev(lon,lat)
        UR.imag = f2.ev(lon,lat)
    else:
        # use scipy regular grid to interpolate values for a given method
        r1 = scipy.interpolate.RegularGridInterpolator((ilon, ilat[::-1]),
            iur[:,::-1], method=METHOD)
        UR = r1.__call__(np.c_[lon,lat])

    # calculate angular coordinates of mean pole at time
    mpx = xmean[0] + xmean[1]*(validation['MJD'] - t0)/365.25
    mpy = ymean[0] + ymean[1]*(validation['MJD'] - t0)/365.25

    # calculate differentials from mean pole positions
    mx = validation['x_p'] - mpx
    my = -(validation['y_p'] - mpy)
    # calculate radial displacement at time
    u_radial = K*atr*np.real((mx*gamma.real + my*gamma.imag)*UR.real +
        (my*gamma.real - mx*gamma.imag)*UR.imag)
    assert np.all((u_radial - validation['u_radial']) < eps)
