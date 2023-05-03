#!/usr/bin/env python
u"""
test_ocean_pole_tide.py
Written by Tyler Sutterley (04/2023)

UPDATE HISTORY:
    Updated 04/2023: using pathlib to define and expand paths
    Updated 12/2022: single implicit import of pyTMD
        use constants class for ellipsoidal parameters
        read header from test file and compare more variables
    Updated 07/2022: define variable formats of test inputs
    Updated 05/2022: change to longcomplex for windows compatibility
    Written 08/2020
"""
import re
import inspect
import pathlib
import pytest
import numpy as np
import scipy.interpolate
import pyTMD

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent

# parameterize interpolation method
@pytest.mark.parametrize("METHOD", ['spline','nearest','linear'])
# test the interpolation of ocean pole tide values
def test_ocean_pole_tide(METHOD):
    # read ocean pole tide test file for header text
    ocean_pole_test_file = filepath.joinpath('opoleloadcmcor.test')
    with ocean_pole_test_file.open(mode='r', encoding='utf8') as f:
        file_contents = f.read().splitlines()

    # extract header parameters
    rx = re.compile(r'(.*?)\s+\=\s+([-+]?\d+\.\d+[e][-+]?\d+)', re.VERBOSE)
    header_contents = [l for l in file_contents if rx.match(l)]
    test_header = {}
    for line in header_contents:
        param,contents = rx.findall(line).pop()
        test_header[param] = np.float64(contents)

    # extract longitude and latitude from header
    lat = re.findall(r'latitude = ([-+]?\d+\.\d+) degrees', file_contents[1]).pop()
    lon = re.findall(r'longitude = ([-+]?\d+\.\d+) degrees', file_contents[1]).pop()
    test_header['latitude'] = np.float64(lat)
    test_header['longitude'] = np.float64(lon)

    # degrees and arcseconds to radians
    dtr = np.pi/180.0
    atr = np.pi/648000.0
    # earth and physical parameters for ellipsoid
    units = pyTMD.constants('IERS')
    # universal constant of gravitation for test [m^3/(kg*s^2)]
    G = 6.673e-11
    # mean equatorial gravitational acceleration [m/s^2]
    ge = 9.7803278
    # density of sea water [kg/m^3]
    rho_w = 1025.0
    # tidal love number differential (1 + kl - hl) for pole tide frequencies
    gamma = 0.6870 + 0.0036j

    # pole tide displacement scale factor
    Hp = np.sqrt(8.0*np.pi/15.0)*(units.omega**2*units.a_axis**4)/units.GM
    K = 4.0*np.pi*G*rho_w*Hp*units.a_axis/(3.0*ge)
    K1 = 4.0*np.pi*G*rho_w*Hp*units.a_axis**3/(3.0*units.GM)

    # determine differences with constants from test data
    eps = np.finfo(np.float16).eps
    assert np.isclose(units.a_axis, test_header['aE'])
    assert np.isclose(units.omega, test_header['Omega'])
    assert np.isclose(units.GM, test_header['GM'])
    assert np.isclose(rho_w, test_header['rhow'])
    assert np.isclose(gamma.real, test_header['gamma2(real)'])
    assert np.isclose(gamma.imag, test_header['gamma2(imag)'])
    assert np.isclose(G, test_header['G'])
    assert np.isclose(Hp, test_header['Hp'])
    assert np.isclose(K, test_header['K'])

    # read test file for values
    names = ('MJD','xbar_p','ybar_p','x_p','y_p','m1','m2','u_radial','u_north','u_east')
    formats = ('i4','f4','f4','f4','f4','f4','f4','f4','f4','f4')
    validation = np.loadtxt(ocean_pole_test_file, skiprows=26,
        dtype=dict(names=names, formats=formats))
    file_lines = len(validation)
    # mean pole coordinates for test
    xmean = np.array([test_header['xmean(t0)'], test_header['xmeandot']])
    ymean = np.array([test_header['ymean(t0)'], test_header['ymeandot']])

    # read ocean pole tide map from Desai (2002)
    ocean_pole_tide_file = pyTMD.utilities.get_data_path(
        ['data','opoleloadcoefcmcor.txt.gz'])
    iur, iun, iue, ilon, ilat = pyTMD.io.ocean_pole_tide(ocean_pole_tide_file)

    # interpolate ocean pole tide map from Desai (2002)
    if (METHOD == 'spline'):
        # use scipy bivariate splines to interpolate to output points
        f1 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
            iur[:,::-1].real, kx=1, ky=1)
        f2 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
            iur[:,::-1].imag, kx=1, ky=1)
        UR = np.zeros((file_lines),dtype=np.longcomplex)
        UR.real = f1.ev(test_header['longitude'], test_header['latitude'])
        UR.imag = f2.ev(test_header['longitude'], test_header['latitude'])
    else:
        # use scipy regular grid to interpolate values for a given method
        r1 = scipy.interpolate.RegularGridInterpolator((ilon, ilat[::-1]),
            iur[:,::-1], method=METHOD)
        UR = r1.__call__(np.c_[test_header['longitude'], test_header['latitude']])

    # calculate differences in time
    dt = (validation['MJD'] - test_header['t0'])/test_header['1 year']
    # calculate angular coordinates of mean pole at time
    mpx = xmean[0] + xmean[1]*dt
    mpy = ymean[0] + ymean[1]*dt

    # calculate differentials from mean pole positions
    mx = validation['x_p'] - mpx
    my = -(validation['y_p'] - mpy)
    # calculate radial displacement at time
    u_radial = K*atr*np.real((mx*gamma.real + my*gamma.imag)*UR.real +
        (my*gamma.real - mx*gamma.imag)*UR.imag)
    assert np.all((u_radial - validation['u_radial']) < eps)
