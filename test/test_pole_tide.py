#!/usr/bin/env python
u"""
test_pole_tide.py
Written by Tyler Sutterley (08/2024)

UPDATE HISTORY:
    Updated 08/2024: add tests for new cartesian pole tides
    Updated 06/2024: use np.clongdouble instead of np.longcomplex
    Updated 02/2024: changed class name for ellipsoid parameters to datum
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
import pyproj
import numpy as np
import scipy.interpolate
import pyTMD.compute
import pyTMD.predict
import timescale.time

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent

# number of days between the Julian day epoch and MJD
_jd_mjd = 2400000.5
# number of days between MJD and the tide epoch (1992-01-01T00:00:00)
_mjd_tide = 48622.0
# number of days between the Julian day epoch and the tide epoch
_jd_tide = _jd_mjd + _mjd_tide

# PURPOSE: compute radial load pole tide displacements
# following IERS Convention (2010) guidelines
@pytest.mark.parametrize("TYPE", ['grid','drift','time series'])
def test_load_pole_tide_displacements(TYPE):

    # create a test dataset for data type
    if (TYPE == 'drift'):
        # number of data points
        n_time = 3000
        y = np.random.randint(-90,90,size=n_time)
        x = np.random.randint(-180,180,size=n_time)
        delta_time = np.random.randint(0,31557600,size=n_time)
    elif (TYPE == 'grid'):
        # number of data points
        n_lat,n_lon,n_time = (181,361,100)
        y = np.linspace(-90,90,n_lat)
        x = np.linspace(-180,180,n_lon)
        delta_time = np.random.randint(0,31557600,size=n_time)
    elif (TYPE == 'time series'):
        n_station,n_time = (300,100)
        y = np.random.randint(-90,90,size=n_station)
        x = np.random.randint(-180,180,size=n_station)
        delta_time = np.random.randint(0,31557600,size=n_time)

    # parameters for ocean pole tide
    EPSG = 4326
    EPOCH = (2018, 1, 1, 0, 0, 0)
    TIME = 'UTC'
    ELLIPSOID = 'WGS84'
    CONVENTION = '2018'
    FILL_VALUE = np.nan

    # compute load pole tides using cartesian functions
    test = pyTMD.compute.LPT_displacements(x, y, delta_time,
        EPSG=EPSG, EPOCH=EPOCH, TYPE=TYPE, TIME=TIME,
        ELLIPSOID=ELLIPSOID, CONVENTION=CONVENTION,
        FILL_VALUE=FILL_VALUE)

    # reform coordinate dimensions for input grids
    # or verify coordinate dimension shapes
    if (TYPE.lower() == 'grid') and (np.size(x) != np.size(y)):
        x,y = np.meshgrid(np.copy(x),np.copy(y))
    elif (TYPE.lower() == 'grid'):
        x = np.atleast_2d(x)
        y = np.atleast_2d(y)
    elif TYPE.lower() in ('time series', 'drift'):
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)

    # converting x,y from EPSG to latitude/longitude
    crs1 = pyTMD.crs().from_input(EPSG)
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon,lat = transformer.transform(x.flatten(), y.flatten())

    # assert delta time is an array
    delta_time = np.atleast_1d(delta_time)
    # convert delta times or datetimes objects to timescale
    if (TIME.lower() == 'datetime'):
        ts = timescale.time.Timescale().from_datetime(
            delta_time.flatten())
    else:
        ts = timescale.time.Timescale().from_deltatime(delta_time,
            epoch=EPOCH, standard=TIME)

    # convert dynamic time to Modified Julian Days (MJD)
    MJD = ts.tt - _jd_mjd
    # convert Julian days to calendar dates
    Y,M,D,h,m,s = timescale.time.convert_julian(ts.tt, format='tuple')
    # calculate time in year-decimal format
    time_decimal = timescale.time.convert_calendar_decimal(Y, M, day=D,
        hour=h, minute=m, second=s)
    # number of time points
    nt = len(time_decimal)

    # degrees and arcseconds to radians
    dtr = np.pi/180.0
    atr = np.pi/648000.0
    # earth and physical parameters for ellipsoid
    units = pyTMD.datum(ellipsoid=ELLIPSOID, units='MKS')
    # tidal love number appropriate for the load tide
    hb2 = 0.6207

    # convert from geodetic latitude to geocentric latitude
    # calculate X, Y and Z from geodetic latitude and longitude
    X,Y,Z = pyTMD.spatial.to_cartesian(lon.flatten(), lat.flatten(),
        a_axis=units.a_axis, flat=units.flat)
    rr = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/dtr
    # geocentric colatitude and longitude in radians
    theta = dtr*(90.0 - latitude_geocentric)
    phi = dtr*lon.flatten()

    # compute normal gravity at spatial location
    gamma_0 = units.gamma_0(theta)
    dfactor = -hb2*atr*(units.omega**2*rr**2)/(2.0*gamma_0)

    # calculate angular coordinates of mean/secular pole at time
    mpx, mpy, fl = timescale.eop.iers_mean_pole(time_decimal, convention=CONVENTION)
    # read and interpolate IERS daily polar motion values
    px, py = timescale.eop.iers_polar_motion(MJD, k=3, s=0)
    # calculate differentials from mean/secular pole positions
    mx = px - mpx
    my = -(py - mpy)

    # calculate radial displacement at time
    if (TYPE == 'grid'):
        ny,nx = np.shape(x)
        Srad = np.ma.zeros((ny,nx,nt), fill_value=FILL_VALUE)
        Srad.mask = np.zeros((ny,nx,nt),dtype=bool)
        approx = np.zeros((ny,nx,nt))
        for i in range(nt):
            SRAD = dfactor*np.sin(2.0*theta)*(mx[i]*np.cos(phi)+my[i]*np.sin(phi))
            # reform grid
            Srad.data[:,:,i] = np.reshape(SRAD, (ny, nx))
            Srad.mask[:,:,i] = np.isnan(Srad.data[:,:,i])
            # approximate values from IERS (2010) conventions
            S = -33.0*np.sin(2.0*theta)*(mx[i]*np.cos(phi) + my[i]*np.sin(phi))
            approx[:,:,i] = np.reshape(S, (ny, nx))/1e3
    elif (TYPE == 'drift'):
        Srad = np.ma.zeros((nt), fill_value=FILL_VALUE)
        Srad.data[:] = dfactor*np.sin(2.0*theta)*(mx*np.cos(phi)+my*np.sin(phi))
        Srad.mask = np.isnan(Srad.data)
        # approximate values from IERS (2010) conventions
        S = -33.0*np.sin(2.0*theta)*(mx*np.cos(phi) + my*np.sin(phi))
        approx = S/1e3
    elif (TYPE == 'time series'):
        nstation = len(x)
        Srad = np.ma.zeros((nstation,nt), fill_value=FILL_VALUE)
        Srad.mask = np.zeros((nstation,nt),dtype=bool)
        approx = np.zeros((nstation,nt))
        for s in range(nstation):
            SRAD = dfactor[s]*np.sin(2.0*theta[s])*(mx*np.cos(phi[s])+my*np.sin(phi[s]))
            Srad.data[s,:] = np.copy(SRAD)
            Srad.mask[s,:] = np.isnan(Srad.data[s,:])
            # approximate values from IERS (2010) conventions
            S = -33.0*np.sin(2.0*theta[s])*(mx*np.cos(phi[s]) + my*np.sin(phi[s]))
            approx[s,:] = np.copy(S)/1e3
    # replace invalid data with fill values
    Srad.data[Srad.mask] = Srad.fill_value
    # compare with functional values
    eps = np.finfo(np.float16).eps
    assert np.all(np.abs(Srad - test) < eps)
    assert np.all(np.abs(approx - test) < eps)

# parameterize interpolation method
@pytest.mark.parametrize("METHOD", ['spline','nearest','linear'])
# test the interpolation of ocean pole tide values
def test_read_ocean_pole(METHOD):
    # read ocean pole tide test file for header text
    ocean_pole_test_file = filepath.joinpath('opoleloadcmcor.test')
    with ocean_pole_test_file.open(mode='r', encoding='utf8') as f:
        file_contents = f.read().splitlines()

    # extract header parameters
    rx = re.compile(r'(.*?)\s+\=\s+([-+]?\d+\.\d+[e][-+]?\d+)', re.VERBOSE)
    header_contents = [l for l in file_contents if rx.match(l)]
    header = {}
    for line in header_contents:
        param,contents = rx.findall(line).pop()
        header[param] = np.float64(contents)

    # extract longitude and latitude from header
    lat = re.findall(r'latitude = ([-+]?\d+\.\d+) degrees', file_contents[1]).pop()
    lon = re.findall(r'longitude = ([-+]?\d+\.\d+) degrees', file_contents[1]).pop()
    header['latitude'] = np.float64(lat)
    # convert longitude from 0:360 to -180:180
    header['longitude'] = np.float64(lon) - 360.0

    # number of points
    npts = len(np.atleast_1d(header['longitude']))
    # read ocean pole tide map from Desai (2002)
    Umap = {}
    Umap['R'], Umap['N'], Umap['E'], ilon, ilat = pyTMD.io.IERS.read_binary_file()
    # interpolate ocean pole tide map from Desai (2002)
    Uint = {}
    if (METHOD == 'spline'):
        # use scipy bivariate splines to interpolate to output points
        for key,val in Umap.items():
            Uint[key] = np.zeros((npts), dtype=np.clongdouble)
            f1 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
                val[:,::-1].real, kx=1, ky=1)
            f2 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
                val[:,::-1].imag, kx=1, ky=1)
            Uint[key].real = f1.ev(header['longitude'], header['latitude'])
            Uint[key].imag = f2.ev(header['longitude'], header['latitude'])
    else:
        # use scipy regular grid to interpolate values for a given method
        for key,val in Umap.items():
            Uint[key] = np.zeros((npts), dtype=np.clongdouble)
            r1 = scipy.interpolate.RegularGridInterpolator((ilon,ilat[::-1]),
                val[:,::-1], bounds_error=False, method=METHOD)
            Uint[key][:] = r1.__call__(np.c_[header['longitude'], header['latitude']])

    # extract coefficients from IERS pole tide map
    U = pyTMD.io.IERS.extract_coefficients(
        header['longitude'], header['latitude'],
        method=METHOD)

    # compare with functional values
    for i, key in enumerate(['R','N','E']):
        assert np.isclose(Uint[key], U[i])

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
    header = {}
    for line in header_contents:
        param,contents = rx.findall(line).pop()
        header[param] = np.float64(contents)

    # extract longitude and latitude from header
    lat = re.findall(r'latitude = ([-+]?\d+\.\d+) degrees', file_contents[1]).pop()
    lon = re.findall(r'longitude = ([-+]?\d+\.\d+) degrees', file_contents[1]).pop()
    header['latitude'] = np.float64(lat)
    # convert longitude from 0:360 to -180:180
    header['longitude'] = np.float64(lon) - 360.0

    # degrees and arcseconds to radians
    dtr = np.pi/180.0
    atr = np.pi/648000.0
    # earth and physical parameters for ellipsoid
    units = pyTMD.datum(ellipsoid='IERS', units='MKS')
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

    # determine differences with constants from test data
    eps = np.finfo(np.float16).eps
    assert np.isclose(units.a_axis, header['aE'])
    assert np.isclose(units.omega, header['Omega'])
    assert np.isclose(units.GM, header['GM'])
    assert np.isclose(rho_w, header['rhow'])
    assert np.isclose(gamma.real, header['gamma2(real)'])
    assert np.isclose(gamma.imag, header['gamma2(imag)'])
    assert np.isclose(G, header['G'])
    assert np.isclose(Hp, header['Hp'])
    assert np.isclose(K, header['K'])

    # read test file for values
    names = ('MJD','xbar_p','ybar_p','x_p','y_p','m1','m2','u_radial','u_north','u_east')
    formats = ('i4','f4','f4','f4','f4','f4','f4','f4','f4','f4')
    validation = np.loadtxt(ocean_pole_test_file, skiprows=26,
        dtype=dict(names=names, formats=formats))
    file_lines = len(validation)
    # mean pole coordinates for test
    xmean = np.array([header['xmean(t0)'], header['xmeandot']])
    ymean = np.array([header['ymean(t0)'], header['ymeandot']])

    # read ocean pole tide map from Desai (2002)
    ocean_pole_tide_file = pyTMD.utilities.get_data_path(
        ['data','opoleloadcoefcmcor.txt.gz'])
    # interpolate ocean pole tide map to coordinates
    iu = {}
    iu['u_radial'], iu['u_north'], iu['u_east'], = \
        pyTMD.io.IERS.extract_coefficients(
            header['longitude'],
            header['latitude'],
            model_file=ocean_pole_tide_file,
            method=METHOD
        )

    # create timescale object from MJD
    ts = timescale.time.Timescale(MJD=validation['MJD'])
    # calculate angular coordinates of mean/secular pole at time
    mpx, mpy, fl = timescale.eop.iers_mean_pole(ts.year,
        convention='2003')
    # read and interpolate IERS daily polar motion values
    px, py = timescale.eop.iers_polar_motion(ts.MJD, k=3, s=0)
    # calculate differentials from mean/secular pole positions
    # using the latest definition from IERS Conventions (2010)
    mx = px - mpx
    my = -(py - mpy)

    # calculate differences in time
    dt = (validation['MJD'] - header['t0'])/header['1 year']
    # calculate angular coordinates of mean pole at time
    mean_pole_x = xmean[0] + xmean[1]*dt
    mean_pole_y = ymean[0] + ymean[1]*dt
    # assert that the mean pole values are close
    assert np.all(np.abs(mpx - mean_pole_x) < eps)
    assert np.all(np.abs(mpy - mean_pole_y) < eps)

    # calculate differentials from mean pole positions
    x_bar = validation['x_p'] - mean_pole_x
    y_bar = -(validation['y_p'] - mean_pole_y)
    # assert that the mean pole differentials are close
    assert np.all(np.abs(mx - x_bar) < eps)
    assert np.all(np.abs(my - y_bar) < eps)

    # calculate ocean pole tide displacements at time
    for key, val in iu.items():
        u = K*atr*np.real((x_bar*gamma.real + y_bar*gamma.imag)*val.real +
            (y_bar*gamma.real - x_bar*gamma.imag)*val.imag)
        # assert that the predicted values are close
        assert np.all(np.abs(u - validation[key]) < eps)

# parameterize interpolation method
@pytest.mark.parametrize("METHOD", ['spline','nearest','linear'])
# test the interpolation of ocean pole tide values
def test_predict_ocean_pole_tide(METHOD):
    # read ocean pole tide test file for header text
    ocean_pole_test_file = filepath.joinpath('opoleloadcmcor.test')
    with ocean_pole_test_file.open(mode='r', encoding='utf8') as f:
        file_contents = f.read().splitlines()

    # extract header parameters
    rx = re.compile(r'(.*?)\s+\=\s+([-+]?\d+\.\d+[e][-+]?\d+)', re.VERBOSE)
    header_contents = [l for l in file_contents if rx.match(l)]
    header = {}
    for line in header_contents:
        param,contents = rx.findall(line).pop()
        header[param] = np.float64(contents)

    # extract longitude and latitude from header
    lat = re.findall(r'latitude = ([-+]?\d+\.\d+) degrees', file_contents[1]).pop()
    lon = re.findall(r'longitude = ([-+]?\d+\.\d+) degrees', file_contents[1]).pop()
    header['latitude'] = np.float64(lat)
    # convert longitude from 0:360 to -180:180
    header['longitude'] = np.float64(lon) - 360.0

    # read test file for values
    names = ('MJD','xbar_p','ybar_p','x_p','y_p','m1','m2','u_radial','u_north','u_east')
    formats = ('i4','f4','f4','f4','f4','f4','f4','f4','f4','f4')
    validation = np.loadtxt(ocean_pole_test_file, skiprows=26,
        dtype=dict(names=names, formats=formats))
    file_lines = len(validation)
    # create timescale object from MJD
    ts = timescale.time.Timescale(MJD=validation['MJD'])

    # degrees to radians
    dtr = np.pi/180.0
    # earth and physical parameters for ellipsoid
    units = pyTMD.datum(ellipsoid='IERS', units='MKS')
    # mean equatorial gravitational acceleration [m/s^2]
    ge = 9.7803278
    # density of sea water [kg/m^3]
    rho_w = 1025.0
    # tidal love number differential (1 + kl - hl) for pole tide frequencies
    gamma = 0.6870 + 0.0036j

    # read ocean pole tide map from Desai (2002)
    # interpolate ocean pole tide map to coordinates
    iu = {}
    iu['u_radial'], iu['u_north'], iu['u_east'], = \
        pyTMD.io.IERS.extract_coefficients(
            header['longitude'],
            header['latitude'],
            method=METHOD
        )

    # convert latitude and longitude to ECEF cartesian coordinates
    X, Y, Z = pyTMD.spatial.to_cartesian(header['longitude'], header['latitude'],
        a_axis=units.a_axis, flat=units.flat)
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/dtr
    # geocentric colatitude and longitude in radians
    theta = dtr*(90.0 - latitude_geocentric)
    phi = dtr*header['longitude']
    # convert pole tide values to cartesian coordinates
    R = np.zeros((file_lines, 3, 3))
    R[:,0,0] = np.cos(phi)*np.cos(theta)
    R[:,0,1] = -np.sin(phi)
    R[:,0,2] = np.cos(phi)*np.sin(theta)
    R[:,1,0] = np.sin(phi)*np.cos(theta)
    R[:,1,1] = np.cos(phi)
    R[:,1,2] = np.sin(phi)*np.sin(theta)
    R[:,2,0] = -np.sin(theta)
    R[:,2,2] = np.cos(theta)
    # calculate pole tide displacements in Cartesian coordinates
    # coefficients reordered to N, E, R to match IERS rotation matrix
    UXYZ = np.einsum('ti...,tji...->tj...',
        np.c_[iu['u_north'], iu['u_east'], iu['u_radial']], R
    )
    # use prediction function to calculate ocean pole tide displacements
    dxi = pyTMD.predict.ocean_pole_tide(ts.tide, np.c_[X, Y, Z], UXYZ,
        gamma_0=ge,
        a_axis=units.a_axis,
        GM=units.GM,
        omega=units.omega,
        rho_w=rho_w,
        g2=gamma,
        convention='2003')
    # calculate components of ocean pole tides
    dui = np.einsum('ti...,tji...->tj...', dxi, np.linalg.inv(R))
    # verify that the predicted values are close
    eps = np.finfo(np.float16).eps
    for i, key in enumerate(['u_north','u_east','u_radial']):
        # assert that the predicted values are close
        assert np.all(np.abs(dui[:,i] - validation[key]) < eps)

# PURPOSE: compute radial load pole tide displacements
# following IERS Convention (2010) guidelines
@pytest.mark.parametrize("TYPE", ['grid','drift','time series'])
@pytest.mark.parametrize("METHOD", ['spline','nearest','linear'])
def test_ocean_pole_tide_displacements(TYPE, METHOD):

    # create a test dataset for data type
    if (TYPE == 'drift'):
        # number of data points
        n_time = 3000
        y = np.random.randint(-90,90,size=n_time)
        x = np.random.randint(-180,180,size=n_time)
        delta_time = np.random.randint(0,31557600,size=n_time)
    elif (TYPE == 'grid'):
        # number of data points
        n_lat,n_lon,n_time = (181,361,100)
        y = np.linspace(-90,90,n_lat)
        x = np.linspace(-180,180,n_lon)
        delta_time = np.random.randint(0,31557600,size=n_time)
    elif (TYPE == 'time series'):
        n_station,n_time = (300,100)
        y = np.random.randint(-90,90,size=n_station)
        x = np.random.randint(-180,180,size=n_station)
        delta_time = np.random.randint(0,31557600,size=n_time)

    # parameters for ocean pole tide
    EPSG = 4326
    EPOCH = (2018, 1, 1, 0, 0, 0)
    TIME = 'UTC'
    ELLIPSOID = 'WGS84'
    CONVENTION = '2018'
    FILL_VALUE = np.nan

    # compute load pole tides using cartesian functions
    test = pyTMD.compute.OPT_displacements(x, y, delta_time,
        EPSG=EPSG, EPOCH=EPOCH, TYPE=TYPE, TIME=TIME,
        ELLIPSOID=ELLIPSOID, CONVENTION=CONVENTION,
        METHOD=METHOD, FILL_VALUE=FILL_VALUE)

    # reform coordinate dimensions for input grids
    # or verify coordinate dimension shapes
    if (TYPE.lower() == 'grid') and (np.size(x) != np.size(y)):
        x,y = np.meshgrid(np.copy(x),np.copy(y))
    elif (TYPE.lower() == 'grid'):
        x = np.atleast_2d(x)
        y = np.atleast_2d(y)
    elif TYPE.lower() in ('time series', 'drift'):
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)

    # converting x,y from EPSG to latitude/longitude
    crs1 = pyTMD.crs().from_input(EPSG)
    crs2 = pyproj.CRS.from_epsg(4326)
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon,lat = transformer.transform(x.flatten(), y.flatten())

    # assert delta time is an array
    delta_time = np.atleast_1d(delta_time)
    # convert delta times or datetimes objects to timescale
    ts = timescale.time.Timescale().from_deltatime(delta_time,
        epoch=EPOCH, standard=TIME)

    # convert dynamic time to Modified Julian Days (MJD)
    MJD = ts.tt - _jd_mjd
    # convert Julian days to calendar dates
    Y,M,D,h,m,s = timescale.time.convert_julian(ts.tt, format='tuple')
    # calculate time in year-decimal format
    time_decimal = timescale.time.convert_calendar_decimal(Y, M, day=D,
        hour=h, minute=m, second=s)
    # number of time points
    nt = len(time_decimal)

    # degrees and arcseconds to radians
    dtr = np.pi/180.0
    atr = np.pi/648000.0
    # earth and physical parameters for ellipsoid
    units = pyTMD.datum(ellipsoid=ELLIPSOID, units='MKS')
    # mean equatorial gravitational acceleration [m/s^2]
    ge = 9.7803278
    # density of sea water [kg/m^3]
    rho_w = 1025.0
    # tidal love number differential (1 + kl - hl) for pole tide frequencies
    gamma = 0.6870 + 0.0036j

    # convert from geodetic latitude to geocentric latitude
    # calculate X, Y and Z from geodetic latitude and longitude
    X,Y,Z = pyTMD.spatial.to_cartesian(lon.flatten(), lat.flatten(),
        a_axis=units.a_axis, flat=units.flat)
    rr = np.sqrt(X**2.0 + Y**2.0 + Z**2.0)
    # calculate geocentric latitude and convert to degrees
    latitude_geocentric = np.arctan(Z / np.sqrt(X**2.0 + Y**2.0))/dtr

    # pole tide displacement scale factor
    Hp = np.sqrt(8.0*np.pi/15.0)*(units.omega**2*units.a_axis**4)/units.GM
    K = 4.0*np.pi*units.G*rho_w*Hp*units.a_axis/(3.0*ge)
    K1 = 4.0*np.pi*units.G*rho_w*Hp*units.a_axis**3/(3.0*units.GM)

    # calculate angular coordinates of mean/secular pole at time
    mpx, mpy, fl = timescale.eop.iers_mean_pole(time_decimal, convention=CONVENTION)
    # read and interpolate IERS daily polar motion values
    px, py = timescale.eop.iers_polar_motion(MJD, k=3, s=0)
    # calculate differentials from mean/secular pole positions
    mx = px - mpx
    my = -(py - mpy)

    # read ocean pole tide map from Desai (2002)
    iur, iun, iue, ilon, ilat = pyTMD.io.IERS.read_binary_file()
    # interpolate ocean pole tide map from Desai (2002)
    if (METHOD == 'spline'):
        # use scipy bivariate splines to interpolate to output points
        f1 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
            iur[:,::-1].real, kx=1, ky=1)
        f2 = scipy.interpolate.RectBivariateSpline(ilon, ilat[::-1],
            iur[:,::-1].imag, kx=1, ky=1)
        UR = np.zeros((len(latitude_geocentric)), dtype=np.clongdouble)
        UR.real = f1.ev(lon.flatten(), latitude_geocentric)
        UR.imag = f2.ev(lon.flatten(), latitude_geocentric)
    else:
        # use scipy regular grid to interpolate values for a given method
        r1 = scipy.interpolate.RegularGridInterpolator((ilon,ilat[::-1]),
            iur[:,::-1], bounds_error=False, method=METHOD)
        UR = r1.__call__(np.c_[lon.flatten(), latitude_geocentric])

    # calculate radial displacement at time
    if (TYPE == 'grid'):
        ny,nx = np.shape(x)
        Urad = np.ma.zeros((ny,nx,nt),fill_value=FILL_VALUE)
        Urad.mask = np.zeros((ny,nx,nt),dtype=bool)
        for i in range(nt):
            URAD = K*atr*np.real((mx[i]*gamma.real + my[i]*gamma.imag)*UR.real +
                (my[i]*gamma.real - mx[i]*gamma.imag)*UR.imag)
            # reform grid
            Urad.data[:,:,i] = np.reshape(URAD, (ny,nx))
            Urad.mask[:,:,i] = np.isnan(Urad.data[:,:,i])
    elif (TYPE == 'drift'):
        Urad = np.ma.zeros((nt),fill_value=FILL_VALUE)
        Urad.data[:] = K*atr*np.real((mx*gamma.real + my*gamma.imag)*UR.real +
            (my*gamma.real - mx*gamma.imag)*UR.imag)
        Urad.mask = np.isnan(Urad.data)
    elif (TYPE == 'time series'):
        nstation = len(x)
        Urad = np.ma.zeros((nstation,nt),fill_value=FILL_VALUE)
        Urad.mask = np.zeros((nstation,nt),dtype=bool)
        for s in range(nstation):
            URAD = K*atr*np.real((mx*gamma.real + my*gamma.imag)*UR.real[s] +
                (my*gamma.real - mx*gamma.imag)*UR.imag[s])
            Urad.data[s,:] = np.copy(URAD)
            Urad.mask[s,:] = np.isnan(Urad.data[s,:])
    # replace invalid data with fill values
    Urad.data[Urad.mask] = Urad.fill_value
    # compare with functional values
    eps = np.finfo(np.float16).eps
    assert np.all(np.abs(Urad - test) < eps)

# PURPOSE: verify inverse of rotation matrix
def test_rotation_matrix():
    # number of data points
    npts = 3000
    lat = np.random.randint(-90,90,size=npts)
    lon = np.random.randint(-180,180,size=npts)
    # colatitude and longitude in radians
    dtr = np.pi/180.0
    theta = dtr*(90.0 - lat)
    phi = dtr*lon
    # convert pole tide values to cartesian coordinates
    R = np.zeros((npts, 3, 3))
    R[:,0,0] = np.cos(phi)*np.cos(theta)
    R[:,0,1] = -np.sin(phi)
    R[:,0,2] = np.cos(phi)*np.sin(theta)
    R[:,1,0] = np.sin(phi)*np.cos(theta)
    R[:,1,1] = np.cos(phi)
    R[:,1,2] = np.sin(phi)*np.sin(theta)
    R[:,2,0] = -np.sin(theta)
    R[:,2,2] = np.cos(theta)
    # rotation matrix for converting from cartesian coordinates
    Rinv = np.zeros((npts, 3, 3))
    Rinv[:,0,0] = np.cos(phi)*np.cos(theta)
    Rinv[:,1,0] = -np.sin(phi)
    Rinv[:,2,0] = np.cos(phi)*np.sin(theta)
    Rinv[:,0,1] = np.sin(phi)*np.cos(theta)
    Rinv[:,1,1] = np.cos(phi)
    Rinv[:,2,1] = np.sin(phi)*np.sin(theta)
    Rinv[:,0,2] = -np.sin(theta)
    Rinv[:,2,2] = np.cos(theta)
    # verify that the rotation matrix is the inverse of the original
    assert np.isclose(np.linalg.inv(R), Rinv).all()
