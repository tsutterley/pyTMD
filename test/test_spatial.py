#!/usr/bin/env python
u"""
test_spatial.py (11/2020)
Verify file read and write with spatial utilities
"""
import boto3
import pytest
import shutil
import inspect
import pathlib
import posixpath
import numpy as np
import pyTMD.spatial
import pyTMD.utilities

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent

# PURPOSE: test the data type function
def test_data_type():
    # test drift type
    exp = 'drift'
    # number of data points
    npts = 30; ntime = 30
    x = np.random.rand(npts)
    y = np.random.rand(npts)
    t = np.random.rand(ntime)
    obs = pyTMD.spatial.data_type(x,y,t)
    assert (obs == exp)
    # test grid type
    exp = 'grid'
    xgrid,ygrid = np.meshgrid(x,y)
    obs = pyTMD.spatial.data_type(xgrid,ygrid,t)
    assert (obs == exp)
    # test grid type with spatial dimensions
    exp = 'grid'
    nx = 30; ny = 20; ntime = 10
    x = np.random.rand(nx)
    y = np.random.rand(ny)
    t = np.random.rand(ntime)
    xgrid,ygrid = np.meshgrid(x,y)
    obs = pyTMD.spatial.data_type(xgrid,ygrid,t)
    assert (obs == exp)
    # test time series type
    exp = 'time series'
    # number of data points
    npts = 1; ntime = 1
    x = np.random.rand(npts)
    y = np.random.rand(npts)
    t = np.random.rand(ntime)
    obs = pyTMD.spatial.data_type(x,y,t)
    assert (obs == exp)
    # test catch for unknown data type
    ermsg = 'Unknown data type'
    # number of data points
    npts = 30; ntime = 10
    x = np.random.rand(npts)
    y = np.random.rand(npts)
    t = np.random.rand(ntime)
    with pytest.raises(ValueError, match=ermsg):
        pyTMD.spatial.data_type(x,y,t)

# PURPOSE: test the read and write of ascii files
def test_ascii():
    # number of data points
    n_time = 30
    # create a test dataset
    output = {}
    output['y'] = np.random.randint(-90,90,size=n_time)
    output['x'] = np.random.randint(-180,180,size=n_time)
    output['data'] = np.random.randn(n_time)
    output['time'] = np.random.randint(0,31557600,size=n_time)

    # output netCDF4 and HDF5 file attributes
    # will be added to YAML header in csv files
    attrib = {}
    # latitude
    attrib['y'] = {}
    attrib['y']['long_name'] = 'Latitude'
    attrib['y']['units'] = 'Degrees_North'
    # longitude
    attrib['x'] = {}
    attrib['x']['long_name'] = 'Longitude'
    attrib['x']['units'] = 'Degrees_East'
    # long-period equilibrium tides
    attrib['data'] = {}
    attrib['data']['long_name'] = 'Height_above_WGS84_ellipsoid'
    attrib['data']['units'] = 'meters'
    # time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'seconds since 2018-01-01T00:00:00'
    attrib['time']['calendar'] = 'standard'

    # create test ascii file
    output_file = filepath.joinpath('test.csv')
    pyTMD.spatial.to_ascii(output, attrib, output_file, delimiter=',',
        columns=['time','y','x','data'], header=True, verbose=True)
    # read test ascii file (change case to test find function)
    input_file = filepath.joinpath('TEST.csv')
    test = pyTMD.spatial.from_ascii(input_file, header='YAML',
        columns=['time','y','x','data'], verbose=True)
    # check that data is valid
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    # read test ascii file as bytes
    fid = output_file.open(mode='r', encoding='utf8')
    test = pyTMD.spatial.from_ascii(fid, compression='bytes', header='YAML',
        columns=['time','y','x','data'])
    # check that data is valid
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    # remove the test file
    output_file.unlink()

# PURPOSE: test the read and write of netCDF4 files
@pytest.mark.parametrize("TYPE", ['drift','grid','time series'])
def test_netCDF4(TYPE):
    # create a test dataset for data type
    output = {}
    if (TYPE == 'drift'):
        # number of data points
        n_time = 3000
        output['y'] = np.random.randint(-90,90,size=n_time)
        output['x'] = np.random.randint(-180,180,size=n_time)
        output['data'] = np.random.randn(n_time)
        output['time'] = np.random.randint(0,31557600,size=n_time)
    elif (TYPE == 'grid'):
        # number of data points
        n_lat,n_lon,n_time = (181,361,100)
        output['y'] = np.linspace(-90,90,n_lat)
        output['x'] = np.linspace(-180,180,n_lon)
        output['data'] = np.random.randn(n_lat,n_lon,n_time)
        output['time'] = np.random.randint(0,31557600,size=n_time)
    elif (TYPE == 'time series'):
        n_station,n_time = (300,100)
        output['y'] = np.random.randint(-90,90,size=n_station)
        output['x'] = np.random.randint(-180,180,size=n_station)
        output['data'] = np.random.randn(n_station,n_time)
        output['time'] = np.random.randint(0,31557600,size=n_time)

    # output netCDF4 and HDF5 file attributes
    # will be added to YAML header in csv files
    attrib = {}
    # latitude
    attrib['y'] = {}
    attrib['y']['long_name'] = 'Latitude'
    attrib['y']['units'] = 'Degrees_North'
    # longitude
    attrib['x'] = {}
    attrib['x']['long_name'] = 'Longitude'
    attrib['x']['units'] = 'Degrees_East'
    # long-period equilibrium tides
    attrib['data'] = {}
    attrib['data']['long_name'] = 'Height_above_WGS84_ellipsoid'
    attrib['data']['units'] = 'meters'
    # time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'seconds since 2018-01-01T00:00:00'
    attrib['time']['calendar'] = 'standard'

    # create test netCDF4 file
    output_file = filepath.joinpath('test.nc')
    pyTMD.spatial.to_netCDF4(output, attrib, output_file,
        data_type=TYPE, verbose=True)
    # read test netCDF4 file (change case to find test function)
    input_file = filepath.joinpath('TEST.nc')
    test = pyTMD.spatial.from_netCDF4(input_file, timename='time',
        xname='x', yname='y', varname='data', verbose=True)
    # check that data is valid
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    # read test netCDF4 file as bytes
    fid = output_file.open(mode='rb')
    test = pyTMD.spatial.from_netCDF4(fid, compression='bytes',
        timename='time', xname='x', yname='y', varname='data')
    # check that data is valid
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    # remove the test file
    output_file.unlink()

# PURPOSE: test the read and write of HDF5 files
@pytest.mark.parametrize("TYPE", ['drift','grid','time series'])
def test_HDF5(TYPE):
    # create a test dataset for data type
    output = {}
    if (TYPE == 'drift'):
        # number of data points
        n_time = 3000
        output['y'] = np.random.randint(-90,90,size=n_time)
        output['x'] = np.random.randint(-180,180,size=n_time)
        output['data'] = np.random.randn(n_time)
        output['time'] = np.random.randint(0,31557600,size=n_time)
    elif (TYPE == 'grid'):
        # number of data points
        n_lat,n_lon,n_time = (181,361,100)
        output['y'] = np.linspace(-90,90,n_lat)
        output['x'] = np.linspace(-180,180,n_lon)
        output['data'] = np.random.randn(n_lat,n_lon,n_time)
        output['time'] = np.random.randint(0,31557600,size=n_time)
    elif (TYPE == 'time series'):
        n_station,n_time = (300,100)
        output['y'] = np.random.randint(-90,90,size=n_station)
        output['x'] = np.random.randint(-180,180,size=n_station)
        output['data'] = np.random.randn(n_station,n_time)
        output['time'] = np.random.randint(0,31557600,size=n_time)

    # output netCDF4 and HDF5 file attributes
    # will be added to YAML header in csv files
    attrib = {}
    # latitude
    attrib['y'] = {}
    attrib['y']['long_name'] = 'Latitude'
    attrib['y']['units'] = 'Degrees_North'
    # longitude
    attrib['x'] = {}
    attrib['x']['long_name'] = 'Longitude'
    attrib['x']['units'] = 'Degrees_East'
    # long-period equilibrium tides
    attrib['data'] = {}
    attrib['data']['long_name'] = 'Height_above_WGS84_ellipsoid'
    attrib['data']['units'] = 'meters'
    # time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'seconds since 2018-01-01T00:00:00'
    attrib['time']['calendar'] = 'standard'

    # create test HDF5 file
    output_file = filepath.joinpath('test.H5')
    pyTMD.spatial.to_HDF5(output, attrib, output_file, verbose=True)
    # read test HDF5 file (change case to test find function)
    input_file = filepath.joinpath('TEST.H5')
    test = pyTMD.spatial.from_HDF5(input_file, timename='time',
        xname='x', yname='y', varname='data', verbose=True)
    # check that data is valid
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    # read test HDF5 file as bytes
    fid = output_file.open(mode='rb')
    test = pyTMD.spatial.from_HDF5(fid, compression='bytes',
        timename='time', xname='x', yname='y', varname='data')
    # check that data is valid
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    # remove the test file
    output_file.unlink()

# PURPOSE: Download IODEM3 from NSIDC
@pytest.fixture(scope="module", autouse=False)
def nsidc_IODEM3(username, password):
    # build urllib2 opener for NSIDC with NASA Earthdata credentials
    pyTMD.utilities.build_opener(username, password,
        password_manager=True, get_ca_certs=False, redirect=False,
        authorization_header=False, urs='https://urs.earthdata.nasa.gov')
    # download NASA Operation IceBridge DMS L3 Photogrammetric DEM
    HOST = ['https://n5eil01u.ecs.nsidc.org','ICEBRIDGE','IODEM3.001',
        '2009.10.25','IODEM3_20091025_212618_02720_DEM.tif']
    granule = filepath.joinpath(HOST[-1])
    # attempt to download file from NSIDC
    remote_buffer = pyTMD.utilities.from_http(HOST, local=granule,
        context=None, verbose=True, mode=0o775)
    # run tests
    yield
    # clean up
    granule.unlink()

# PURPOSE: Download IODEM3 from AWS S3 bucket
@pytest.fixture(scope="module", autouse=True)
def AWS_IODEM3(aws_access_key_id, aws_secret_access_key, aws_region_name):
    # get aws session object
    session = boto3.Session(
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=aws_region_name)
    # get s3 object and bucket object for pytmd data
    s3 = session.resource('s3')
    bucket = s3.Bucket('pytmd')
    # retrieve IODEM3 files
    HOST = ['ICEBRIDGE','IODEM3.001','2009','10','25',
        'IODEM3_20091025_212618_02720_DEM.tif']
    obj = bucket.Object(key=posixpath.join(*HOST))
    response = obj.get()
    granule = filepath.joinpath(HOST[-1])
    with granule.open(mode='wb') as destination:
        shutil.copyfileobj(response['Body'], destination)
    assert granule.exists()
    # run tests
    yield
    # clean up
    granule.unlink()

# PURPOSE: test the read and write of geotiff files
def test_geotiff():
    # read IODEM3 geotiff file
    granule = filepath.joinpath('IODEM3_20091025_212618_02720_DEM.tif')
    dinput = pyTMD.spatial.from_geotiff(granule, verbose=True)
    # copy global geotiff attributes for projection and grid parameters
    attrib = {a:dinput['attributes'][a] for a in ['wkt','spacing','extent']}
    # copy variable attributes for data
    attrib['data'] = {}
    for key,val in dinput['attributes']['data'].items():
        if isinstance(val,np.float32):
            attrib['data'][key] = np.float64(val)
        else:
            attrib['data'][key] = np.copy(val)
    # create test geotiff file
    output_file = filepath.joinpath('test.tif')
    output = {'data':dinput['data'].astype(np.float64)}
    pyTMD.spatial.to_geotiff(output, attrib, output_file,
        driver='GTiff', verbose=True)
    # check that data is valid
    test = pyTMD.spatial.from_geotiff(output_file, verbose=True)
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in dinput.items())
    # check that data is valid from in-memory object
    with granule.open(mode='rb') as fid:
        test = pyTMD.spatial.from_geotiff(fid, compression='bytes')
    eps = np.finfo(np.float32).eps
    assert np.all((np.abs(v-test[k]) < eps) for k,v in dinput.items())
    # remove the test files
    output_file.unlink()

# PURPOSE: test the default field mapping function
def test_field_mapping():
    # test without data variable
    variables = ['delta_time','latitude','longitude']
    field_mapping = pyTMD.spatial.default_field_mapping(variables)
    validation = dict(time='delta_time', y='latitude', x='longitude')
    assert field_mapping == validation
    # test with data variable
    variables = ['time','y','x','height']
    field_mapping = pyTMD.spatial.default_field_mapping(variables)
    validation = dict(time='time', y='y', x='x', data='height')
    assert field_mapping == validation

# PURPOSE: test ellipsoidal parameters within spatial module
def test_ellipsoid_parameters():
    assert pyTMD.spatial._wgs84.a_axis == 6378137.0
    assert pyTMD.spatial._wgs84.flat == (1.0/298.257223563)

# PURPOSE: test ellipsoid conversion
def test_convert_ellipsoid():
    # semimajor axis (a) and flattening (f) for TP and WGS84 ellipsoids
    atop,ftop = (6378136.3,1.0/298.257)
    awgs,fwgs = (6378137.0,1.0/298.257223563)
    # create latitude and height array in WGS84
    lat_WGS84 = 90.0 - np.arange(181,dtype=np.float64)
    elev_WGS84 = 3000.0 + np.zeros((181),dtype=np.float64)
    # convert from WGS84 to Topex/Poseidon Ellipsoids
    lat_TPX,elev_TPX = pyTMD.spatial.convert_ellipsoid(lat_WGS84, elev_WGS84,
        awgs, fwgs, atop, ftop, eps=1e-12, itmax=10)
    # check minimum and maximum are expected from NSIDC IDL program
    # ftp://ftp.nsidc.org/DATASETS/icesat/tools/idl/ellipsoid/test_ce.pro
    minlat = np.min(lat_TPX-lat_WGS84)
    maxlat = np.max(lat_TPX-lat_WGS84)
    explat = [-1.2305653e-7,1.2305653e-7]
    minelev = 100.0*np.min(elev_TPX-elev_WGS84)
    maxelev = 100.0*np.max(elev_TPX-elev_WGS84)
    expelev = [70.000000,71.368200]
    assert np.isclose([minlat,maxlat],explat).all()
    assert np.isclose([minelev,maxelev],expelev).all()
    # convert back from Topex/Poseidon to WGS84 Ellipsoids
    phi_WGS84,h_WGS84 = pyTMD.spatial.convert_ellipsoid(lat_TPX, elev_TPX,
        atop, ftop, awgs, fwgs, eps=1e-12, itmax=10)
    # check minimum and maximum are expected from NSIDC IDL program
    # ftp://ftp.nsidc.org/DATASETS/icesat/tools/idl/ellipsoid/test_ce.pro
    minlatdel = np.min(phi_WGS84-lat_WGS84)
    maxlatdel = np.max(phi_WGS84-lat_WGS84)
    explatdel = [-2.1316282e-14,2.1316282e-14]
    minelevdel = 100.0*np.min(h_WGS84-elev_WGS84)
    maxelevdel = 100.0*np.max(h_WGS84-elev_WGS84)
    expelevdel = [-1.3287718e-7,1.6830199e-7]
    assert np.isclose([minlatdel,maxlatdel],explatdel).all()
    assert np.isclose([minelevdel,maxelevdel],expelevdel,atol=1e-5).all()

# PURPOSE: verify cartesian to geodetic conversions
def test_convert_geodetic():
    # choose a random set of locations
    latitude = -90.0 + 180.0*np.random.rand(100)
    longitude = -180.0 + 360.0*np.random.rand(100)
    height = 2000.0*np.random.rand(100)
    # ellipsoidal parameters
    a_axis = pyTMD.spatial._wgs84.a_axis
    flat = pyTMD.spatial._wgs84.flat
    # convert to cartesian coordinates
    x, y, z = pyTMD.spatial.to_cartesian(longitude, latitude, h=height,
        a_axis=a_axis, flat=flat)
    # convert back to geodetic coordinates
    ln1, lt1, h1 = pyTMD.spatial.to_geodetic(x, y, z,
        a_axis=a_axis, flat=flat, method='moritz')
    ln2, lt2, h2 = pyTMD.spatial.to_geodetic(x, y, z,
        a_axis=a_axis, flat=flat, method='bowring')
    ln3, lt3, h3 = pyTMD.spatial.to_geodetic(x, y, z,
        a_axis=a_axis, flat=flat, method='zhu')
    # validate outputs for Moritz iterative method
    assert np.isclose(longitude, ln1).all()
    assert np.isclose(latitude, lt1).all()
    assert np.isclose(height, h1).all()
     # validate outputs for Bowring iterative method
    assert np.isclose(longitude, ln2).all()
    assert np.isclose(latitude, lt2).all()
    assert np.isclose(height, h2).all()
    # validate outputs for Zhu closed-form method
    assert np.isclose(longitude, ln3).all()
    assert np.isclose(latitude, lt3).all()
    assert np.isclose(height, h3).all()

# PURPOSE: test wrap longitudes
def test_wrap_longitudes():
    # number of data points
    lon = np.arange(360)
    obs = pyTMD.spatial.wrap_longitudes(lon)
    # expected wrapped longitudes
    exp = np.zeros((360))
    exp[:181] = np.arange(181)
    exp[181:] = np.arange(-179,0)
    assert np.isclose(obs,exp).all()
