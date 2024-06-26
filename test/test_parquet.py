#!/usr/bin/env python
u"""
test_parquet.py (06/2024)
Verify (geo)parquet file read and write with spatial utilities
"""
import inspect
import pathlib
import numpy as np
import pyTMD.spatial
import pyTMD.utilities
# attempt imports
geopandas = pyTMD.utilities.import_dependency('geopandas')
pyproj = pyTMD.utilities.import_dependency('pyproj')

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent

# PURPOSE: test the read and write of parquet files
def test_parquet():
    # number of data points
    n_time = 30000
    # create a test dataset
    output = {}
    # random number generator
    rng = np.random.default_rng()
    # randomly generated data and times
    output['data'] = 2.0*rng.standard_normal(n_time)
    output['time'] = rng.uniform(0.0, high=31557600.0, size=n_time)
    # use a global range of coordinates
    output['y'] = rng.uniform(-90.0, high=90.0, size=n_time)
    output['x'] = rng.uniform(-180.0, high=180.0, size=n_time)
    # coordinate reference system
    crs = 4326
    crs2 = pyTMD.crs().from_input(crs)
    # validation tolerance
    eps = np.finfo(np.float64).eps

    # output file attributes
    attrib = {}
    # latitude
    attrib['y'] = {}
    attrib['y']['long_name'] = 'Latitude'
    attrib['y']['units'] = 'Degrees_North'
    # longitude
    attrib['x'] = {}
    attrib['x']['long_name'] = 'Longitude'
    attrib['x']['units'] = 'Degrees_East'
    # data
    attrib['data'] = {}
    attrib['data']['long_name'] = 'Height_above_WGS84_ellipsoid'
    attrib['data']['units'] = 'meters'
    # time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'seconds since 2018-01-01T00:00:00'
    attrib['time']['calendar'] = 'standard'

    # create test parquet file
    output_file = filepath.joinpath('test.parquet')
    pyTMD.spatial.to_parquet(output, attrib, output_file, crs=crs)
    # read test parquet file
    test = pyTMD.spatial.from_file(output_file, format='parquet',
        columns=['time','y','x','data'])
    # check that the crs is valid (retrieved from pyTMD metadata)
    crs1 = pyTMD.crs().from_input(test.attrs['wkt'])
    assert crs1.equals(crs2)
    # check that data is valid
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    # test that parquet file can be read by pandas
    df = geopandas.pd.read_parquet(output_file)
    # check that data is valid
    assert np.all((np.abs(v-df[k].values) < eps) for k,v in output.items())
    # remove the test file
    output_file.unlink()

# PURPOSE: test the read and write of geoparquet files
def test_geoparquet():
    # number of data points
    n_time = 30000
    # create a test dataset
    output = {}
    # random number generator
    rng = np.random.default_rng()
    # randomly generated data and times
    output['data'] = 1000.0*(3.0 - rng.standard_normal(n_time))
    output['time'] = rng.uniform(0.0, high=31557600.0, size=n_time)
    # use range of Bamber 1km Antarctic DEM
    output['y'] = rng.uniform(-560.*5e3, high=560.*5e3, size=n_time)
    output['x'] = rng.uniform(-560.*5e3, high=560.*5e3, size=n_time)
    # coordinate reference system
    crs = 3031
    crs2 = pyTMD.crs().from_input(crs)
    # dictionary of coordinate reference system variables
    cs_to_cf = crs2.cs_to_cf()
    # validation tolerance
    eps = np.finfo(np.float64).eps

    # output file attributes
    attrib = {}
    # x and y
    attrib['x'],attrib['y'] = ({},{})
    for att_name in ['long_name','standard_name','units']:
        attrib['x'][att_name] = cs_to_cf[0][att_name]
        attrib['y'][att_name] = cs_to_cf[1][att_name]
    # data
    attrib['data'] = {}
    attrib['data']['long_name'] = 'Height_above_WGS84_ellipsoid'
    attrib['data']['units'] = 'meters'
    # time
    attrib['time'] = {}
    attrib['time']['long_name'] = 'Time'
    attrib['time']['units'] = 'seconds since 2018-01-01T00:00:00'
    attrib['time']['calendar'] = 'standard'

    # create test geoparquet file
    output_file = filepath.joinpath('test.parquet')
    pyTMD.spatial.to_parquet(output, attrib, output_file,
        geoparquet=True, geometry_encoding='WKB', crs=crs)
    # read test geoparquet file
    test = pyTMD.spatial.from_file(output_file, format='parquet')
    # check that data is valid
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    assert 'geometry' in test.columns
    # check that the crs is valid (retrieved from geo metadata)
    crs1 = pyTMD.crs().from_input(test.attrs['wkt'])
    assert crs1.equals(crs2)
    # test that geoparquet file can be read by geopandas
    gdf = geopandas.read_parquet(output_file)
    # check that the crs is valid
    crs1 = pyTMD.crs().from_input(gdf.crs.to_wkt())
    assert crs1.equals(crs2)
    # check that data is valid
    assert np.all((np.abs(v-gdf[k].values) < eps) for k,v in output.items())
    # remove the test file
    output_file.unlink()
