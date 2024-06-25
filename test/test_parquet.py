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

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent

# PURPOSE: test the read and write of (geo)parquet files
def test_parquet():
    # number of data points
    n_time = 30000
    # create a test dataset
    output = {}
    rng = np.random.default_rng()
    output['y'] = rng.uniform(-90.0, high=90.0, size=n_time)
    output['x'] = rng.uniform(-180.0, high=180.0, size=n_time)
    output['data'] = rng.standard_normal(n_time)
    output['time'] = rng.uniform(0.0, high=31557600.0, size=n_time)
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
    # long-period equilibrium tides
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
    pyTMD.spatial.to_parquet(output, attrib, output_file, crs=4326)
    # read test parquet file
    test = pyTMD.spatial.from_file(output_file, format='parquet',
        columns=['time','y','x','data'])
    # check that data is valid
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    # test that parquet file can be read by pandas
    df = geopandas.pd.read_parquet(output_file)
    # check that data is valid
    assert np.all((np.abs(v-df[k].values) < eps) for k,v in output.items())
    # remove the test file
    output_file.unlink()

    # create test geoparquet file
    pyTMD.spatial.to_parquet(output, attrib, output_file,
        geoparquet=True, geometry_encoding='WKB', crs=4326)
    # read test geoparquet file
    test = pyTMD.spatial.from_file(output_file, format='parquet',
        columns=['time','y','x','data'])
    # check that data is valid
    assert np.all((np.abs(v-test[k]) < eps) for k,v in output.items())
    # test that geoparquet file can be read by geopandas
    gdf = geopandas.read_parquet(output_file)
    # check that data is valid
    assert np.all((np.abs(v-gdf[k].values) < eps) for k,v in output.items())
    # remove the test file
    output_file.unlink()
