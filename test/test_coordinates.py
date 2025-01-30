#!/usr/bin/env python
u"""
test_coordinates.py (09/2024)
Verify forward and backwards coordinate conversions

UPDATE HISTORY:
    Updated 09/2024: add test for Arctic regions with new projection
        using new JSON dictionary format for model projections
    Updated 07/2024: add check for if projections are geographic
    Updated 12/2023: use new crs class for coordinate reprojection
    Written 08/2020
"""
import pytest
import numpy as np
import pyTMD.crs

# parameterize projections
models = ['AOTIM-5-2018','Gr1km-v2','CATS2008','TPXO9-atlas-v5']
@pytest.mark.parametrize("MODEL", models)
# PURPOSE: verify forward and backwards coordinate conversions
def test_coordinates(MODEL):
    startlat = {'Gr1km-v2':60, 'CATS2008':-60,
        'AOTIM-5-2018':60,'TPXO9-atlas-v5':90}
    endlat = {'Gr1km-v2':70, 'CATS2008':-70,
        'AOTIM-5-2018':70,'TPXO9-atlas-v5':-90}
    is_geographic = {'Gr1km-v2':False, 'CATS2008':False,
        'AOTIM-5-2018':False, 'TPXO9-atlas-v5':True}
    i1 = np.arange(-180,180+1,1)
    i2 = np.linspace(startlat[MODEL],endlat[MODEL],len(i1))
    # convert latitude and longitude to and from projection
    model = pyTMD.models.elevation[MODEL]
    crs = pyTMD.crs().get(model['projection'])
    o1, o2 = crs.transform(i1, i2, direction='FORWARD')
    lon, lat = crs.transform(o1, o2, direction='INVERSE')
    # calculate great circle distance between inputs and outputs
    cdist = np.arccos(np.sin(i2*np.pi/180.0)*np.sin(lat*np.pi/180.0) +
        np.cos(i2*np.pi/180.0)*np.cos(lat*np.pi/180.0)*
        np.cos((lon-i1)*np.pi/180.0),dtype=np.float32)
    # test that forward and backwards conversions are within tolerance
    eps = np.finfo(np.float32).eps
    assert np.all(cdist < eps)
    # convert latitude and longitude to and from projection
    # using the convert function
    o1, o2 = pyTMD.crs().convert(i1, i2, model['projection'], 'F')
    lon, lat = pyTMD.crs().convert(o1, o2, model['projection'], 'B')
    # calculate great circle distance between inputs and outputs
    cdist = np.arccos(np.sin(i2*np.pi/180.0)*np.sin(lat*np.pi/180.0) +
        np.cos(i2*np.pi/180.0)*np.cos(lat*np.pi/180.0)*
        np.cos((lon-i1)*np.pi/180.0),dtype=np.float32)
    # test that forward and backwards conversions are within tolerance
    eps = np.finfo(np.float32).eps
    assert np.all(cdist < eps)
    assert crs.is_geographic == is_geographic[MODEL]

# PURPOSE: verify coordinate conversions are close for Arctic regions
def test_arctic_projection():
    # generate random latitude and longitude coordinates
    N = 10000
    i1 = -180.0 + 360.0*np.random.rand(N)
    i2 = 60.0 + 30.0*np.random.rand(N)
    # convert latitude and longitude to and from projection
    model = pyTMD.models.elevation['AOTIM-5-2018']
    crs = pyTMD.crs().get(model['projection'])
    o1, o2 = crs.transform(i1, i2, direction='FORWARD')
    lon, lat = crs.transform(o1, o2, direction='INVERSE')
    # calculate great circle distance between inputs and outputs
    cdist = np.arccos(np.sin(i2*np.pi/180.0)*np.sin(lat*np.pi/180.0) +
        np.cos(i2*np.pi/180.0)*np.cos(lat*np.pi/180.0)*
        np.cos((lon-i1)*np.pi/180.0),dtype=np.float32)
    # test that forward and backwards conversions are within tolerance
    eps = np.finfo(np.float32).eps
    assert np.all(cdist < eps)
    # convert projected coordinates from latitude and longitude
    x = (90.0 - i2)*111.7*np.cos(i1/180.0*np.pi)
    y = (90.0 - i2)*111.7*np.sin(i1/180.0*np.pi)
    assert np.isclose(o1, x).all()
    assert np.isclose(o2, y).all()
    # convert latitude and longitude from projected coordinates
    ln = np.arctan2(y, x)*180.0/np.pi
    lt = 90.0 - np.sqrt(x**2 + y**2)/111.7
    # adjust longitudes to be -180:180
    ii, = np.nonzero(ln < 0)
    ln[ii] += 360.0
    # calculate great circle distance between inputs and outputs
    cdist = np.arccos(np.sin(lat*np.pi/180.0)*np.sin(lt*np.pi/180.0) +
        np.cos(lat*np.pi/180.0)*np.cos(lt*np.pi/180.0)*
        np.cos((lon-ln)*np.pi/180.0),dtype=np.float32)
    # test that forward and backwards conversions are within tolerance
    eps = np.finfo(np.float32).eps
    assert np.all(cdist < eps)
