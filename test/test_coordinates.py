#!/usr/bin/env python
u"""
test_coordinates.py (09/2024)
Verify forward and backwards coordinate conversions

UPDATE HISTORY:
    Updated 09/2024: add test for Arctic regions with new projection
    Updated 07/2024: add check for if projections are geographic
    Updated 12/2023: use new crs class for coordinate reprojection
    Written 08/2020
"""
import pytest
import numpy as np
import pyTMD.crs

# parameterize projections
@pytest.mark.parametrize("PROJ", ['3031','CATS2008','3976','AEDNorth','4326'])
# PURPOSE: verify forward and backwards coordinate conversions
def test_coordinates(PROJ):
    startlat = {'3031':-60,'CATS2008':-60,'3976':-60,'AEDNorth':60,'4326':90}
    endlat = {'3031':-70,'CATS2008':-70,'3976':-70,'AEDNorth':70,'4326':-90}
    is_geographic = {'3031':False,'CATS2008':False,'3976':False,
        'AEDNorth':False,'4326':True}
    i1 = np.arange(-180,180+1,1)
    i2 = np.linspace(startlat[PROJ],endlat[PROJ],len(i1))
    # convert latitude and longitude to and from projection
    transform = pyTMD.crs().get(PROJ)
    o1, o2 = pyTMD.crs().convert(i1,i2,PROJ,'F')
    lon, lat = pyTMD.crs().convert(o1,o2,PROJ,'B')
    # calculate great circle distance between inputs and outputs
    cdist = np.arccos(np.sin(i2*np.pi/180.0)*np.sin(lat*np.pi/180.0) +
        np.cos(i2*np.pi/180.0)*np.cos(lat*np.pi/180.0)*
        np.cos((lon-i1)*np.pi/180.0),dtype=np.float32)
    # test that forward and backwards conversions are within tolerance
    eps = np.finfo(np.float32).eps
    assert np.all(cdist < eps)
    assert transform.is_geographic == is_geographic[PROJ]

# PURPOSE: verify coordinate conversions are close for Arctic regions
def test_arctic_projection():
    # generate random latitude and longitude coordinates
    N = 10000
    i1 = -180.0 + 360.0*np.random.rand(N)
    i2 = 60.0 + 30.0*np.random.rand(N)
    # convert latitude and longitude to and from projection
    o1, o2 = pyTMD.crs().convert(i1,i2,'AEDNorth','F')
    lon, lat = pyTMD.crs().convert(o1,o2,'AEDNorth','B')
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
