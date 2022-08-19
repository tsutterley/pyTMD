#!/usr/bin/env python
u"""
test_coordinates.py (08/2020)
Verify forward and backwards coordinate conversions
"""
import warnings
import pytest
import numpy as np
import pyTMD.convert_ll_xy

# parameterize projections
@pytest.mark.parametrize("PROJ", ['3031','CATS2008','3976','PSNorth','4326'])
# PURPOSE: verify forward and backwards coordinate conversions
def test_coordinates(PROJ):
    startlat = {'3031':-60,'CATS2008':-60,'3976':-60,'PSNorth':60,'4326':90}
    endlat = {'3031':-70,'CATS2008':-70,'3976':-70,'PSNorth':70,'4326':-90}
    i1 = np.arange(-180,180+1,1)
    i2 = np.linspace(startlat[PROJ],endlat[PROJ],len(i1))
    # convert latitude and longitude to and from projection
    o1,o2 = pyTMD.convert_ll_xy(i1,i2,PROJ,'F')
    lon,lat = pyTMD.convert_ll_xy(o1,o2,PROJ,'B')
    # calculate great circle distance between inputs and outputs
    cdist = np.arccos(np.sin(i2*np.pi/180.0)*np.sin(lat*np.pi/180.0) +
        np.cos(i2*np.pi/180.0)*np.cos(lat*np.pi/180.0)*
        np.cos((lon-i1)*np.pi/180.0),dtype=np.float32)
    # test that forward and backwards conversions are within tolerance
    eps = np.finfo(np.float32).eps
    assert np.all(cdist < eps)
