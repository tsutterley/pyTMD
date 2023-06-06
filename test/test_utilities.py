#!/usr/bin/env python
u"""
test_utilities.py (12/2020)
Verify file utility functions
"""
import io
import gzip
import pytest
import posixpath
import pyTMD.utilities

def test_hash():
    # get hash of compressed file
    ocean_pole_tide_file = pyTMD.utilities.get_data_path(['data',
        'opoleloadcoefcmcor.txt.gz'])
    TEST = pyTMD.utilities.get_hash(ocean_pole_tide_file)
    assert (TEST == '9c66edc2d0fbf627e7ae1cb923a9f0e5')
    # get hash of uncompressed file
    with gzip.open(ocean_pole_tide_file) as fid:
        TEST = pyTMD.utilities.get_hash(io.BytesIO(fid.read()))
        assert (TEST == 'cea08f83d613ed8e1a81f3b3a9453721')

def test_urlsplit():
    HOST = ['https://cddis.nasa.gov','archive','products','iers','deltat.data']
    url = posixpath.join(*HOST)
    TEST = pyTMD.utilities.url_split(url)
    assert (HOST == list(TEST))

def test_roman():
    for s,i in zip(['MMXVIII','MCMXVI','DXCIV'],[2018,1916,594]):
        TEST = pyTMD.utilities.roman_to_int(s)
        assert (TEST == i)

def test_even():
    for s,i in zip([2015,1916,591,99,10,3],[2014,1916,590,98,10,2]):
        TEST = pyTMD.utilities.even(s)
        assert (TEST == i)

def test_ceil():
    for s,i in zip([-2.5, 0.0, 2.5], [-2, 0, 3]):
        TEST = pyTMD.utilities.ceil(s)
        assert (TEST == i)

def test_token(username, password):
    # attempt to login to NASA Earthdata
    urs = 'urs.earthdata.nasa.gov'
    opener = pyTMD.utilities.attempt_login(urs,
        username=username,
        password=password,
        password_manager=False,
        get_ca_certs=False,
        redirect=False,
        authorization_header=True)
    # get data access token
    token = pyTMD.utilities.get_token(build=False)
    # get list of tokens
    token_list = pyTMD.utilities.list_tokens(build=False)
    # assert that token is in list of tokens
    access_tokens = [t['access_token'] for t in token_list]
    assert (token['access_token'] in access_tokens)
    # revoke the access token
    pyTMD.utilities.revoke_token(token['access_token'], build=False)
