#!/usr/bin/env python
u"""
test_math.py (11/2024)
"""
import pytest
import numpy as np
import pyTMD.math

def test_normalize_angle():
    """
    Tests the normalization of angles to between 0 and 360 degrees
    """
    # test angles
    angles = np.array([-180, -90, 0, 90, 180, 270, 360, 450])
    # expected values
    exp = np.array([180, 270, 0, 90, 180, 270, 0, 90])
    # test normalization of angles
    test = pyTMD.math.normalize_angle(angles)
    assert np.all(exp == test)

def test_aliasing():
    """
    Tests the calculation of an aliasing frequency
    """
    # test frequencies
    frequency = 39000
    sampling = np.array([37500, 38000, 38500, 44000, 44500, 45000])
    # expected values
    exp = np.array([1500, 1000, 500, 5000, 5500, 6000])
    # test aliasing frequencies
    test = pyTMD.math.aliasing(frequency, sampling)
    assert np.all(exp == test)

@pytest.mark.parametrize("l", [1, 2, 3])
def test_legendre(l, x=[-1.0, -0.9, -0.8]):
    """test the calculation of unnormalized Legendre polynomials
    """
    # calculate Legendre polynomials
    nx = len(x)
    obs = np.zeros((l+1, nx))
    for m in range(l+1):
        obs[m,:] = pyTMD.math.legendre(l, x, m=m)
    # expected values for each spherical harmonic degree
    if (l == 1):
        expected = np.array([
            [-1.00000, -0.90000, -0.80000],
            [ 0.00000, -0.43589, -0.60000]
        ])
    elif (l == 2):
        expected = np.array([
            [1.00000, 0.71500, 0.46000],
            [0.00000, 1.17690, 1.44000],
            [0.00000, 0.57000, 1.08000]
        ])
    elif (l == 3):
        expected = np.array([
            [-1.00000, -0.47250, -0.08000],
            [0.00000, -1.99420, -1.98000],
            [0.00000, -2.56500, -4.32000],
            [0.00000, -1.24229, -3.24000]
        ])
    # check with expected values
    assert np.isclose(obs, expected, atol=1e-05).all()

# PURPOSE: test fully-normalized Legendre polynomials
def test_normalized(lmax=3, x=[-1.0, -0.9, -0.8]):
    """test the calculation of normalized Legendre polynomials
    """
    obs = pyTMD.math.assoc_legendre(lmax, x)
    expected = np.array([
        [-2.64575, -1.25012, -0.21166],
        [-0.00000, 2.15398, 2.13864],
        [0.00000, -0.87611, -1.47556],
        [-0.00000, 0.17323, 0.45180]
    ])
    assert np.isclose(obs[lmax:,:], expected, atol=1e-05).all()

# PURPOSE: test fully-normalized zonal Legendre polynomials
def test_zonal(lmax=3, x=[-1.0, -0.9, -0.8]):
    """test the calculation of zonal Legendre polynomials
    """
    obs = pyTMD.math.assoc_legendre(lmax, x)
    expected = np.array([
        [1.00000, 1.00000, 1.00000],
        [-1.73205, -1.55885, -1.38564],
        [2.23607, 1.59879, 1.02859],
        [-2.64575, -1.25012, -0.21166],
    ])
    assert np.isclose(obs[:lmax+1,0,:], expected, atol=1e-05).all()
