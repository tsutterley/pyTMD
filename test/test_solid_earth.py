"""
test_solid_earth.py (04/2023)
Tests the steps for calculating the solid earth tides

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Written 04/2023
"""
import pytest
import numpy as np
import pyTMD.astro
import pyTMD.predict

def test_out_of_phase_diurnal():
    """Test out-of-phase diurnal corrections with IERS outputs
    """
    # station locations and planetary ephemerides
    XYZ = np.array([[4075578.385, 931852.890, 4801570.154]])
    SXYZ = np.array([[137859926952.015, 54228127881.4350, 23509422341.6960]])
    LXYZ = np.array([[-179996231.920342, -312468450.131567, -169288918.592160]])
    # factors for sun and moon
    F2_solar = 0.163271964478954
    F2_lunar = 0.321989090026845
    # expected results
    dx_expected = -0.2836337012840008001e-3
    dy_expected = 0.1125342324347507444e-3
    dz_expected = -0.2471186224343683169e-3
    # calculate expected displacements
    dx, dy, dz = pyTMD.predict._out_of_phase_diurnal(XYZ, SXYZ, LXYZ,
        F2_solar, F2_lunar)
    # assert matching
    assert np.isclose(np.c_[dx_expected, dy_expected, dz_expected],
        np.c_[dx, dy, dz]).all()

def test_out_of_phase_semidiurnal():
    """Test out-of-phase semidiurnal corrections with IERS outputs
    """
    # station locations and planetary ephemerides
    XYZ = np.array([[4075578.385, 931852.890, 4801570.154]])
    SXYZ = np.array([[137859926952.015, 54228127881.4350, 23509422341.6960]])
    LXYZ = np.array([[-179996231.920342, -312468450.131567, -169288918.592160]])
    # factors for sun and moon
    F2_solar = 0.163271964478954
    F2_lunar = 0.321989090026845
    # expected results
    dx_expected = -0.2801334805106874015e-3
    dy_expected = 0.2939522229284325029e-4
    dz_expected = -0.6051677912316721561e-4
    # calculate expected displacements
    dx, dy, dz = pyTMD.predict._out_of_phase_semidiurnal(XYZ, SXYZ, LXYZ,
        F2_solar, F2_lunar)
    # assert matching
    assert np.isclose(np.c_[dx_expected, dy_expected, dz_expected],
        np.c_[dx, dy, dz]).all()

def test_latitude_dependence():
    """Test latitude dependence corrections with IERS outputs
    """
    # station locations and planetary ephemerides
    XYZ = np.array([[4075578.385, 931852.890, 4801570.154]])
    SXYZ = np.array([[137859926952.015, 54228127881.4350, 23509422341.6960]])
    LXYZ = np.array([[-179996231.920342, -312468450.131567, -169288918.592160]])
    # factors for sun and moon
    F2_solar = 0.163271964478954
    F2_lunar = 0.321989090026845
    # expected results
    dx_expected = 0.2367189532359759044e-3
    dy_expected = 0.5181609907284959182e-3
    dz_expected = -0.3014881422940427977e-3
    # calculate expected displacements
    dx, dy, dz = pyTMD.predict._latitude_dependence(XYZ, SXYZ, LXYZ,
        F2_solar, F2_lunar)
    # assert matching
    assert np.isclose(np.c_[dx_expected, dy_expected, dz_expected],
        np.c_[dx, dy, dz]).all()

def test_frequency_dependence_diurnal():
    """Test diurnal band frequency dependence corrections
    with IERS outputs
    """
    # station locations
    XYZ = np.array([[4075578.385, 931852.890, 4801570.154]])
    MJD = 55414.0
    # convert from MJD to centuries relative to 2000-01-01T12:00:00
    T = (MJD - 51544.5)/36525.0
    T_expected = 0.1059411362080767
    assert np.isclose(T_expected, T)
    # expected results
    dx_expected = 0.4193085327321284701e-2
    dy_expected = 0.1456681241014607395e-2
    dz_expected = 0.5123366597450316508e-2
    # calculate expected displacements
    dx, dy, dz = pyTMD.predict._frequency_dependence_diurnal(XYZ, MJD)
    # assert matching
    assert np.isclose(np.c_[dx_expected, dy_expected, dz_expected],
        np.c_[dx, dy, dz]).all()

def test_frequency_dependence_long_period():
    """Test long period band frequency dependence corrections
    with IERS outputs
    """
    # station locations
    XYZ = np.array([[4075578.385, 931852.890, 4801570.154]])
    MJD = 55414.0
    # convert from MJD to centuries relative to 2000-01-01T12:00:00
    T = (MJD - 51544.5)/36525.0
    T_expected = 0.1059411362080767
    assert np.isclose(T_expected, T)
    # expected results
    dx_expected = -0.9780962849562107762e-4
    dy_expected = -0.2236349699932734273e-4
    dz_expected = 0.3561945821351565926e-3
    # calculate expected displacements
    dx, dy, dz = pyTMD.predict._frequency_dependence_long_period(XYZ, MJD)
    # assert matching
    assert np.isclose(np.c_[dx_expected, dy_expected, dz_expected],
        np.c_[dx, dy, dz]).all()

def test_phase_angles():
    """Test that longitudes and phase angles match between functions
    """
    MJD = 55414.0
    dtr = np.pi/180.0
    # convert from MJD to centuries relative to 2000-01-01T12:00:00
    T = (MJD - 51544.5)/36525.0
    s, h, p, N, PP = pyTMD.astro.mean_longitudes(MJD, ASTRO5=True)
    PR = dtr* pyTMD.astro.polynomial_sum(np.array([0.0, 1.396971278,
        3.08889e-4, 2.1e-8, 7.0e-9]), T)
    S, H, P, TAU, ZNS, PS = pyTMD.astro.phase_angles(MJD)
    assert np.isclose(dtr*s + PR, S)
    assert np.isclose(dtr*h, H)
    assert np.isclose(dtr*p, P)
    assert np.isclose(2.0*np.pi - N*dtr, ZNS)
