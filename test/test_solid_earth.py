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
import pyTMD.time
from pyTMD.compute_tide_corrections import compute_SET_corrections

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
    # calculate displacements
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
    # calculate displacements
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
    # calculate displacements
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
    # calculate displacements
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
    # calculate displacements
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
    PR = dtr*pyTMD.astro.polynomial_sum(np.array([0.0, 1.396971278,
        3.08889e-4, 2.1e-8, 7.0e-9]), T)
    S, H, P, TAU, ZNS, PS = pyTMD.astro.phase_angles(MJD)
    assert np.isclose(dtr*s + PR, S)
    assert np.isclose(dtr*h, H)
    assert np.isclose(dtr*p, P)
    assert np.isclose(2.0*np.pi - N*dtr, ZNS)

def test_solid_earth_tide():
    """Test solid earth tides with IERS outputs
    """
    # case 1 from IERS
    XYZ = np.array([[4075578.385, 931852.890, 4801570.154]])
    SXYZ = np.array([[137859926952.015, 54228127881.4350, 23509422341.6960]])
    LXYZ = np.array([[-179996231.920342, -312468450.131567, -169288918.592160]])
    tide_time = pyTMD.time.convert_calendar_dates(2009, 4, 13,
        hour=0, minute=0, second=0,
        epoch=pyTMD.time._tide_epoch)
    # expected results
    dx_expected = 0.7700420357108125891e-01
    dy_expected = 0.6304056321824967613e-01
    dz_expected = 0.5516568152597246810e-01
    # calculate solid earth tides
    dxt = pyTMD.predict.solid_earth_tide(tide_time, XYZ, SXYZ, LXYZ)
    # assert matching
    assert np.isclose(np.c_[dx_expected, dy_expected, dz_expected],dxt).all()
    # case 2 from IERS
    XYZ = np.array([[1112200.5696, -4842957.8511, 3985345.9122]])
    SXYZ = np.array([[100210282451.6279, 103055630398.316, 56855096480.4475]])
    LXYZ = np.array([[369817604.4348, 1897917.5258, 120804980.8284]])
    tide_time = pyTMD.time.convert_calendar_dates(2015, 7, 15,
        hour=0, minute=0, second=0,
        epoch=pyTMD.time._tide_epoch)
    # expected results
    dx_expected = 0.00509570869172363845
    dy_expected = 0.0828663025983528700
    dz_expected = -0.0636634925404189617
    # calculate solid earth tides
    dxt = pyTMD.predict.solid_earth_tide(tide_time, XYZ, SXYZ, LXYZ)
    # assert matching
    assert np.isclose(np.c_[dx_expected, dy_expected, dz_expected],dxt).all()

def test_solid_earth_radial():
    """Test radial solid tides with predictions from ICESat-2
    """
    times = np.array(['2018-10-14 00:21:48','2018-10-14 00:21:48',
        '2018-10-14 00:21:48','2018-10-14 00:21:48',
        '2022-07-23 13:53:08','2022-07-23 13:53:08',
        '2022-07-23 13:53:08','2022-07-23 13:53:08'], dtype=np.datetime64)
    longitudes = np.array([-136.79534534,-136.79545175,
        -136.79548250,-136.79549453,-71.77356870,-71.77374742,
        -71.77392700,-71.77410705])
    latitudes = np.array([68.95910366,68.95941755,
        68.95950895,68.95954490,-79.00591611,-79.00609103,
        -79.00626593,-79.00644081])
    # expected results (tide-free)
    tide_earth = np.array([-0.14320290,-0.14320324,
        -0.14320339,-0.14320345,-0.11887791,-0.11887763,
        -0.11887724,-0.11887683])
    # tide_mean = tide_free + tide_earth_free2mean
    tide_earth_free2mean = np.array([-0.09726650,-0.09726728,
        -0.09726749,-0.09726755,-0.11400376,-0.11400391,
        -0.11400412,-0.11400434])
    # predict radial solid earth tides
    tide_free = compute_SET_corrections(longitudes, latitudes, times,
        EPSG=4326, TYPE='drift', TIME='datetime', ELLIPSOID='WGS84')
    tide_mean = compute_SET_corrections(longitudes, latitudes, times,
        EPSG=4326, TYPE='drift', TIME='datetime', ELLIPSOID='WGS84',
        TIDE_SYSTEM='mean_tide')
    # as using estimated ephemerides, assert within 1 mm
    assert np.isclose(tide_earth, tide_free, atol=1e-3).all()
    # check permanent tide offsets (additive correction in ICESat-2)
    # expected results (mean-tide)
    tide_expected = tide_earth + tide_earth_free2mean
    # # sign differences with ATLAS product?
    # assert np.isclose(tide_expected, tide_mean, atol=1e-3).all()
