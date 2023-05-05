"""
test_solid_earth.py (04/2023)
Tests the steps for calculating the solid earth tides

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 04/2023: added test for using JPL ephemerides for positions
    Written 04/2023
"""
import pytest
import numpy as np
import pyTMD.astro
import pyTMD.predict
import pyTMD.time
import pyTMD.utilities
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
    dXYZ = pyTMD.predict._out_of_phase_diurnal(XYZ, SXYZ, LXYZ,
        F2_solar, F2_lunar)
    # assert matching
    assert np.isclose(np.c_[dx_expected, dy_expected, dz_expected], dXYZ).all()

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
    dXYZ = pyTMD.predict._out_of_phase_semidiurnal(XYZ, SXYZ, LXYZ,
        F2_solar, F2_lunar)
    # assert matching
    assert np.isclose(np.c_[dx_expected, dy_expected, dz_expected], dXYZ).all()

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
    dXYZ = pyTMD.predict._latitude_dependence(XYZ, SXYZ, LXYZ,
        F2_solar, F2_lunar)
    # assert matching
    assert np.isclose(np.c_[dx_expected, dy_expected, dz_expected], dXYZ).all()

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
    dXYZ = pyTMD.predict._frequency_dependence_diurnal(XYZ, MJD)
    # assert matching
    assert np.isclose(np.c_[dx_expected, dy_expected, dz_expected], dXYZ).all()

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
    dXYZ = pyTMD.predict._frequency_dependence_long_period(XYZ, MJD)
    # assert matching
    assert np.isclose(np.c_[dx_expected, dy_expected, dz_expected], dXYZ).all()

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

def test_fundamental_arguments():
    """Test fundamental (Delaunay) arguments with IERS outputs
    """
    T = 0.07995893223819302
    # convert to MJD from centuries relative to 2000-01-01T12:00:00
    MJD = T*36525.0 + 51544.5
    assert np.isclose(MJD, 54465)
    L_expected = 2.291187512612069099
    LP_expected = 6.212931111003726414
    F_expected = 3.658025792050572989
    D_expected = 4.554139562402433228
    OM_expected = -0.5167379217231804489 + 2.0*np.pi
    # test outputs from function
    l, lp, F, D, Om = pyTMD.astro.delaunay_arguments(MJD)
    # assert matching
    assert np.isclose(L_expected, l)
    assert np.isclose(LP_expected, lp)
    assert np.isclose(F_expected, F)
    assert np.isclose(D_expected, D)
    assert np.isclose(OM_expected, Om)

def test_mean_obliquity():
    """Test that the mean obliquity values matches expected outputs
    """
    MJD = 54465.0
    expected = 0.40907444424006084
    mean_obliquity = pyTMD.astro.mean_obliquity(MJD)
    assert np.isclose(expected, mean_obliquity)

def test_precession_matrix():
    """Test that the precession matrix matches expected outputs
    """
    MJD = 54465.0
    # convert from MJD to centuries relative to 2000-01-01T12:00:00
    T = (MJD - 51544.5)/36525.0
    expected = np.array([
        [ 9.99998100e-01, -1.78795448e-03, -7.76914888e-04],
        [ 1.78795449e-03,  9.99998402e-01, -6.84570121e-07],
        [ 7.76914871e-04, -7.04519640e-07,  9.99999698e-01]
    ])
    P = pyTMD.astro._precession_matrix(T)
    assert np.isclose(expected, P[:,:,0]).all()

def test_nutation_matrix():
    """Test that the nutation matrix matches expected outputs
    """
    MJD = 54465.0
    # convert from MJD to centuries relative to 2000-01-01T12:00:00
    T = (MJD - 51544.5)/36525.0
    expected = np.array([
        [ 9.99998100e-01, -1.78795448e-03, -7.76914888e-04],
        [ 1.78795449e-03,  9.99998402e-01, -6.84570121e-07],
        [ 7.76914871e-04, -7.04519640e-07,  9.99999698e-01]
    ])
    P = pyTMD.astro._precession_matrix(T)
    assert np.isclose(expected, P[:,:,0]).all()

def test_frame_bias_matrix():
    """Test that the frame bias matrix matches expected outputs
    """
    expected = np.array([
        [ 1.00000000e+00, -7.07827974e-08,  8.05614894e-08],
        [ 7.07827974e-08,  1.00000000e+00,  3.30604145e-08],
        [-8.05614894e-08, -3.30604145e-08,  1.00000000e+00]
    ])
    B = pyTMD.astro._frame_bias_matrix()
    assert np.isclose(expected, B).all()

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

# parameterize ephemerides
@pytest.mark.parametrize("EPHEMERIDES", ['approximate','JPL'])
def test_solid_earth_radial(EPHEMERIDES):
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
        EPSG=4326, TYPE='drift', TIME='datetime', ELLIPSOID='WGS84',
        EPHEMERIDES=EPHEMERIDES)
    tide_mean = compute_SET_corrections(longitudes, latitudes, times,
        EPSG=4326, TYPE='drift', TIME='datetime', ELLIPSOID='WGS84',
        TIDE_SYSTEM='mean_tide', EPHEMERIDES=EPHEMERIDES)
    # as using estimated ephemerides, assert within 1/2 mm
    assert np.isclose(tide_earth, tide_free, atol=5e-4).all()
    # check permanent tide offsets (additive correction in ICESat-2)
    tide_expected = tide_earth + tide_earth_free2mean
    predicted = 0.06029 - 0.180873*np.sin(latitudes*np.pi/180.0)**2
    assert np.isclose(tide_expected, tide_mean, atol=5e-4).all()
    assert np.isclose(tide_earth_free2mean, predicted, atol=5e-4).all()
    assert np.isclose(tide_mean-tide_free, predicted, atol=5e-4).all()

# PURPOSE: Download JPL ephemerides from Solar System Dynamics server
@pytest.fixture(scope="module", autouse=True)
def download_jpl_ephemerides():
    """Download JPL ephemerides from Solar System Dynamics server
    """
    # get path to default ephemerides
    de440s = pyTMD.astro._default_kernel
    # download JPL ephemerides if not existing
    if not de440s.exists():
        pyTMD.utilities.from_jpl_ssd(de440s.name)
        # run tests
        yield
        # clean up
        de440s.unlink()
    else:
        # run tests
        yield

def test_solar_ecef():
    """Test solar ECEF coordinates with ephemeride predictions
    """
    # calculate solar ephemerides
    MJD = 55414.0
    x1, y1, z1 = pyTMD.astro.solar_ecef(MJD)
    r1 = np.sqrt(x1**2 + y1**2 + z1**2)
    # predict solar ephemerides
    x2, y2, z2 = pyTMD.astro.solar_ephemerides(MJD)
    r2 = np.sqrt(x2**2 + y2**2 + z2**2)
    # test distances
    assert np.isclose(np.c_[x1,y1,z1], np.c_[x2,y2,z2], atol=1e9).all()
    # test absolute distance
    assert np.isclose(r1, r2, atol=1e9).all()

def test_lunar_ecef():
    """Test lunar ECEF coordinates with ephemeride predictions
    """
    # calculate solar ephemerides
    MJD = 55414.0
    x1, y1, z1 = pyTMD.astro.lunar_ecef(MJD)
    r1 = np.sqrt(x1**2 + y1**2 + z1**2)
    # predict solar ephemerides
    x2, y2, z2 = pyTMD.astro.lunar_ephemerides(MJD)
    r2 = np.sqrt(x2**2 + y2**2 + z2**2)
    # test distances
    assert np.isclose(np.c_[x1,y1,z1], np.c_[x2,y2,z2], atol=5e6).all()
    # test absolute distance
    assert np.isclose(r1, r2, atol=5e6).all()

def test_earth_rotation_angle():
    """Test that the Earth rotation angle (ERA) matches expected outputs
    """
    # create timescale from modified Julian dates
    ts = pyTMD.time.timescale(MJD=55414.0)
    # expected earth rotation angle as fraction of a turn
    expected = 0.8730204642501604
    assert np.isclose(360.0*expected, ts.era).all()

def test_greenwich():
    """Test approximations of Greenwich Hour Angle in degrees
    using Meeus approximation and calculation within pyTMD
    """
    # create timescale from modified Julian dates
    ts = pyTMD.time.timescale(MJD=55414.0)
    # Meeus approximation
    GHA = np.mod(280.46061837504 + 360.9856473662862*(ts.T*36525.0), 360.0)
    # compare with pyTMD calculation
    assert np.isclose(GHA, ts.gha)

def test_sidereal():
    # create timescale from modified Julian dates
    ts = pyTMD.time.timescale(MJD=55414.0)
    # expected side real time in hours
    expected = 20.96154017401333
    assert np.isclose(expected, 24.0*ts.st).all()
