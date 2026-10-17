#!/usr/bin/env python
"""
test_math.py (11/2024)
"""

import pytest
import numpy as np
import pyTMD.ellipse
import pyTMD.math
from scipy.special import factorial, gammaln


def test_asec2rad():
    """
    Tests the conversion of arcseconds to radians
    """
    # test angles in arcseconds
    angles = np.array([-180, -90, 0, 90, 180, 270, 360]) * 3600.0
    # expected values in radians
    exp = np.array(
        [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi]
    )
    # test conversion to radians
    test = pyTMD.math.asec2rad(angles)
    assert np.allclose(exp, test)
    # arcseconds to radians
    atr = np.pi / 648000.0
    # test conversion to radians
    test = angles * atr
    assert np.allclose(exp, test)
    # test reverse conversion
    arcseconds = pyTMD.math.rad2asec(test)
    assert np.allclose(angles, arcseconds)


def test_masec2rad():
    """
    Tests the conversion of microarcseconds to radians
    """
    # test angles in microarcseconds
    angles = np.array([-180, -90, 0, 90, 180, 270, 360]) * 3600.0 * 1e6
    # expected values in radians
    exp = np.array(
        [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi]
    )
    # test conversion to radians
    test = pyTMD.math.masec2rad(angles)
    assert np.allclose(exp, test)
    # microarcseconds to radians
    atr = np.pi / 648e9
    # test conversion to radians
    test = angles * atr
    assert np.allclose(exp, test)
    # test reverse conversion
    microarcseconds = pyTMD.math.rad2masec(test)
    assert np.allclose(angles, microarcseconds)


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


def test_condon_shortley():
    """
    Tests the calculation of the Condon-Shortley phase
    """
    # test values for m
    m = np.array([-2, -1, 0, 1, 2])
    # expected values
    exp = np.array([1, -1, 1, -1, 1])
    # test Condon-Shortley phase
    test = pyTMD.math._condon_shortley(m)
    assert np.all(exp == test)


def test_kronecker_delta():
    """
    Tests the calculation of the Kronecker delta function
    """
    # test values for m
    m = np.array([-2, -1, 0, 1, 2])
    # expected values
    exp = np.array([0, 0, 1, 0, 0])
    # test Kronecker delta function
    test = pyTMD.math._kronecker_delta(m, 0)
    assert np.all(exp == test)


@pytest.mark.parametrize("l", [1, 2, 3, 4, 5, 10, 20, 30, 40, 50])
def test_factorial(l):
    """
    Tests the calculation of factorials
    """
    m = 0
    r = np.floor_divide(l - m, 2)
    for k in range(r + 1):
        # calculate factorial terms for Legendre polynomials
        Fexact = (
            factorial(2 * l - 2 * k)
            / factorial(k)
            / factorial(l - k)
            / factorial(l - m - 2 * k)
        )
        # calculate approximate factorial terms using gamma functions
        Fgamma = np.exp(
            gammaln(2.0 * l - 2.0 * k + 1)
            - gammaln(k + 1)
            - gammaln(l - k + 1)
            - gammaln(l - m - 2.0 * k + 1)
        )
        assert np.allclose(Fexact, Fgamma)
    # for each order
    for m in range(l + 1):
        # calculate factorial terms for the normalization
        # using integer arithmetic for exact calculations
        Fexact = factorial(l - m) / factorial(l + m)
        # using gamma functions for cases with potential overflow
        Fgamma = np.exp(gammaln(l - m + 1) - gammaln(l + m + 1))
        assert np.allclose(Fexact, Fgamma)


@pytest.mark.parametrize("l", [1, 2, 3])
def test_assoc_legendre(l, x=[-1.0, -0.9, -0.8]):
    """test the calculation of unnormalized Legendre polynomials"""
    # calculate Legendre polynomials
    nx = len(x)
    obs = np.zeros((l + 1, nx))
    test = np.zeros((l + 1, nx))
    for m in range(l + 1):
        obs[m, :], _ = pyTMD.math.legendre(l, x, m=m)
        test[m, :] = pyTMD.math._assoc_legendre(l, m, x)
    # expected values for each spherical harmonic degree
    if l == 1:
        expected = np.array(
            [
                [-1.00000, -0.90000, -0.80000],
                [0.00000, -0.43589, -0.60000],
            ]
        )
    elif l == 2:
        expected = np.array(
            [
                [1.00000, 0.71500, 0.46000],
                [0.00000, 1.17690, 1.44000],
                [0.00000, 0.57000, 1.08000],
            ]
        )
    elif l == 3:
        expected = np.array(
            [
                [-1.00000, -0.47250, -0.08000],
                [0.00000, -1.99420, -1.98000],
                [0.00000, -2.56500, -4.32000],
                [0.00000, -1.24229, -3.24000],
            ]
        )
    # check with expected values
    assert np.allclose(obs, expected, atol=1e-05)
    # check that the two methods give the same values
    assert np.allclose(obs, test)


@pytest.mark.parametrize("l", [1, 2, 3, 4])
def test_legendre(l):
    """test the calculation of unnormalized Legendre polynomials"""
    # test over the range of latitudes
    lat = np.arange(-90, 91, 1)
    th = np.radians(90.0 - lat)
    # test values for x
    x = np.cos(th)
    u = np.sqrt(1.0 - x**2)
    for m in range(l + 1):
        obs, dobs = pyTMD.math.legendre(l, x, m=m)
        # since tides only use low-degree harmonics:
        # functions are hard coded rather than using a recursion relation
        if (l == 0) and (m == 0):
            Plm = 1.0
            dPlm = 0.0
        elif (l == 1) and (m == 0):
            Plm = x
            dPlm = -u
        elif (l == 1) and (m == 1):
            Plm = u
            dPlm = x
        elif (l == 2) and (m == 0):
            Plm = 0.5 * (3.0 * x**2 - 1.0)
            dPlm = -3.0 * u * x
        elif (l == 2) and (m == 1):
            Plm = 3.0 * x * u
            dPlm = 3.0 * (1.0 - 2.0 * u**2)
        elif (l == 2) and (m == 2):
            Plm = 3.0 * u**2
            dPlm = 6.0 * u * x
        elif (l == 3) and (m == 0):
            Plm = 0.5 * (5.0 * x**2 - 3.0) * x
            dPlm = u * (1.5 - 7.5 * x**2)
        elif (l == 3) and (m == 1):
            Plm = 1.5 * (5.0 * x**2 - 1.0) * u
            dPlm = -1.5 * x * (10.0 * u**2 - 5.0 * x**2 + 1.0)
        elif (l == 3) and (m == 2):
            Plm = 15.0 * x * u**2
            dPlm = 15.0 * u * (3.0 * x**2 - 1.0)
        elif (l == 3) and (m == 3):
            Plm = 15.0 * u**3
            dPlm = 45.0 * x * u**2
        elif (l == 4) and (m == 0):
            Plm = 0.125 * (35.0 * x**4 - 30.0 * x**2 + 3.0)
            dPlm = -2.5 * (7.0 * x**2 - 3.0) * u * x
        elif (l == 4) and (m == 1):
            Plm = 2.5 * (7.0 * x**2 - 3.0) * u * x
            dPlm = 2.5 * (28.0 * x**4 - 27.0 * x**2 + 3.0)
        elif (l == 4) and (m == 2):
            Plm = 7.5 * (7.0 * x**2 - 1.0) * u**2
            dPlm = (105 * x**2 - 105 * u**2 - 15.0) * u * x
        elif (l == 4) and (m == 3):
            Plm = 105.0 * x * u**3
            dPlm = (420.0 * x**2 - 105.0) * u**2
        elif (l == 4) and (m == 4):
            Plm = 105.0 * u**4
            dPlm = 420.0 * x * u**3
        # apply Condon-Shortley phase
        Plm *= np.power(-1.0, m)
        dPlm *= np.power(-1.0, m)
        # check with expected values
        assert np.allclose(obs, Plm, atol=1e-05)
        assert np.allclose(dobs, dPlm, atol=1e-05)


def test_legendre_hw95():
    """test the calculation of fully-normalized Legendre polynomials
    and their derivative versus values from HW95
    """
    # colatitude for test values
    theta = np.radians(np.array([30.0, 147.86]))
    # validation values from HW95 (l, m, Plm, dPlm)
    validation = np.array([
        [0, 0, 1.0, 0.0, 1.0, -0.0],
        [1, 0, 1.5, -0.86602540378444, -1.46661528359239, -0.92143345388216],
        [1, 1, 0.86602540378444, 1.5, 0.92143345388216, -1.46661528359239],
        [2, 0, 1.39754248593737, -2.90473750965556, 1.28681283579950, 3.0217962957189],
        [2, 1, 1.67705098312484, 1.93649167310371, -1.7446349047695, 1.68077249984623],
        [2, 2, 0.48412291827593, 1.67705098312484, 0.5480527115903, -1.74463490476952],
        [3, 0, 0.85923294280422, -5.45686207907072, -0.65521141818753, -5.45749169751143],
        [3, 1, 2.2277546150777, 0.35078038001005, 2.22801165572977, 0.336385482397345],
        [3, 2, 1.10926495933118, 3.20217211436237, -1.22779802359212, 3.13710012570862],
        [3, 3, 0.2614562582919, 1.35856656995526, 0.31491915259677, -1.50373933249918],
        [4, 0, 0.0703125, -7.3070893444312, -0.1939318382809, 6.82082063071094],
        [4, 1, 2.31070453947492, -3.55756236768943, -2.15693287045128, -4.6596516566514],
        [4, 2, 1.78186666957014, 3.63092188706945, 1.9074843406585, -3.07893114132216],
        [4, 3, 0.67928328497763, 3.13747509950278, -0.799974065904, 3.31727574631542],
        [4, 4, 0.1386581199164, 0.96065163430871, 0.1776964222809, -1.13133417354818],
        [5, 0, -0.74051002865529, -7.19033890096581, 0.98385838968044, -5.81480164751363],
        [5, 1, 1.85653752113519, -8.95158333012718, 1.50137532948798, 10.01062398467376],
        [5, 2, 2.29938478949397, 1.85857059805883, -2.34343900747049, 0.4845821062241],
        [5, 3, 1.24653144252643, 4.78747153809058, 1.42384088377768, -4.68161769034657],
        [5, 4, 0.39826512815546, 2.52932326844337, -0.49903453974047, 2.86365789536277],
        [5, 5, 0.07271293151948, 0.62971245879506, 0.09914672538807, -0.78904288833685],
        [6, 0, -1.34856068213155, -4.35442243247701, -1.46613935548827, 1.93862045690768],
        [6, 1, 0.95021287641141, -14.00557979016896, -0.42304166607473, -14.11073063990035],
        [6, 2, 2.47470311782905, -2.56294916449777, 2.33756705084144, 4.76570656668609],
        [6, 3, 1.85592870532597, 5.20453026842398, -2.03449392554487, 4.31069136306377],
        [6, 4, 0.81047568870385, 4.55019988574613, 0.98663447485822, -4.86180775544823],
        [6, 5, 0.22704605589841, 1.83519142087945, -0.30269526691992, 2.21877545058208],
        [6, 6, 0.0378410093164, 0.39325530447417, 0.05489879051610, -0.52428358151593],
    ])
    # test for two points provided with HW95 catalog
    PLM = np.zeros((2))
    DPLM = np.zeros((2))
    # calculate legendre polynomials with Mohlenkamp recursion relation
    P1, D1 = pyTMD.math.legendreP(6, np.cos(theta))
    # check each row of values
    for l, m, PLM[0], DPLM[0], PLM[1], DPLM[1] in validation:
        # verify that degree and order are integers
        l = int(l)
        m = int(m)
        # compare HW95 against recursion relation
        assert np.allclose(P1[l, m], PLM, atol=1e-12)
        assert np.allclose(D1[l, m], DPLM, atol=1e-12)
        # HW95 normalization of degree l and order m:
        # fully-normalized without Condon-Shortley phase
        cs = pyTMD.math._condon_shortley(m)
        # normalization factors (unapply Condon-Shortley phase)
        norm = cs * pyTMD.math._legendre_norm(l, m, normalization="full")
        # Legendre polynomials and their first derivative
        P2, D2 = pyTMD.math.legendre(l, np.cos(theta), m=m, norm=norm)
        assert np.allclose(P2, PLM, atol=1e-12)
        assert np.allclose(D2, DPLM, atol=1e-12)


# PURPOSE: test the calculation of ellipse coordinates
def test_ellipse_xy():
    """test the calculation of ellipse coordinates"""
    # number of points
    npts = 30
    # define ellipse parameters
    umajor = 4.0 + 2.0 * np.random.rand(npts)
    uminor = 2.0 + np.random.rand(npts)
    uincl = 180.0 * np.random.rand(npts)
    # center of the ellipse
    xy = (10.0 - 20.0 * np.random.rand(1), 10.0 - 20.0 * np.random.rand(1))
    # calculate coordinates
    x, y = pyTMD.ellipse._xy(umajor, uminor, uincl, phase=0.0, xy=xy)
    # verify that the coordinates match the ellipse equation
    phi = uincl * np.pi / 180.0
    X = (x - xy[0]) * np.cos(phi) + (y - xy[1]) * np.sin(phi)
    Y = -(x - xy[0]) * np.sin(phi) + (y - xy[1]) * np.cos(phi)
    test = (uminor * X) ** 2 + (umajor * Y) ** 2
    validation = (umajor * uminor) ** 2
    assert np.allclose(test, validation)
