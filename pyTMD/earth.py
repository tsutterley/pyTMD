#!/usr/bin/env python
"""
earth.py
Written by Tyler Sutterley (09/2026)
Calculates Earth reference and deformation parameters

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 09/2026: add reader function for load Love number tables
    Written 05/2026
"""

import re
import pathlib
import numpy as np
from pyTMD.constituents import frequency
from pyTMD.utilities import get_data_path

__all__ = [
    "_ellipsoids",
    "_units",
    "datum",
    "love_numbers",
    "complex_love_numbers",
    "degree_love_numbers",
    "load_love_numbers",
    "adjust_load_love_numbers",
    "greens_function",
    "_melchior_table_52",
    "_infconv",
]

_ellipsoids = [
    "CLK66",
    "CLK80",
    "GRS67",
    "GRS80",
    "WGS60",
    "WGS66",
    "WGS72",
    "WGS84",
    "ATS77",
    "NAD27",
    "NAD83",
    "INTER",
    "KRASS",
    "SGS90",
    "AIRY",
    "MAIRY",
    "HLMRT",
    "HOUGH",
    "HGH80",
    "MERIT",
    "TOPEX",
    "EGM96",
    "IAG75",
    "IAU64",
    "IAU76",
    "IERS89",
    "IERS",
]
_units = ["MKS", "CGS"]


class datum:
    """
    Class for gravitational and ellipsoidal parameters
    :cite:p:`HofmannWellenhof:2006hy,Urban:2013vl`

    Parameters
    ----------
    ellipsoid: str, default 'WGS84'
        Reference ellipsoid name

            - ``'CLK66'``: Clarke 1866
            - ``'CLK80'``: Clarke 1880
            - ``'GRS67'``: Geodetic Reference System 1967
            - ``'GRS80'``: Geodetic Reference System 1980
            - ``'HGH80'``: Hughes 1980 Ellipsoid
            - ``'WGS60'``: World Geodetic System 1960
            - ``'WGS66'``: World Geodetic System 1966
            - ``'WGS72'``: World Geodetic System 1972
            - ``'WGS84'``: World Geodetic System 1984
            - ``'ATS77'``: Quasi-earth centred ellipsoid for ATS77
            - ``'NAD27'``: North American Datum 1927
            - ``'NAD83'``: North American Datum 1983
            - ``'INTER'``: International
            - ``'KRASS'``: Krassovsky (USSR)
            - ``'SGS90'``: Soviet Geodetic System 1990
            - ``'HLMRT'``: Helmert 1906 Ellipsoid
            - ``'HOUGH'``: Hough 1960 Ellipsoid
            - ``'AIRY'``: Airy (1830)
            - ``'MAIRY'``: Modified Airy (1849)
            - ``'MERIT'``: MERIT 1983 ellipsoid
            - ``'TOPEX'``: TOPEX/POSEIDON ellipsoid
            - ``'EGM96'``: EGM 1996 gravity model
            - ``'IAG75'``: International Association of Geodesy (1975)
            - ``'IAU64'``: International Astronomical Union (1964)
            - ``'IAU76'``: International Astronomical Union (1976)
            - ``'IERS89'``: IERS Numerical Standards (1989)
            - ``'IERS'``: IERS Numerical Standards (2010)
    units: str, default `MKS`
        Output units

            - ``'MKS'``: meters, kilograms, seconds
            - ``'CGS'``: centimeters, grams, seconds

    Attributes
    ----------
    a_axis: float
        Semi-major axis of the ellipsoid
    flat: float
        Flattening of the ellipsoid
    omega: float
        Angular velocity of the Earth
    GM: float
        Geocentric gravitational constant
    """

    np.seterr(invalid="ignore")

    def __init__(self, **kwargs):
        # set default keyword arguments
        kwargs.setdefault("ellipsoid", "WGS84")
        kwargs.setdefault("units", "MKS")
        kwargs.setdefault("a_axis", None)
        kwargs.setdefault("flat", None)
        kwargs.setdefault("GM", None)
        kwargs.setdefault("omega", None)
        # set ellipsoid name and units
        self.units = kwargs["units"].upper()
        if (kwargs["a_axis"] is not None) and (kwargs["flat"] is not None):
            self.name = "user_defined"
        else:
            self.name = kwargs["ellipsoid"].upper()
        # validate ellipsoid and units
        if self.name not in _ellipsoids + ["user_defined"]:
            raise ValueError("Invalid ellipsoid")
        if self.units not in _units:
            raise ValueError("Invalid output units")

        # set parameters for ellipsoid
        if self.name in ("CLK66", "NAD27"):
            # Clarke 1866
            self.a_axis = 6378206.4  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 294.9786982  # flattening of the ellipsoid

        elif self.name == "CLK80":
            # Clarke 1880
            self.a_axis = 6378249.145  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 293.465  # flattening of the ellipsoid

        elif self.name in ("GRS80", "NAD83"):
            # Geodetic Reference System 1980
            # North American Datum 1983
            self.a_axis = 6378135.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.26  # flattening of the ellipsoid
            self.GM = 3.986005e14  # [m^3/s^2] Geocentric Gravitational Constant

        elif self.name == "GRS67":
            # Geodetic Reference System 1967
            self.a_axis = 6378160.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.247167427  # flattening of the ellipsoid
            self.GM = 3.98603e14  # [m^3/s^2] Geocentric Gravitational Constant
            self.omega = (
                7292115.1467e-11  # angular velocity of the Earth [rad/s]
            )

        elif self.name == "WGS60":
            # World Geodetic System 1960
            self.a_axis = 6378165.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.30  # flattening of the ellipsoid

        elif self.name == "WGS66":
            # World Geodetic System 1966
            self.a_axis = 6378145.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.25  # flattening of the ellipsoid

        elif self.name == "WGS72":
            # World Geodetic System 1972
            self.a_axis = 6378135.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.26  # flattening of the ellipsoid

        elif self.name == "WGS84":
            # World Geodetic System 1984
            self.a_axis = 6378137.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.257223563  # flattening of the ellipsoid

        elif self.name == "ATS77":
            # Quasi-earth centred ellipsoid for ATS77
            self.a_axis = 6378135.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.257  # flattening of the ellipsoid

        elif self.name == "KRASS":
            # Krassovsky (USSR)
            self.a_axis = 6378245.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.3  # flattening of the ellipsoid

        elif self.name == "SGS90":
            # Soviet Geodetic System 1990
            self.a_axis = 6378136.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.2578393  # flattening of the ellipsoid

        elif self.name == "INTER":
            # International 1924
            self.a_axis = 6378388.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1 / 297.0  # flattening of the ellipsoid

        elif self.name == "AIRY":
            # Airy 1830 Ellipsoid
            self.a_axis = 6377563.396  # [m] semi-major axis of the ellipsoid
            self.flat = 1 / 299.3249646  # flattening of the ellipsoid

        elif self.name == "MAIRY":
            # Modified Airy 1849 Ellipsoid
            self.a_axis = 6377340.189  # [m] semi-major axis of the ellipsoid
            self.flat = 1 / 299.3249646  # flattening of the ellipsoid

        elif self.name == "HLMRT":
            # Helmert 1906 Ellipsoid
            self.a_axis = 6378200.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.3  # flattening of the ellipsoid

        elif self.name == "HOUGH":
            # Hough 1960 Ellipsoid
            self.a_axis = 6378270.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 297.0  # flattening of the ellipsoid

        elif self.name == "HGH80":
            # Hughes 1980 Ellipsoid
            self.a_axis = 6378273.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.279411123064  # flattening of the ellipsoid

        elif self.name == "MERIT":
            # MERIT 1983 ellipsoid
            self.a_axis = 6378137.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.257  # flattening of the ellipsoid

        elif self.name == "TOPEX":
            # TOPEX/POSEIDON ellipsoid
            self.a_axis = 6378136.3  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.257  # flattening of the ellipsoid
            self.GM = 3.986004415e14  # [m^3/s^2]

        elif self.name == "EGM96":
            # EGM 1996 gravity model
            self.a_axis = 6378136.3  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.256415099  # flattening of the ellipsoid
            self.GM = 3.986004415e14  # [m^3/s^2]

        elif self.name == "IAG75":
            # International Association of Geodesy (IAG 1975)
            self.a_axis = 6378140.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.257  # flattening of the ellipsoid

        elif self.name == "IAU64":
            # International Astronomical Union (IAU 1964)
            self.a_axis = 6378160.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.25  # flattening of the ellipsoid

        elif self.name == "IAU76":
            # International Astronomical Union (IAU 1964)
            self.a_axis = 6378140.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.257  # flattening of the ellipsoid

        elif self.name == "IERS89":
            # IERS Numerical Standards (1989)
            self.a_axis = 6378136.0  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.257  # flattening of the ellipsoid

        elif self.name == "IERS":
            # IERS Numerical Standards (2010)
            self.a_axis = 6378136.6  # [m] semi-major axis of the ellipsoid
            self.flat = 1.0 / 298.25642  # flattening of the ellipsoid

        elif self.name == "user_defined":
            # custom datum
            self.a_axis = np.float64(kwargs["a_axis"])
            self.flat = np.float64(kwargs["flat"])

        # set default parameters if not listed as part of ellipsoid
        # Geocentric Gravitational Constant
        if kwargs["GM"] is not None:
            # user defined Geocentric Gravitational Constant
            self.GM = np.float64(kwargs["GM"])
        elif self.name not in ("GRS80", "GRS67", "NAD83", "TOPEX", "EGM96"):
            # for ellipsoids not listing the Geocentric Gravitational Constant
            self.GM = 3.986004418e14  # [m^3/s^2]

        # angular velocity of the Earth
        if kwargs["omega"] is not None:
            # user defined angular velocity of the Earth
            self.omega = np.float64(kwargs["omega"])
        elif self.name not in ("GRS67"):
            # for ellipsoids not listing the angular velocity of the Earth
            self.omega = 7292115e-11  # [rad/s]

        # universal gravitational constant [N*m^2/kg^2]
        self.G = 6.67430e-11

        # standard gravitational acceleration [m/s^2]
        # (World Meteorological Organization)
        self.gamma = 9.80665

        # convert units to CGS
        if self.units == "CGS":
            self.a_axis *= 100.0
            self.GM *= 1e6
            self.G *= 1000.0  # [dyn*cm^2/g^2]
            self.gamma *= 100.0

    # mean radius of the Earth having the same volume
    # (4pi/3)R^3 = (4pi/3)(a^2)b = (4pi/3)(a^3)(1 - f)
    @property
    def rad_e(self) -> float:
        """Average radius of the Earth with same volume as ellipsoid"""
        return self.a_axis * (1.0 - self.flat) ** (1.0 / 3.0)

    # semiminor axis of the ellipsoid
    @property
    def b_axis(self) -> float:
        """Semi-minor axis of the ellipsoid"""
        return (1.0 - self.flat) * self.a_axis

    # Ratio between ellipsoidal axes
    @property
    def ratio(self) -> float:
        """Ratio between ellipsoidal axes"""
        return 1.0 - self.flat

    # Polar radius of curvature
    @property
    def rad_p(self) -> float:
        """Polar radius of curvature"""
        return self.a_axis / (1.0 - self.flat)

    # Linear eccentricity
    @property
    def ecc(self) -> float:
        """Linear eccentricity"""
        return np.sqrt((2.0 * self.flat - self.flat**2) * self.a_axis**2)

    # first numerical eccentricity
    @property
    def ecc1(self) -> float:
        """First numerical eccentricity"""
        return self.ecc / self.a_axis

    # second numerical eccentricity
    @property
    def ecc2(self) -> float:
        """Second numerical eccentricity"""
        return self.ecc / self.b_axis

    # m parameter [omega^2*a^2*b/(GM)]
    # p. 70, Eqn.(2-137)
    @property
    def m(self) -> float:
        """:math:`m` Parameter"""
        return self.omega**2 * ((1 - self.flat) * self.a_axis**3) / self.GM

    # flattening f2 component
    # p. 80, Eqn.(2-200)
    @property
    def f2(self) -> float:
        """:math:`f_2` component"""
        return (
            -self.flat
            + (5.0 / 2.0) * self.m
            + (1.0 / 2.0) * self.flat**2.0
            - (26.0 / 7.0) * self.flat * self.m
            + (15.0 / 4.0) * self.m**2.0
        )

    # flattening f4 component
    # p. 80, Eqn.(2-200)
    @property
    def f4(self) -> float:
        """:math:`f_4` component"""
        return -(1.0 / 2.0) * self.flat**2.0 + (5.0 / 2.0) * self.flat * self.m

    # q
    # p. 67, Eqn.(2-113)
    @property
    def q(self) -> float:
        """:math:`q` Parameter"""
        return 0.5 * (
            (1.0 + 3.0 / (self.ecc2**2)) * np.arctan(self.ecc2)
            - 3.0 / self.ecc2
        )

    # q_0
    # p. 67, Eqn.(2-113)
    @property
    def q0(self) -> float:
        r""":math:`q_0` Parameter"""
        return (
            3
            * (1.0 + 1.0 / (self.ecc2**2))
            * (1.0 - 1.0 / self.ecc2 * np.arctan(self.ecc2))
            - 1.0
        )

    # J_2 p. 75 Eqn.(2-167), p. 76 Eqn.(2-172)
    @property
    def J2(self) -> float:
        """Oblateness :math:`J_2` coefficient"""
        return (
            (self.ecc1**2)
            * (1.0 - 2.0 * self.m * self.ecc2 / (15.0 * self.q))
            / 3.0
        )

    # Normalized C20 harmonic
    # p. 60, Eqn.(2-80)
    @property
    def C20(self) -> float:
        r"""Normalized :math:`C_{20}` harmonic"""
        return -self.J2 / np.sqrt(5.0)

    # Normal gravity at the equator
    # p. 79, Eqn.(2-286)
    @property
    def gamma_a(self) -> float:
        """Normal gravity at the equator"""
        return (self.GM / (self.a_axis * self.b_axis)) * (
            1.0 - (3.0 / 2.0) * self.m - (3.0 / 14.0) * self.ecc2**2.0 * self.m
        )

    # Normal gravity at the pole
    # p. 79, Eqn.(2-286)
    @property
    def gamma_b(self) -> float:
        """Normal gravity at the pole"""
        return (self.GM / (self.a_axis**2)) * (
            1.0 + self.m + (3.0 / 7.0) * self.ecc2**2.0 * self.m
        )

    # Normal gravity at location
    # p. 80, Eqn.(2-199)
    def gamma_0(self, theta: float | np.ndarray) -> float:
        """Normal gravity at colatitudes

        Parameters
        ----------
        theta: float
            Colatitudes (radians)
        """
        return self.gamma_a * (
            1.0
            + self.f2 * np.cos(theta) ** 2.0
            + self.f4 * np.cos(theta) ** 4.0
        )

    # Normal gravity at location
    # p. 82, Eqn.(2-215)
    def gamma_h(
        self,
        theta: float | np.ndarray,
        height: float | np.ndarray,
    ) -> float:
        """Normal gravity at colatitudes and heights

        Parameters
        ----------
        theta: float
            Colatitudes (radians)
        height: float
            Height above ellipsoid (same as ``units``)
        """
        return self.gamma_0(theta) * (
            1.0
            - (2.0 / self.a_axis)
            * (
                1.0
                + self.flat
                + self.m
                - 2.0 * self.flat * np.cos(theta) ** 2.0
            )
            * height
            + (3.0 / self.a_axis**2.0) * height**2.0
        )

    # ratio between gravity at pole versus gravity at equator
    @property
    def dk(self) -> float:
        """Ratio between gravity at pole versus gravity at equator"""
        return self.b_axis * self.gamma_b / (self.a_axis * self.gamma_b) - 1.0

    # Normal potential at the ellipsoid
    # p. 68, Eqn.(2-123)
    @property
    def U0(self) -> float:
        """Normal potential at the ellipsoid"""
        return (
            self.GM / self.ecc * np.arctan(self.ecc2)
            + (1.0 / 3.0) * self.omega**2 * self.a_axis**2
        )

    # Surface area of the reference ellipsoid
    @property
    def area(self) -> float:
        """Surface area of the ellipsoid"""
        return (
            np.pi
            * self.a_axis**2.0
            * (
                2.0
                + ((1.0 - self.ecc1**2) / self.ecc1)
                * np.log((1.0 + self.ecc1) / (1.0 - self.ecc1))
            )
        )

    # Volume of the reference ellipsoid
    @property
    def volume(self) -> float:
        """Volume of the ellipsoid"""
        return (
            (4.0 * np.pi / 3.0)
            * (self.a_axis**3.0)
            * (1.0 - self.ecc1**2.0) ** 0.5
        )

    # Average density
    @property
    def rho_e(self) -> float:
        """Average density"""
        return self.GM / (self.G * self.volume)

    def __str__(self):
        """String representation of the ``datum`` object"""
        properties = ["pyTMD.datum"]
        properties.append(f"    name: {self.name}")
        properties.append(f"    units: {self.units}")
        return "\n".join(properties)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)


# helper variables with common ellipsoid parameters
# using meters-kilogram-seconds standard
_wgs84 = datum(ellipsoid="WGS84", units="MKS")
_iers = datum(ellipsoid="IERS", units="MKS")


def love_numbers(
    omega: np.ndarray,
    model: str = "PREM",
    **kwargs,
):
    """
    Compute the body tide Love/Shida numbers for a given frequency
    :cite:p:`Wahr:1979vx,Wahr:1981ea,Wahr:1981if,Mathews:1995go`

    Parameters
    ----------
    omega: np.ndarray
        Angular frequency (radians per second)
    model: str, default 'PREM'
        Earth model to use for Love numbers

            - '1066A'
            - '1066A-N' (neutral)
            - '1066A-S' (stable)
            - 'PEM-C'
            - 'C2'
            - 'PREM'
    astype: np.dtype, default np.float64
        Data type for the output Love numbers

    Returns
    -------
    h2: float
        Degree-2 Love number of vertical displacement
    k2: float
        Degree-2 Love number of gravitational potential
    l2: float
        Degree-2 Love (Shida) number of horizontal displacement
    """
    # set default keyword arguments
    kwargs.setdefault("astype", np.float64)
    # free core nutation frequencies (cycles per sidereal day) and
    # Love number parameters from Wahr (1981) table 6
    # and Mathews et al. (1995) table 3
    if model == "1066A":
        lambda_fcn = 1.0021714
        h0, h1 = np.array([6.03e-1, -2.46e-3])
        k0, k1 = np.array([2.98e-1, -1.23e-3])
        l0, l1 = np.array([8.42e-2, 7.81e-5])
    elif model == "1066A-N":
        lambda_fcn = 1.0021716
        h0, h1 = np.array([6.03e-1, -2.46e-3])
        k0, k1 = np.array([2.98e-1, -1.24e-3])
        l0, l1 = np.array([8.42e-2, 7.82e-5])
    elif model == "1066A-S":
        lambda_fcn = 1.0021708
        h0, h1 = np.array([6.03e-1, -2.46e-3])
        k0, k1 = np.array([2.98e-1, -1.24e-3])
        l0, l1 = np.array([8.42e-2, 7.83e-5])
    elif model == "PEM-C":
        lambda_fcn = 1.0021771
        h0, h1 = np.array([6.02e-1, -2.46e-3])
        k0, k1 = np.array([2.98e-1, -1.24e-3])
        l0, l1 = np.array([8.39e-2, 7.69e-5])
    elif model == "C2":
        lambda_fcn = 1.0021844
        h0, h1 = np.array([6.02e-1, -2.45e-3])
        k0, k1 = np.array([2.98e-1, -1.23e-3])
        l0, l1 = np.array([8.46e-2, 7.58e-5])
    elif model == "PREM":
        lambda_fcn = 1.0023214
        h0, h1 = np.array([5.994e-1, -2.532e-3])
        k0, k1 = np.array([2.962e-1, -1.271e-3])
        l0, l1 = np.array([8.378e-2, 7.932e-5])
    else:
        raise ValueError(f"Unknown Earth model: {model}")
    # Love numbers for different frequency bands
    if omega > 1e-4:
        # tides in the semi-diurnal band
        h2 = 0.609
        k2 = 0.302
        l2 = 0.0852
    elif omega < 2e-5:
        # tides in the long period band
        h2 = 0.606
        k2 = 0.299
        l2 = 0.0840
    else:
        # calculate love numbers following J. Wahr (1979)
        # use resonance formula for tides in the diurnal band
        # frequency of the o1 tides (radians/second)
        (omega_o1,) = frequency("o1", **kwargs)
        # convert frequency from cycles per sidereal day
        # frequency of free core nutation (radians/second)
        omega_fcn = lambda_fcn * 7292115e-11
        # Love numbers for frequency using equation 4.18 of Wahr (1981)
        # (simplification to use only the free core nutation term)
        ratio = (omega - omega_o1) / (omega_fcn - omega)
        h2 = h0 + h1 * ratio
        k2 = k0 + k1 * ratio
        l2 = l0 + l1 * ratio
    # return the Love numbers for frequency
    return np.array([h2, k2, l2], dtype=kwargs["astype"])


def complex_love_numbers(omega: np.ndarray, **kwargs):
    """
    Compute the complex body tide Love/Shida numbers with in-phase
    and out-of-phase components for a given frequency
    :cite:p:`Mathews:1997js,Mathews:2002cr,Petit:2010tp`

    Parameters
    ----------
    omega: np.ndarray
        Angular frequency (radians per second)
    kwargs: dict
        Keyword arguments for :py:func:`pyTMD.earth.love_numbers`

    Returns
    -------
    h2: complex
        Degree-2 Love number of vertical displacement
    k2: complex
        Degree-2 Love number of gravitational potential
    l2: complex
        Degree-2 Love (Shida) number of horizontal displacement
    """
    # number of sidereal days per solar day
    sidereal_ratio = 1.002737909
    # number of seconds in a sidereal day (approximately 86164.1)
    sidereal_day = 86400.0 / sidereal_ratio
    # frequency in cycles per sidereal day
    f = omega * sidereal_day / (2.0 * np.pi)
    # Love numbers for different frequency bands
    if omega == 0.0:
        # use real-valued body tide Love numbers for the permanent tide
        # to prevent singularities in frequency-dependent model
        h2, k2, l2 = love_numbers(omega, **kwargs)
    elif omega > 1e-4:
        # in-phase and out-of-phase components for the semi-diurnal band
        # table 7.3a (IERS conventions 2010)
        h2 = 0.6078 - 0.0022j
        # table 6.5c (IERS conventions 2010)
        k2 = 0.30102 - 0.0013j
        # table 7.3a (IERS conventions 2010)
        l2 = 0.0847 - 0.0007j
    elif omega < 2e-5:
        # compute in-phase and out-of-phase components for the long period band
        # variation largely due to mantle anelasticity
        alpha = 0.15
        # frequency equivalent to 200s
        fm = sidereal_day / 200.0
        factor = np.tan(alpha * np.pi / 2.0) ** (-1)
        anelasticity_model = (
            factor * (1.0 - (fm / f) ** alpha) + 1j * (fm / f) ** alpha
        )
        # model for the variation of Love numbers across the zonal tide band
        # equation 7.4a (IERS conventions 2010)
        h2 = 0.5998 - 9.96e-4 * anelasticity_model
        # equation 6.12 (IERS conventions 2010)
        k2 = 0.29525 - 5.796e-4 * anelasticity_model
        # equation 7.4b (IERS conventions 2010)
        l2 = 0.0831 - 3.01e-4 * anelasticity_model
    else:
        # in-phase and out-of-phase components for the diurnal band
        # following IERS conventions and Mathews et al. (2002)
        # values from equation 6.10 of IERS conventions 2010
        # and from Mathews et al. (2002)
        sigma = np.zeros((4), dtype=np.complex128)
        # factor for calculating L0
        sigma[0] = f - 1.0
        # Chandler wobble
        sigma[1] = -0.0026010 - 0.0001361j
        # retrograde free core nutation
        sigma[2] = 1.0023181 + 0.000025j
        # prograde free core nutation
        sigma[3] = 0.999026 + 0.000780j
        # frequency dependence of Love number h2 (vertical)
        # table 7.1 (IERS conventions 2010)
        H2 = np.zeros((4), dtype=np.complex128)
        H2[0] = 0.60671 - 0.242e-2j
        H2[1] = -0.15777e-2 - 0.7630e-4j
        H2[2] = 0.18053e-3 - 0.6292e-5j
        H2[3] = -0.18616e-5 + 0.1379e-6j
        # frequency dependence of Love number k2 (potential)
        # table 6.4 (IERS conventions 2010)
        K2 = np.zeros((4), dtype=np.complex128)
        K2[0] = 0.29954 - 0.1412e-2j
        K2[1] = -0.77896e-3 - 0.3711e-4j
        K2[2] = 0.90963e-4 - 0.2963e-5j
        K2[3] = -0.11416e-5 + 0.5325e-7j
        # frequency dependence of Love number l2 (horizontal)
        # table 7.1 (IERS conventions 2010)
        L2 = np.zeros((4), dtype=np.complex128)
        L2[0] = 0.84963e-1 - 0.7395e-3j
        L2[1] = -0.22107e-3 - 0.9646e-5j
        L2[2] = -0.54710e-5 - 0.2990e-6j
        L2[3] = -0.29904e-7 - 0.7717e-8j
        # estimate the complex Love numbers for diurnal tides
        # equation 6.9 (IERS conventions 2010)
        h2 = np.sum(H2 / (f - sigma))
        k2 = np.sum(K2 / (f - sigma))
        l2 = np.sum(L2 / (f - sigma))

    # return the Love numbers as a complex number
    # to include the in-phase and out-of-phase components
    return np.array([h2, k2, l2], dtype=np.complex128)


def degree_love_numbers(
    l: int,
    model: str = "Longman",
    **kwargs,
):
    """
    Extracts body tide Love/Shida numbers for a given degree
    :cite:p:`Melchior:1983wd`

    Parameters
    ----------
    l: int
        Degree of the spherical harmonics
    model: str, default "Longman"
        Earth model

        - ``'Longman'``: :cite:t:`Longman:1959hw`
        - ``'Kaula'``: :cite:t:`Kaula:1966un`
        - ``'Takeuchi'``: :cite:t:`Takeuchi:1950hi`

    Returns
    -------
    hl: float
        Body tide Love number of vertical displacement
    kl: float
        Body tide Love number of gravitational potential
    ll: float
        Body tide Love (Shida) number of horizontal displacement
    """
    # get the table of body tide love numbers
    table = _melchior_table_52(model)
    # provide zero values for degrees not provided in the table
    lmin, lmax = np.array([table[0, 0], table[-1, 0]], dtype="i")
    if (l < lmin) | (l > lmax):
        return (0.0, 0.0, 0.0)
    # column 1: Spherical harmonic degree
    n = table[l - lmin, 0]
    # verify that the rows match
    if n != l:
        raise ValueError(f"Mismatched row for degree {l}")
    # column 2: Love number of vertical displacement
    hl = table[l - lmin, 1]
    # column 3: Love number of gravitational potential
    kl = table[l - lmin, 2]
    # column 4: Shida number of horizontal displacement
    ll = table[l - lmin, 3]
    # return the Love numbers for degree l
    return (hl, kl, ll)


# tables of load Love/Shida numbers
# PREM outputs from Han and Wahr (1995)
# https://doi.org/10.1111/j.1365-246X.1995.tb01819.x
_han_wahr_lln_table = get_data_path(["data", "han-wahr-lln.txt"])
# PREM outputs from Gegout (2005)
# http://gemini.gsfc.nasa.gov/aplo/
_gegout_lln_table = get_data_path(["data", "gegout-lln.txt"])
# PREM outputs from Wang et al. (2012)
# https://doi.org/10.1016/j.cageo.2012.06.022
_wang_prem_lln_table = get_data_path(["data", "wang-prem-lln.txt"])
# default maximum degree and order in case of infinite
_default_max_degree = 100000
# dictionary for tables of load Love/Shida numbers
_lln_table = {}
_lln_table["han-wahr"] = _han_wahr_lln_table
_lln_table["gegout"] = _gegout_lln_table
_lln_table["wang-prem"] = _wang_prem_lln_table


# PURPOSE: read tables of load Love/Shida numbers
def load_love_numbers(
    filename: str | pathlib.Path = _han_wahr_lln_table,
    reference: str = "CE",
    **kwargs,
):
    """
    Read tables of Load Love/Shida numbers and applies isomorphic
    parameters :cite:p:`Blewitt:2003bz,Farrell:1972cm,Wahr:1998hy`

    Parameters
    ----------
    filename: str | pathlib.Path
        Elastic load Love/Shida numbers file
    reference: str, default 'CE'
        Reference frame of degree 1 Love/Shida numbers

            - ``'CF'``: Center of Surface Figure
            - ``'CL'``: Center of Surface Lateral Figure
            - ``'CH'``: Center of Surface Height Figure
            - ``'CM'``: Center of Mass of Earth System
            - ``'CE'``: Center of Mass of Solid Earth
    lmax: int, default from file
        Upper bound of spherical harmonic degrees

    Returns
    -------
    hl: np.ndarray
        Load Love number of vertical displacement
    kl: np.ndarray
        Load Love number of gravitational potential
    ll: np.ndarray
        Load Love (Shida) number of horizontal displacement
    """
    # read load Love/Shida number data file
    columns = ("l", "hl", "kl", "ll")
    formats = ("i4", "f8", "f8", "f8")
    love = np.loadtxt(
        filename,
        dtype={"names": columns, "formats": formats},
        converters={0: _infconv},
    )

    # calculate isomorphic parameters for different reference frames
    # following Blewitt (2003) and Farrell (1972)
    if reference.upper() == "CF":
        # Center of Surface Figure
        alpha = (love["hl"][1] + 2.0 * love["ll"][1]) / 3.0
    elif reference.upper() == "CL":
        # Center of Surface Lateral Figure
        alpha = love["ll"][1].copy()
    elif reference.upper() == "CH":
        # Center of Surface Height Figure
        alpha = love["hl"][1].copy()
    elif reference.upper() == "CM":
        # Center of Mass of Earth System
        alpha = 1.0
    elif reference.upper() == "CE":
        # Center of Mass of Solid Earth
        alpha = 0.0
    else:
        raise ValueError(f"Invalid Reference Frame {reference}")

    # maximum degree of the load Love/Shida numbers
    lmax = kwargs.get("lmax", np.max(love["l"]))
    # array of spherical harmonic degrees up to lmax
    l = np.arange(lmax + 1)
    # dictionary of output Love/Shida numbers
    # interpolated over all degrees up to lmax
    tmp = {}
    for col in ("hl", "kl", "ll"):
        # apply isomorphic parameters
        love[col][1] -= alpha
        # linearly interpolate missing values
        tmp[col] = np.interp(l, love["l"], love[col])

    # return the load love numbers
    return (tmp["hl"], tmp["kl"], tmp["ll"])


def adjust_load_love_numbers(omega: np.ndarray, **kwargs):
    """
    Compute the adjustments for load Love/Shida numbers
    to account for free core nutation and other resonances
    :cite:p:`Mathews:1997js,Mathews:2002cr,Petit:2010tp`

    Parameters
    ----------
    omega: np.ndarray
        Angular frequency (radians per second)

    Returns
    -------
    dh21: float
        Load Love number adjustments for vertical displacement
    dk21: float
        Load Love number adjustments for gravitational potential
    dl21: float
        Load Love/Shida number adjustments for horizontal displacement
    """
    # number of sidereal days per solar day
    sidereal_ratio = 1.002737909
    # number of seconds in a sidereal day (approximately 86164.1)
    sidereal_day = 86400.0 / sidereal_ratio
    # frequency in cycles per sidereal day
    f = omega * sidereal_day / (2.0 * np.pi)
    # load Love numbers adjustments for different frequency bands
    if (omega > 1e-4) or (omega < 2e-5):
        # semi-diurnal and long-period bands
        dh21 = 0.0
        dk21 = 0.0
        dl21 = 0.0
    else:
        # in-phase adjustments for the diurnal band
        # following IERS conventions and Mathews et al. (2002)
        # values from equation 6.10 of IERS conventions 2010
        # and from Mathews et al. (2002)
        sigma = np.zeros((3))
        # Chandler wobble
        sigma[0] = -0.0026010
        # retrograde free core nutation
        sigma[1] = 1.0023181
        # prograde free core nutation
        sigma[2] = 0.999026
        # frequency dependence of load Love number h2,1 (vertical)
        # table 7.1 (IERS conventions 2010)
        H21 = np.zeros((3))
        H21[0] = 1.6583e-3
        H21[1] = 2.8018e-4
        H21[2] = 5.5852e-7
        # frequency dependence of load Love number k2,1 (potential)
        # table 6.4 (IERS conventions 2010)
        K21 = np.zeros((3))
        K21[0] = 8.1874e-4
        K21[1] = 1.4116e-4
        K21[2] = 3.4618e-7
        # frequency dependence of load Love number l2,1 (horizontal)
        # table 7.1 (IERS conventions 2010)
        L21 = np.zeros((3))
        L21[0] = 2.3232e-4
        L21[1] = -8.4659e-6
        L21[2] = 1.0724e-8
        # estimate the load Love numbers for diurnal tides
        # equation 6.9 (IERS conventions 2010)
        dh21 = np.sum(H21 / (f - sigma))
        dk21 = np.sum(K21 / (f - sigma))
        dl21 = np.sum(L21 / (f - sigma))

    # return the load Love numbers adjustments
    return (dh21, dk21, dl21)


def greens_function(
    alpha: float | np.ndarray,
    lmax: int | None = None,
    lln: str = "han-wahr",
    reference: str = "CE",
    hinf: float | None = None,
):
    """
    Computes Green's functions for vertical displacement following
    :cite:t:`Farrell:1972cm` using the Legendre recurrence relation
    from :cite:t:`HofmannWellenhof:2006hy`

    Parameters
    ----------
    alpha: float or np.ndarray,
        Angular distance (radians)
    lmax: int or None, default None
        Upper bound of spherical harmonic degrees
    lln: str, default 'han-wahr'
        Name of the Load Love number dataset to use

            - ``'han-wahr'``: :cite:t:`Han:1995go`
            - ``'gegout'``: :cite:t:`Gegout:2010gc`
            - ``'wang-prem'``: :cite:t:`Wang:2012gc`
    reference: str, default 'CE'
        Reference frame of degree 1 load Love numbers

            - ``'CF'``: Center of Surface Figure
            - ``'CL'``: Center of Surface Lateral Figure
            - ``'CH'``: Center of Surface Height Figure
            - ``'CM'``: Center of Mass of Earth System
            - ``'CE'``: Center of Mass of Solid Earth
    hinf: float | None, default None
        Asymptotic value of the Love number for vertical displacement

    Returns
    -------
    GF: numpy.ndarray
        Green's function for vertical displacement
    """
    # read load Love numbers from table
    hl, kl, ll = load_love_numbers(_lln_table[lln], reference=reference)
    # maximum degree for spherical harmonics
    if lmax is None:
        lmax = len(hl) - 1
    # Love number asymtoptic value
    if hinf is None:
        hinf = np.copy(hl[-1])
    # calculate the cosine of the angular distance
    x = np.cos(alpha)
    # allocate for Legendre polynomial series
    # P0 and P1
    Pm2 = np.ones_like(x)
    Pm1 = np.copy(x)
    P = np.zeros_like(x)
    # calculate Legendre asymptotic value using Kummer's transformation
    Pinf = np.divide(1.0, 2.0 * np.sin(alpha / 2.0))
    # calculate Green's function at the asymptote
    GF = hinf * Pinf
    # add first two degrees to the output Green's function
    GF += (hl[0] - hinf) * Pm2
    GF += (hl[1] - hinf) * Pm1
    # for each spherical harmonic degree
    for l in range(2, lmax + 1):
        # calculate unnormalized Legendre polynomial for degree l
        P[:] = np.divide(2 * l - 1, l) * x * Pm1 - np.divide(l - 1, l) * Pm2
        # add to output Green's function
        GF += (hl[l] - hinf) * P
        # update Legendre polynomial recursion series
        Pm2[:] = np.copy(Pm1)
        Pm1[:] = np.copy(P)
    # return the Green's function
    return GF


def _melchior_table_52(model: str):
    """
    Body tide Love numbers for an Earth model provided in
    Table 5.2 of :cite:t:`Melchior:1983wd`

    Parameters
    ----------
    model: str
        Earth model

        - ``'Longman'``: :cite:t:`Longman:1959hw`
        - ``'Kaula'``: :cite:t:`Kaula:1966un`
        - ``'Takeuchi'``: :cite:t:`Takeuchi:1950hi`
    """
    # table 5.2 from Melchior (1983)
    # column 1: Spherical harmonic degree
    # column 2: Love number of vertical displacement
    # column 3: Love number of gravitational potential
    # column 4: Shida number of horizontal displacement
    if model == "Longman":
        # values from Longman (1959)
        table_52 = np.array(
            [
                [2, 0.612, 0.302, 0.083],
                [3, 0.290, 0.093, 0.014],
                [4, 0.175, 0.042, 0.010],
                [5, 0.129, 0.025, 0.008],
                [6, 0.107, 0.017, 0.007],
                [7, 0.095, 0.013, 0.005],
                [8, 0.087, 0.010, 0.004],
                [9, 0.081, 0.008, 0.004],
                [10, 0.076, 0.007, 0.003],
                [11, 0.072, 0.006, 0.002],
                [12, 0.069, 0.005, 0.002],
                [13, 0.066, 0.005, 0.002],
                [14, 0.064, 0.004, 0.001],
                [15, 0.062, 0.004, 0.001],
                [16, 0.060, 0.003, 0.001],
                [17, 0.058, 0.003, 0.001],
                [18, 0.056, 0.003, 0.001],
                [19, 0.055, 0.003, 0.001],
                [20, 0.053, 0.002, 0.001],
                [21, 0.052, 0.002, 0.001],
                [22, 0.051, 0.002, 0.000],
                [23, 0.050, 0.002, 0.000],
                [24, 0.048, 0.002, 0.000],
                [25, 0.047, 0.002, 0.000],
            ]
        )
    elif model == "Kaula":
        # values from Kaula (1966)
        table_52 = np.array(
            [
                [2, 0.624, 0.317, 0.085],
                [3, 0.293, 0.095, 0.014],
                [4, 0.177, 0.042, 0.010],
                [5, 0.130, 0.025, 0.008],
                [6, 0.107, 0.017, 0.007],
                [7, 0.095, 0.013, 0.005],
                [8, 0.087, 0.010, 0.004],
                [9, 0.081, 0.008, 0.004],
                [10, 0.076, 0.007, 0.003],
                [11, 0.072, 0.006, 0.002],
                [12, 0.069, 0.005, 0.002],
                [13, 0.066, 0.005, 0.002],
                [14, 0.064, 0.004, 0.001],
                [15, 0.062, 0.004, 0.001],
                [16, 0.060, 0.003, 0.001],
            ]
        )
    elif model == "Takeuchi":
        # values from Takeuchi (1950)
        table_52 = np.array(
            [
                [2, 0.592, 0.280, 0.076],
                [3, 0.274, 0.083, 0.010],
                [4, 0.161, 0.035, 0.007],
                [5, 0.116, 0.020, 0.006],
                [6, 0.094, 0.013, 0.005],
                [7, 0.081, 0.009, 0.004],
                [8, 0.073, 0.007, 0.003],
                [9, 0.068, 0.006, 0.0025],
                [10, 0.063, 0.005, 0.002],
                [11, 0.059, 0.004, 0.002],
                [12, 0.055, 0.003, 0.002],
                [13, 0.053, 0.003, 0.0015],
                [14, 0.051, 0.003, 0.001],
                [15, 0.0495, 0.0025, 0.001],
                [16, 0.048, 0.002, 0.001],
            ]
        )
    else:
        raise ValueError(f"Unknown Earth model: {model}")
    # return the table of love numbers
    return table_52


def _infconv(fld, lmax: int = _default_max_degree) -> int:
    """
    Convert degree to integer and check if degree is infinite

    Parameters
    ----------
    fld: str
        Degree field from load Love/Shida numbers file
    lmax: int, default 100000
        Truncation in case of infinite degree
    """
    if re.search(f"inf", fld, re.IGNORECASE):
        return lmax
    else:
        return int(fld)
