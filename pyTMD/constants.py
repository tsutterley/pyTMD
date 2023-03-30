#!/usr/bin/env python
u"""
constants.py
Written by Tyler Sutterley (03/2023)

Gravitational and ellipsoidal parameters

CALLING SEQUENCE
    wgs84 = constants('WGS84', units='MKS')

INPUT:
    ellipsoid - reference ellipsoid name
        CLK66 = Clarke 1866
        GRS67 = Geodetic Reference System 1967 (IAU ellipsoid)
        GRS80 = Geodetic Reference System 1980
        WGS72 = World Geodetic System 1972
        WGS84 = World Geodetic System 1984
        ATS77 = Quasi-earth centred ellipsoid for ATS77
        NAD27 = North American Datum 1927
        NAD83 = North American Datum 1983
        INTER = International
        KRASS = Krassovsky (USSR)
        MAIRY = Modified Airy (Ireland 1965/1975)
        HGH80 = Hughes 1980 Ellipsoid used in some NSIDC data
        TOPEX = TOPEX/POSEIDON ellipsoid
        EGM96 = EGM 1996 gravity model
        IERS = IERS Numerical Standards (2010)

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

REFERENCE:
    Hofmann-Wellenhof and Moritz, "Physical Geodesy" (2006)
    IERS Numerical Standards
        https://iers-conventions.obspm.fr/content/tn36.pdf

UPDATE HISTORY:
    Updated 03/2023: add basic variable typing
        set ellipsoid name and output units as constants attributes
    Updated 01/2023: include main ellipsoid attributes in docstring
    Written 12/2022
"""
import numpy as np

_ellipsoids = ['CLK66', 'GRS67', 'GRS80', 'WGS72', 'WGS84', 'ATS77',
    'NAD27', 'NAD83', 'INTER', 'KRASS', 'MAIRY', 'HGH80', 'TOPEX',
    'EGM96', 'IERS']
_units = ['MKS', 'CGS']

class constants(object):
    """
    Class for gravitational and ellipsoidal parameters

    Parameters
    ----------
    ellipsoid: str, default 'WGS84'
        Reference ellipsoid name

            - ``'CLK66'``: Clarke 1866
            - ``'GRS67'``: Geodetic Reference System 1967
            - ``'GRS80'``: Geodetic Reference System 1980
            - ``'HGH80'``: Hughes 1980 Ellipsoid
            - ``'WGS72'``: World Geodetic System 1972
            - ``'WGS84'``: World Geodetic System 1984
            - ``'ATS77'``: Quasi-earth centred ellipsoid for ATS77
            - ``'NAD27'``: North American Datum 1927
            - ``'NAD83'``: North American Datum 1983
            - ``'INTER'``: International
            - ``'KRASS'``: Krassovsky (USSR)
            - ``'MAIRY'``: Modified Airy (Ireland 1965/1975)
            - ``'TOPEX'``: TOPEX/POSEIDON ellipsoid
            - ``'EGM96'``: EGM 1996 gravity model
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
    np.seterr(invalid='ignore')
    def __init__(self, ellipsoid: str = 'WGS84', units: str = 'MKS'):
        # set ellipsoid name and units
        self.name = ellipsoid.upper()
        self.units = units.upper()
        # validate ellipsoid and units
        assert self.name in _ellipsoids
        assert self.units in _units

        # set parameters for ellipsoid
        if self.name in ('CLK66','NAD27'):
            # Clarke 1866
            self.a_axis = 6378206.4# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/294.9786982# flattening of the ellipsoid

        elif self.name in ('GRS80','NAD83'):
            # Geodetic Reference System 1980
            # North American Datum 1983
            self.a_axis = 6378135.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.26# flattening of the ellipsoid
            self.GM = 3.986005e14# [m^3/s^2] Geocentric Gravitational Constant

        elif (self.name == 'GRS67'):
            # Geodetic Reference System 1967
            # International Astronomical Union (IAU ellipsoid)
            self.a_axis = 6378160.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.247167427# flattening of the ellipsoid
            self.GM = 3.98603e14# [m^3/s^2] Geocentric Gravitational Constant
            self.omega = 7292115.1467e-11# angular velocity of the Earth [rad/s]

        elif (self.name == 'WGS72'):
            # World Geodetic System 1972
            self.a_axis = 6378135.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.26# flattening of the ellipsoid

        elif (self.name == 'WGS84'):
            # World Geodetic System 1984
            self.a_axis = 6378137.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.257223563# flattening of the ellipsoid

        elif (self.name == 'ATS77'):
            # Quasi-earth centred ellipsoid for ATS77
            self.a_axis = 6378135.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.257# flattening of the ellipsoid

        elif (self.name == 'KRASS'):
            # Krassovsky (USSR)
            self.a_axis = 6378245.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.3# flattening of the ellipsoid

        elif (self.name == 'INTER'):
            # International
            self.a_axis = 6378388.0# [m] semimajor axis of the ellipsoid
            self.flat = 1/297.0# flattening of the ellipsoid

        elif (self.name == 'MAIRY'):
            # Modified Airy (Ireland 1965/1975)
            self.a_axis = 6377340.189# [m] semimajor axis of the ellipsoid
            self.flat = 1/299.3249646# flattening of the ellipsoid

        elif (self.name == 'HGH80'):
            # Hughes 1980 Ellipsoid used in some NSIDC data
            self.a_axis = 6378273.0# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.279411123064# flattening of the ellipsoid

        elif (self.name == 'TOPEX'):
            # TOPEX/POSEIDON ellipsoid
            self.a_axis = 6378136.3# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.257# flattening of the ellipsoid
            self.GM = 3.986004415e14# [m^3/s^2]

        elif (self.name == 'EGM96'):
            # EGM 1996 gravity model
            self.a_axis = 6378136.3# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.256415099# flattening of the ellipsoid
            self.GM = 3.986004415e14# [m^3/s^2]

        elif (self.name == 'IERS'):
            # IERS Numerical Standards
            self.a_axis = 6378136.6# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/298.25642# flattening of the ellipsoid

        # set default parameters if not listed as part of ellipsoid
        if self.name not in ('GRS80','GRS67','NAD83','TOPEX','EGM96'):
            # for ellipsoids not listing the Geocentric Gravitational Constant
            self.GM = 3.986004418e14# [m^3/s^2]

        if self.name not in ('GRS67'):
            # for ellipsoids not listing the angular velocity of the Earth
            self.omega = 7292115e-11# [rad/s]

        # universal gravitational constant [N*m^2/kg^2]
        self.G = 6.67430e-11

        # convert units to CGS
        if (self.units == 'CGS'):
            self.a_axis *= 100.0
            self.GM *= 1e6
            self.G *= 1000.0 # [dyn*cm^2/g^2]

    # mean radius of the Earth having the same volume
    # (4pi/3)R^3 = (4pi/3)(a^2)b = (4pi/3)(a^3)(1 - f)
    @property
    def rad_e(self) -> float:
        """Average radius of the Earth with same volume as ellipsoid
        """
        return self.a_axis*(1.0 - self.flat)**(1.0/3.0)

    # semiminor axis of the ellipsoid
    @property
    def b_axis(self) -> float:
        """Semi-minor axis of the ellipsoid
        """
        return (1.0 - self.flat)*self.a_axis

    # Ratio between ellipsoidal axes
    @property
    def ratio(self) -> float:
        """Ratio between ellipsoidal axes
        """
        return (1.0 - self.flat)

    # Polar radius of curvature
    @property
    def rad_p(self) -> float:
        """Polar radius of curvature
        """
        return self.a_axis/(1.0 - self.flat)

    # Linear eccentricity
    @property
    def ecc(self) -> float:
        """Linear eccentricity
        """
        return np.sqrt((2.0*self.flat - self.flat**2)*self.a_axis**2)

    # first numerical eccentricity
    @property
    def ecc1(self) -> float:
        """First numerical eccentricity
        """
        return self.ecc/self.a_axis

    # second numerical eccentricity
    @property
    def ecc2(self) -> float:
        """Second numerical eccentricity
        """
        return self.ecc/self.b_axis

    # m parameter [omega^2*a^2*b/(GM)]
    # p. 70, Eqn.(2-137)
    @property
    def m(self) -> float:
        """m Parameter
        """
        return self.omega**2*((1 - self.flat)*self.a_axis**3)/self.GM

    # flattening f2 component
    # p. 80, Eqn.(2-200)
    @property
    def f2(self) -> float:
        """f2 component
        """
        return -self.flat + (5.0/2.0)*self.m + (1.0/2.0)*self.flat**2.0 - \
            (26.0/7.0)*self.flat*self.m + (15.0/4.0)*self.m**2.0

    # flattening f4 component
    # p. 80, Eqn.(2-200)
    @property
    def f4(self) -> float:
        """f4 component
        """
        return -(1.0/2.0)*self.flat**2.0 + (5.0/2.0)*self.flat*self.m

    # q
    # p. 67, Eqn.(2-113)
    @property
    def q(self) -> float:
        """q Parameter
        """
        return 0.5*((1.0 + 3.0/(self.ecc2**2))*np.arctan(self.ecc2)-3.0/self.ecc2)

    # q_0
    # p. 67, Eqn.(2-113)
    @property
    def q0(self) -> float:
        """q\ :sub:`0` Parameter
        """
        return 3*(1.0 + 1.0/(self.ecc2**2)) * \
            (1.0 -1.0/self.ecc2*np.arctan(self.ecc2)) - 1.0

    # J_2 p. 75 Eqn.(2-167), p. 76 Eqn.(2-172)
    @property
    def J2(self) -> float:
        """Oblateness coefficient
        """
        return (self.ecc1**2)*(1.0 - 2.0*self.m*self.ecc2/(15.0*self.q))/3.0

    # Normalized C20 harmonic
    # p. 60, Eqn.(2-80)
    @property
    def C20(self) -> float:
        """Normalized C\ :sub:`20` harmonic
        """
        return -self.J2/np.sqrt(5.0)

    # Normal gravity at the equator
    # p. 79, Eqn.(2-286)
    @property
    def gamma_a(self) -> float:
        """Normal gravity at the equator
        """
        return (self.GM/(self.a_axis*self.b_axis)) * \
            (1.0 - (3.0/2.0)*self.m - (3.0/14.0)*self.ecc2**2.0*self.m)

    # Normal gravity at the pole
    # p. 79, Eqn.(2-286)
    @property
    def gamma_b(self) -> float:
        """Normal gravity at the pole
        """
        return (self.GM/(self.a_axis**2)) * \
            (1.0 + self.m + (3.0/7.0)*self.ecc2**2.0*self.m)

    # Normal gravity at location
    # p. 80, Eqn.(2-199)
    def gamma_0(self, theta) -> float:
        """Normal gravity at colatitudes

        Parameters
        ----------
        theta: float
            Colatitudes in radians
        """
        return self.gamma_a*(1.0 + self.f2*np.cos(theta)**2.0 + self.f4*np.cos(theta)**4.0)

    # Normal gravity at location
    # p. 82, Eqn.(2-215)
    def gamma_h(self, theta, height) -> float:
        """Normal gravity at colatitudes and heights

        Parameters
        ----------
        theta: float
            Colatitudes in radians
        height: float
            Height above ellipsoid
        """
        return self.gamma_0(theta) * \
            (1.0 - (2.0/self.a_axis) * \
            (1.0 + self.flat + self.m - 2.0*self.flat*np.cos(theta)**2.0)*height +
            (3.0/self.a_axis**2.0)*height**2.0)

    # ratio between gravity at pole versus gravity at equator
    @property
    def dk(self) -> float:
        """Ratio between gravity at pole versus gravity at equator
        """
        return self.b_axis*self.gamma_b/(self.a_axis*self.gamma_b) - 1.0

    # Normal potential at the ellipsoid
    # p. 68, Eqn.(2-123)
    @property
    def U0(self) -> float:
        """Normal potential at the ellipsoid
        """
        return self.GM/self.ecc*np.arctan(self.ecc2) + \
            (1.0/3.0)*self.omega**2*self.a_axis**2

    # Surface area of the reference ellipsoid
    @property
    def area(self) -> float:
        """Surface area of the ellipsoid
        """
        return np.pi*self.a_axis**2.0 * \
            (2.0 + ((1.0 - self.ecc1**2)/self.ecc1) *
                np.log((1.0 + self.ecc1)/(1.0 - self.ecc1)))

    # Volume of the reference ellipsoid
    @property
    def volume(self) -> float:
        """Volume of the ellipsoid
        """
        return (4.0*np.pi/3.0)*(self.a_axis**3.0)*(1.0 - self.ecc1**2.0)**0.5

    # Average density
    @property
    def rho_e(self) -> float:
        """Average density
        """
        return self.GM/(self.G*self.volume)
