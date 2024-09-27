#!/usr/bin/env python
u"""
crs.py
Written by Tyler Sutterley (09/2024)
Coordinates Reference System (CRS) routines

CALLING SEQUENCE:
    x, y = pyTMD.crs().convert(lon, lat, PROJ, 'F')
    lon, lat = pyTMD.crs().convert(x, y, PROJ, 'B')

INPUTS:
    i1: longitude ('F') or projection easting x ('B')
    i2: latitude ('F') or projection northing y ('B')
    PROJ: spatial reference system for coordinate transformations
    BF: backwards ('B') or forward ('F') translations

OPTIONS:
    EPSG: spatial reference system code for input (F) and output (B) coordinates

OUTPUTS:
    o1: projection easting x ('F') or longitude ('B')
    o2: projection northing y ('F') or latitude ('B')

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
        https://pyproj4.github.io/pyproj/

UPDATE HISTORY:
    Updated 09/2024: added function for idealized Arctic Azimuthal projection
        complete refactor to use JSON dictionary format for model projections
    Updated 07/2024: added function to get the CRS transform
    Updated 05/2024: make subscriptable and allow item assignment
    Updated 04/2024: use wrapper to importlib for optional dependencies
    Updated 02/2024: changed class name for ellipsoid parameters to datum
    Updated 12/2023: converted conversion functions to class
    Updated 03/2023: add basic variable typing to function inputs
        renamed coordinate reference system conversion functions
    Updated 02/2023: use named exception before passing to custom
    Updated 11/2022: place some imports within try/except statements
        use f-strings for formatting verbose or ascii output
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 09/2021: added function for using custom projections
    Updated 06/2021: added 3413 for new 1km Greenland model from ESR
    Updated 08/2020: using conversion protocols following pyproj-2 updates
        https://pyproj4.github.io/pyproj/stable/gotchas.html
    Updated 07/2020: added function docstrings. changed function name
    Updated 03/2020: remove commented coordinate conversion functions
    Updated 11/2019: using pyproj for coordinate conversions
    Written 09/2017
"""

from __future__ import annotations

import logging
import numpy as np
from pyTMD.utilities import import_dependency
# attempt imports
pyproj = import_dependency('pyproj')

__all__ = [
    'crs',
    'datum'
]

class crs:
    """Coordinate Reference System transformations for tide models

    Attributes
    ----------
    name: str
        Projection name
    transformer: obj
        ``pyproj`` transformer for changing coordinate reference system
    """
    def __init__(self):
        self.name = None
        self.transformer = None
        self._direction = None

    def convert(self,
            i1: np.ndarray,
            i2: np.ndarray,
            PROJ: str | dict,
            BF: str,
            EPSG: int | str = 4326
        ):
        """
        Converts points to and from Coordinates Reference Systems (CRS)

        Parameters
        ----------
        i1: np.ndarray
            Input x-coordinates
        i2: np.ndarray
            Input y-coordinates
        PROJ: str or dict
            Spatial reference system for coordinate transformations
        BF: str
            Direction of transformation

                - ``'B'``: backwards
                - ``'F'``: forwards
        EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
            input (``'F'``) or output (``'B'``) coordinate system

        Returns
        -------
        o1: np.ndarray
            Output transformed x-coordinates
        o2: np.ndarray
            Output transformed y-coordinates
        """
        # name of the projection
        self.name = PROJ
        # get the CRS and transform direction
        self.get(PROJ)
        self._direction = BF[0].upper()
        # run conversion program and return values
        return self.transform(i1, i2, EPSG=EPSG)

    # PURPOSE: try to get the projection information
    def get(self, PROJ: str | dict):
        """
        Tries to get the coordinate reference system

        Parameters
        ----------
        PROJ: str or dict
            Spatial reference system for coordinate transformations
        """
        # get the coordinate reference system
        try:
            self.crs = self.from_input(PROJ)
            self.name = self.crs.name
        except Exception as exc:
            pass
        else:
            return self
        # projection not found or available
        raise pyproj.exceptions.CRSError

    def transform(self,
            i1: np.ndarray,
            i2: np.ndarray,
            EPSG: int | str = 4326,
            **kwargs):
        """
        Performs Coordinates Reference System (CRS) transformations

        Parameters
        ----------
        i1: np.ndarray
            Input x-coordinates
        i2: np.ndarray
            Input y-coordinates
        EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
            input (``'F'``) or output (``'B'``) coordinate system
        kwargs: dict
            Keyword arguments for the transformation

        Returns
        -------
        o1: np.ndarray
            Output transformed x-coordinates
        o2: np.ndarray
            Output transformed y-coordinates
        """
        # set the direction of the transformation
        kwargs.setdefault('direction', self.direction)
        # get the coordinate reference system and transform
        source_crs = self.from_input(EPSG)
        self.transformer = pyproj.Transformer.from_crs(
            source_crs, self.crs, always_xy=True)
        # convert coordinate reference system
        o1, o2 = self.transformer.transform(i1, i2, **kwargs)
        # return the transformed coordinates
        return (o1, o2)

    # PURPOSE: try to get the projection information
    def from_input(self, PROJECTION: int | str | dict):
        """
        Attempt to retrieve the Coordinate Reference System

        Parameters
        ----------
        PROJECTION: int, str or dict
            Coordinate Reference System
        """
        # coordinate reference system dictoinary
        try:
            CRS = pyproj.CRS.from_user_input(PROJECTION)
        except (ValueError, pyproj.exceptions.CRSError):
            pass
        else:
            return CRS
        # EPSG projection code
        try:
            CRS = pyproj.CRS.from_epsg(int(PROJECTION))
        except (ValueError, pyproj.exceptions.CRSError):
            pass
        else:
            return CRS
        # coordinate reference system string
        try:
            CRS = pyproj.CRS.from_string(PROJECTION)
        except (ValueError, pyproj.exceptions.CRSError):
            pass
        else:
            return CRS
        # no projection can be made
        raise pyproj.exceptions.CRSError

    @property
    def direction(self):
        """
        ``pyproj`` direction of the coordinate transform
        """
        # convert from input coordinates to model coordinates
        if (self._direction is None) or (self._direction.upper() == 'F'):
            return pyproj.enums.TransformDirection.FORWARD
        # convert from model coordinates to coordinates
        elif (self._direction.upper() == 'B'):
            return pyproj.enums.TransformDirection.INVERSE

    @property
    def is_geographic(self):
        """
        Check if the coordinate reference system is geographic
        """
        return self.crs.is_geographic

    def __str__(self):
        """String representation of the ``crs`` object
        """
        properties = ['pyTMD.crs']
        properties.append(f"    name: {self.name}")
        properties.append(f"    direction: {self.direction.name}")
        return '\n'.join(properties)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

_ellipsoids = ['CLK66', 'GRS67', 'GRS80', 'WGS72', 'WGS84', 'ATS77',
    'NAD27', 'NAD83', 'INTER', 'KRASS', 'MAIRY', 'HGH80', 'TOPEX',
    'EGM96', 'IERS']
_units = ['MKS', 'CGS']

class datum:
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
    def __init__(self, **kwargs):
        # set default keyword arguments
        kwargs.setdefault('ellipsoid', 'WGS84')
        kwargs.setdefault('units', 'MKS')
        kwargs.setdefault('a_axis', None)
        kwargs.setdefault('flat', None)
        kwargs.setdefault('GM', None)
        kwargs.setdefault('omega', None)
        # set ellipsoid name and units
        self.units = kwargs['units'].upper()
        if ((kwargs['a_axis'] is not None) and (kwargs['flat'] is not None)):
            self.name = 'user_defined'
        else:
            self.name = kwargs['ellipsoid'].upper()
        # validate ellipsoid and units
        assert self.name in _ellipsoids + ['user_defined']
        assert self.units in _units

        # set parameters for ellipsoid
        if self.name in ('CLK66', 'NAD27'):
            # Clarke 1866
            self.a_axis = 6378206.4# [m] semimajor axis of the ellipsoid
            self.flat = 1.0/294.9786982# flattening of the ellipsoid

        elif self.name in ('GRS80', 'NAD83'):
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

        elif (self.name == 'user_defined'):
            # custom datum
            self.a_axis = np.float64(kwargs['a_axis'])
            self.flat = np.float64(kwargs['flat'])

        # set default parameters if not listed as part of ellipsoid
        # Geocentric Gravitational Constant
        if kwargs['GM'] is not None:
            # user defined Geocentric Gravitational Constant
            self.GM = np.float64(kwargs['GM'])
        elif self.name not in ('GRS80', 'GRS67', 'NAD83', 'TOPEX', 'EGM96'):
            # for ellipsoids not listing the Geocentric Gravitational Constant
            self.GM = 3.986004418e14# [m^3/s^2]

        # angular velocity of the Earth
        if kwargs['omega'] is not None:
            # user defined angular velocity of the Earth
            self.omega = np.float64(kwargs['omega'])
        elif self.name not in ('GRS67'):
            # for ellipsoids not listing the angular velocity of the Earth
            self.omega = 7292115e-11# [rad/s]

        # universal gravitational constant [N*m^2/kg^2]
        self.G = 6.67430e-11

        # standard gravitational acceleration [m/s^2]
        # (World Meteorological Organization)
        self.gamma = 9.80665

        # convert units to CGS
        if (self.units == 'CGS'):
            self.a_axis *= 100.0
            self.GM *= 1e6
            self.G *= 1000.0 # [dyn*cm^2/g^2]
            self.gamma *= 100.0

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
        r"""q\ :sub:`0` Parameter
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
        r"""Normalized C\ :sub:`20` harmonic
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

    def __str__(self):
        """String representation of the ``datum`` object
        """
        properties = ['pyTMD.datum']
        properties.append(f"    name: {self.name}")
        properties.append(f"    units: {self.units}")
        return '\n'.join(properties)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)
