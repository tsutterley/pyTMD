#!/usr/bin/env python
u"""
response.py
Written by Tyler Sutterley (11/2024)
Routines for estimating tidal constituents using the response method

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://scipy.org

PROGRAM DEPENDENCIES:
    arguments.py: loads nodal corrections for tidal constituents
    astro.py: computes the basic astronomical mean longitudes
    math.py: special mathematical functions
    spatial.py: utilities for spatial operations

UPDATE HISTORY:
    Written 11/2024
"""

from __future__ import annotations

import numpy as np
import scipy.linalg
import scipy.optimize
import pyTMD.arguments
import pyTMD.astro
import pyTMD.math
import pyTMD.spatial

__all__ = [
    'response',
    '_gravitational',
    '_radiational',
    '_kappa'
]

def response(
        t: np.ndarray,
        lon: np.ndarray,
        lat: np.ndarray,
        ht: np.ndarray,
        constituents: str | list | np.ndarray,
        deltat: float | np.ndarray = 0.0,
        corrections: str = 'OTIS',
        **kwargs: dict
    ):
    """
    Estimate tidal constituents using the response method [1]_ [2]_ [3]_ [4]_

    Parameters
    ----------
    t: np.ndarray
        Time (days relative to January 1, 1992)
    lon: np.ndarray
        longitude coordinates
    lat: np.ndarray
        latitude coordinates
    ht: np.ndarray
        elevation time series (meters)
    constituents: str, list or np.ndarray
        tidal constituent ID(s)
    deltat: float or np.ndarray, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS/ATLAS or GOT/FES models

    Returns
    -------
    amp: np.ndarray
        amplitude of each harmonic constant (meters)
    phase: np.ndarray
        phase of each harmonic constant (degrees)

    References
    ----------
    .. [1] W. H. Munk, D. E. Cartwright, and E. C. Bullard, "Tidal
        spectroscopy and prediction," *Philosophical Transactions of the
        Royal Society of London. Series A, Mathematical and Physical
        Sciences*, 259(1105), 533--581, (1966).
        `doi: 10.1098/rsta.1966.0024 <https://doi.org/10.1098/rsta.1966.0024>`_
    .. [2] G. W. Groves and R. W. Reynolds, "An orthogonalized convolution
        method of tide prediction," *Journal of Geophysical Research*,
        80(30), 4131--4138, (1975).
        `doi: 10.1029/jc080i030p04131 <https://doi.org/10.1029/jc080i030p04131>`_
    .. [3] B. D. Zetler and W. T. Munk, "The optimum wiggliness of tidal
        admittances," *Journal of Marine Research*, 33(S), (1975).
    .. [4] D. E. Cartwright and R. D. Ray, "Oceanic tides from Geosat altimetry,"
        *Journal of Geophysical Research: Oceans*, 95(C3), 3069--3090, (1990).
        `doi: 10.1029/JC095iC03p03069 <https://doi.org/10.1029/JC095iC03p03069>`_
    """
    # default keyword arguments
    kwargs.setdefault('ephemerides', 'approximate')
    # check if input constituents is a string
    if isinstance(constituents, str):
        constituents = [constituents]
    # number of time values and constituents
    nt = len(np.atleast_1d(t))
    nc = len(constituents)
    # verify dimensions
    singular_values = (np.ndim(t) == 0)
    t = np.atleast_1d(t).flatten()
    lon = np.atleast_1d(lon).flatten()
    lat = np.atleast_1d(lat).flatten()
    # assert dimensions
    assert len(t) == len(lat), 'coordinates must have the same dimensions'
    assert len(lon) == len(lat), 'coordinates must have the same dimensions'
    # check that the number of time values matches the number of height values
    if (nt != len(np.atleast_1d(ht))):
        raise ValueError('Dimension mismatch between input variables')
    # degrees to radians
    dtr = np.pi/180.0
    # solar ephemerides (convert time to Modified Julian Days)
    SXYZ = pyTMD.astro.solar_ecef(t + 48622.0, **kwargs)
    # lunar ephemerides (convert time to Modified Julian Days)
    LXYZ = pyTMD.astro.lunar_ecef(t + 48622.0, **kwargs)

def _gravitational(
        t: float | np.ndarray,
        lon: np.ndarray,
        lat: np.ndarray,
        ht: np.ndarray,
        SXYZ: np.ndarray,
        LXYZ: np.ndarray,
        l: int,
        **kwargs: dict
    ):
    """
    Estimate gravitational tides using the response method [1]_

    Parameters
    ----------
    t: np.ndarray
        Time (days relative to January 1, 1992)
    lon: np.ndarray
        longitude coordinates
    lat: np.ndarray
        latitude coordinates
    ht: np.ndarray
        elevation time series (meters)
    SXYZ: np.ndarray
        Earth-centered Earth-fixed coordinates of the sun (meters)
    LXYZ: np.ndarray
        Earth-centered Earth-fixed coordinates of the moon (meters)
    l: int
        spherical harmonic degree

    Returns
    -------
    amp: np.ndarray
        amplitude of each harmonic constant (meters)
    phase: np.ndarray
        phase of each harmonic constant (degrees)

    References
    ----------
    .. [1] W. H. Munk, D. E. Cartwright, and E. C. Bullard, "Tidal
        spectroscopy and prediction," *Philosophical Transactions of the
        Royal Society of London. Series A, Mathematical and Physical
        Sciences*, 259(1105), 533--581, (1966).
        `doi: 10.1098/rsta.1966.0024 <https://doi.org/10.1098/rsta.1966.0024>`_
    """
    # default keyword arguments
    kwargs.setdefault('a_axis', 6378137.0)
    kwargs.setdefault('flat', 1.0/298.257223563)
    kwargs.setdefault('omega', 7.2921151467e-5)
    kwargs.setdefault('AU', 1.495978707e11)
    kwargs.setdefault('LD', 3.84399e8)
    # mass ratios between earth and sun/moon
    kwargs.setdefault('mass_ratio_solar', 332946.0482)
    kwargs.setdefault('mass_ratio_lunar', 0.0123000371)
    # degrees to radians
    dtr = np.pi/180.0
    # average radius of the earth
    rad_e = kwargs['a_axis']*(1.0 - kwargs['flat'])**(1.0/3.0)
    # solar and lunar radii from ephemerides
    solar_radius = np.sqrt(SXYZ[:,0]**2 + SXYZ[:,1]**2 + SXYZ[:,2]**2)
    lunar_radius = np.sqrt(LXYZ[:,0]**2 + LXYZ[:,1]**2 + LXYZ[:,2]**2)
    # convert solar and lunar ephemerides from ECEF to zenith angle
    solar_zenith = pyTMD.spatial.to_zenith(SXYZ[:,0], SXYZ[:,1], SXYZ[:,2],
        lon, lat, **kwargs)
    lunar_zenith = pyTMD.spatial.to_zenith(LXYZ[:,0], LXYZ[:,1], LXYZ[:,2],
        lon, lat, **kwargs)
    # associated Legendre functions of zenith angle for degree l
    solar_legendre = pyTMD.math.legendre(l, np.cos(dtr*solar_zenith))
    lunar_legendre = pyTMD.math.legendre(l, np.cos(dtr*lunar_zenith))
    # k values from Munk and Cartwright (1966)
    solar_k = rad_e*kwargs['mass_ratio_solar']*(rad_e/solar_radius)**(l+1)
    lunar_k = rad_e*kwargs['mass_ratio_lunar']*(rad_e/lunar_radius)**(l+1)
    # gravitational potential of the sun and moon
    solar_grav = solar_k*solar_legendre*(kwargs['AU']/solar_radius)**(l+1)
    lunar_grav = lunar_k*lunar_legendre*(kwargs['LD']/lunar_radius)**(l+1)

def _radiational(
        t: float | np.ndarray,
        lon: np.ndarray,
        lat: np.ndarray,
        ht: np.ndarray,
        SXYZ: np.ndarray,
        l: int,
        **kwargs: dict
    ):
    """
    Estimate radiational tides using the response method [1]_

    Parameters
    ----------
    t: np.ndarray
        Time (days relative to January 1, 1992)
    lon: np.ndarray
        longitude coordinates
    lat: np.ndarray
        latitude coordinates
    ht: np.ndarray
        elevation time series (meters)
    SXYZ: np.ndarray
        Earth-centered Earth-fixed coordinates of the sun (meters)
    l: int
        spherical harmonic degree

    Returns
    -------
    amp: np.ndarray
        amplitude of each harmonic constant (meters)
    phase: np.ndarray
        phase of each harmonic constant (degrees)

    References
    ----------
    .. [1] W. H. Munk, D. E. Cartwright, and E. C. Bullard, "Tidal
        spectroscopy and prediction," *Philosophical Transactions of the
        Royal Society of London. Series A, Mathematical and Physical
        Sciences*, 259(1105), 533--581, (1966).
        `doi: 10.1098/rsta.1966.0024 <https://doi.org/10.1098/rsta.1966.0024>`_
    """
    # default keyword arguments
    kwargs.setdefault('a_axis', 6378137.0)
    kwargs.setdefault('solar_constant', 1380.0)
    kwargs.setdefault('AU', 1.495978707e11)
    # degrees to radians
    dtr = np.pi/180.0
    # 1/AU in terms of Earth equatorial radius (~23455)
    # defined as parallax in Munk and Cartwright (1966)
    xi = kwargs['a_axis']/kwargs['AU']
    # kappa value (Munk and Cartwright, 1966)
    kappa = _kappa(l, xi)
    # solar radius from ephemerides
    solar_radius = np.sqrt(SXYZ[:,0]**2 + SXYZ[:,1]**2 + SXYZ[:,2]**2)
    # convert solar ephemerides from ECEF to zenith angle
    solar_zenith = pyTMD.spatial.to_zenith(SXYZ[:,0], SXYZ[:,1], SXYZ[:,2],
        lon, lat, **kwargs)
    # create radiation function
    S = np.copy(kwargs['solar_constant'])
    R = S*(kwargs['AU']/solar_radius)*kappa*np.cos(dtr*solar_zenith)**l
    # radiation is only valid during the day
    R = np.where(solar_zenith < 0.5*np.pi, R, 0.0)

def _kappa(l: int, xi: float) -> float:
    """
    Kappa values as a function of degree [1]_

    Parameters
    ----------
    l: int
        spherical harmonic degree
    xi: float
        parallax

    Returns
    -------
    kappa: float
        kappa value for the given degree

    References
    ----------
    .. [1] W. H. Munk, D. E. Cartwright, and E. C. Bullard, "Tidal
        spectroscopy and prediction," *Philosophical Transactions of the
        Royal Society of London. Series A, Mathematical and Physical
        Sciences*, 259(1105), 533--581, (1966).
        `doi: 10.1098/rsta.1966.0024 <https://doi.org/10.1098/rsta.1966.0024>`_
    """
    kappa = np.zeros((4))
    kappa[0] = 1.0/4.0 + xi/6.0
    kappa[1] = 1.0/2.0 + 3.0*xi/8.0
    kappa[2] = 5.0/16.0 + xi/3.0
    return kappa[l]
