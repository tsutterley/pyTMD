#!/usr/bin/env python
u"""
solve.py
Written by Tyler Sutterley (01/2024)
Routines for estimating the harmonic constants for ocean tides

REFERENCES:
    G. D. Egbert and S. Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://scipy.org

PROGRAM DEPENDENCIES:
    arguments.py: loads nodal corrections for tidal constituents
    astro.py: computes the basic astronomical mean longitudes
    constants.py: calculate reference parameters for common ellipsoids

UPDATE HISTORY:
    Updated 01/2024: add functions for tide generating forces and potentials
    Written 12/2023
"""

from __future__ import annotations

import numpy as np
import scipy.linalg
import pyTMD.arguments
import pyTMD.constants

# get WGS84 ellipsoid parameters in meters-kg-seconds
_wgs84 = pyTMD.constants.constants(ellipsoid='WGS84', units='MKS')

def constants(t: float | np.ndarray,
        ht: np.ndarray,
        constituents: str | list | np.ndarray,
        deltat: float | np.ndarray = 0.0,
        corrections: str = 'OTIS',
        solver: str = 'lstsq'
    ):
    """
    Estimate the harmonic constants for an elevation time series [1]_

    Parameters
    ----------
    t: float or np.ndarray
        days relative to 1992-01-01T00:00:00
    ht: np.ndarray
        elevation time series (meters)
    constituents: str, list or np.ndarray
        tidal constituent ID(s)
    deltat: float or np.ndarray, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS/ATLAS or GOT/FES models
    solver: str, default 'lstsq'
        least squares solver to use

        - ``'lstsq'``: least squares solution
        - ``'gelsy'``: complete orthogonal factorization
        - ``'gelss'``: singular value decomposition (SVD)
        - ``'gelsd'``: SVD with divide and conquer method

    Returns
    -------
    amp: np.ndarray
        amplitude of each harmonic constant (meters)
    phase: np.ndarray
        phase of each harmonic constant (degrees)

    References
    ----------
    .. [1] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of
        Barotropic Ocean Tides," *Journal of Atmospheric and Oceanic
        Technology*, 19(2), 183--204, (2002).
        `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

    .. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
    """
    # check if input constituents is a string
    if isinstance(constituents, str):
        constituents = [constituents]
    # check that there are enough values for a time series fit
    nt = len(np.atleast_1d(t))
    nc = len(constituents)
    if (nt <= 2*nc):
        raise ValueError('Not enough values for fit')
    # check that the number of time values matches the number of height values
    if (nt != len(np.atleast_1d(ht))):
        raise ValueError('Dimension mismatch between input variables')

    # load the nodal corrections
    # convert time to Modified Julian Days (MJD)
    pu, pf, G = pyTMD.arguments.arguments(t + 48622.0, constituents,
        deltat=deltat, corrections=corrections)

    # create design matrix
    M = []
    # add constant term for mean
    M.append(np.ones_like(t))
    # add constituent terms
    for k,c in enumerate(constituents):
        if corrections in ('OTIS', 'ATLAS', 'TMD3', 'netcdf'):
            amp, ph, omega, alpha, species = \
                pyTMD.arguments._constituent_parameters(c)
            th = omega*t*86400.0 + ph + pu[:,k]
        elif corrections in ('GOT', 'FES'):
            th = G[:,k]*np.pi/180.0 + pu[:,k]
        # add constituent to design matrix
        M.append(pf[:,k]*np.cos(th))
        M.append(-pf[:,k]*np.sin(th))
    # take the transpose of the design matrix
    M = np.transpose(M)

    # use a least-squares fit to solve for parameters
    if (solver == 'lstsq'):
        p, res, rnk, s = np.linalg.lstsq(M, ht, rcond=-1)
    elif solver in ('gelsd', 'gelsy', 'gelss'):
        p, res, rnk, s = scipy.linalg.lstsq(M, ht,
            lapack_driver=solver)

    # calculate amplitude and phase for each constituent
    amp = np.zeros((nc))
    ph = np.zeros((nc))
    # skip over the first indice in the fit (constant term)
    for k,c in enumerate(constituents):
        amp[k] = np.abs(1j*p[2*k+2] + p[2*k+1])
        ph[k] = np.arctan2(-p[2*k+2], p[2*k+1])
    # convert phase to degrees
    phase = ph*180.0/np.pi
    phase[phase < 0] += 360.0

    # return the amplitude and phase
    return (amp, phase)

# PURPOSE: calculate the astronomical tide generating force
def _generating_force(
        lon: np.ndarray,
        lat: np.ndarray,
        c: str,
    ):
    """
    Computes the astronomical tide generating force (the horizontal gradient
    of the tide generating potential) for a tidal constituent [1]_

    Parameters
    ----------
    lon: np.ndarray
        longitudes of grid centers (degrees east)
    lat: np.ndarray
        latitudes of grid centers (degrees north)
    c: str
        tidal constituent ID

    References
    ----------
    .. [1] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of
        Barotropic Ocean Tides," *Journal of Atmospheric and Oceanic
        Technology*, 19(2), 183--204, (2002).
        `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

    .. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
    """

    # grid spacing in radians
    dph = (lon[1] - lon[0])*np.pi/180.0
    dth = (lat[1] - lat[0])*np.pi/180.0
    # calculate meshgrid from latitude and longitude
    gridlon, gridlat = np.meshgrid(lon, lat)

    # longitude (in radians) of u and v nodes
    phi_v = gridlon*np.pi/180.0
    # x-coordinates for u transports
    phi_u = phi_v - dph/2.0

    # colatitudes (in radians) of u and v nodes
    th_u = (90.0 - gridlat)*np.pi/180.0
    # y-coordinates for v transports
    th_v = th_u + dth/2.0

    # load parameters for each constituent
    amp, ph, omega, alpha, species = \
        pyTMD.arguments._constituent_parameters(c)
    # uniform zonal dependence of astronomical forcing
    ph = 0.0

    # loading effects
    c = (alpha*_wgs84.gamma*amp/_wgs84.rad_e)
    # calculate forcing for constituent
    Fu = c*np.exp(1j*phi_u*species + 1j*ph)
    Fv = c*np.exp(1j*phi_v*species + 1j*ph)
    # calculate latitudinal dependence of forcing
    # for a given spherical harmonic dependence
    if (species == 1):
        # diurnal species
        Fu *= 2j*np.cos(th_u)
        Fv *= 2.0*(2.0*np.sin(th_v)**2 - 1.0)
    elif (species == 2):
        # semidiurnal species
        Fu *= 2j*np.sin(th_u)
        Fv *= -2.0*(np.sin(th_v)*np.cos(th_v))
    else:
        # long-period species
        Fu *= 0.0 + 0j
        Fv *= -3.0*(np.sin(th_v)*np.cos(th_v))
    # return the generating forces
    return (Fu, Fv)

# PURPOSE: calculate the astronomical tide generating potential
# converted to an equilibrium tide elevation
def _generating_potential(
        lon: np.ndarray,
        lat: np.ndarray,
        c: str,
    ):
    """
    Computes the astronomical tide generating potential for a tidal
    constituent converted to an equilibrium tide elevation [1]_

    Parameters
    ----------
    lon: np.ndarray
        longitudes of grid centers (degrees east)
    lat: np.ndarray
        latitudes of grid centers (degrees north)
    c: str
        tidal constituent ID

    References
    ----------
    .. [1] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of
        Barotropic Ocean Tides," *Journal of Atmospheric and Oceanic
        Technology*, 19(2), 183--204, (2002).
        `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

    .. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
    """

    # calculate meshgrid from latitude and longitude
    gridlon, gridlat = np.meshgrid(lon, lat)

    # longitude (in radians)
    phi_h = gridlon*np.pi/180.0
    # colatitudes (in radians)
    th_h = (90.0 - gridlat)*np.pi/180.0

    # load parameters for each constituent
    amp, ph, omega, alpha, species = \
        pyTMD.arguments._constituent_parameters(c)
    # uniform zonal dependence of astronomical potential
    ph = 0.0

    # calculate potential for constituent
    c = (alpha*amp)
    zeta = c*np.exp(1j*phi_h*species + 1j*ph)

    # calculate latitudinal dependence of potential
    # for a given spherical harmonic dependence
    if (species == 1):
        # diurnal species
        zeta *= 2*0*np.sin(th_h)*np.cos(th_h)
    elif (species == 2):
        # semidiurnal species
        zeta *= np.sin(th_h)**2
    else:
        # long-period species
        zeta *= -1.0*(3.0/2.0*np.cos(th_h)**2 - 1.0/2.0)
    # return the generating potential
    return zeta
