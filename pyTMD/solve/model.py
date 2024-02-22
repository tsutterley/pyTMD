#!/usr/bin/env python
u"""
model.py
Written by Tyler Sutterley (02/2024)
Class for estimating tidal constituents using hydrodynamic modeling

REFERENCES:
    G. D. Egbert and S. Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/

PROGRAM DEPENDENCIES:
    arguments.py: loads nodal corrections for tidal constituents
    astro.py: computes the basic astronomical mean longitudes
    grid.py: sets up finite difference grids for tidal modeling

UPDATE HISTORY:
    Written 02/2024
"""
from __future__ import annotations

import numpy as np
import scipy.sparse
import scipy.ndimage
import pyTMD.arguments

class _kernels:
    """
    Class for defining convolution kernels for finite differences
    """
    # convolution arrays for differentials
    dx = np.array([[1, -1]]) # dz/dx
    dy = np.array([[1], [-1]]) # dz/dy
    dx2 = np.array([[1, -2, 1]]) # d2z/dx2
    dy2 = np.array([[1], [-2], [1]]) # d2z/dy2

class model:
    """
    Class for estimating tidal constituents using hydrodynamic modeling

    Parameters
    ----------
    grid: pyTMD.solve.grid
        Finite difference grid
    beta: float, default 0.90
        load tide and self attraction correction
    rho_w: float, default 1035.0
        average density of sea water [kg/m^3]
    c_drag: float, default 0.003
        non-dimensional quadratic bottom drag coefficient
    h2: float, default 0.609
        Love number of vertical displacement
    k2: float, default 0.302
        Love number of gravitational potential
    N_0: float, default 5.24e-3
        Brunt-Vaisala frequency
    b_len: float, default 1300
        Brunt-Vaisala decay length [m]

    References
    ----------
    .. [1] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of
        Barotropic Ocean Tides," *Journal of Atmospheric and Oceanic
        Technology*, 19(2), 183--204, (2002).
        `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

    .. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
    """
    # seconds per day
    day = 86400.0
    # 360 degrees
    turndeg = 360.0
    # degrees to radians
    deg2rad = np.pi/180.0
    # convolution kernels
    _kernel = _kernels()

    def __init__(self, grid, **kwargs):
        # set initial attributes
        # set the finite difference grid
        assert isinstance(grid, pyTMD.solve.grid)
        self.grid = grid
        # model parameters and constants
        # load tide and self attraction correction
        self.beta = 0.90
        # average density of sea water [kg/m^3]
        self.rho_w = 1035.0
        # non-dimensional quadratic bottom drag coefficient
        # ranges from 0.0025 to 0.003 [Bowden and Fairburn, 1952]
        self.c_drag = 0.003
        # Love numbers appropriate for body tides [Wahr et al., (1985)]
        self.h2 = 0.609
        self.k2 = 0.302
        # Brunt-Vaisala frequency
        self.N_0 = 5.24e-3
        # Brunt-Vaisala decay length [m]
        self.b_len = 1300
        # set optional fields or redefine constants
        for key, val in kwargs.items():
            setattr(self, key, val)
        # validate inputs
        self.__validate__()

    # PURPOSE: calculate the astronomical tide generating force
    def generating_force(self, c: str):
        """
        Computes the astronomical tide generating force (the horizontal gradient
        of the tide generating potential) for a tidal constituent

        Parameters
        ----------
        c: str
            tidal constituent ID
        """
        # longitudes and latitudes for u and v transports
        lon_u, lat_u = self.grid.get_latlon(self.grid.x_u, self.grid.y_u)
        lon_v, lat_v = self.grid.get_latlon(self.grid.x_v, self.grid.y_v)
        # longitudes and colatitudes for u transports in radians
        phi_u = self.deg2rad*lon_u
        th_u = self.deg2rad*(90.0 - lat_u)
        # longitudes and colatitudes for v transports in radians
        phi_v = self.deg2rad*lon_v
        th_v = self.deg2rad*(90.0 - lat_v)
        # gravitational acceleration at each colatitude [m/s^2]
        gamma_u = self.grid.gamma_0(th_u)
        gamma_v = self.grid.gamma_0(th_v)

        # load parameters for each constituent
        amp, ph, omega, alpha, species = \
            pyTMD.arguments._constituent_parameters(c)
        # uniform zonal dependence of astronomical forcing
        ph = 0.0

        # calculate forcing for constituent
        Fu = alpha*amp*gamma_u*np.exp(1j*phi_u*species + 1j*ph)/self.grid.rad_e
        Fv = alpha*amp*gamma_v*np.exp(1j*phi_v*species + 1j*ph)/self.grid.rad_e
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
    def generating_potential(self, c: str):
        """
        Computes the astronomical tide generating potential for a tidal
        constituent converted to an equilibrium tide elevation

        Parameters
        ----------
        c: str
            tidal constituent ID
        """
        # longitudes and latitudes for zeta nodes
        lon_z, lat_z = self.grid.get_latlon(self.grid.x_z, self.grid.y_z)
        # longitudes and colatitudes for zeta nodes in radians
        phi_z = self.deg2rad*lon_z
        th_z = self.deg2rad*(90.0 - lat_z)

        # load parameters for each constituent
        amp, ph, omega, alpha, species = \
            pyTMD.arguments._constituent_parameters(c)
        # uniform zonal dependence of astronomical potential
        ph = 0.0

        # calculate potential for constituent
        zeta = alpha*amp*np.exp(1j*phi_z*species + 1j*ph)

        # calculate latitudinal dependence of potential
        # for a given spherical harmonic dependence
        if (species == 1):
            # diurnal species
            zeta *= 2.0*np.sin(th_z)*np.cos(th_z)
        elif (species == 2):
            # semidiurnal species
            zeta *= np.sin(th_z)**2
        else:
            # long-period species
            zeta *= -1.0*(3.0/2.0*np.cos(th_z)**2 - 1.0/2.0)
        # return the generating potential and the angular frequency
        return (zeta, omega)

    def convolve(self,
            array: np.ndarray,
            kernel: np.ndarray,
            **kwargs
        ):
        """Convolve a two-dimensional array with a given kernel

        Parameters
        ----------
        array: np.ndarray
            Input array
        kernel: np.ndarray
            Weights kernel
        kwargs: dict
            Keyword arguments for scipy.ndimage.convolve
        """
        # define default for how array will be extended at boundaries
        if self.grid.is_global:
            kwargs.setdefault('mode', 'wrap')
        else:
            kwargs.setdefault('mode', 'nearest')
        # calculate the convolution
        return scipy.ndimage.convolve(array, kernel, **kwargs)

    def __validate__(self):
        """Check if class inputs are appropriate
        """
        assert isinstance(self.beta, (int, float))
        assert isinstance(self.rho_w, (int, float))
        assert isinstance(self.c_drag, (int, float))
        assert isinstance(self.h2, (int, float))
        assert isinstance(self.k2, (int, float))
        assert isinstance(self.N_0, (int, float))
        assert isinstance(self.b_len, (int, float))

    def __str__(self):
        """String representation of the ``model`` object
        """
        properties = ['pyTMD.solve.model']
        return '\n'.join(properties)
