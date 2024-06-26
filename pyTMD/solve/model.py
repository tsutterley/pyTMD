#!/usr/bin/env python
u"""
model.py
Written by Tyler Sutterley (05/2024)
Base class for estimating tidal constituents using hydrodynamic modeling

REFERENCES:
    G. D. Egbert and S. Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/docmatrix/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
        https://pyproj4.github.io/pyproj/

PROGRAM DEPENDENCIES:
    arguments.py: loads nodal corrections for tidal constituents
    astro.py: computes the basic astronomical mean longitudes
    grid.py: sets up finite difference grids for tidal modeling
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 05/2024: make subscriptable and allow item assignment
    Written 02/2024
"""
from __future__ import annotations

import numpy as np
import pyTMD.arguments
import pyTMD.solve.grid as fdgrid

class model:
    """
    Base class for estimating tidal constituents using
    hydrodynamic modeling

    Parameters
    ----------
    grid: pyTMD.solve.grid
        Finite difference grid
    beta: float, default 0.10
        scaling factor for load tides and self attraction effects
    rho_w: float, default 1035.0
        average density of sea water [kg/m^3]
    c_drag: float, default 0.003
        non-dimensional quadratic bottom drag coefficient
    nu: float, default 100.0
        Eddy viscosity [m^2/s^2]
    Rd: float or NoneType, default None
        Rossby radius of deformation
    h2: float, default 0.609
        Love number of vertical displacement
    k2: float, default 0.302
        Love number of gravitational potential
    dt: float, default 30.0
        time step in seconds

    References
    ----------
    .. [1] G. D. Egbert and S. Y. Erofeeva, "Efficient Inverse Modeling of
        Barotropic Ocean Tides," *Journal of Atmospheric and Oceanic
        Technology*, 19(2), 183--204, (2002).
        `doi: 10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2`__

    .. __: https://doi.org/10.1175/1520-0426(2002)019<0183:EIMOBO>2.0.CO;2
    """
    # degrees to radians
    deg2rad = np.pi/180.0

    def __init__(self, grid: fdgrid, constituents: list = [], **kwargs):
        # set initial attributes
        # set the finite difference grid
        assert isinstance(grid, fdgrid)
        self.grid = grid
        # set the tidal constituents
        assert isinstance(constituents, (list, np.ndarray))
        self.constituents = constituents
        # model parameters and constants
        # load tide and self attraction scaling factor
        self.beta = 0.10
        # average density of sea water [kg/m^3]
        self.rho_w = 1035.0
        # non-dimensional quadratic bottom drag coefficient
        # ranges from 0.0025 to 0.003 [Bowden and Fairburn, 1952]
        self.c_drag = 0.003
        # eddy viscosity [m^2/s^2]
        self.nu = 100.0
        # Rossby radius of deformation
        self.Rd = None
        # Love numbers appropriate for body tides [Wahr et al., (1985)]
        self.h2 = 0.609
        self.k2 = 0.302
        # time step
        self.dt = 30.0
        # set optional fields or redefine constants
        for key, val in kwargs.items():
            setattr(self, key, val)
        # validate inputs
        self.__validate__()

    # PURPOSE: load constituent parameters
    def constituent_parameters(self, c: str):
        """
        Loads parameters for a given tidal constituent

        Parameters
        ----------
        c: str
            tidal constituent ID

        Returns
        -------
        amplitude: float
            amplitude of equilibrium tide for tidal constituent (meters)
        phase: float
            phase of tidal constituent (radians)
        omega: float
            angular frequency of constituent (radians)
        alpha: float
            load love number of tidal constituent
        species: float
            spherical harmonic dependence of quadrupole potential
        """
        return pyTMD.arguments._constituent_parameters(c)

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
        amp, ph, omega, alpha, species = self.constituent_parameters(c)
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
        amp, ph, omega, alpha, species = self.constituent_parameters(c)
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

    @property
    def shape(self) -> tuple:
        """Shape of the grid
        """
        return self.grid.shape

    @property
    def size(self) -> int:
        """Total number of nodes
        """
        return self.grid.size

    @property
    def nx(self) -> int:
        """Number of nodes in the x-direction
        """
        return self.grid.shape[1]

    @property
    def ny(self) -> int:
        """Number of nodes in the y-direction
        """
        return self.grid.shape[0]

    def __validate__(self):
        """Check if class inputs are appropriate
        """
        NoneType = type(None)
        assert isinstance(self.beta, (int, float))
        assert isinstance(self.rho_w, (int, float))
        assert isinstance(self.c_drag, (int, float))
        assert isinstance(self.nu, (int, float))
        assert isinstance(self.Rd, (int, float, NoneType))
        assert isinstance(self.h2, (int, float))
        assert isinstance(self.k2, (int, float))
        assert isinstance(self.dt, (int, float))

    def __str__(self):
        """String representation of the ``model`` object
        """
        properties = ['pyTMD.solve.model']
        shape = ', '.join(map(str, self.shape))
        properties.append(f"    shape: {shape}")
        return '\n'.join(properties)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)
