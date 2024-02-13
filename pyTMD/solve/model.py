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

UPDATE HISTORY:
    Written 02/2024
"""
from __future__ import annotations

import numpy as np
import scipy.sparse
import scipy.ndimage
import pyTMD.arguments
from pyTMD.constants import constants

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
    extent: list, tuple or np.ndarray
        spatial extent of the grid in degrees
    spacing: list, tuple or np.ndarray
        grid spacing in degrees
    a_axis: float, default 6378137.0
        Semi-major axis of the Earth's ellipsoid [m]
    flat: float, default 1.0/298.257223563
        Flattening of the Earth's ellipsoid
    GM: float, default 3.986004418e14
        Geocentric Gravitational Constant [m^3/s^2]
    gamma: float, default 9.80665
        gravitational acceleration [m/s^2]
    omega: float, default 7.292115e-5
        average angular rotation rate of the Earth [rad/s]
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

    def __init__(self, **kwargs):
        # set initial attributes
        self.extent = [None, None, None, None]
        self.spacing = [None, None]
        # model parameters and constants
        # semi-major axis of the ellipsoid [m] (WGS84 definition)
        self.a_axis = 6378137.0
        # flattening of the ellipsoid (WGS84 definition)
        self.flat = 1.0/298.257223563
        # Geocentric Gravitational Constant [m^3/s^2]
        self.GM = 3.986004418e14
        # standard gravitational acceleration [m/s^2]
        # (World Meteorological Organization)
        self.gamma = 9.80665
        # average angular rotation rate of the Earth [rad/s]
        self.omega = 7292115e-11
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

    # PURPOSE: calculate the astronomical tide generating force
    def generating_force(self, c: str):
        """
        Computes the astronomical tide generating force (the horizontal gradient
        of the tide generating potential) for a tidal constituent [1]_

        Parameters
        ----------
        c: str
            tidal constituent ID
        """
        # x and y coordinates for u transports in radians
        phi_u, th_u = np.meshgrid(
            self.deg2rad*self.lon_u,
            self.deg2rad*(90.0 - self.lat_u)
        )
        # x and y coordinates for v transports in radians
        phi_v, th_v = np.meshgrid(
            self.deg2rad*self.lon_v,
            self.deg2rad*(90.0 - self.lat_v)
        )

        # load parameters for each constituent
        amp, ph, omega, alpha, species = \
            pyTMD.arguments._constituent_parameters(c)
        # uniform zonal dependence of astronomical forcing
        ph = 0.0

        # calculate forcing for constituent
        coeff = (alpha*self.gamma*amp/self.rad_e)
        Fu = coeff*np.exp(1j*phi_u*species + 1j*ph)
        Fv = coeff*np.exp(1j*phi_v*species + 1j*ph)
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
        constituent converted to an equilibrium tide elevation [1]_

        Parameters
        ----------
        c: str
            tidal constituent ID
        """
        # x and y coordinates for heights in radians
        phi_z, th_z = np.meshgrid(
            self.deg2rad*self.lon_z,
            self.deg2rad*(90.0 - self.lat_z)
        )

        # load parameters for each constituent
        amp, ph, omega, alpha, species = \
            pyTMD.arguments._constituent_parameters(c)
        # uniform zonal dependence of astronomical potential
        ph = 0.0

        # calculate potential for constituent
        coeff = (alpha*amp)
        zeta = coeff*np.exp(1j*phi_z*species + 1j*ph)

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

    def f(self, lat: np.ndarray):
        """Coriolis parameter at latitudes

        Parameters
        ----------
        lat: np.ndarray
            Latitude [degrees]
        """
        theta = self.deg2rad*(90.0 - lat)
        return 2.0*self.omega*np.cos(theta)

    def gamma_0(self, lat: np.ndarray):
        """Normal gravity at latitudes

        Parameters
        ----------
        lat: np.ndarray
            Latitude [degrees]
        """
        theta = self.deg2rad*(90.0 - lat)
        return self.ellipsoid.gamma_0(theta)

    # PURPOSE: construct masks for u, v and zeta nodes
    def _mask_nodes(self, hz: np.ndarray):
        """
        Construct masks for u, v and zeta nodes on a C-grid

        Parameters
        ----------
        hz: np.ndarray
            bathymetry at grid centers
        """
        # shape of input bathymetry
        ny, nx = self.shape
        # for grid center mask: find where bathymetry is greater than 0
        self.mask_z = (hz > 0).astype(int)
        # initialize integer masks for u and v grids
        self.mask_u = np.zeros((ny, nx), dtype=int)
        self.mask_v = np.zeros((ny, nx), dtype=int)
        # wrap mask if global
        if self.global_grid:
            # x-indices
            indx = np.zeros((nx), dtype=int)
            indx[:-1] = np.arange(1, nx)
            indx[-1] = 0
            # y-indices
            indy = np.zeros((ny), dtype=int)
            indy[:-1] = np.arange(1, ny)
            indy[-1] = 0
            # calculate masks on u and v grids
            self.mask_u[indy,:] = self.mask_z*self.mask_z[indy,:]
            self.mask_v[:,indx] = self.mask_z*self.mask_z[:,indx]
        else:
            # x-indices
            indx = np.zeros((nx), dtype=int)
            indx[0] = 0
            indx[1:] = np.arange(nx-1)
            # y-indices
            indy = np.zeros((ny), dtype=int)
            indy[0] = 0
            indy[1:] = np.arange(ny-1)
            # calculate masks on u and v grids
            self.mask_u[:,:] = self.mask_z*self.mask_z[indy,:]
            self.mask_v[:,:] = self.mask_z*self.mask_z[:,indx]
        # return the masks
        return self

    # PURPOSE: interpolate data to u and v nodes
    def _interpolate_to_nodes(self, data_z: np.ndarray):
        """
        Interpolate data from zeta nodes to u and v nodes on a C-grid

        Parameters
        ----------
        data_z: np.ndarray
            data at grid centers
        """
        # shape of input data
        ny, nx = self.shape
        # initialize data for u and v grids
        data_u = np.zeros((ny, nx), dtype=data_z.dtype)
        data_v = np.zeros((ny, nx), dtype=data_z.dtype)
        # wrap data if global
        if self.global_grid:
            # x-indices
            indx = np.zeros((nx), dtype=int)
            indx[:-1] = np.arange(1, nx)
            indx[-1] = 0
            # y-indices
            indy = np.zeros((ny), dtype=int)
            indy[:-1] = np.arange(1, ny)
            indy[-1] = 0
            # calculate data at u and v nodes
            data_u[indy,:] = self.mask_u*(data_z + data_z[indy,:])/2.0
            data_v[:,indx] = self.mask_v*(data_z + data_z[:,indx])/2.0
        else:
            # x-indices
            indx = np.zeros((nx), dtype=int)
            indx[0] = 0
            indx[1:] = np.arange(nx-1)
            # y-indices
            indy = np.zeros((ny), dtype=int)
            indy[0] = 0
            indy[1:] = np.arange(ny-1)
            # calculate data at u and v nodes
            data_u[:,:] = self.mask_u*(data_z + data_z[indy,:])/2.0
            data_v[:,:] = self.mask_v*(data_z + data_z[:,indx])/2.0
        # return the interpolated data values
        return (data_u, data_v)

    def _convolve(self,
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
        if self.global_grid:
            kwargs.setdefault('mode', 'wrap')
        else:
            kwargs.setdefault('mode', 'nearest')
        # calculate the convolution
        return scipy.ndimage.convolve(array, kernel, **kwargs)

    # grid dimensions
    @property
    def dimensions(self) -> list:
        """Dimensions of the grid
        """
        dims = [None, None]
        # calculate y dimensions with extents
        dims[0] = np.int64((self.extent[3] - self.extent[2])/self.spacing[1])
        # calculate x dimensions with extents
        dims[1] = np.int64((self.extent[1] - self.extent[0])/self.spacing[0])
        return dims

    @property
    def shape(self) -> tuple:
        """Shape of the grid
        """
        return (self.dimensions[0], self.dimensions[1], )

    @property
    def ndim(self) -> int:
        """Number of dimensions in the grid
        """
        return len(self.shape)

    @property
    def size(self) -> int:
        """Total number of nodes in the grid
        """
        return np.prod(self.shape)

    @property
    def global_grid(self) -> bool:
        """Defines if the grid is global in terms of longitude
        """
        return np.isclose(
            self.lon_z[-1] - self.lon_z[0],
            self.turndeg - self.spacing[0]
        )

    @property
    def ellipsoid(self):
        """Ellipsoidal parameters based on user definitions
        """
        return constants(
            ellipsoid=None,
            a_axis=self.a_axis,
            flat=self.flat,
            GM=self.GM,
            omega=self.omega
        )

    @property
    def rad_e(self) -> float:
        """Radius of the Earth with same volume as ellipsoid [m]
        """
        return self.ellipsoid.rad_e

    @property
    def lon_z(self):
        """Longitude coordinates for zeta nodes on an Arakawa C-grid
        """
        return self.extent[0] + self.spacing[0]*np.arange(self.dimensions[1]) + \
            self.spacing[0]/2.0

    @property
    def lat_z(self):
        """Latitude coordinates for zeta nodes on an Arakawa C-grid
        """
        return self.extent[2] + self.spacing[1]*np.arange(self.dimensions[0]) + \
            self.spacing[1]/2.0

    @property
    def lon_u(self):
        """Longitude coordinates for U nodes on an Arakawa C-grid
        """
        return self.extent[0] + self.spacing[0]*np.arange(self.dimensions[1])

    @property
    def lat_u(self):
        """Latitude coordinates for U nodes on an Arakawa C-grid
        """
        return self.extent[2] + self.spacing[1]*np.arange(self.dimensions[0]) + \
            self.spacing[1]/2.0

    @property
    def lon_v(self):
        """Longitude coordinates for V nodes on an Arakawa C-grid
        """
        return self.extent[0] + self.spacing[0]*np.arange(self.dimensions[1]) + \
            self.spacing[0]/2.0

    @property
    def lat_v(self):
        """Latitude coordinates for V nodes on an Arakawa C-grid
        """
        return self.extent[2] + self.spacing[1]*np.arange(self.dimensions[0])

    def __str__(self):
        """String representation of the ``model`` object
        """
        properties = ['pyTMD.solve.model']
        extent = ', '.join(map(str, self.extent))
        properties.append(f"    extent: {extent}")
        shape = ', '.join(map(str, self.shape))
        properties.append(f"    shape: {shape}")
        return '\n'.join(properties)
