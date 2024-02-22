#!/usr/bin/env python
u"""
grid.py
Written by Tyler Sutterley (02/2024)
Class for setting up finite difference grids for tidal modeling

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
        https://pyproj4.github.io/pyproj/

PROGRAM DEPENDENCIES:
    crs.py: Coordinate Reference System (CRS) routines
    spatial.py: utilities for working with geospatial data

UPDATE HISTORY:
    Written 02/2024
"""
from __future__ import annotations

import logging
import numpy as np
from pyTMD.crs import datum
from pyTMD.spatial import scale_areas

# attempt imports
try:
    import pyproj
except (AttributeError, ImportError, ModuleNotFoundError) as exc:
    logging.critical("pyproj not available")

class grid:
    """
    Class for setting up finite difference grids for tidal modeling

    Parameters
    ----------
    extent: list, tuple or np.ndarray
        spatial extent of the grid
    spacing: list, tuple or np.ndarray
        grid spacing
    crs: int, str, pyproj.CRS or NoneType, default None
        Coordinate Reference System definition
    a_axis: float, default 6378137.0
        Semi-major axis of the Earth's ellipsoid [m]
    flat: float, default 1.0/298.257223563
        Flattening of the Earth's ellipsoid
    GM: float, default 3.986004418e14
        Geocentric Gravitational Constant [m^3/s^2]
    gamma: float, default 9.80665
        standard gravitational acceleration [m/s^2]
    omega: float, default 7.292115e-5
        average angular rotation rate of the Earth [rad/s]
    """
    # 360 degrees
    turndeg = 360.0
    # degrees to radians
    deg2rad = np.pi/180.0

    def __init__(self, **kwargs):
        # set initial attributes
        self.extent = [None, None, None, None]
        self.spacing = [None, None]
        # set the spatial projection reference information
        if ('crs' in kwargs):
            self._crs(kwargs.pop('crs'))
            self._get_geod()
        # model parameters and constants
        # Geocentric Gravitational Constant [m^3/s^2]
        self.GM = 3.986004418e14
        # standard gravitational acceleration [m/s^2]
        # (World Meteorological Organization)
        self.gamma = 9.80665
        # average angular rotation rate of the Earth [rad/s]
        self.omega = 7292115e-11
        # set optional fields or redefine constants
        for key, val in kwargs.items():
            setattr(self, key, val)
        # validate inputs
        self.__validate__()

    def get_latlon(self,
            x: np.ndarray,
            y: np.ndarray
        ):
        """
        Get the latitude and longitude of grid cells

        Parameters
        ----------
        x: np.ndarray
            x-coordinates of grid cells
        y: np.ndarray
            y-coordinates of grid cells

        Returns
        -------
        longitude: np.ndarray
            longitude coordinates of grid cells
        latitude: np.ndarray
            latitude coordinates of grid cells
        """
        # target spatial reference (WGS84 latitude and longitude)
        target_crs = pyproj.CRS.from_epsg(4326)
        # create transformation
        transformer = pyproj.Transformer.from_crs(self.crs, target_crs,
            always_xy=True)
        # create meshgrid of points in original projection
        xgrid, ygrid = np.meshgrid(x, y)
        # convert coordinates to latitude and longitude
        longitude, latitude = transformer.transform(xgrid, ygrid)
        # return the latitude and longitude of the grid cells
        return longitude, latitude

    def get_step_size(self,
            x: np.ndarray,
            y: np.ndarray,
            scale: float = 1.0
        ):
        """Calculate the step size of the grid cells in meters

        Parameters
        ----------
        x: np.ndarray
            x-coordinates of grid cells
        y: np.ndarray
            y-coordinates of grid cells
        scale: float, default 1.0
            scaling factor to convert the step size to meters

        Returns
        -------
        dx: np.ndarray
            step size in the x-direction [m]
        dy: np.ndarray
            step size in the y-direction [m]
        """
        # Get the latitude and longitude of grid cells
        longitude, latitude = self.get_latlon(x, y)
        if self.is_geographic:
            # calculate the step sizes in radians
            dphi = self.deg2rad*np.abs(x[1] - x[0])
            dth = self.deg2rad*np.abs(y[1] - y[0])
            # convert the step sizes to meters
            dx = self.rad_e*dphi*np.cos(self.deg2rad*latitude)
            dy = self.rad_e*dth*np.ones_like(latitude)
            # assume a spherical Earth
        elif self.is_stereographic:
            # stereographic scaling factors for areal distortion
            reference_latitude = self._cf['standard_parallel']
            ps_scale = scale_areas(latitude, flat=self.flat,
                ref=reference_latitude)
            # calculate the step sizes in meters
            # taking into account the mapping distortion
            dx = scale*np.abs(x[1] - x[0])*np.sqrt(ps_scale)
            dy = scale*np.abs(y[1] - y[0])*np.sqrt(ps_scale)
        else:
            # calculate the step sizes in meters
            dx = scale*np.abs(x[1] - x[0])*np.ones_like(latitude)
            dy = scale*np.abs(y[1] - y[0])*np.ones_like(latitude)
        # return the grid cell step sizes
        return (dx, dy)

    # PURPOSE: try to get the projection information
    def _crs(self, projection: int | str | pyproj.CRS):
        """
        Attempt to retrieve the Coordinate Reference System

        Parameters
        ----------
        projection: int, str or pyproj.CRS
            Coordinate Reference System definition
        """
        if isinstance(projection, pyproj.CRS):
            self.crs = projection
            return self
        # EPSG projection code
        try:
            self.crs = pyproj.CRS.from_epsg(int(projection))
        except (ValueError, pyproj.exceptions.CRSError):
            pass
        else:
            return self
        # coordinate reference system string
        try:
            self.crs = pyproj.CRS.from_string(projection)
        except (ValueError, pyproj.exceptions.CRSError):
            pass
        else:
            return self
        # no projection can be made
        raise pyproj.exceptions.CRSError

    def _get_geod(self):
        """Get the geodetic parameters describing the ellipsoid
        """
        self.a_axis = self.crs.get_geod().a
        self.flat = self.crs.get_geod().f
        return self

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
        if self.is_global:
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
        if self.is_global:
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
    def is_geographic(self) -> bool:
        """Defines if the grid is geographic
        """
        return self.crs.is_geographic

    @property
    def is_projected(self) -> bool:
        """Defines if the grid is projected
        """
        return self.crs.is_projected

    @property
    def is_stereographic(self) -> bool:
        """Defines if the grid is in a stereographic projection
        """
        return self.is_projected and ('stereo' in self._cf['grid_mapping_name'])

    @property
    def is_global(self) -> bool:
        """Defines if the grid is global in terms of longitude
        """
        return self.crs.is_geographic and \
            np.isclose(self.x_z[-1] - self.x_z[0],
                self.turndeg - self.spacing[0]
            )

    def f(self, theta: np.ndarray):
        """Coriolis parameter at colatitudes

        Parameters
        ----------
        theta: np.ndarray
            Colatitude [radians]
        """
        return 2.0*self.omega*np.cos(theta)

    def gamma_0(self, theta: np.ndarray):
        """Normal gravity at colatitudes

        Parameters
        ----------
        theta: np.ndarray
            Colatitude [radians]
        """
        return self._ellipsoid.gamma_0(theta)

    @property
    def rad_e(self) -> float:
        """Radius of the Earth with same volume as ellipsoid [m]
        """
        return self._ellipsoid.rad_e

    @property
    def _ellipsoid(self):
        """Ellipsoidal parameters based on user definitions
        """
        return datum(
            ellipsoid=None,
            a_axis=self.a_axis,
            flat=self.flat,
            GM=self.GM,
            omega=self.omega
        )

    @property
    def _cf(self) -> dict:
        """CF Convention metadata for the coordinate reference system
        """
        return self.crs.to_cf()

    @property
    def x_z(self):
        """x-coordinates for zeta nodes on an Arakawa C-grid
        """
        return self.extent[0] + self.spacing[0]*np.arange(self.dimensions[1]) + \
            self.spacing[0]/2.0

    @property
    def y_z(self):
        """y-coordinates for zeta nodes on an Arakawa C-grid
        """
        return self.extent[2] + self.spacing[1]*np.arange(self.dimensions[0]) + \
            self.spacing[1]/2.0

    @property
    def x_u(self):
        """x-coordinates for U nodes on an Arakawa C-grid
        """
        return self.extent[0] + self.spacing[0]*np.arange(self.dimensions[1])

    @property
    def y_u(self):
        """y-coordinates for U nodes on an Arakawa C-grid
        """
        return self.extent[2] + self.spacing[1]*np.arange(self.dimensions[0]) + \
            self.spacing[1]/2.0

    @property
    def x_v(self):
        """x-coordinates for V nodes on an Arakawa C-grid
        """
        return self.extent[0] + self.spacing[0]*np.arange(self.dimensions[1]) + \
            self.spacing[0]/2.0

    @property
    def y_v(self):
        """y-coordinates for V nodes on an Arakawa C-grid
        """
        return self.extent[2] + self.spacing[1]*np.arange(self.dimensions[0])

    def __validate__(self):
        """Check if class inputs are appropriate
        """
        assert isinstance(self.extent, (list, tuple, np.ndarray))
        assert len(self.extent) == 4
        assert isinstance(self.spacing, (list, tuple, np.ndarray))
        assert len(self.spacing) == 2
        assert isinstance(self.crs, pyproj.CRS)
        assert isinstance(self.a_axis, (int, float))
        assert isinstance(self.flat, (int, float))
        assert isinstance(self.GM, (int, float))
        assert isinstance(self.omega, (int, float))

    def __str__(self):
        """String representation of the ``grid`` object
        """
        properties = ['pyTMD.solve.grid']
        properties.append(f"    crs: {self.crs.name}")
        extent = ', '.join(map(str, self.extent))
        properties.append(f"    extent: {extent}")
        spacing = ', '.join(map(str, self.spacing))
        properties.append(f"    spacing: {spacing}")
        shape = ', '.join(map(str, self.shape))
        properties.append(f"    shape: {shape}")
        return '\n'.join(properties)
