#!/usr/bin/env python
u"""
crs.py
Written by Tyler Sutterley (12/2023)
Coordinates Reference System (CRS) routines

CALLING SEQUENCE:
    x, y = pyTMD.crs.convert(lon, lat, PROJ, 'F')
    lon, lat = pyTMD.crs.convert(x, y, PROJ, 'B')

INPUTS:
    i1: longitude ('F') or projection easting x ('B')
    i2: latitude ('F') or projection northing y ('B')
    PROJ: spatial reference system code for coordinate transformations
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

# attempt imports
try:
    import pyproj
except (AttributeError, ImportError, ModuleNotFoundError) as exc:
    logging.critical("pyproj not available")

class crs:
    """Coordinate Reference System transformations for tide models

    Attributes
    ----------
    transformer: obj
        pyproj transformer for changing coordinate reference system
    direction: obj
        pyproj transform direction
    name: str
        Projection name
    """
    def __init__(self):
        self.transformer = None
        self.direction = pyproj.enums.TransformDirection.FORWARD
        self.name = None
    
    def convert(self,
                i1: np.ndarray,
                i2: np.ndarray,
                PROJ: str,
                BF: str,
                EPSG: int | str = 4326):
        """
        Converts points to and from Coordinates Reference Systems (CRS)

        Parameters
        ----------
        i1: np.ndarray
            Longitude (``'F'``) or projected x-coordinates (``'B'``)
        i2: np.ndarray
            Latitude (``'F'``) or projected y-coordinates (``'B'``)
        PROJ: str
            Spatial reference system code for coordinate transformations
        BF: str
            Direction of translation

                - ``'B'``: backwards
                - ``'F'``: forwards
        EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
            input (``'F'``) or output (``'B'``) coordinate system

        Returns
        -------
        o1: np.ndarray
            Projected x-coordinates (``'F'``) or longitude (``'B'``)
        o2: np.ndarray
            Projected y-coordinates (``'F'``) or latitude (``'B``')
        """
        self.name = PROJ
        # python dictionary with named conversion functions
        transforms = {}
        transforms['3031'] = self._EPSG3031
        transforms['3413'] = self._EPSG3413
        transforms['CATS2008'] = self._CATS2008
        transforms['3976'] = self._EPSG3976
        transforms['PSNorth'] = self._PSNorth
        transforms['4326'] = self._EPSG4326
        # convert lat/lon to Polar-Stereographic x/y
        if (BF.upper() == 'F'):
            self.direction = pyproj.enums.TransformDirection.FORWARD
        # convert Polar-Stereographic x/y to lat/lon
        elif (BF.upper() == 'B'):
            self.direction = pyproj.enums.TransformDirection.INVERSE
        # check that PROJ for conversion was entered correctly
        # run named conversion program and return values
        try:
            transforms[PROJ](EPSG)
        except KeyError as exc:
            pass
        else:
            # return the output variables
            return self.transform(i1, i2)
        # try changing the projection using a custom projection
        # run custom conversion program and return values
        try:
            self._custom(PROJ, EPSG=EPSG)
        except Exception as exc:
            pass
        else:
            return self.transform(i1, i2)
        # projection not found or available
        raise Exception(f'PROJ: {PROJ} conversion function not found')
    
    def transform(self, i1, i2):
        """
        Performs Coordinates Reference System (CRS) transformations

        Parameters
        ----------
        i1: np.ndarray
            Input x-coordinates
        i2: np.ndarray
            Input y-coordinates
        
        Returns
        -------
        o1: np.ndarray
            Output transformed x-coordinates
        o2: np.ndarray
            Output transformed y-coordinates
        """
        if (self.name == 'PSNorth') and (self.direction.name == 'FORWARD'):
            # convert input coordinate reference system to lat/lon
            lon, lat = self.transformer.transform(i1, i2,
                direction=self.direction)
            # convert lat/lon to (idealized) Polar-Stereographic x/y
            o1 = (90.0 - lat)*111.7*np.cos(lon/180.0*np.pi)
            o2 = (90.0 - lat)*111.7*np.sin(lon/180.0*np.pi)
        elif (self.name == 'PSNorth') and (self.direction.name == 'INVERSE'):
            # convert (idealized) Polar-Stereographic x/y to lat/lon
            lon = np.arctan2(i2, i1)*180.0/np.pi
            lat = 90.0 - np.sqrt(i1**2 + i2**2)/111.7
            # adjust longitudes to be -180:180
            ii, = np.nonzero(lon < 0)
            lon[ii] += 360.0
            # convert to output coordinate reference system
            o1, o2 = self.transformer.transform(lon, lat,
                direction=self.direction)
        else:
            # convert coordinate reference system
            o1, o2 = self.transformer.transform(i1, i2,
                direction=self.direction)
        # return the transformed coordinates
        return (o1, o2)

    # PURPOSE: try to get the projection information
    def from_input(self, PROJECTION: int | str):
        """
        Attempt to get the Coordinate Reference System

        Parameters
        ----------
        PROJECTION: int or str
            Coordinate Reference System code
        """
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
    
    def _EPSG3031(self, EPSG: int | str = 4326):
        """
        Transform for models in EPSG:3031 (Antarctic Polar Stereographic)

        Parameters
        ----------
        EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
            input (``'F'``) or output (``'B'``) coordinate system
        """
        # coordinate reference system information
        crs1 = self.from_input(EPSG)
        crs2 = pyproj.CRS.from_user_input({'proj':'stere',
            'lat_0':-90, 'lat_ts':-71, 'lon_0':0, 'x_0':0., 'y_0':0.,
            'ellps':'WGS84', 'datum':'WGS84', 'units':'km'})
        self.transformer = pyproj.Transformer.from_crs(crs1, crs2,
            always_xy=True)
        return self

    # wrapper function for models in EPSG 3413 (Sea Ice Polar Stereographic North)
    def _EPSG3413(self, EPSG: int | str = 4326):
        """
        Transform for models in EPSG:3413 (Sea Ice Polar Stereographic North)

        Parameters
        ----------
        EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
            input (``'F'``) or output (``'B'``) coordinate system
        """
        # coordinate reference system information
        crs1 = self.from_input(EPSG)
        crs2 = pyproj.CRS.from_user_input({'proj':'stere',
            'lat_0':90, 'lat_ts':70, 'lon_0':-45, 'x_0':0., 'y_0':0.,
            'ellps':'WGS84', 'datum':'WGS84', 'units':'km'})
        self.transformer = pyproj.Transformer.from_crs(crs1, crs2,
            always_xy=True)
        return self

    # wrapper function for CATS2008 tide models
    def _CATS2008(self, EPSG: int | str = 4326):
        """
        Transform for Circum-Antarctic Tidal Simulation models

        Parameters
        ----------
        EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
            input (``'F'``) or output (``'B'``) coordinate system
        """
        # coordinate reference system information
        crs1 = self.from_input(EPSG)
        crs2 = pyproj.CRS.from_user_input({'proj':'stere',
            'lat_0':-90, 'lat_ts':-71, 'lon_0':-70, 'x_0':0., 'y_0':0.,
            'ellps':'WGS84', 'datum':'WGS84', 'units':'km'})
        self.transformer = pyproj.Transformer.from_crs(crs1, crs2,
            always_xy=True)
        return self

    # wrapper function for models in EPSG 3976 (NSIDC Sea Ice Stereographic South)
    def _EPSG3976(self, EPSG: int | str = 4326):
        """
        Transform for models in EPSG:3976 (Sea Ice Polar Stereographic South)

        Parameters
        ----------
        EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
            input (``'F'``) or output (``'B'``) coordinate system
        """
        # coordinate reference system information
        crs1 = self.from_input(EPSG)
        crs2 = pyproj.CRS.from_user_input({'proj':'stere',
            'lat_0':-90, 'lat_ts':-70, 'lon_0':0, 'x_0':0., 'y_0':0.,
            'ellps':'WGS84', 'datum':'WGS84', 'units':'km'})
        self.transformer = pyproj.Transformer.from_crs(crs1, crs2,
            always_xy=True)
        return self
    
    # function for models in (idealized) PSNorth projection
    def _PSNorth(self, EPSG: int | str = 4326):
        """
        Transform for idealized Arctic Polar Stereographic models

        Parameters
        ----------
        EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
            input (``'F'``) or output (``'B'``) coordinate system
        """
        # projections for converting to and from input EPSG
        crs1 = self.from_input(EPSG)
        crs2 = pyproj.CRS.from_epsg(4326)
        self.transformer = pyproj.Transformer.from_crs(crs1, crs2,
            always_xy=True)
        return self

    # wrapper function to pass lat/lon values or convert if EPSG
    def _EPSG4326(self, EPSG: int | str = 4326):
        """
        Transform for models in EPSG:4326 (WGS84 Latitude/Longitude)

        Parameters
        ----------
        EPSG: int, default 4326 (WGS84 Latitude/Longitude)
            input (``'F'``) or output (``'B'``) coordinate system
        """
        crs1 = self.from_input(EPSG)
        crs2 = pyproj.CRS.from_epsg(4326)
        self.transformer = pyproj.Transformer.from_crs(crs1, crs2,
            always_xy=True)
        return self

    # wrapper function for using custom projections
    def _custom(self, PROJ: int | str, EPSG: int | str = 4326):
        """
        Transform for models in a custom projection

        Parameters
        ----------
        PROJ: int or str
            Spatial reference system code for coordinate transformations
        EPSG: int or str, default 4326 (WGS84 Latitude/Longitude)
            input (``'F'``) or output (``'B'``) coordinate system
        """
        # coordinate reference system information
        crs1 = self.from_input(EPSG)
        crs2 = self.from_input(PROJ)
        self.transformer = pyproj.Transformer.from_crs(crs1, crs2,
            always_xy=True)
        return self

# simplify calling crs class
crs = crs()
