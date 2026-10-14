#!/usr/bin/env python
"""
dataset.py
Written by Tyler Sutterley (10/2026)
An xarray.Dataset extension for tidal model data

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    pint: Python package to define, operate and manipulate physical quantities
        https://pypi.org/project/Pint/
        https://pint.readthedocs.io/en/stable
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
        https://pyproj4.github.io/pyproj/
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    xarray: N-D labeled arrays and datasets in Python
        https://docs.xarray.dev/en/stable/

UPDATE HISTORY:
    Updated 10/2026: added function to convert coordinates to geographic
        added grid cell area calculators for geographic and projected models
    Updated 08/2026: add isel_bounds; gridded crop uses index hyperslabs
        Pacific/dateline crop via dual lon hyperslabs (no half-globe pad)
        use np.hypot to calculate harmonic amplitudes
    Updated 06/2026: moved peak finding algorithm to prediction module
        drift type renamed to trajectory. drift still accepted as an alias
        added function to infer minor constituents and add to dataset
    Updated 05/2026: added parameters to allow for extrapolation with
        inverse distance weighting (IDW) in addition to nearest-neighbors (NN)
    Updated 04/2026: add barycentric interpolation for unstructured grids
        add support for unstructured (e.g. finite element) grids
        added function to calculate the high and low peaks of a prediction
        added function to try to convert units into a pint-friendly format
        added combine_attrs to merge conflicts into a list
    Updated 03/2026: allow caching of the kd-tree for extrapolation
    Updated 02/2026: create subaccessor registration functions
        add functions to test if units are compatible with known groups
    Updated 01/2026: handle scalar inputs for coordinate transformations
    Updated 12/2025: add coords functions to transform coordinates
        set units attribute for amplitude and phase data arrays
        add functions for assigning coordinates to datasets
    Updated 11/2025: get crs directly using pyproj.CRS.from_user_input
        set variable name to constituent for to_dataarray method
        added is_global property for models covering a global domain
        added pad function to pad global datasets along boundaries
        added inpaint function to fill missing data in datasets
    Updated 09/2025: added argument to limit the list of constituents
        when converting to an xarray DataArray
    Written 08/2025
"""

import re
import pint
import pyproj
import warnings
import numpy as np
import xarray as xr
from typing import Any
from xarray.core.utils import equivalent

# suppress warnings
warnings.filterwarnings("ignore", category=UserWarning)

__all__ = [
    "DataTree",
    "Dataset",
    "DataArray",
    "combine_attrs",
    "equivalent_attrs",
    "isel_bounds",
    "register_datatree_subaccessor",
    "register_dataset_subaccessor",
    "register_dataarray_subaccessor",
    "_transform",
    "_coords",
]

# pint unit registry
__ureg__ = pint.UnitRegistry()
# default units for pyTMD outputs
_default_units = {
    "elevation": "m",
    "current": "cm/s",
    "transport": "m^2/s",
}


@xr.register_datatree_accessor("tmd")
class DataTree:
    """Accessor for extending an ``xarray.DataTree`` for tidal model data"""

    def __init__(self, dtree):
        # initialize DataTree
        self._dtree = dtree

    def assign_coords(
        self,
        x: np.ndarray,
        y: np.ndarray,
        crs: str | int | dict = 4326,
        **kwargs,
    ):
        """
        Assign new coordinates to the ``DataTree``

        Parameters
        ----------
        x: np.ndarray
            Updated x-coordinates
        y: np.ndarray
            Updated y-coordinates
        crs: str, int, or dict, default 4326 (WGS84 Latitude/Longitude)
            Coordinate reference system of coordinates
        kwargs: dict
            Keyword arguments for ``xarray.Dataset.assign_coords``

        Returns
        -------
        dtree: xarray.DataTree
            ``DataTree`` with updated coordinates
        """
        # assign new coordinates to each dataset
        dtree = self._dtree.copy()
        for key, ds in self._dtree.items():
            ds = ds.to_dataset().assign_coords(dict(x=x, y=y), **kwargs)
            ds.attrs["crs"] = crs
            dtree[key] = ds
        # return the datatree
        return dtree

    def coords_as(
        self,
        x: np.ndarray,
        y: np.ndarray,
        crs: str | int | dict = 4326,
        **kwargs,
    ):
        """
        Transform coordinates into ``DataArrays`` in the ``DataTree``
        coordinate reference system

        Parameters
        ----------
        x: np.ndarray
            Input x-coordinates
        y: np.ndarray
            Input y-coordinates
        crs: str, int, or dict, default 4326 (WGS84 Latitude/Longitude)
            Coordinate reference system of input coordinates

        Returns
        -------
        X: xarray.DataArray
            Transformed x-coordinates
        Y: xarray.DataArray
            Transformed y-coordinates
        """
        # convert coordinate reference system to that of the datatree
        # and format as xarray DataArray with appropriate dimensions
        X, Y = _coords(x, y, source_crs=crs, target_crs=self.crs, **kwargs)
        # return the transformed coordinates
        return X, Y

    def crop(self, *args, **kwargs):
        """
        Crop ``DataTree`` to input bounding box

        Maps :meth:`Dataset.crop` over each group. For regular grids this
        is an index hyperslab (lazy when children are chunked).
        """
        # create copy of datatree
        dtree = self._dtree.copy()
        # crop each dataset in the datatree
        for key, ds in dtree.items():
            ds = ds.to_dataset()
            dtree[key] = ds.tmd.crop(*args, **kwargs)
        # return the datatree
        return dtree

    def inpaint(self, **kwargs):
        """
        Inpaint over missing data in ``DataTree``
        """
        # create copy of datatree
        dtree = self._dtree.copy()
        # inpaint each dataset in the datatree
        for key, ds in dtree.items():
            ds = ds.to_dataset()
            dtree[key] = ds.tmd.inpaint(**kwargs)
        # return the datatree
        return dtree

    def interp(
        self,
        x: np.ndarray,
        y: np.ndarray,
        **kwargs,
    ):
        """
        Interpolate ``DataTree`` to new coordinates

        Parameters
        ----------
        x: np.ndarray
            Interpolation x-coordinates
        y: np.ndarray
            Interpolation y-coordinates
        """
        # create copy of datatree
        dtree = self._dtree.copy()
        # interpolate each dataset in the datatree
        for key, ds in dtree.items():
            ds = ds.to_dataset()
            dtree[key] = ds.tmd.interp(x, y, **kwargs)
        # return the datatree
        return dtree

    def subset(self, c: str | list):
        """
        Reduce to a subset of constituents

        Parameters
        ----------
        c: str or list
            List of constituents names
        """
        # create copy of datatree
        dtree = self._dtree.copy()
        # subset each dataset in the datatree
        for key, ds in dtree.items():
            ds = ds.to_dataset()
            dtree[key] = ds.tmd.subset(c)
        # return the datatree
        return dtree

    def transform_as(
        self,
        x: np.ndarray,
        y: np.ndarray,
        crs: str | int | dict = 4326,
        **kwargs,
    ):
        """
        Transform coordinates to/from the ``DataTree`` coordinate reference system

        Parameters
        ----------
        x: np.ndarray
            Input x-coordinates
        y: np.ndarray
            Input y-coordinates
        crs: str, int, or dict, default 4326 (WGS84 Latitude/Longitude)
            Coordinate reference system of input coordinates
        direction: str, default 'FORWARD'
            Direction of transformation

            - ``'FORWARD'``: from input crs to model crs
            - ``'INVERSE'``: from model crs to input crs

        Returns
        -------
        X: np.ndarray
            Transformed x-coordinates
        Y: np.ndarray
            Transformed y-coordinates
        """
        # convert coordinate reference system to that of the datatree
        X, Y = _transform(x, y, source_crs=crs, target_crs=self.crs, **kwargs)
        # return the transformed coordinates
        return (X, Y)

    def to_ellipse(self, **kwargs):
        """
        Expresses tidal currents in terms of four ellipse parameters

        Returns
        -------
        dtree: xr.DataTree
            ``DataTree`` containing:

            - ``major``: amplitude of the semi-major axis
            - ``minor``: amplitude of the semi-minor axis
            - ``incl``: angle of inclination of the northern semi-major axis
            - ``phase``: phase lag of the current behind the tidal potential
        """
        from pyTMD.ellipse import ellipse

        # get u and v components from datatree
        dsu = (
            self._dtree.get("u", None) or self._dtree.get("U", None)
        ).to_dataset()
        dsv = (
            self._dtree.get("v", None) or self._dtree.get("V", None)
        ).to_dataset()
        # calculate ellipse parameters for each constituent
        dmajor = xr.Dataset()
        dminor = xr.Dataset()
        dincl = xr.Dataset()
        dphase = xr.Dataset()
        # for each constituent in the u-component
        for c in dsu.tmd.constituents:
            # assert units between datasets are the same
            if dsu[c].attrs.get("units", "") != dsv[c].attrs.get("units", ""):
                raise ValueError(
                    f"Incompatible units for {c} in u and v datasets"
                )
            # calculate ellipse parameters
            major, minor, incl, phase = ellipse(dsu[c].values, dsv[c].values)
            # create xarray DataArray for ellipse parameters
            dmajor[c] = xr.DataArray(major, dims=dsu[c].dims, coords=dsu.coords)
            dminor[c] = xr.DataArray(minor, dims=dsu[c].dims, coords=dsu.coords)
            dincl[c] = xr.DataArray(incl, dims=dsu[c].dims, coords=dsu.coords)
            dphase[c] = xr.DataArray(phase, dims=dsu[c].dims, coords=dsu.coords)
            # add attributes to each variable
            dmajor[c].attrs["units"] = dsu[c].attrs.get("units", "")
            dminor[c].attrs["units"] = dsu[c].attrs.get("units", "")
            dincl[c].attrs["units"] = "degrees"
            dphase[c].attrs["units"] = "degrees"
        # create output datatree
        dtree = xr.DataTree()
        # add datasets to output datatree
        dtree["major"] = dmajor
        dtree["minor"] = dminor
        dtree["incl"] = dincl
        dtree["phase"] = dphase
        # return datatree
        return dtree

    def from_ellipse(self, **kwargs):
        """
        Calculates tidal currents from the four ellipse parameters

        - ``major``: amplitude of the semi-major axis
        - ``minor``: amplitude of the semi-minor axis
        - ``incl``: angle of inclination of the northern semi-major axis
        - ``phase``: phase lag of the current behind the tidal potential

        Returns
        -------
        dtree: xr.DataTree
            ``DataTree`` containing transports or currents
        """
        from pyTMD.ellipse import inverse

        # get ellipse parameters from datatree
        dmajor = self._dtree["major"].to_dataset()
        dminor = self._dtree["minor"].to_dataset()
        dincl = self._dtree["incl"].to_dataset()
        dphase = self._dtree["phase"].to_dataset()
        # calculate currents for each constituent
        dsu = xr.Dataset()
        dsv = xr.Dataset()
        # for each constituent in the major parameter
        for c in dmajor.tmd.constituents:
            # calculate ellipse parameters
            u, v = inverse(
                dmajor[c].values,
                dminor[c].values,
                dincl[c].values,
                dphase[c].values,
            )
            # create xarray DataArray for ellipse parameters
            dsu[c] = xr.DataArray(u, dims=dmajor[c].dims, coords=dmajor.coords)
            dsv[c] = xr.DataArray(v, dims=dmajor[c].dims, coords=dmajor.coords)
            # add attributes to each variable
            dsu[c].attrs["units"] = dmajor[c].attrs.get("units", "")
            dsv[c].attrs["units"] = dmajor[c].attrs.get("units", "")
            if dmajor[c].tmd.group == "current":
                ukey, vkey = "u", "v"
            elif dmajor[c].tmd.group == "transport":
                ukey, vkey = "U", "V"
        # create output datatree
        dtree = xr.DataTree()
        # add datasets to output datatree
        dtree[ukey] = dsu
        dtree[vkey] = dsv
        # return the datatree
        return dtree

    @property
    def crs(self):
        """Coordinate reference system of the ``DataTree``"""
        # inherit CRS from one of the datasets
        for key, ds in self._dtree.items():
            ds = ds.to_dataset()
            return ds.tmd.crs


@xr.register_dataset_accessor("tmd")
class Dataset:
    """Accessor for extending an ``xarray.Dataset`` for tidal model data"""

    def __init__(self, ds):
        # initialize Dataset
        self._ds = ds

    def to_dataarray(self, **kwargs):
        """
        Converts ``Dataset`` to a ``DataArray`` with constituents as a dimension
        """
        kwargs.setdefault("constituents", self.constituents)
        # reduce dataset to constituents and convert to dataarray
        da = self._ds[kwargs["constituents"]].to_dataarray(dim="constituent")
        # stack constituents as the last dimension
        da = da.transpose(*da.dims[1:], da.dims[0])
        da = da.assign_coords(constituent=kwargs["constituents"])
        return da

    def assign_coords(
        self,
        x: np.ndarray,
        y: np.ndarray,
        crs: str | int | dict = 4326,
        **kwargs,
    ):
        """
        Assign new coordinates to the ``Dataset``

        Parameters
        ----------
        x: np.ndarray
            Updated x-coordinates
        y: np.ndarray
            Updated y-coordinates
        crs: str, int, or dict, default 4326 (WGS84 Latitude/Longitude)
            Coordinate reference system of coordinates
        kwargs: dict
            Keyword arguments for ``xarray.Dataset.assign_coords``

        Returns
        -------
        ds: xarray.Dataset
            ``Dataset`` with updated coordinates
        """
        # assign new coordinates to dataset
        ds = self._ds.assign_coords(dict(x=x, y=y), **kwargs)
        ds.attrs["crs"] = crs
        # return the dataset
        return ds

    def barycentric_interp(
        self,
        x: np.ndarray,
        y: np.ndarray,
        **kwargs,
    ):
        """
        Interpolate unstructured ``Datasets`` using a barycentric
        method with first or second order triangular finite elements

        Parameters
        ----------
        x: np.ndarray
            Interpolation x-coordinates
        y: np.ndarray
            Interpolation y-coordinates
        order: int
            Polynomial order of the triangular elements

            - ``1``: linear
            - ``2``: quadratic
        cutoff: int or float, default np.inf
            Maximum distance to check for elements

        Returns
        -------
        other: xarray.Dataset
            Interpolated ``Dataset``
        """
        # import barycentric interpolation functions
        from pyTMD.interpolate import (
            _to_barycentric,
            _inside_triangle,
            _shape_functions,
            _winding_number,
        )

        # get the polynomial order of the finite elements
        order = self._ds["element"].attrs.get("order", 1)
        # default order is same as the tide model
        kwargs.setdefault("order", order)

        # get cutoff distance to crop elements to bounding box
        cutoff = kwargs.get("cutoff", np.inf)
        # crop dataset to bounding box of other dataset plus buffer
        if np.isfinite(cutoff) and self.crs.is_geographic:
            # use twice the cutoff distance as a buffer
            cutoff_km = cutoff * __ureg__.parse_units("km")
            a_axis = 6378.137 * __ureg__.parse_units("km")
            buffer = 2.0 * (cutoff_km / a_axis).to(self.axis_units).magnitude
            # bounds of interpolation coordinates
            bounds = [np.min(x), np.max(x), np.min(y), np.max(y)]
            # crop dataset to bounding box of other dataset plus buffer
            ds = self.crop(bounds=bounds, buffer=buffer)
        elif np.isfinite(cutoff):
            # use twice the cutoff distance as a buffer
            cutoff_km = cutoff * __ureg__.parse_units("km")
            buffer = 2.0 * cutoff_km.to(self.axis_units).magnitude
            # bounds of interpolation coordinates
            bounds = [np.min(x), np.max(x), np.min(y), np.max(y)]
            # crop dataset to bounding box of other dataset plus buffer
            ds = self.crop(bounds=bounds, buffer=buffer)
        else:
            # copy dataset without cropping
            ds = self._ds.copy()

        # allocate for barycentric coordinates
        xi = xr.full_like(x, np.nan)
        eta = xr.full_like(x, np.nan)
        null_points = xi.isnull()
        # allocate for indices of valid elements
        element = xr.zeros_like(x, dtype="i")
        # find the valid elements and barycentric coordinates
        for i, elem in enumerate(ds.element):
            # x and y coordinates of element vertices
            x_elem = ds.x.isel(element=i).drop_vars("element")
            y_elem = ds.y.isel(element=i).drop_vars("element")
            # copy x-coordinates to not affect outside array
            xtmp = x.copy(deep=False)
            # if model is geographic:
            # check if element crosses a meridian
            if self.crs.is_geographic:
                # calculate winding number of triangle element
                # negative winding numbers are clockwise
                wind = _winding_number(x_elem, y_elem)
                # shift coordinates for meridian crossings
                if (wind < 0) & (x_elem < 0.0).any():
                    # adjust points to be 0:360
                    x_elem = x_elem.where(x_elem >= 0.0, x_elem + 360.0)
                    xtmp = xtmp.where(xtmp >= 0, xtmp + 360.0)
                elif (wind < 0) & (x_elem > 180.0).any():
                    # adjust points to be -180:180
                    x_elem = x_elem.where(x_elem <= 180.0, x_elem - 360.0)
                    xtmp = xtmp.where(xtmp <= 180.0, xtmp - 360.0)
            # convert model coordinates to barycentric
            xi_elem, eta_elem = _to_barycentric(x_elem, y_elem, xtmp, y)
            # drop dimensions
            xi_elem = xi_elem.drop_vars("vertex", errors="ignore")
            eta_elem = eta_elem.drop_vars("vertex", errors="ignore")
            # determine if points are within element and need values
            inside_element = _inside_triangle(xi_elem, eta_elem)
            # skip if nothing is inside the element
            if not np.any(inside_element & null_points):
                continue
            # save barycentric coordinates and indices
            update_element = np.logical_not(inside_element & null_points)
            xi = xi.where(update_element, xi_elem, drop=False)
            eta = eta.where(update_element, eta_elem, drop=False)
            element = element.where(update_element, i, drop=False)
            # can quit search if all interpolation points have values
            null_points = xi.isnull()
            if not null_points.any():
                break
        # get shape functions and convert to DataArray
        N = _shape_functions(xi, eta, kwargs["order"])
        beta = xr.concat(N, dim="node")
        # allocate for output dataset
        other = xr.Dataset()
        # copy attributes
        for att_name, att_val in self._ds.attrs.items():
            other.attrs[att_name] = att_val
        # iterate over variables in dataset
        for i, v in enumerate(ds.data_vars.keys()):
            # tide model variable for valid elements
            var = ds[v].isel(element=element)
            # calculate dot product over elements and nodes
            other[v] = var.dot(beta, dim="node")
            # copy variable attributes
            for att_name, att_val in self._ds[v].attrs.items():
                other[v].attrs[att_name] = att_val
        # add coordinates to output dataset
        other.coords["x"] = x
        other.coords["y"] = y
        # return the interpolated dataset
        # drop empty vertex coordinates
        return other.drop_vars("vertex", errors="ignore").compute()

    def cell_area(self, assume_spherical: bool = False):
        """
        Calculate the area of each grid cell in the ``Dataset``

        Parameters
        ----------
        assume_spherical: bool, default False
            use a spherical Earth model for the area calculation

        Returns
        -------
        area: xarray.DataArray
            Area of each grid cell
        """
        # import area functions
        from pyTMD.spatial import scale_factors

        # get PROJ4 parameters for dataset projection
        crs = self.crs.to_dict()
        # get geodetic parameters
        geod = self.crs.get_geod()
        # ellipsoid semi-major and semi-minor axes
        a_axis = geod.a
        b_axis = geod.b
        # ellipsoidal flattening
        flat = geod.f
        # first numerical eccentricity and its square
        e12 = geod.es
        ecc = np.sqrt(e12)
        # authalic radius (same area as ellipsoid)
        rad_e = np.sqrt(0.5 * (a_axis**2 + b_axis**2 * np.arctanh(ecc) / ecc))
        # coordinates and attributes for output DataArray
        coords = dict(y=self._ds.y, x=self._ds.x)
        attrs = dict(units="m^2", long_name="Grid Cell Area")
        # calculate areas based on the coordinate reference system
        if self.crs.is_geographic and assume_spherical:
            # geographic coordinates (spherical Earth model)
            gridy, _ = xr.broadcast(np.radians(self._ds.y), self._ds.x)
            # grid spacing in the x and y directions
            dx = np.abs(np.radians(self._x[1] - self._x[0]))
            dy = np.abs(np.radians(self._y[1] - self._y[0]))
            # calculate area of each grid cell
            area = (rad_e * dy) * (rad_e * dx * np.cos(gridy))
        elif self.crs.is_geographic:
            # geographic coordinates (ellipsoidal Earth model)
            gridy, _ = xr.broadcast(np.radians(self._ds.y), self._ds.x)
            # grid spacing in the x and y directions
            dx = np.abs(np.radians(self._x[1] - self._x[0]))
            dy = np.abs(np.radians(self._y[1] - self._y[0]))
            # radius of curvature in prime vertical direction (east-west)
            N = a_axis / np.sqrt(1.0 - e12 * np.sin(gridy) ** 2)
            # radius of curvature in meridional direction (north-south)
            M = a_axis * (1.0 - e12) / (1.0 - e12 * np.sin(gridy) ** 2) ** 1.5
            # calculate area of each grid cell
            area = (M * dy) * (N * np.cos(gridy) * dx)
        elif self.crs.is_projected and crs.get("proj") == "stere":
            # stereographic projection
            geodetic_crs = getattr(self.crs, "geodetic_crs", 4326)
            # get latitude and true-scale latitude
            _, lat = self.to_geographic(crs=geodetic_crs)
            lat_ts = crs.get("lat_ts", 90.0)
            # calculate scaling factors for area distortions
            ps_scale = scale_factors(lat, flat=flat, reference_latitude=lat_ts)
            # calculate scaling factors to convert from axis units to meters
            axis_units = 1.0 * __ureg__.parse_units(self.axis_units)
            axis_scale = axis_units.to(__ureg__.meter).magnitude
            # grid spacing in the x and y directions
            dx = axis_scale * np.abs(self._x[1] - self._x[0])
            dy = axis_scale * np.abs(self._y[1] - self._y[0])
            # calculate area of each grid cell
            area = ps_scale * dx * dy
        else:
            # projected coordinates (assume Cartesian)
            ny, nx = len(self._y), len(self._x)
            # calculate scaling factors for area distortions
            gridx, gridy = np.meshgrid(self._x, self._y)
            factors = pyproj.Proj(self.crs).get_factors(gridx, gridy)
            areal_scale = factors.areal_scale
            # calculate scaling factors to convert from axis units to meters
            axis_units = 1.0 * __ureg__.parse_units(self.axis_units)
            axis_scale = axis_units.to(__ureg__.meter).magnitude
            # grid spacing in the x and y directions
            dx = axis_scale * np.abs(self._x[1] - self._x[0])
            dy = axis_scale * np.abs(self._y[1] - self._y[0])
            # calculate area of each grid cell
            area = dx * dy * np.ones((ny, nx)) / areal_scale
        # return area as xarray DataArray
        return xr.DataArray(area, coords=coords, dims=["y", "x"], attrs=attrs)

    def coords_as(
        self,
        x: np.ndarray,
        y: np.ndarray,
        crs: str | int | dict = 4326,
        **kwargs,
    ):
        """
        Transform coordinates into ``DataArrays`` in the ``Dataset``
        coordinate reference system

        Parameters
        ----------
        x: np.ndarray
            Input x-coordinates
        y: np.ndarray
            Input y-coordinates
        crs: str, int, or dict, default 4326 (WGS84 Latitude/Longitude)
            Coordinate reference system of input coordinates

        Returns
        -------
        X: xarray.DataArray
            Transformed x-coordinates
        Y: xarray.DataArray
            Transformed y-coordinates
        """
        # convert coordinate reference system to that of the dataset
        # and format as xarray DataArray with appropriate dimensions
        X, Y = _coords(x, y, source_crs=crs, target_crs=self.crs, **kwargs)
        # return the transformed coordinates
        return X, Y

    def crop(
        self,
        bounds: list | tuple,
        buffer: int | float = 0,
    ):
        """
        Crop ``Dataset`` to input bounding box

        For regular grids, uses :func:`isel_bounds` so chunked / netCDF-backed
        arrays can hyperslab without materializing the full domain. For
        unstructured FE meshes, builds an element AABB mask from vertex
        coordinates and ``isel`` along ``element`` (no full-domain
        ``where``). Prefer a chunked ``open_*`` then ``tmd.crop`` (or
        ``DataTree.tmd.crop``) over driver-specific open-time windows.

        Global geographic grids handle Pacific / dateline longitude wrap by
        one or two index hyperslabs (and a continuous ``x`` rebuild), not by
        padding half the globe. ``xmin > xmax`` means an eastward wrap
        (e.g. ``170`` to ``-170``).

        Parameters
        ----------
        bounds: list, tuple
            Bounding box ``[min_x, max_x, min_y, max_y]``
        buffer: int or float, default 0
            Buffer to add to bounds for cropping
        """
        # unpack bounds and buffer
        xmin = bounds[0] - buffer
        xmax = bounds[1] + buffer
        ymin = bounds[2] - buffer
        ymax = bounds[3] + buffer
        # crop dataset to bounding box
        if self.grid_type == "unstructured":
            # Element AABB test on vertex coords, then ragged isel along
            # ``element``. Avoids ``where`` over every data variable and the
            # eager ``chunk(-1).compute()`` that forced a full FE materialise.
            # (KDTree / point-location search is separate; see pyTMD#553.)
            x = self._ds["x"]
            y = self._ds["y"]
            mask = (
                (x.max(dim="vertex") >= xmin)
                & (x.min(dim="vertex") <= xmax)
                & (y.max(dim="vertex") >= ymin)
                & (y.min(dim="vertex") <= ymax)
            )
            if hasattr(mask.data, "compute"):
                mask = mask.compute()
            idx = np.flatnonzero(np.asarray(mask.values))
            return self._ds.isel(element=idx)
        # global geographic: dateline-aware lon windows
        if self.crs.is_geographic and self.is_global:
            lon_wrap = self.crs.to_dict().get("lon_wrap", 0)
            return _crop_global_longitude(
                self._ds, xmin, xmax, ymin, ymax, lon_wrap=lon_wrap
            )
        # regional / projected regular grids
        return isel_bounds(self._ds, "x", "y", bounds=[xmin, xmax, ymin, ymax])

    def extrap_like(self, other: xr.Dataset, **kwargs):
        """
        Extrapolate missing values in ``Dataset`` using nearest-neighbors
        or inverse distance weighting

        Parameters
        ----------
        other: xarray.Dataset
            ``Dataset`` with missing values to be extrapolated
        kwargs: dict
            Keyword arguments for :py:func:`pyTMD.interpolate._query_tree`

        Returns
        -------
        other: xarray.Dataset
            ``Dataset`` with extrapolated values
        """
        # import extrapolate functions
        from pyTMD.interpolate import (
            _to_cartesian,
            _build_tree,
            _query_tree,
        )

        # get extrapolation cutoff distance
        cutoff = kwargs.get("cutoff", np.inf)
        # check if chunks are present
        if hasattr(other, "chunks") and bool(other.chunks):
            other = other.chunk(-1).compute()
        # bounds of other dataset
        bounds = [
            other.x.values.min(),
            other.x.values.max(),
            other.y.values.min(),
            other.y.values.max(),
        ]
        # crop dataset to bounding box of other dataset plus buffer
        if np.isfinite(cutoff) and self.crs.is_geographic:
            # use twice the cutoff distance as a buffer
            cutoff_km = cutoff * __ureg__.parse_units("km")
            a_axis = 6378.137 * __ureg__.parse_units("km")
            buffer = 2.0 * (cutoff_km / a_axis).to(self.axis_units).magnitude
            # crop dataset to bounding box of other dataset plus buffer
            ds = self.crop(bounds=bounds, buffer=buffer)
        elif np.isfinite(cutoff):
            # use twice the cutoff distance as a buffer
            cutoff_km = cutoff * __ureg__.parse_units("km")
            buffer = 2.0 * cutoff_km.to(self.axis_units).magnitude
            # crop dataset to bounding box of other dataset plus buffer
            ds = self.crop(bounds=bounds, buffer=buffer)
        else:
            # copy dataset without cropping
            ds = self._ds.copy()
        # check if extrapolating from grid or mesh
        if self.grid_type == "unstructured":
            # get the polynomial order of the finite elements
            order = self._ds["element"].attrs.get("order", 1)
            # reduce to data at the vertices
            nodes = [0, 1, 2] if (order == 1) else [0, 2, 4]
            ds = ds.isel(node=nodes)
        # calculate meshgrid of cropped model coordinates
        # or extract unstructured mesh x and y values
        gridy, gridx = xr.broadcast(ds.y, ds.x)
        # initialize valid mask for building tree
        valid_mask = np.zeros_like(gridx, dtype=bool)
        tree = None
        # iterate over variables in dataset
        for i, v in enumerate(other.data_vars.keys()):
            # check for missing values
            invalid = other[v].isnull()
            if not invalid.any():
                # no missing values
                continue
            # find valid values
            mask = ds[v].notnull().values
            # build tree if on the first iteration
            # or if the valid mask has changed
            if (tree is None) or (mask != valid_mask).any():
                # get indices of valid points
                valid_indices = np.nonzero(mask)
                # reduce to valid original values
                p_in = _to_cartesian(
                    gridx.values[valid_indices],
                    gridy.values[valid_indices],
                    is_geographic=self.crs.is_geographic,
                )
                # build kd-tree for valid points
                tree = _build_tree(p_in)
                # copy valid mask for next iteration
                valid_mask = np.copy(mask)
            # reduce model to valid original values
            flattened = ds[v].values[valid_indices]
            # extrapolate missing values using NN or IDW
            if other[v].ndim == 0:
                # single point extrapolation
                p_out = _to_cartesian(
                    other.x.values,
                    other.y.values,
                    is_geographic=self.crs.is_geographic,
                )
                (other[v].values,) = _query_tree(
                    tree, p_out, flattened, **kwargs
                )
            else:
                # only extrapolate invalid points
                p_out = _to_cartesian(
                    other.x.values[invalid],
                    other.y.values[invalid],
                    is_geographic=self.crs.is_geographic,
                )
                other[v].values[invalid] = _query_tree(
                    tree, p_out, flattened, **kwargs
                )
        # return xarray dataset
        return other

    def grid_interp(
        self,
        x: np.ndarray,
        y: np.ndarray,
        method="linear",
        **kwargs,
    ):
        """
        Interpolate a regular or rectilinear ``Dataset`` to new coordinates

        Parameters
        ----------
        x: np.ndarray
            Interpolation x-coordinates
        y: np.ndarray
            Interpolation y-coordinates
        method: str, default 'linear'
            Interpolation method

        Returns
        -------
        other: xarray.Dataset
            Interpolated ``Dataset``
        """
        # pad global grids along x-dimension (if necessary)
        if self.is_global:
            self._ds = self.pad(n=1)
        # verify longitudinal convention for geographic models
        if self.crs.is_geographic:
            # grid spacing in x-direction
            dx = self._x[1] - self._x[0]
            # adjust input longitudes to be consistent with model
            if (np.min(x) < 0.0) & (self._x.max() > (180.0 + dx)):
                # input points convention (-180:180)
                # tide model convention (0:360)
                x = xr.where(x < 0.0, x + 360.0, x)
            elif (np.max(x) > 180.0) & (self._x.min() < (0.0 - dx)):
                # input points convention (0:360)
                # tide model convention (-180:180)
                x = xr.where(x > 180.0, x - 360.0, x)
        # interpolate dataset using built-in xarray methods
        other = self._ds.interp(x=x, y=y, method=method)
        # return xarray dataset
        return other

    def infer(self, t: float | np.ndarray, **kwargs):
        """
        Infer minor tides from ``Dataset`` at times

        Parameters
        ----------
        t: float or np.ndarray
            Days relative to 1992-01-01T00:00:00 UTC
        kwargs: dict
            Keyword arguments for :py:func:`pyTMD.predict.infer_minor`

        Returns
        -------
        darr: xarray.DataArray
            Predicted tides
        """
        from pyTMD.predict import infer_minor

        # infer minor tides at times
        darr = infer_minor(t, self._ds, **kwargs)
        # return the inferred tides
        return darr

    def inpaint(self, **kwargs):
        """
        Inpaint over missing data in ``Dataset``

        Parameters
        ----------
        kwargs: dict
            Keyword arguments for :py:func:`pyTMD.interpolate.inpaint`

        Returns
        -------
        ds: xarray.Dataset
            Interpolated ``Dataset``
        """
        # import inpaint function
        from pyTMD.interpolate import inpaint

        # create copy of dataset
        ds = self._ds.copy()
        # inpaint each variable in the dataset
        for v in ds.data_vars.keys():
            ds[v].values = inpaint(
                self._x, self._y, self._ds[v].values, **kwargs
            )
        # return the dataset
        return ds

    def interp(
        self,
        x: np.ndarray,
        y: np.ndarray,
        **kwargs,
    ):
        """
        Interpolate ``Dataset`` to new coordinates

        Parameters
        ----------
        x: np.ndarray
            Interpolation x-coordinates
        y: np.ndarray
            Interpolation y-coordinates
        extrapolate: bool, default False
            Spatially extrapolate values beyond model domain
        cutoff: int or float, default np.inf
            Maximum distance for extrapolation
        k: int, default 1
            Number of nearest neighbors to use for extrapolation

            - ``1``: nearest neighbor (NN)
            - ``>1``: inverse distance weighting (IDW)
        power: int or float, default 2
            Power parameter for inverse distance weighting extrapolation
        workers: int, default 1
            Number of workers to use for parallel extrapolation
        kwargs: dict
            Keyword arguments for interpolation functions

        Returns
        -------
        other: xarray.Dataset
            Interpolated ``Dataset``
        """
        # set default keyword arguments
        kwargs.setdefault("method", "linear")
        kwargs.setdefault("extrapolate", False)
        # check if interpolating from a grid or mesh
        if self.grid_type == "unstructured":
            # use barycentric interpolation if data is unstructured
            other = self.barycentric_interp(x, y, **kwargs)
        else:
            # use built-in xarray interpolation methods
            other = self.grid_interp(x, y, **kwargs)
        # extrapolate using nearest-neighbors or inverse distance weighting
        if kwargs["extrapolate"]:
            # get extrapolation parameters
            # maximum distance for extrapolation
            cutoff = kwargs.get("cutoff", np.inf)
            # number of nearest neighbors to use for extrapolation
            k = kwargs.get("k", 1)
            # power parameter for inverse distance weighting extrapolation
            power = kwargs.get("power", 2)
            # number of workers to use for parallel extrapolation
            workers = kwargs.get("workers", 1)
            # use extrapolation (NN or IDW) to fill in missing values
            other = self.extrap_like(
                other, k=k, cutoff=cutoff, power=power, workers=workers
            )
        # return xarray dataset
        return other

    def minor_admittances(self, **kwargs):
        """
        Interpolate admittances for minor constituents from ``Dataset``
        :cite:p:`Doodson:1941td,Schureman:1958ty,Foreman:1989dt`

        Parameters
        ----------
        kwargs: dict
            Keyword arguments for admittance interpolation functions

        Returns
        -------
        ds: xarray.Dataset
            Merged ``Dataset`` with interpolated minor constituents
        """
        from pyTMD.predict.ocean_load import minor_admittance

        # copy dataset
        ds = self._ds.copy()
        # interpolate admittances for each species
        other = minor_admittance(self._ds, **kwargs)
        # merge datasets with minor constituents into original dataset
        ds = xr.merge(
            [ds, other],
            combine_attrs=combine_attrs,
            join="outer",
            compat="override",
        )
        # return xarray dataset
        return ds

    def node_equilibrium(self, crs: str | int | dict = 4326):
        """
        Compute the equilibrium amplitude and phase of the 18.6 year
        node tide :cite:p:`Cartwright:1971iz,Cartwright:1973em`

        Parameters
        ----------
        crs: str, int, or dict, default 4326 (WGS84 Latitude/Longitude)
            Coordinate reference system for geographic coordinates
        """
        # get name of first listed constituent
        c = self._ds.tmd.constituents[0]
        # Cartwright and Edden potential amplitude
        amajor = 0.027929  # node
        # Love numbers for long-period tides (Wahr, 1981)
        k2 = 0.299
        h2 = 0.606
        # tilt factor: response with respect to the solid earth
        gamma_2 = 1.0 + k2 - h2
        # transform model coordinates to lat/lon coordinates
        _, lat = self.to_geographic(crs=crs)
        # colatitude in radians
        th = np.radians(90.0 - lat)
        # 2nd degree Legendre polynomials
        P20 = 0.5 * (3.0 * np.cos(th) ** 2 - 1.0)
        # normalization for spherical harmonics
        dfactor = np.sqrt((4.0 + 1.0) / (4.0 * np.pi))
        # copy dataset
        ds = self._ds.copy()
        # allocate for output node equilibrium tide
        ds["node"] = xr.zeros_like(ds[c])
        # calculate equilibrium node constants
        ds["node"].values = (
            dfactor * P20 * gamma_2 * amajor * np.exp(-1j * np.pi)
        )
        ds["node"].attrs["units"] = "m"
        # return xarray dataset
        return ds

    def pad(
        self,
        n: int = 1,
        chunks=None,
    ):
        """
        Pad ``Dataset`` by repeating edge values in the x-direction

        Parameters
        ----------
        n: int, default 1
            Number of padding values to add on each side

        Returns
        -------
        ds: xarray.Dataset
            Padded ``Dataset``
        """
        # (possibly) unchunk x-coordinates and pad to wrap at meridian
        x = xr.DataArray(self._x, dims="x").pad(
            x=n, mode="reflect", reflect_type="odd"
        )
        # pad dataset and re-assign x-coordinates
        ds = self._ds.copy()
        ds = ds.pad(x=n, mode="wrap").assign_coords(x=x)
        # rechunk dataset (if specified)
        if chunks is not None:
            ds = ds.chunk(chunks)
        # return the dataset
        return ds

    def predict(self, t: float | np.ndarray, **kwargs):
        """
        Predict tides from ``Dataset`` at times

        Parameters
        ----------
        t: float or np.ndarray
            Days relative to 1992-01-01T00:00:00 UTC
        kwargs: dict
            Keyword arguments for :py:func:`pyTMD.predict.time_series`

        Returns
        -------
        darr: xarray.DataArray
            Predicted tides
        """
        from pyTMD.predict import time_series

        # predict tides at times
        darr = time_series(t, self._ds, **kwargs)
        # return the predicted tides
        return darr

    def subset(self, c: str | list):
        """
        Reduce to a subset of constituents

        Parameters
        ----------
        c: str or list
            List of constituents names
        """
        # create copy of dataset
        ds = self._ds.copy()
        # if no constituents are specified, return self
        # else return reduced dataset
        if c is None:
            return ds
        elif isinstance(c, str):
            return ds[[c]]
        else:
            return ds[c]

    def to_geographic(self, crs: str | int | dict = 4326):
        """
        Get latitude and longitude coordinates for the ``Dataset``

        Parameters
        ----------
        crs: str, int, or dict, default 4326 (WGS84 Latitude/Longitude)
            Coordinate reference system for geographic coordinates

        Returns
        -------
        lon: xarray.DataArray
            Longitude coordinates for the dataset
        lat: xarray.DataArray
            Latitude coordinates for the dataset
        """
        # check dimensions
        y, x = xr.broadcast(self._ds.y, self._ds.x)
        # allocate for output latitude and longitude
        lon = xr.zeros_like(x)
        lat = xr.zeros_like(y)
        # convert coordinates to latitude and longitude
        lon.values, lat.values = _transform(
            x, y, source_crs=self.crs, target_crs=crs, direction="FORWARD"
        )
        # return the geographic coordinates
        return lon, lat

    def transform_as(
        self,
        x: np.ndarray,
        y: np.ndarray,
        crs: str | int | dict = 4326,
        **kwargs,
    ):
        """
        Transform coordinates to/from the ``Dataset`` coordinate reference system

        Parameters
        ----------
        x: np.ndarray
            Input x-coordinates
        y: np.ndarray
            Input y-coordinates
        crs: str, int, or dict, default 4326 (WGS84 Latitude/Longitude)
            Coordinate reference system of input coordinates
        direction: str, default 'FORWARD'
            Direction of transformation

            - ``'FORWARD'``: from input crs to model crs
            - ``'INVERSE'``: from model crs to input crs

        Returns
        -------
        X: np.ndarray
            Transformed x-coordinates
        Y: np.ndarray
            Transformed y-coordinates
        """
        # convert coordinate reference system to that of the dataset
        X, Y = _transform(x, y, source_crs=crs, target_crs=self.crs, **kwargs)
        # return the transformed coordinates
        return (X, Y)

    def to_units(
        self,
        units: str,
        value: float = 1.0,
    ):
        """Convert ``Dataset`` to specified tide units

        Parameters
        ----------
        units: str
            Output units
        value: float, default 1.0
            Scaling factor to apply
        """
        # create copy of dataset
        ds = self._ds.copy()
        # convert each constituent in the dataset
        for c in self.constituents:
            ds[c] = ds[c].tmd.to_units(units, value=value)
        # return the dataset
        return ds

    def to_base_units(self):
        """Convert ``Dataset`` to base units"""
        # create copy of dataset
        ds = self._ds.copy()
        # convert each constituent in the dataset
        for c in self.constituents:
            ds[c] = ds[c].tmd.to_base_units()
        # return the dataset
        return ds

    def to_default_units(self):
        """Convert ``Dataset`` to default tide units"""
        # create copy of dataset
        ds = self._ds.copy()
        # convert each constituent in the dataset
        for c in self.constituents:
            ds[c] = ds[c].tmd.to_default_units()
        # return the dataset
        return ds

    @property
    def constituents(self):
        """List of tidal constituent names in the ``Dataset``"""
        # import constituents parser
        from pyTMD.constituents import _parse_name

        # output list of tidal constituents
        cons = []
        # parse list of model constituents
        for i, c in enumerate(self._ds.data_vars.keys()):
            try:
                cons.append(_parse_name(c))
            except ValueError:
                pass
        # return list of constituents
        return cons

    @property
    def crs(self):
        """Coordinate reference system of the ``Dataset``"""
        # return the CRS of the dataset
        # default is EPSG:4326 (WGS84)
        CRS = self._ds.attrs.get("crs", 4326)
        return pyproj.CRS.from_user_input(CRS)

    @property
    def is_global(self) -> bool:
        """Determine if ``Dataset`` covers a global domain"""
        # grid spacing in x-direction
        dx = self._x[1] - self._x[0]
        # cyclic: cell centers span 360-dx, or endpoints include both 0 and 360
        # (FES-style grids often store a duplicated meridian)
        span = self._x[-1] - self._x[0]
        cyclic = np.isclose(span, 360.0 - dx) or np.isclose(span, 360.0)
        return self.crs.is_geographic and cyclic

    @property
    def area_of_use(self) -> str | None:
        """Area of use from the ``Dataset`` CRS"""
        if self.crs.area_of_use is not None:
            return self.crs.area_of_use.name.replace(".", "").lower()

    @property
    def axis_units(self) -> str:
        """Units of the coordinate axes from the ``Dataset`` CRS"""
        return self.crs.axis_info[0].unit_name

    @property
    def grid_type(self) -> str:
        """Spatial structure of the ``Dataset``"""
        return self._ds.attrs.get("grid_type", "grid")

    @property
    def _x(self):
        """x-coordinates of the ``Dataset``"""
        return self._ds.x.values

    @property
    def _y(self):
        """y-coordinates of the ``Dataset``"""
        return self._ds.y.values


@xr.register_dataarray_accessor("tmd")
class DataArray:
    """Accessor for extending an ``xarray.DataArray`` for tidal model data"""

    def __init__(self, da):
        # initialize DataArray
        self._da = da

    @property
    def amplitude(self):
        """
        Calculate the amplitude of a tide model constituent

        Returns
        -------
        amp: xarray.DataArray
            Tide model constituent amplitude
        """
        # calculate constituent amplitude
        amp = np.hypot(self._da.real, self._da.imag)
        amp.attrs["units"] = self._da.attrs.get("units", "")
        return amp

    @property
    def phase(self):
        """
        Calculate the phase of a tide model constituent

        Returns
        -------
        ph: xarray.DataArray
            Tide model constituent phase (degrees)
        """
        # calculate constituent phase and convert to degrees
        ph = np.degrees(np.arctan2(-self._da.imag, self._da.real))
        ph = ph.where(ph >= 0, ph + 360.0, drop=False)
        ph.attrs["units"] = "degrees"
        return ph

    def find_peaks(self, **kwargs):
        """
        Find peaks in the ``DataArray``

        Parameters
        ----------
        kwargs: dict
            Keyword arguments for ``xarray.DataArray.differentiate``

        Returns
        -------
        high_peaks: xarray.DataArray
            Boolean array indicating locations of high tide peaks
        low_peaks: xarray.DataArray
            Boolean array indicating locations of low tide peaks
        """
        from pyTMD.predict import find_peaks as _find_peaks

        # find the high and low tides using a first order differentiation
        high_peaks, low_peaks = _find_peaks(self._da, dim="time", **kwargs)
        # return the peaks
        return (high_peaks, low_peaks)

    def to_units(
        self,
        units: str,
        value: float = 1.0,
    ):
        """Convert ``DataArray`` to specified tide units

        Parameters
        ----------
        units: str
            Output units
        value: float, default 1.0
            Scaling factor to apply
        """
        # convert to specified units
        conversion = value * self.quantity.to(units)
        da = self._da * conversion.magnitude
        da.attrs["units"] = str(conversion.units)
        return da

    def to_base_units(self, value=1.0):
        """Convert ``DataArray`` to base units

        Parameters
        ----------
        value: float, default 1.0
            Scaling factor to apply
        """
        # convert to base units
        conversion = value * self.quantity.to_base_units()
        da = self._da * conversion.magnitude
        da.attrs["units"] = str(conversion.units)
        return da

    def to_default_units(self, value=1.0):
        """Convert ``DataArray`` to default tide units

        Parameters
        ----------
        value: float, default 1.0
            Scaling factor to apply
        """
        # convert to default units
        default_units = _default_units.get(self.group, self.units)
        da = self.to_units(default_units, value=value)
        return da

    @property
    def units(self):
        """Units of the ``DataArray``"""
        try:
            return self._parse_units(self._units)
        except TypeError as exc:
            raise ValueError(f"Unknown units: {self._units}") from exc
        except AttributeError as exc:
            raise AttributeError("DataArray has no attribute 'units'") from exc

    @property
    def quantity(self):
        """``Pint`` Quantity of the ``DataArray``"""
        return 1.0 * self.units

    @property
    def group(self):
        """Variable group of the ``DataArray``"""
        if self.units.is_compatible_with("m"):
            return "elevation"
        elif self.units.is_compatible_with("m/s"):
            return "current"
        elif self.units.is_compatible_with("m^2/s"):
            return "transport"
        elif self.units.is_compatible_with("m/s^2"):
            return "acceleration"
        elif self.units.is_compatible_with("degrees"):
            return "angle"
        else:
            raise ValueError(f"Unknown unit group: {self._units}")

    @staticmethod
    def _parse_units(units: str):
        """
        Convert units attributes to ``pint`` units
        """
        # fix the exponent notation in units string
        units = re.sub(
            r"(\w)([-]?\d+)",
            lambda m: m.group(1) + r"^" + m.group(2),
            units,
            flags=re.IGNORECASE,
        )
        # parse units string using pint
        return __ureg__.parse_units(units.lower())

    @property
    def _units(self):
        """Units attribute of the ``DataArray`` as a string"""
        return self._da.attrs.get("units")

    @property
    def _has_compatible_units(self):
        """Tests that units are compatible with known groups"""
        try:
            unit_group = self.group
        except (TypeError, ValueError, AttributeError) as exc:
            return False
        else:
            return True


def _longitude_model_range(
    xvals: np.ndarray, lon_wrap: int | float = 0
) -> tuple[float, float]:
    """
    Infer the primary longitude range of a global grid.

    Prefers the coordinate extent; falls back to ``lon_wrap``
    (``180`` → 0–360, else −180–180).
    """
    xvals = np.asarray(xvals, dtype=float)
    dx = float(xvals[1] - xvals[0]) if xvals.size > 1 else 1.0
    xmin = float(np.min(xvals))
    xmax = float(np.max(xvals))
    # 0–360 style (FES / GOT / ATLAS)
    if (xmin >= -dx) and (xmax > (180.0 + dx)):
        return 0.0, 360.0
    # −180–180 style
    if xmin < -dx:
        return -180.0, 180.0
    if lon_wrap == 180:
        return 0.0, 360.0
    return -180.0, 180.0


def _shift_longitude_into(
    lon: float, model_lo: float, model_hi: float
) -> float:
    """Shift ``lon`` by multiples of 360 into ``[model_lo, model_hi)``."""
    period = model_hi - model_lo
    while lon < model_lo:
        lon += period
    while lon >= model_hi:
        lon -= period
    return lon


def _longitude_model_segments(
    xmin: float,
    xmax: float,
    model_lo: float,
    model_hi: float,
    grid_xmin: float | None = None,
    grid_xmax: float | None = None,
) -> list[tuple[float, float]]:
    """
    Map a longitude window to one or two segments in model coordinates.

    Travel is eastward from ``xmin`` to ``xmax``. If ``xmin > xmax``, the
    window wraps (Pacific dateline short arc).

    For non-wrapping windows whose span is ≥ one period (e.g. a huge
    ``buffer``), intersect with the native grid extent instead of forcing a
    modular full-globe read — matching historical ``where`` crop behavior.
    """
    period = model_hi - model_lo
    if xmin <= xmax:
        span = xmax - xmin
        # oversized linear box: clip to native grid (not modular full globe)
        if span >= period - 1e-12:
            x0 = model_lo if grid_xmin is None else grid_xmin
            x1 = model_hi if grid_xmax is None else grid_xmax
            lo = max(xmin, x0)
            hi = min(xmax, x1)
            if lo > hi:
                return []
            return [(lo, hi)]
        start = xmin
    else:
        # eastward wrap: e.g. 170 → -170 is a 20° arc, not the long way
        span = (xmax + period) - xmin
        start = xmin
        if span >= period - 1e-12:
            return [(model_lo, model_hi)]
    start_m = _shift_longitude_into(start, model_lo, model_hi)
    end_m = start_m + span
    if end_m <= model_hi + 1e-12:
        return [(start_m, end_m)]
    return [(start_m, model_hi), (model_lo, end_m - period)]


def _crop_global_longitude(
    ds: xr.Dataset,
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
    lon_wrap: int | float = 0,
) -> xr.Dataset:
    """
    Crop a global geographic grid with dateline-safe longitude windows.

    Uses at most two :func:`isel_bounds` hyperslabs and rebuilds a
    continuous ``x`` coordinate for the eastward path (no half-globe pad).
    """
    xvals = np.asarray(ds.x.values)
    model_lo, model_hi = _longitude_model_range(xvals, lon_wrap=lon_wrap)
    period = model_hi - model_lo
    segments = _longitude_model_segments(
        xmin,
        xmax,
        model_lo,
        model_hi,
        grid_xmin=float(np.min(xvals)),
        grid_xmax=float(np.max(xvals)),
    )
    # Oversized non-wrapping buffer: native-grid intersect only (no ±360 rebuild)
    if xmin <= xmax and (xmax - xmin) >= period - 1e-12 and len(segments) == 1:
        lo, hi = segments[0]
        return isel_bounds(ds, "x", "y", bounds=[lo, hi, ymin, ymax])
    pieces: list[xr.Dataset] = []
    x_cursor: float | None = None
    # FES-style grids store both period endpoints (0 and 360, or −180 and 180).
    # A wrap split then selects model_hi on the west piece and model_lo on the
    # east piece — the same meridian twice. Drop the closing sample so the
    # east piece stays one cell east of x_cursor (no extra +360 shift).
    drop_closing_meridian = (
        len(segments) > 1
        and np.isclose(float(np.min(xvals)), model_lo)
        and np.isclose(float(np.max(xvals)), model_hi)
    )
    for lo, hi in segments:
        part = isel_bounds(ds, "x", "y", bounds=[lo, hi, ymin, ymax])
        if part.sizes.get("x", 0) == 0:
            continue
        if (
            drop_closing_meridian
            and x_cursor is None
            and np.isclose(float(part.x.values[-1]), model_hi)
        ):
            part = part.isel(x=slice(0, -1))
            if part.sizes.get("x", 0) == 0:
                continue
        x = np.asarray(part.x.values, dtype=float).copy()
        if x_cursor is None:
            # place the first segment near the requested start longitude
            while (np.min(x) - xmin) > 180.0:
                x -= 360.0
            while (xmin - np.min(x)) > 180.0:
                x += 360.0
        else:
            # keep subsequent segments east of the previous piece
            while np.min(x) <= x_cursor + 1e-12:
                x += 360.0
        pieces.append(part.assign_coords(x=("x", x)))
        x_cursor = float(np.max(x))
    if not pieces:
        return isel_bounds(ds, "x", "y", bounds=[xmin, xmax, ymin, ymax])
    if len(pieces) == 1:
        return pieces[0]
    return xr.concat(pieces, dim="x")


def isel_bounds(
    ds: xr.Dataset,
    x_name: str,
    y_name: str,
    bounds: list | tuple,
) -> xr.Dataset:
    """
    Index-select a gridded ``Dataset`` to ``bounds``.

    Uses 1-D coordinate values only (cheap), then ``isel`` so netCDF/dask
    backends can hyperslab amplitude/phase or complex fields instead of
    reading the globe. Callers that need a halo should expand ``bounds``
    themselves (``Dataset.tmd.crop`` applies ``buffer`` before calling).

    Parameters
    ----------
    ds: xarray.Dataset
        Dataset with coordinates ``x_name`` / ``y_name``
    x_name, y_name: str
        Coordinate names (e.g. ``x``/``y`` or file ``lon``/``lat``)
    bounds: list or tuple
        Bounding box ``[xmin, xmax, ymin, ymax]`` in coordinate units

    Returns
    -------
    ds: xarray.Dataset
        Windowed dataset. If no cells fall in bounds, returns empty
        slices along ``x_name`` / ``y_name`` (same spirit as
        ``where(..., drop=True)``).
    """
    xmin, xmax, ymin, ymax = map(float, bounds)
    # load coordinate vectors only
    xvals = np.asarray(ds[x_name].values)
    yvals = np.asarray(ds[y_name].values)
    xind = np.flatnonzero((xvals >= xmin) & (xvals <= xmax))
    yind = np.flatnonzero((yvals >= ymin) & (yvals <= ymax))
    if xind.size == 0 or yind.size == 0:
        return ds.isel(
            {
                x_name: slice(0, 0),
                y_name: slice(0, 0),
            }
        )
    return ds.isel(
        {
            x_name: slice(int(xind[0]), int(xind[-1]) + 1),
            y_name: slice(int(yind[0]), int(yind[-1]) + 1),
        }
    )


def combine_attrs(
    attrs_list: list[dict],
    context: str | None,
    **kwargs,
) -> dict:
    """
    Combine attributes from multiple datasets into a single dictionary
    merging conflicting values into a list

    Parameters
    ----------
    attrs_list: list of dict
        List of attribute dictionaries from multiple datasets
    context: str
        Context for the attributes being combined
    skip_keys: list of str, default ["units"]
        List of attribute keys to skip from comparison

    Returns
    -------
    result: dict
        Combined attributes dictionary
    """
    # set default keyword arguments
    skip_keys = kwargs.get("skip_keys", ["units"])
    # return an empty dictionary when no attributes are provided
    if not attrs_list:
        return {}
    # initialize combined attributes with the first dictionary in the list
    result = attrs_list[0].copy()
    append_keys = set()
    # for each attribute key, check if values are equivalent
    for attrs in attrs_list:
        for key, value in attrs.items():
            # skip keys that have already been identified as conflicts
            # and keys that should be skipped from comparison
            if key in append_keys or key in skip_keys:
                continue
            # check if the attribute values are equivalent
            if not equivalent_attrs(result.get(key), value):
                append_keys.add(key)
    # combine conflicting attributes into lists
    for key in append_keys:
        # build list of values for this key across all datasets
        combined_values = []
        for attrs in attrs_list:
            # check if the key is present
            # if a list or tuple: extend the combined values
            # if a single value: append to the combined values
            if key in attrs and isinstance(attrs[key], (list, tuple)):
                combined_values.extend(attrs[key])
            elif key in attrs:
                combined_values.append(attrs[key])
        # clean up combined results: removing duplicates and null values
        result[key] = sorted(set(filter(None, combined_values)))
        # if only one unique value remains, simplify back to a single value
        if len(result[key]) == 1:
            result[key] = result[key].pop()
    # return the combined attributes
    return result


def equivalent_attrs(a: Any, b: Any) -> bool:
    """
    Check if two attribute values are equivalent (ignoring case for strings)

    Adapted from ``xarray.structure.merge.equivalent_attrs``

    Parameters
    ----------
    a: Any
        First attribute value
    b: Any
        Second attribute value
    """
    # if both attributes are strings, compare them case-insensitively
    if isinstance(a, str) and isinstance(b, str):
        return equivalent(a.casefold(), b.casefold())
    # otherwise, compare the attributes directly
    # exceptions would indicate comparison is ambiguous
    try:
        return equivalent(a, b)
    except (TypeError, ValueError):
        return False


def register_datatree_subaccessor(name):
    """Register a custom subaccessor on ``DataTree`` objects

    Parameters
    ----------
    name: str
        Name of the subaccessor
    """
    return xr.core.extensions._register_accessor(name, DataTree)


def register_dataset_subaccessor(name):
    """Register a custom subaccessor on ``Dataset`` objects

    Parameters
    ----------
    name: str
        Name of the subaccessor
    """
    return xr.core.extensions._register_accessor(name, Dataset)


def register_dataarray_subaccessor(name):
    """Register a custom subaccessor on ``DataArray`` objects

    Parameters
    ----------
    name: str
        Name of the subaccessor
    """
    return xr.core.extensions._register_accessor(name, DataArray)


def _transform(
    i1: np.ndarray,
    i2: np.ndarray,
    source_crs: str | int | dict = 4326,
    target_crs: str | int | dict = None,
    **kwargs,
):
    """
    Transform coordinates to/from the dataset coordinate reference system

    Parameters
    ----------
    i1: np.ndarray
        Input x-coordinates
    i2: np.ndarray
        Input y-coordinates
    source_crs: str, int, or dict, default 4326 (WGS84 Latitude/Longitude)
        Coordinate reference system of input coordinates
    target_crs: str, int, or dict, default None
        Coordinate reference system of output coordinates
    direction: str, default 'FORWARD'
        Direction of transformation

        - ``'FORWARD'``: from input crs to model crs
        - ``'INVERSE'``: from model crs to input crs

    Returns
    -------
    o1: np.ndarray
        Transformed x-coordinates
    o2: np.ndarray
        Transformed y-coordinates
    """
    # set the direction of the transformation
    kwargs.setdefault("direction", "FORWARD")
    if kwargs["direction"] not in ("FORWARD", "INVERSE", "IDENT"):
        raise ValueError("Invalid transformation direction")
    # get the coordinate reference system and transform
    source_crs = pyproj.CRS.from_user_input(source_crs)
    target_crs = pyproj.CRS.from_user_input(target_crs)
    transformer = pyproj.Transformer.from_crs(
        source_crs, target_crs, always_xy=True
    )
    # convert coordinate reference system
    o1, o2 = transformer.transform(i1, i2, **kwargs)
    # return the transformed coordinates
    return (o1, o2)


def _coords(
    x: np.ndarray,
    y: np.ndarray,
    source_crs: str | int | dict = 4326,
    target_crs: str | int | dict = None,
    **kwargs,
):
    """
    Transform coordinates into DataArrays in a new
    coordinate reference system

    Parameters
    ----------
    x: np.ndarray
        Input x-coordinates
    y: np.ndarray
        Input y-coordinates
    source_crs: str, int, or dict, default 4326 (WGS84 Latitude/Longitude)
        Coordinate reference system of input coordinates
    target_crs: str, int, or dict, default None
        Coordinate reference system of output coordinates
    type: str or None, default None
        Coordinate data type

        If not provided: must specify ``time`` parameter to auto-detect

        - ``None``: determined from input variable dimensions
        - ``'trajectory'``: drift buoys or satellite/airborne altimetry
        - ``'grid'``: spatial grids or images
        - ``'time series'``: time series at a single point
    time: np.ndarray or None, default None
        Time variable for determining coordinate data type

    Returns
    -------
    X: xarray.DataArray
        Transformed x-coordinates
    Y: xarray.DataArray
        Transformed y-coordinates
    """
    from pyTMD.spatial import data_type

    # set default keyword arguments
    kwargs.setdefault("type", None)
    kwargs.setdefault("time", None)
    # determine coordinate data type if possible
    if (np.ndim(x) == 0) and (np.ndim(y) == 0):
        coord_type = "time series"
    elif kwargs["type"] is None and kwargs["time"] is None:
        # must provide time variable to determine data type
        raise ValueError("Must provide time parameter to determine type")
    elif kwargs["type"] is None:
        # auto-detect coordinate data type from input variable dimensions
        coord_type = data_type(x, y, np.ravel(kwargs["time"]))
    elif isinstance(kwargs["type"], str) and kwargs["type"].lower() == "drift":
        # handle legacy 'drift' type for trajectory data
        coord_type = "trajectory"
    elif isinstance(kwargs["type"], str):
        # use provided coordinate data type
        # and verify that it is lowercase
        coord_type = kwargs.get("type").lower()
    else:
        raise ValueError("Coordinate data type must be a string")
    # convert coordinates to a new coordinate reference system
    if (coord_type == "grid") and (np.size(x) != np.size(y)):
        gridx, gridy = np.meshgrid(x, y)
        mx, my = _transform(
            gridx,
            gridy,
            source_crs=source_crs,
            target_crs=target_crs,
            direction="FORWARD",
        )
    else:
        mx, my = _transform(
            x,
            y,
            source_crs=source_crs,
            target_crs=target_crs,
            direction="FORWARD",
        )
    # convert to xarray DataArray with appropriate dimensions
    if (np.ndim(x) == 0) and (np.ndim(y) == 0):
        X = xr.DataArray(mx)
        Y = xr.DataArray(my)
    elif coord_type == "grid":
        X = xr.DataArray(mx, dims=("y", "x"))
        Y = xr.DataArray(my, dims=("y", "x"))
    elif coord_type == "trajectory":
        X = xr.DataArray(mx, dims=("time"))
        Y = xr.DataArray(my, dims=("time"))
    elif coord_type == "time series":
        X = xr.DataArray(mx, dims=("station"))
        Y = xr.DataArray(my, dims=("station"))
    else:
        raise ValueError(f"Unknown coordinate data type: {coord_type}")
    # set attributes for coordinate type
    X.attrs["type"] = coord_type
    Y.attrs["type"] = coord_type
    # return the transformed coordinates
    return (X, Y)
