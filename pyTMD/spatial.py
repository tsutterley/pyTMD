#!/usr/bin/env python
"""
spatial.py
Written by Tyler Sutterley (08/2026)

Spatial transformation routines

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 08/2026: use np.hypot to calculate euclidean distances
    Updated 06/2026: use item() to extract scalars from 0-dimensional arrays
        standardize use of lambda (lmda) to denote longitudes
        drift type renamed to trajectory to better fit CF conventions
    Updated 05/2026: moved datum ellipsoidal parameters to earth module
    Updated 04/2026: updated scale factors to add case where reference
        latitude is at the pole
    Updated 02/2026: add function to compute geocentric latitudes
    Updated 12/2025: add units to input and output variables in docstrings
    Updated 08/2025: convert angles with numpy radians and degrees functions
    Updated 03/2025: add more ellipsoidal parameters to datum class
    Updated 02/2025: major refactor to move io routines out of this module
    Updated 12/2024: add latitude and longitude as potential dimension names
    Updated 11/2024: added function to calculate the altitude and azimuth
    Updated 09/2024: deprecation fix case where an array is output to scalars
    Updated 08/2024: changed from 'geotiff' to 'GTiff' and 'cog' formats
        added functions to convert to and from East-North-Up coordinates
    Updated 07/2024: added functions to convert to and from DMS
    Updated 06/2024: added function to write parquet files with metadata
    Updated 05/2024: added function to read from parquet files
        allowing for decoding of the geometry column from WKB
        deprecation update to use exceptions with osgeo osr
    Updated 04/2024: use timescale for temporal operations
        use wrapper to importlib for optional dependencies
    Updated 03/2024: can calculate polar stereographic distortion for distances
    Updated 02/2024: changed class name for ellipsoid parameters to datum
    Updated 10/2023: can read from netCDF4 or HDF5 variable groups
        apply no formatting to columns in ascii file output
    Updated 09/2023: add function to invert field mapping keys and values
        use datetime64[ns] for parsing dates from ascii files
    Updated 08/2023: remove possible crs variables from output fields list
        place PyYAML behind try/except statement to reduce build size
    Updated 05/2023: use datetime parser within pyTMD.time module
    Updated 04/2023: copy inputs in cartesian to not modify original arrays
        added iterative methods for converting from cartesian to geodetic
        allow netCDF4 and HDF5 outputs to be appended to existing files
        using pathlib to define and expand paths
    Updated 03/2023: add basic variable typing to function inputs
    Updated 02/2023: use outputs from constants class for WGS84 parameters
        include more possible dimension names for gridded and drift outputs
    Updated 01/2023: added default field mapping for reading from netCDF4/HDF5
        split netCDF4 output into separate functions for grid and drift types
    Updated 12/2022: add software information to output HDF5 and netCDF4
    Updated 11/2022: place some imports within try/except statements
        added encoding for writing ascii files
        use f-strings for formatting verbose or ascii output
    Updated 10/2022: added datetime parser for ascii time columns
    Updated 06/2022: added field_mapping options to netCDF4 and HDF5 reads
        added from_file wrapper function to read from particular formats
    Updated 04/2022: add option to reduce input GDAL raster datasets
        updated docstrings to numpy documentation format
        use gzip virtual filesystem for reading compressed geotiffs
        include utf-8 encoding in reads to be windows compliant
    Updated 03/2022: add option to specify output GDAL driver
    Updated 01/2022: use iteration breaks in convert ellipsoid function
        remove fill_value attribute after creating netCDF4 and HDF5 variables
    Updated 11/2021: added empty cases to netCDF4 and HDF5 output for crs
        try to get grid mapping attributes from netCDF4 and HDF5
    Updated 10/2021: add pole case in stereographic area scale calculation
        using python logging for handling verbose output
    Updated 09/2021: can calculate height differences between ellipsoids
    Updated 07/2021: added function for determining input variable type
    Updated 03/2021: added polar stereographic area scale calculation
        add routines for converting to and from cartesian coordinates
        replaced numpy bool/int to prevent deprecation warnings
    Updated 01/2021: add streaming from bytes for ascii, netCDF4, HDF5, geotiff
        set default time for geotiff files to 0
    Updated 12/2020: added module for converting ellipsoids
    Updated 11/2020: output data as masked arrays if containing fill values
        add functions to read from and write to geotiff image formats
    Written 09/2020
"""

from __future__ import annotations

import warnings
import numpy as np

# get ellipsoid parameters
from pyTMD.earth import _wgs84

__all__ = [
    "data_type",
    "convert_ellipsoid",
    "compute_delta_h",
    "wrap_longitudes",
    "to_dms",
    "from_dms",
    "to_cartesian",
    "to_sphere",
    "to_geodetic",
    "_moritz_iterative",
    "_bowring_iterative",
    "_zhu_closed_form",
    "to_ENU",
    "from_ENU",
    "to_horizontal",
    "to_zenith",
    "geocentric_latitude",
    "scale_factors",
]


def data_type(
    x: np.ndarray,
    y: np.ndarray,
    t: np.ndarray,
) -> str:
    """
    Determines input data type based on variable dimensions

    Parameters
    ----------
    x: np.ndarray
        x-dimension coordinates
    y: np.ndarray
        y-dimension coordinates
    t: np.ndarray
        Time coordinates

    Returns
    -------
    String denoting input data type

        - ``'time series'``
        - ``'trajectory'``
        - ``'grid'``
    """
    xsize = np.size(x)
    ysize = np.size(y)
    tsize = np.size(t)
    if (xsize == 1) and (ysize == 1) and (tsize >= 1):
        return "time series"
    elif (xsize == ysize) & (xsize == tsize):
        return "trajectory"
    elif (np.ndim(x) > 1) & (xsize == ysize):
        return "grid"
    elif xsize != ysize:
        return "grid"
    else:
        raise ValueError("Unknown data type")


def convert_ellipsoid(
    lat1: np.ndarray,
    h1: np.ndarray,
    a1: float,
    f1: float,
    a2: float,
    f2: float,
    eps: float = 1e-12,
    itmax: int = 10,
):
    """
    Convert latitudes and heights to a different ellipsoid using
    Newton-Raphson :cite:p:`Meeus:1991vh`

    Parameters
    ----------
    lat1: np.ndarray
        Latitude of input ellipsoid (degrees)
    h1: np.ndarray
        Height above input ellipsoid (meters)
    a1: float
        Semi-major axis of input ellipsoid
    f1: float
        Flattening of input ellipsoid
    a2: float
        Semi-major axis of output ellipsoid
    f2: float
        Flattening of output ellipsoid
    eps: float, default 1e-12
        Tolerance to prevent division by small numbers and
        to determine convergence
    itmax: int, default 10
        Maximum number of iterations to use in Newton-Raphson

    Returns
    -------
    lat2: np.ndarray
        Latitude of output ellipsoid (degrees)
    h2: np.ndarray
        Height above output ellipsoid (meters)
    """
    if len(lat1) != len(h1):
        raise ValueError("lat and h have incompatible dimensions")
    # semiminor axis of input and output ellipsoid
    b1 = (1.0 - f1) * a1
    b2 = (1.0 - f2) * a2
    # initialize output arrays
    npts = len(lat1)
    lat2 = np.zeros(npts)
    h2 = np.zeros(npts)
    # for each point
    for N in range(npts):
        # force lat1 into range -90 <= lat1 <= 90
        if np.abs(lat1[N]) > 90.0:
            lat1[N] = np.sign(lat1[N]) * 90.0
        # handle special case near the equator
        # lat2 = lat1 (latitudes congruent)
        # h2 = h1 + a1 - a2
        if np.abs(lat1[N]) < eps:
            lat2[N] = np.copy(lat1[N])
            h2[N] = h1[N] + a1 - a2
        # handle special case near the poles
        # lat2 = lat1 (latitudes congruent)
        # h2 = h1 + b1 - b2
        elif (90.0 - np.abs(lat1[N])) < eps:
            lat2[N] = np.copy(lat1[N])
            h2[N] = h1[N] + b1 - b2
        # handle case if latitude is within 45 degrees of equator
        elif np.abs(lat1[N]) <= 45:
            # convert lat1 to radians
            lat1r = np.radians(lat1[N])
            sinlat1 = np.sin(lat1r)
            coslat1 = np.cos(lat1r)
            # prevent division by very small numbers
            coslat1 = np.copy(eps) if (coslat1 < eps) else coslat1
            # calculate tangent
            tanlat1 = sinlat1 / coslat1
            u1 = np.arctan(b1 / a1 * tanlat1)
            hpr1sin = b1 * np.sin(u1) + h1[N] * sinlat1
            hpr1cos = a1 * np.cos(u1) + h1[N] * coslat1
            # set initial value for u2
            u2 = np.copy(u1)
            # setup constants
            k0 = b2 * b2 - a2 * a2
            k1 = a2 * hpr1cos
            k2 = b2 * hpr1sin
            # perform newton-raphson iteration to solve for u2
            # cos(u2) will not be close to zero since abs(lat1) <= 45
            for i in range(0, itmax + 1):
                cosu2 = np.cos(u2)
                fu2 = k0 * np.sin(u2) + k1 * np.tan(u2) - k2
                fu2p = k0 * cosu2 + k1 / (cosu2 * cosu2)
                if np.abs(fu2p) < eps:
                    break
                else:
                    delta = fu2 / fu2p
                    u2 -= delta
                    if np.abs(delta) < eps:
                        break
            # convert latitude to degrees and verify values between +/- 90
            lat2r = np.arctan(a2 / b2 * np.tan(u2))
            lat2[N] = np.clip(np.degrees(lat2r), -90.0, 90.0)
            # calculate height
            h2[N] = (hpr1cos - a2 * np.cos(u2)) / np.cos(lat2r)
        # handle final case where latitudes are between 45 degrees and pole
        else:
            # convert lat1 to radians
            lat1r = np.radians(lat1[N])
            sinlat1 = np.sin(lat1r)
            coslat1 = np.cos(lat1r)
            # prevent division by very small numbers
            coslat1 = np.copy(eps) if (coslat1 < eps) else coslat1
            # calculate tangent
            tanlat1 = sinlat1 / coslat1
            u1 = np.arctan(b1 / a1 * tanlat1)
            hpr1sin = b1 * np.sin(u1) + h1[N] * sinlat1
            hpr1cos = a1 * np.cos(u1) + h1[N] * coslat1
            # set initial value for u2
            u2 = np.copy(u1)
            # setup constants
            k0 = a2 * a2 - b2 * b2
            k1 = b2 * hpr1sin
            k2 = a2 * hpr1cos
            # perform newton-raphson iteration to solve for u2
            # sin(u2) will not be close to zero since abs(lat1) > 45
            for i in range(0, itmax + 1):
                sinu2 = np.sin(u2)
                fu2 = k0 * np.cos(u2) + k1 / np.tan(u2) - k2
                fu2p = -1 * (k0 * sinu2 + k1 / (sinu2 * sinu2))
                if np.abs(fu2p) < eps:
                    break
                else:
                    delta = fu2 / fu2p
                    u2 -= delta
                    if np.abs(delta) < eps:
                        break
            # convert latitude to degrees and verify values between +/- 90
            lat2r = np.arctan(a2 / b2 * np.tan(u2))
            lat2[N] = np.clip(np.degrees(lat2r), -90.0, 90.0)
            # calculate height
            h2[N] = (hpr1sin - b2 * np.sin(u2)) / np.sin(lat2r)

    # return the latitude and height
    return (lat2, h2)


def compute_delta_h(
    lat: np.ndarray,
    a1: float,
    f1: float,
    a2: float,
    f2: float,
):
    """
    Compute difference in elevation for two ellipsoids at a given
    latitude using a simplified empirical relation :cite:p:`Meeus:1991vh`

    Parameters
    ----------
    lat: np.ndarray
        Latitudes (degrees north)
    a1: float
        Semi-major axis of input ellipsoid
    f1: float
        Flattening of input ellipsoid
    a2: float
        Semi-major axis of output ellipsoid
    f2: float
        Flattening of output ellipsoid

    Returns
    -------
    delta_h: np.ndarray
        Difference in elevation for two ellipsoids
    """
    # force latitudes to be within -90 to 90 and convert to radians
    phi = np.radians(np.clip(lat, -90.0, 90.0))
    # semi-minor axis of input and output ellipsoid
    b1 = (1.0 - f1) * a1
    b2 = (1.0 - f2) * a2
    # compute differences in semi-major and semi-minor axes
    delta_a = a2 - a1
    delta_b = b2 - b1
    # compute differences between ellipsoids
    # delta_h = -(delta_a * cos(phi)^2 + delta_b * sin(phi)^2)
    delta_h = -(delta_a * np.cos(phi) ** 2 + delta_b * np.sin(phi) ** 2)
    return delta_h


def wrap_longitudes(lon: float | np.ndarray):
    """
    Wraps longitudes to range from -180 to +180

    Parameters
    ----------
    lon: float or np.ndarray
        Longitude (degrees east)
    """
    lmda = np.arctan2(np.sin(np.radians(lon)), np.cos(np.radians(lon)))
    # convert longitudes from radians to degrees
    return np.degrees(lmda)


def to_dms(d: np.ndarray):
    """
    Convert decimal degrees to degrees, minutes and seconds

    Parameters
    ----------
    d: np.ndarray
        Angle (decimal degrees)

    Returns
    -------
    degree: np.ndarray
        Degrees
    minute: np.ndarray
        Minutes (arcminutes)
    second: np.ndarray
        Seconds (arcseconds)
    """
    sign = np.sign(d)
    minute, second = np.divmod(np.abs(d) * 3600.0, 60.0)
    degree, minute = np.divmod(minute, 60.0)
    return (sign * degree, minute, second)


def from_dms(
    degree: np.ndarray,
    minute: np.ndarray,
    second: np.ndarray,
):
    """
    Convert degrees, minutes and seconds to decimal degrees

    Parameters
    ----------
    degree: np.ndarray
        Degrees
    minute: np.ndarray
        Minutes (arcminutes)
    second: np.ndarray
        Seconds (arcseconds)

    Returns
    -------
    d: np.ndarray
        Angle (decimal degrees)
    """
    sign = np.sign(degree)
    d = np.abs(degree) + minute / 60.0 + second / 3600.0
    return sign * d


def to_cartesian(
    lon: np.ndarray,
    lat: np.ndarray,
    h: float | np.ndarray = 0.0,
    a_axis: float = _wgs84.a_axis,
    flat: float = _wgs84.flat,
):
    """
    Converts geodetic coordinates to Cartesian coordinates

    Parameters
    ----------
    lon: np.ndarray
        Longitude (degrees east)
    lat: np.ndarray
        Latitude (degrees north)
    h: float or np.ndarray, default 0.0
        Height above ellipsoid (or sphere)
    a_axis: float, default 6378137.0
        Semi-major axis of the ellipsoid

        For spherical coordinates set to radius of the Earth
    flat: float, default 1.0/298.257223563
        Ellipsoidal flattening

        For spherical coordinates set to 0

    Returns
    -------
    x: np.ndarray
        Cartesian x-coordinates (meters)
    y: np.ndarray
        Cartesian y-coordinates (meters)
    z: np.ndarray
        Cartesian z-coordinates (meters)
    """
    # verify axes and copy to not modify inputs
    singular_values = np.ndim(lon) == 0
    lon = np.atleast_1d(np.copy(lon)).astype(np.float64)
    lat = np.atleast_1d(np.copy(lat)).astype(np.float64)
    # fix coordinates to be 0:360
    lon = np.where(lon < 0, lon + 360, lon)
    # Linear eccentricity and first numerical eccentricity
    lin_ecc = np.sqrt((2.0 * flat - flat**2) * a_axis**2)
    ecc1 = lin_ecc / a_axis
    # convert from geodetic latitude to geocentric latitude
    # geodetic latitude in radians
    latitude_geodetic_rad = np.radians(lat)
    # prime vertical radius of curvature
    N = a_axis / np.sqrt(1.0 - ecc1**2.0 * np.sin(latitude_geodetic_rad) ** 2.0)
    # calculate X, Y and Z from geodetic latitude and longitude
    x = (N + h) * np.cos(latitude_geodetic_rad) * np.cos(np.radians(lon))
    y = (N + h) * np.cos(latitude_geodetic_rad) * np.sin(np.radians(lon))
    z = (N * (1.0 - ecc1**2.0) + h) * np.sin(latitude_geodetic_rad)
    # return the cartesian coordinates
    # flattened to singular values if necessary
    if singular_values:
        return (x.item(), y.item(), z.item())
    else:
        return (x, y, z)


def to_sphere(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
):
    """
    Convert from cartesian coordinates to spherical coordinates

    Parameters
    ----------
    x, np.ndarray
        Cartesian x-coordinates (meters)
    y, np.ndarray
        Cartesian y-coordinates (meters)
    z, np.ndarray
        Cartesian z-coordinates (meters)

    Returns
    -------
    lon: np.ndarray
        Longitude (degrees east)
    lat: np.ndarray
        Latitude (degrees north)
    rad: np.ndarray
        Radius (meters)
    """
    # verify axes and copy to not modify inputs
    singular_values = np.ndim(x) == 0
    x = np.atleast_1d(np.copy(x)).astype(np.float64)
    y = np.atleast_1d(np.copy(y)).astype(np.float64)
    z = np.atleast_1d(np.copy(z)).astype(np.float64)
    # calculate radius
    rad = np.sqrt(x**2.0 + y**2.0 + z**2.0)
    # calculate angular coordinates
    # lmda: azimuthal angle
    lmda = np.arctan2(y, x)
    # th: polar angle
    th = np.arccos(z / rad)
    # convert to degrees and fix to 0:360
    lon = np.degrees(lmda)
    lon = np.where(lon < 0.0, lon + 360.0, lon)
    # convert to degrees and fix to -90:90
    lat = 90.0 - np.degrees(th)
    np.clip(lat, -90, 90, out=lat)
    # return longitude, latitude and radius
    # flattened to singular values if necessary
    if singular_values:
        return (lon.item(), lat.item(), rad.item())
    else:
        return (lon, lat, rad)


def to_geodetic(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    a_axis: float = _wgs84.a_axis,
    flat: float = _wgs84.flat,
    method: str = "bowring",
    eps: float = np.finfo(np.float64).eps,
    iterations: int = 10,
):
    """
    Convert from cartesian coordinates to geodetic coordinates
    using either iterative or closed-form methods

    Parameters
    ----------
    x, np.ndarray
        Cartesian x-coordinates (meters)
    y, np.ndarray
        Cartesian y-coordinates (meters)
    z, np.ndarray
        Cartesian z-coordinates (meters)
    a_axis: float, default 6378137.0
        Semi-major axis of the ellipsoid
    flat: float, default 1.0/298.257223563
        Ellipsoidal flattening
    method: str, default 'bowring'
        Method to use for conversion

            - ``'moritz'``: iterative solution
            - ``'bowring'``: iterative solution
            - ``'zhu'``: closed-form solution
    eps: float, default np.finfo(np.float64).eps
        Tolerance for iterative methods
    iterations: int, default 10
        Maximum number of iterations

    Returns
    -------
    lon: np.ndarray
        Longitude (degrees east)
    lat: np.ndarray
        Latitude (degrees north)
    h: np.ndarray
        Height above ellipsoid (meters)
    """
    # verify axes and copy to not modify inputs
    singular_values = np.ndim(x) == 0
    x = np.atleast_1d(np.copy(x)).astype(np.float64)
    y = np.atleast_1d(np.copy(y)).astype(np.float64)
    z = np.atleast_1d(np.copy(z)).astype(np.float64)
    # calculate the geodetic coordinates using the specified method
    if method.lower() == "moritz":
        lon, lat, h = _moritz_iterative(
            x, y, z, a_axis=a_axis, flat=flat, eps=eps, iterations=iterations
        )
    elif method.lower() == "bowring":
        lon, lat, h = _bowring_iterative(
            x, y, z, a_axis=a_axis, flat=flat, eps=eps, iterations=iterations
        )
    elif method.lower() == "zhu":
        lon, lat, h = _zhu_closed_form(x, y, z, a_axis=a_axis, flat=flat)
    else:
        raise ValueError(f"Unknown conversion method: {method}")
    # return longitude, latitude and height
    # flattened to singular values if necessary
    if singular_values:
        return (lon.item(), lat.item(), h.item())
    else:
        return (lon, lat, h)


def _moritz_iterative(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    a_axis: float = _wgs84.a_axis,
    flat: float = _wgs84.flat,
    eps: float = np.finfo(np.float64).eps,
    iterations: int = 10,
):
    """
    Convert from cartesian coordinates to geodetic coordinates
    using the iterative solution of :cite:t:`HofmannWellenhof:2006hy`

    Parameters
    ----------
    x, np.ndarray
        Cartesian x-coordinates (meters)
    y, np.ndarray
        Cartesian y-coordinates (meters)
    z, np.ndarray
        Cartesian z-coordinates (meters)
    a_axis: float, default 6378137.0
        Semi-major axis of the ellipsoid
    flat: float, default 1.0/298.257223563
        Ellipsoidal flattening
    eps: float, default np.finfo(np.float64).eps
        Tolerance for iterative method
    iterations: int, default 10
        Maximum number of iterations
    """
    # Linear eccentricity and first numerical eccentricity
    lin_ecc = np.sqrt((2.0 * flat - flat**2) * a_axis**2)
    ecc1 = lin_ecc / a_axis
    # calculate longitude
    lon = np.degrees(np.arctan2(y, x))
    # set initial estimate of height to 0
    h = np.zeros_like(lon)
    h0 = np.inf * np.ones_like(lon)
    # calculate radius of parallel
    p = np.hypot(x, y)
    # initial estimated value for phi using h=0
    phi = np.arctan(z / (p * (1.0 - ecc1**2)))
    # iterate to tolerance or to maximum number of iterations
    i = 0
    while np.any(np.abs(h - h0) > eps) and (i <= iterations):
        # copy previous iteration of height
        h0 = np.copy(h)
        # calculate radius of curvature
        N = a_axis / np.sqrt(1.0 - ecc1**2 * np.sin(phi) ** 2)
        # estimate new value of height
        h = p / np.cos(phi) - N
        # estimate new value for latitude using heights
        phi = np.arctan(z / (p * (1.0 - ecc1**2 * N / (N + h))))
        # add to iterator
        i += 1
    # convert latitudes and fix values
    lat = np.clip(np.degrees(phi), -90.0, 90.0)
    # return longitude, latitude and height
    return (lon, lat, h)


def _bowring_iterative(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    a_axis: float = _wgs84.a_axis,
    flat: float = _wgs84.flat,
    eps: float = np.finfo(np.float64).eps,
    iterations: int = 10,
):
    """
    Convert from cartesian coordinates to geodetic coordinates using
    the iterative solution of :cite:t:`Bowring:1976jh,Bowring:1985du`

    Parameters
    ----------
    x, np.ndarray
        Cartesian x-coordinates (meters)
    y, np.ndarray
        Cartesian y-coordinates (meters)
    z, np.ndarray
        Cartesian z-coordinates (meters)
    a_axis: float, default 6378137.0
        Semi-major axis of the ellipsoid
    flat: float, default 1.0/298.257223563
        Ellipsoidal flattening
    eps: float, default np.finfo(np.float64).eps
        Tolerance for iterative method
    iterations: int, default 10
        Maximum number of iterations
    """
    # semiminor axis of the WGS84 ellipsoid [m]
    b_axis = (1.0 - flat) * a_axis
    # Linear eccentricity
    lin_ecc = np.sqrt((2.0 * flat - flat**2) * a_axis**2)
    # square of first and second numerical eccentricity
    e12 = lin_ecc**2 / a_axis**2
    e22 = lin_ecc**2 / b_axis**2
    # calculate longitude
    lon = np.degrees(np.arctan2(y, x))
    # calculate radius of parallel
    p = np.hypot(x, y)
    # initial estimated value for reduced parametric latitude
    u = np.arctan(a_axis * z / (b_axis * p))
    # initial estimated value for latitude
    phi = np.arctan(
        (z + e22 * b_axis * np.sin(u) ** 3)
        / (p - e12 * a_axis * np.cos(u) ** 3)
    )
    phi0 = np.inf * np.ones_like(lon)
    # iterate to tolerance or to maximum number of iterations
    i = 0
    while np.any(np.abs(phi - phi0) > eps) and (i <= iterations):
        # copy previous iteration of phi
        phi0 = np.copy(phi)
        # calculate reduced parametric latitude
        u = np.arctan(b_axis * np.tan(phi) / a_axis)
        # estimate new value of latitude
        phi = np.arctan(
            (z + e22 * b_axis * np.sin(u) ** 3)
            / (p - e12 * a_axis * np.cos(u) ** 3)
        )
        # add to iterator
        i += 1
    # calculate final radius of curvature
    N = a_axis / np.sqrt(1.0 - e12 * np.sin(phi) ** 2)
    # estimate final height (Bowring, 1985)
    h = p * np.cos(phi) + z * np.sin(phi) - a_axis**2 / N
    # convert latitudes and fix values
    lat = np.clip(np.degrees(phi), -90.0, 90.0)
    # return longitude, latitude and height
    return (lon, lat, h)


def _zhu_closed_form(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    a_axis: float = _wgs84.a_axis,
    flat: float = _wgs84.flat,
):
    """
    Convert from cartesian coordinates to geodetic coordinates
    using the closed-form solution of :cite:t:`Zhu:1993ja`

    Parameters
    ----------
    x, np.ndarray
        Cartesian x-coordinates (meters)
    y, np.ndarray
        Cartesian y-coordinates (meters)
    z, np.ndarray
        Cartesian z-coordinates (meters)
    a_axis: float, default 6378137.0
        Semi-major axis of the ellipsoid
    flat: float, default 1.0/298.257223563
        Ellipsoidal flattening
    """
    # semiminor axis of the WGS84 ellipsoid [m]
    b_axis = (1.0 - flat) * a_axis
    # Linear eccentricity
    lin_ecc = np.sqrt((2.0 * flat - flat**2) * a_axis**2)
    # square of first numerical eccentricity
    e12 = lin_ecc**2 / a_axis**2
    # calculate longitude
    lon = np.degrees(np.arctan2(y, x))
    # calculate radius of parallel
    w = np.hypot(x, y)
    # allocate for output latitude and height
    lat = np.zeros_like(lon)
    h = np.zeros_like(lon)
    if np.any(w == 0):
        # special case where w == 0 (exact polar solution)
        (ind,) = np.nonzero(w == 0)
        h[ind] = np.sign(z[ind]) * z[ind] - b_axis
        lat[ind] = 90.0 * np.sign(z[ind])
    else:
        # all other cases
        (ind,) = np.nonzero(w != 0)
        l = e12 / 2.0
        m = (w[ind] / a_axis) ** 2.0
        n = ((1.0 - e12) * z[ind] / b_axis) ** 2.0
        i = -(2.0 * l**2 + m + n) / 2.0
        k = (l**2.0 - m - n) * l**2.0
        q = (1.0 / 216.0) * (m + n - 4.0 * l**2) ** 3.0 + m * n * l**2.0
        D = np.sqrt((2.0 * q - m * n * l**2) * m * n * l**2)
        B = i / 3.0 - (q + D) ** (1.0 / 3.0) - (q - D) ** (1.0 / 3.0)
        t = np.sqrt(np.sqrt(B**2 - k) - (B + i) / 2.0) - np.sign(
            m - n
        ) * np.sqrt((B - i) / 2.0)
        wi = w / (t + l)
        zi = (1.0 - e12) * z[ind] / (t - l)
        # calculate latitude and height
        lat[ind] = np.degrees(np.arctan2(zi, ((1.0 - e12) * wi)))
        h[ind] = np.sign(t - 1.0 + l) * np.sqrt(
            (w - wi) ** 2.0 + (z[ind] - zi) ** 2.0
        )
    # return longitude, latitude and height
    return (lon, lat, h)


def to_ENU(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    lon0: float | np.ndarray = 0.0,
    lat0: float | np.ndarray = 0.0,
    h0: float | np.ndarray = 0.0,
    a_axis: float = _wgs84.a_axis,
    flat: float = _wgs84.flat,
):
    """
    Convert from Earth-Centered Earth-Fixed (ECEF) cartesian coordinates
    to East-North-Up coordinates (ENU)

    Parameters
    ----------
    x, np.ndarray
        Cartesian x-coordinates (meters)
    y, np.ndarray
        Cartesian y-coordinates (meters)
    z, np.ndarray
        Cartesian z-coordinates (meters)
    lon0: float or np.ndarray, default 0.0
        Reference longitude (degrees east)
    lat0: float or np.ndarray, default 0.0
        Reference latitude (degrees north)
    h0: float or np.ndarray, default 0.0
        Reference height (meters)
    a_axis: float, default 6378137.0
        Semi-major axis of the ellipsoid
    flat: float, default 1.0/298.257223563
        Ellipsoidal flattening

    Returns
    -------
    E: np.ndarray
        East coordinates (meters)
    N: np.ndarray
        North coordinates (meters)
    U: np.ndarray
        Up coordinates (meters)
    """
    # verify axes and copy to not modify inputs
    singular_values = np.ndim(x) == 0
    x = np.atleast_1d(np.copy(x)).astype(np.float64)
    y = np.atleast_1d(np.copy(y)).astype(np.float64)
    z = np.atleast_1d(np.copy(z)).astype(np.float64)
    # convert latitude and longitude to ECEF
    X0, Y0, Z0 = to_cartesian(lon0, lat0, h=h0, a_axis=a_axis, flat=flat)
    # latitude and longitude in radians
    theta0 = np.radians(lat0)
    lambda0 = np.radians(lon0)
    # calculate the rotation matrix
    R = np.zeros((3, 3))
    R[0, 0] = -np.sin(lambda0)
    R[0, 1] = np.cos(lambda0)
    R[0, 2] = 0.0
    R[1, 0] = -np.sin(theta0) * np.cos(lambda0)
    R[1, 1] = -np.sin(theta0) * np.sin(lambda0)
    R[1, 2] = np.cos(theta0)
    R[2, 0] = np.cos(theta0) * np.cos(lambda0)
    R[2, 1] = np.cos(theta0) * np.sin(lambda0)
    R[2, 2] = np.sin(theta0)
    # calculate the ENU coordinates
    E, N, U = np.dot(R, np.vstack((x - X0, y - Y0, z - Z0)))
    # return the ENU coordinates
    # flattened to singular values if necessary
    if singular_values:
        return (E.item(), N.item(), U.item())
    else:
        return (E, N, U)


def from_ENU(
    E: np.ndarray,
    N: np.ndarray,
    U: np.ndarray,
    lon0: float | np.ndarray = 0.0,
    lat0: float | np.ndarray = 0.0,
    h0: float | np.ndarray = 0.0,
    a_axis: float = _wgs84.a_axis,
    flat: float = _wgs84.flat,
):
    """
    Convert from East-North-Up coordinates (ENU) to
    Earth-Centered Earth-Fixed (ECEF) cartesian coordinates

    Parameters
    ----------
    E, np.ndarray
        East coordinates (meters)
    N, np.ndarray
        North coordinates (meters)
    U, np.ndarray
        Up coordinates (meters)
    lon0: float or np.ndarray, default 0.0
        Reference longitude (degrees east)
    lat0: float or np.ndarray, default 0.0
        Reference latitude (degrees north)
    h0: float or np.ndarray, default 0.0
        Reference height (meters)
    a_axis: float, default 6378137.0
        Semi-major axis of the ellipsoid
    flat: float, default 1.0/298.257223563
        Ellipsoidal flattening

    Returns
    -------
    x, float
        Cartesian x-coordinates (meters)
    y, float
        Cartesian y-coordinates (meters)
    z, float
        Cartesian z-coordinates (meters)
    """
    # verify axes and copy to not modify inputs
    singular_values = np.ndim(E) == 0
    E = np.atleast_1d(np.copy(E)).astype(np.float64)
    N = np.atleast_1d(np.copy(N)).astype(np.float64)
    U = np.atleast_1d(np.copy(U)).astype(np.float64)
    # convert latitude and longitude to ECEF
    X0, Y0, Z0 = to_cartesian(lon0, lat0, h=h0, a_axis=a_axis, flat=flat)
    # latitude and longitude in radians
    theta0 = np.radians(lat0)
    lambda0 = np.radians(lon0)
    # calculate the rotation matrix
    R = np.zeros((3, 3))
    R[0, 0] = -np.sin(lambda0)
    R[1, 0] = np.cos(lambda0)
    R[2, 0] = 0.0
    R[0, 1] = -np.sin(theta0) * np.cos(lambda0)
    R[1, 1] = -np.sin(theta0) * np.sin(lambda0)
    R[2, 1] = np.cos(theta0)
    R[0, 2] = np.cos(theta0) * np.cos(lambda0)
    R[1, 2] = np.cos(theta0) * np.sin(lambda0)
    R[2, 2] = np.sin(theta0)
    # calculate the ECEF coordinates
    x, y, z = np.dot(R, np.vstack((E, N, U)))
    # add reference coordinates
    x += X0
    y += Y0
    z += Z0
    # return the ECEF coordinates
    # flattened to singular values if necessary
    if singular_values:
        return (x.item(), y.item(), z.item())
    else:
        return (x, y, z)


def to_horizontal(
    E: np.ndarray,
    N: np.ndarray,
    U: np.ndarray,
):
    """
    Convert from East-North-Up coordinates (ENU) to a
    celestial horizontal coordinate system (alt-az)

    Parameters
    ----------
    E: np.ndarray
        East coordinates (meters)
    N: np.ndarray
        North coordinates (meters)
    U: np.ndarray
        Up coordinates (meters)

    Returns
    -------
    alpha: np.ndarray
        Altitude (elevation) angle (degrees)
    phi: np.ndarray
        Azimuth angle (degrees)
    D: np.ndarray
        Distance from observer to object (meters)
    """
    # calculate distance to object
    # convert coordinates to unit vectors
    D = np.sqrt(E**2 + N**2 + U**2)
    # altitude (elevation) angle in degrees
    alpha = np.degrees(np.arcsin(U / D))
    # azimuth angle in degrees (fixed to 0 to 360)
    phi = np.mod(np.degrees(np.arctan2(E / D, N / D)), 360.0)
    return (alpha, phi, D)


def to_zenith(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    lon0: float | np.ndarray = 0.0,
    lat0: float | np.ndarray = 0.0,
    h0: float | np.ndarray = 0.0,
    a_axis: float = _wgs84.a_axis,
    flat: float = _wgs84.flat,
):
    """
    Calculate zenith angle of an object from Earth-Centered
    Earth-Fixed (ECEF) cartesian coordinates

    Parameters
    ----------
    x, np.ndarray
        Cartesian x-coordinates (meters)
    y, np.ndarray
        Cartesian y-coordinates (meters)
    z, np.ndarray
        Cartesian z-coordinates (meters)
    lon0: float or np.ndarray, default 0.0
        Reference longitude (degrees east)
    lat0: float or np.ndarray, default 0.0
        Reference latitude (degrees north)
    h0: float or np.ndarray, default 0.0
        Reference height (meters)
    a_axis: float, default 6378137.0
        Semi-major axis of the ellipsoid
    flat: float, default 1.0/298.257223563
        Ellipsoidal flattening

    Returns
    -------
    zenith: np.ndarray
        Zenith angle of object (degrees)
    """
    # convert from ECEF to ENU
    E, N, U = to_ENU(
        x, y, z, lon0=lon0, lat0=lat0, h0=h0, a_axis=a_axis, flat=flat
    )
    # convert from ENU to horizontal coordinates
    alpha, phi, D = to_horizontal(E, N, U)
    # calculate zenith angle in degrees
    zenith = 90.0 - alpha
    # return zenith angle
    return zenith


def geocentric_latitude(
    lat: np.ndarray,
    flat: float = _wgs84.flat,
):
    """
    Compute the geocentric latitude from a geodetic latitude
    using a simplified empirical relation :cite:p:`Snyder:1982gf`

    Parameters
    ----------
    lat: np.ndarray
        Geodetic latitudes (degrees north)
    flat: float, default 1.0/298.257223563
        Ellipsoidal flattening

    Returns
    -------
    geolat: np.ndarray
        Geocentric latitude (degrees)
    """
    # convert latitude to radians
    phi = np.radians(lat)
    # compute difference between geodetic and geocentric latitudes
    d2 = (flat + flat**2 / 2.0) * np.sin(2.0 * phi)
    d4 = (flat**2 / 2.0) * np.sin(4.0 * phi)
    geolat = lat - np.degrees(d2 - d4)
    return geolat


def scale_areas(*args, **kwargs):
    warnings.warn(
        "Deprecated. Please use pyTMD.spatial.scale_factors instead",
        DeprecationWarning,
    )
    return scale_factors(*args, **kwargs)


def scale_factors(
    lat: np.ndarray,
    flat: float = _wgs84.flat,
    reference_latitude: float = 70.0,
    metric: str = "area",
):
    """
    Calculates scaling factors to account for polar stereographic
    distortion including special case of at the exact pole
    :cite:p:`Snyder:1982gf`

    Parameters
    ----------
    lat: np.ndarray
        Latitude (degrees north)
    flat: float, default 1.0/298.257223563
        Ellipsoidal flattening
    reference_latitude: float, default 70.0
        Reference latitude (true scale latitude)
    metric: str, default 'area'
        Metric to calculate scaling factors

            - ``'distance'``: scale factors for distance
            - ``'area'``: scale factors for area

    Returns
    -------
    scale: np.ndarray
        Scaling factors at input latitudes
    """
    if metric.lower() not in ["distance", "area"]:
        raise ValueError("Unknown metric")
    # power for scaling factors
    power = 1.0 if metric.lower() == "distance" else 2.0
    # convert latitude to positive radians
    phi = np.radians(np.abs(lat))
    # convert reference latitude to positive radians
    phi_ref = np.radians(np.abs(reference_latitude))
    # square of the eccentricity of the ellipsoid
    # ecc2 = (1-b**2/a**2) = 2.0*flat - flat^2
    ecc2 = 2.0 * flat - flat**2
    # eccentricity of the ellipsoid
    ecc = np.sqrt(ecc2)
    # get p values following equations 17.33 and 17.35
    p = np.sqrt(np.power(1.0 + ecc, 1.0 + ecc) * np.power(1.0 - ecc, 1.0 - ecc))
    # calculate m factors using equation 12.15
    m = np.cos(phi) / np.sqrt(1.0 - ecc2 * np.sin(phi) ** 2)
    m_ref = np.cos(phi_ref) / np.sqrt(1.0 - ecc2 * np.sin(phi_ref) ** 2)
    # calculate t factors using equation 13.9
    t = np.tan(np.pi / 4.0 - phi / 2.0) / np.power(
        (1.0 - ecc * np.sin(phi)) / (1.0 + ecc * np.sin(phi)), ecc / 2.0
    )
    t_ref = np.tan(np.pi / 4.0 - phi_ref / 2.0) / np.power(
        (1.0 - ecc * np.sin(phi_ref)) / (1.0 + ecc * np.sin(phi_ref)), ecc / 2.0
    )
    # calculate scaling factors following Snyder (1982)
    # ignore divide by zero and invalid value warnings
    with np.errstate(divide="ignore", invalid="ignore"):
        # check if reference latitude is at the pole
        if np.isclose(phi_ref, np.pi / 2.0):
            # equations 17.32 and 17.33
            k = 2.0 * t / (p * m)
            # at the pole (true scale)
            k_pole = 1.0
        else:
            # equations 17.32 and 17.34
            k = (m_ref / m) * (t / t_ref)
            # at the pole from equation 17.35
            k_pole = 0.5 * m_ref * p / t_ref
        # distance and area scaling factors with special case at the pole
        scale = np.where(
            np.isclose(phi, np.pi / 2.0),
            np.power(1.0 / k_pole, power),
            np.power(1.0 / k, power),
        )
    # return the scaling factors
    return scale
