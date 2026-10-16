#!/usr/bin/env python
"""
potential.py
Written by Tyler Sutterley (10/2026)
Spherical harmonic expansions and summations

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    timescale: Python tools for time and astronomical calculations
        https://pypi.org/project/timescale/
    xarray: N-D labeled arrays and datasets in Python
        https://docs.xarray.dev/en/stable/

PROGRAM DEPENDENCIES:
    constituents.py: calculates constituent parameters and nodal arguments
    earth.py: calculates Earth parameters and Body Tide Love numbers
    math.py: Special functions of mathematical physics
    spatial.py: utilities for working with geospatial data

UPDATE HISTORY:
    Updated 10/2026: added function to estimate crustal deformation for
        non-geographic coordinates using a Green's function approach
    Updated 09/2026: added function to estimate sea water densities
    Written 09/2026
"""

from __future__ import annotations

import numpy as np
import xarray as xr
import pyTMD.constituents
import pyTMD.earth
import pyTMD.math
import pyTMD.spatial

# get ellipsoid parameters and load love number tables
from pyTMD.earth import _lln_table, _wgs84

__all__ = [
    "ocean_harmonics",
    "crustal_loading",
    "_harmonic_summation",
    "_greens_function",
    "_seawater_density",
]


# PURPOSE: convert tidal constituents into spherical harmonics
def ocean_harmonics(
    ds: xr.Dataset,
    lmax: int = 696,
    a_axis: float = _wgs84.a_axis,
    flat: float = _wgs84.flat,
    GM: float = _wgs84.GM,
    rho_w: float | xr.DataArray = 1025.0,
    lln: str = "han-wahr",
    reference: str = "CE",
    **kwargs,
):
    r"""
    Converts ocean tide model constituents into spherical
    harmonic coefficients :cite:p:`Petit:2010tp,Ray:1989wf,Wahr:1998hy`

    Parameters
    ----------
    ds: xarray.Dataset
        Dataset containing tidal harmonic constants
    lmax: int, default 696
        Upper bound of spherical harmonic degrees
    a_axis: float, default 6378136.3
        Semi-major axis of the Earth (meters)
    flat: float, default 1.0/298.257223563
        Ellipsoidal flattening
    GM: float, default 3.986004418e14
        Geocentric gravitational constant (m\ :sup:`3` s\ :sup:`-2`)
    rho_w: float or xarray.DataArray, default 1025.0
        Density of sea water (kg m\ :sup:`-3`)

        Can be spatially uniform or a grid matching the tidal dataset
    lln: str, default 'han-wahr'
        Name of the Load Love number dataset to use

            - ``'han-wahr'``: :cite:t:`Han:1995go`
            - ``'gegout'``: :cite:t:`Gegout:2010gc`
            - ``'wang-prem'``: :cite:t:`Wang:2012gc`
    reference: str, default 'CE'
        Reference frame of degree 1 load Love numbers

            - ``'CF'``: Center of Surface Figure
            - ``'CL'``: Center of Surface Lateral Figure
            - ``'CH'``: Center of Surface Height Figure
            - ``'CM'``: Center of Mass of Earth System
            - ``'CE'``: Center of Mass of Solid Earth
    kwargs: dict
        Keyword arguments for :py:func:`pyTMD.constituents.frequency`

    Returns
    -------
    Ylms: xr.Dataset
        Fully-normalized spherical harmonic coefficients

            - ``clm``: cosine spherical harmonics (complex)
            - ``slm``: sine spherical harmonics (complex)
    """
    # verify units of input data are in meters
    ds = ds.tmd.to_units("meters")
    # verify coordinate reference system (global and geographic)
    if not ds.tmd.is_global:
        raise ValueError("Unsupported coordinate reference system")
    # verify load Love number table name
    if lln not in _lln_table.keys():
        raise ValueError(f"Unknown load Love number dataset {lln}")
    # verify maximum degree and order
    lmax = int(lmax)
    if lmax < 2:
        raise ValueError("Degree of truncation should be at least 2")
    # list of tidal constituents
    constituents = ds.tmd.constituents
    nc = len(constituents)
    # angular frequencies for constituents
    omega = pyTMD.constituents.frequency(constituents, **kwargs)

    # Earth parameters and physical constants
    # universal gravitational constant [N*m^2/kg^2]
    G = 6.67430e-11
    # average radius of the Earth with same volume as ellipsoid [m]
    rad_e = a_axis * np.power(1.0 - flat, 1.0 / 3.0)
    # average density of the Earth [kg / m^3]
    rho_e = 0.75 * GM / (G * np.pi * rad_e**3)

    # convert from geodetic latitude to geocentric latitude
    geolat = pyTMD.spatial.geocentric_latitude(ds.y, flat=flat)
    # calculate colatitude and longitude (radians)
    theta = np.radians(90.0 - geolat)
    lmda = np.radians(ds.x)
    # convert longitudes to range 0:360 (if previously -180:180)
    lmda = lmda.where(lmda >= 0, lmda + 2.0 * np.pi, drop=False)
    # multiply sin(th) with differentials of theta and lambda
    # to calculate the integration factor at each latitude
    dlam = np.abs(lmda[1] - lmda[0])
    dth = np.abs(theta[1] - theta[0])
    int_fact = np.sin(theta) * dlam * dth

    # spherical harmonic degree and order
    l = np.arange(lmax + 1)
    m = np.arange(lmax + 1)
    # calculate polynomials using Martin Mohlenkamp's relation
    Plm, dPlm = pyTMD.math.legendreP(lmax, np.cos(theta))
    # read load Love numbers from table
    hl, kl, ll = pyTMD.earth.load_love_numbers(
        _lln_table[lln], reference=reference, lmax=lmax
    )
    # convert to xarray data arrays
    Plm = xr.DataArray(
        Plm,
        dims=("l", "m", "y"),
        coords={"l": l, "m": m, "y": ds.y},
    )
    kl = xr.DataArray(
        kl,
        dims=("l",),
        coords={"l": l},
    )
    # allocate for frequency-dependent load Love numbers adjustments
    dk = xr.DataArray(
        np.zeros((lmax + 1, lmax + 1)),
        dims=("l", "m"),
        coords={"l": l, "m": m},
    )

    # calculate cos/sin of lambda arrays
    m_cos = np.cos(Plm.m.dot(lmda))
    m_sin = np.sin(Plm.m.dot(lmda))
    # integration coefficients for converting to spherical harmonics
    int_coeff = int_fact * Plm

    # allocate for output spherical harmonics
    clm = np.zeros((lmax + 1, lmax + 1, nc), dtype=np.complex128)
    slm = np.zeros((lmax + 1, lmax + 1, nc), dtype=np.complex128)
    # for each constituent
    for i, c in enumerate(constituents):
        # get constituent and replace nans with 0
        data = ds[c].fillna(0.0)
        # multiply gridded data with sin/cos of m#lambda
        # sum through all lambdas in the dot product
        # multiply heights by sea water density
        # sea water density can be a scalar value for uniform
        # or a map of the column averages
        d_cos = m_cos.dot(rho_w * data)
        d_sin = m_sin.dot(rho_w * data)
        # adjust load Love numbers for frequency dependence
        dh, dk[2, 1], dl = pyTMD.earth.adjust_load_love_numbers(omega[i])
        # degree dependent factors for converting from water equivalent
        # taking into account frequency dependence of load Love numbers
        # modified from Wahr et al., (2018)
        dfactor = (
            3.0
            * (1.0 + kl + dk)
            / (1.0 + 2.0 * Plm.l)
            / (4.0 * np.pi * rad_e * rho_e)
        )
        # integrate over all latitudes
        # fully-normalize output spherical harmonics
        clm[:, :, i] = dfactor * int_coeff.dot(d_cos, dim="y")
        slm[:, :, i] = dfactor * int_coeff.dot(d_sin, dim="y")
    # convert to xarray dataset
    Ylms = xr.Dataset(
        data_vars=dict(
            clm=(["l", "m", "constituent"], clm),
            slm=(["l", "m", "constituent"], slm),
        ),
        coords={"l": l, "m": m, "constituent": constituents},
    )
    # add attributes for dimensions
    Ylms.l.attrs["long_name"] = "spherical harmonic degree"
    Ylms.m.attrs["long_name"] = "spherical harmonic order"
    Ylms.l.attrs["standard_name"] = "degree"
    Ylms.m.attrs["standard_name"] = "order"
    Ylms.l.attrs["units"] = "wavenumber"
    Ylms.m.attrs["units"] = "wavenumber"
    # add attributes for spherical harmonics
    Ylms.clm.attrs["long_name"] = "cosine spherical harmonics (complex)"
    Ylms.slm.attrs["long_name"] = "sine spherical harmonics (complex)"
    Ylms.clm.attrs["description"] = (
        "cosine spherical harmonic coefficients containing the in-phase (real) "
        "and out-of-phase (imag) parts of the tidal constituents"
    )
    Ylms.slm.attrs["description"] = (
        "sine spherical harmonic coefficients containing the in-phase (real) "
        "and out-of-phase (imag) parts of the tidal constituents"
    )
    # copy attributes from original dataset
    Ylms.attrs.update(ds.attrs)
    Ylms.attrs["product_type"] = "gravity_field"
    Ylms.attrs["normalization"] = "fully-normalized"
    # add attributes for degree of truncation
    Ylms.attrs["max_degree"] = lmax
    Ylms.attrs["max_order"] = lmax
    # add attributes for earth model and love numbers
    Ylms.attrs["earth_love_numbers"] = lln
    Ylms.attrs["reference_frame"] = reference
    # add attributes for earth and model parameters
    Ylms.attrs["earth_radius"] = f"{rad_e:0.3f} m"
    Ylms.attrs["earth_density"] = f"{rho_e:0.3f} kg/m^3"
    Ylms.attrs["earth_inverse_flattening"] = f"{1.0 / flat:0.3f}"
    Ylms.attrs["earth_gravity_constant"] = f"{GM:0.3f} m^3/s^2"
    # add attribute for seawater density (uniform or gridded)
    if np.ndim(rho_w) == 0:
        Ylms.attrs["seawater_density"] = f"{float(rho_w):0.3f} kg/m^3"
    else:
        Ylms.attrs["seawater_density"] = "gridded"
    # check if chunks were present in original dataset
    if hasattr(ds, "chunks") and bool(ds.chunks):
        Ylms = Ylms.chunk("auto")
    # return the spherical harmonic dataset
    return Ylms


# PURPOSE: convert spherical harmonics into load tide maps
def crustal_loading(
    ds: xr.Dataset,
    other: xr.Dataset,
    **kwargs,
):
    r"""
    Calculates the crustal deformation induced by ocean tide loading
    :cite:p:`Desai:2014ee,Farrell:1972cm,Ray:2025uo`

    Parameters
    ----------
    ds: xarray.Dataset
        Dataset with spatial coordinates
    other: xarray.Dataset
        Auxiliary dataset for the crustal loading calculation
    method: str, default "SHA"
        Method to use for the crustal loading calculation

        - `'SHA'`: spherical harmonic analysis
        - `'GF'`: Green's function approach
    kwargs: dict
        Keyword arguments for the crustal loading function

    Returns
    -------
    ds: xr.Dataset
        Dataset containing tidal harmonic constants
    """
    # get the calculation method
    method = kwargs.get("method", "SHA")
    if method.upper() == "SHA":
        return _harmonic_summation(ds, other, **kwargs)
    elif method.upper() == "GF":
        return _greens_function(ds, other, **kwargs)
    else:
        raise ValueError(f"Unknown method {method}")


def _harmonic_summation(
    ds: xr.Dataset,
    Ylms: xr.Dataset,
    lmax: int | None = None,
    a_axis: float = _wgs84.a_axis,
    flat: float = _wgs84.flat,
    lln: str = "han-wahr",
    reference: str = "CE",
    **kwargs,
):
    r"""
    Calculates the crustal deformation induced by ocean tide loading
    via spherical harmonic summation :cite:p:`Desai:2014ee,Ray:2025uo`

    Parameters
    ----------
    ds: xarray.Dataset
        Dataset with spatial coordinates
    Ylms: xarray.Dataset
        Dataset with spherical harmonic coefficients

            - ``clm``: cosine spherical harmonics (complex)
            - ``slm``: sine spherical harmonics (complex)
    lmax: int or None, default None
        Upper bound of spherical harmonic degrees
    a_axis: float, default 6378136.3
        Semi-major axis of the Earth (meters)
    flat: float, default 1.0/298.257223563
        Ellipsoidal flattening
    lln: str, default 'han-wahr'
        Name of the Load Love number dataset to use

            - ``'han-wahr'``: :cite:t:`Han:1995go`
            - ``'gegout'``: :cite:t:`Gegout:2010gc`
            - ``'wang-prem'``: :cite:t:`Wang:2012gc`
    reference: str, default 'CE'
        Reference frame of degree 1 load Love numbers

            - ``'CF'``: Center of Surface Figure
            - ``'CL'``: Center of Surface Lateral Figure
            - ``'CH'``: Center of Surface Height Figure
            - ``'CM'``: Center of Mass of Earth System
            - ``'CE'``: Center of Mass of Solid Earth
    kwargs: dict
        Keyword arguments for :py:func:`pyTMD.constituents.frequency`

    Returns
    -------
    ds: xr.Dataset
        Dataset containing tidal harmonic constants
    """
    # verify load Love number table name
    if lln not in _lln_table.keys():
        raise ValueError(f"Unknown load Love number dataset {lln}")
    # verify maximum degree and order
    if lmax is None:
        lmax = int(Ylms.l.max().values)
    if lmax < 2:
        raise ValueError("Degree of truncation should be at least 2")
    # truncate to degree and order lmax
    if lmax < Ylms.l.max():
        Ylms = Ylms.where((Ylms.l <= lmax) & (Ylms.m <= lmax), drop=True)
    # extract harmonics and convert to datasets
    clm = Ylms.clm.to_dataset(dim="constituent")
    slm = Ylms.slm.to_dataset(dim="constituent")
    # list of tidal constituents
    constituents = clm.tmd.constituents
    # angular frequencies for constituents
    omega = pyTMD.constituents.frequency(constituents, **kwargs)

    # Earth parameters and physical constants
    # average radius of the Earth with same volume as ellipsoid [m]
    rad_e = a_axis * np.power(1.0 - flat, 1.0 / 3.0)
    # convert from geodetic latitude to geocentric latitude
    geolat = pyTMD.spatial.geocentric_latitude(ds.y, flat=flat)
    # calculate colatitude and longitude (radians)
    theta = np.radians(90.0 - geolat)
    lmda = np.radians(ds.x)
    # convert longitudes to range 0:360 (if previously -180:180)
    lmda = lmda.where(lmda >= 0, lmda + 2.0 * np.pi, drop=False)

    # spherical harmonic degree and order
    l = np.arange(lmax + 1)
    m = np.arange(lmax + 1)
    # calculate polynomials using Martin Mohlenkamp's relation
    Plm, dPlm = pyTMD.math.legendreP(lmax, np.cos(theta))
    # read load Love numbers from table
    hl, kl, ll = pyTMD.earth.load_love_numbers(
        _lln_table[lln], reference=reference, lmax=lmax
    )
    # convert to xarray data arrays
    Plm = xr.DataArray(
        Plm,
        dims=("l", "m", "y"),
        coords={"l": l, "m": m, "y": ds.y},
    )
    hl = xr.DataArray(
        hl,
        dims=("l",),
        coords={"l": l},
    )
    # allocate for frequency-dependent load Love numbers adjustments
    dh = xr.DataArray(
        np.zeros((lmax + 1, lmax + 1)),
        dims=("l", "m"),
        coords={"l": l, "m": m},
    )

    # calculate cos/sin of lambda arrays
    m_cos = np.cos(Plm.m.dot(lmda))
    m_sin = np.sin(Plm.m.dot(lmda))

    # create output dataset
    tmp = xr.Dataset(coords=ds.coords)
    # for each constituent
    for i, c in enumerate(constituents):
        # adjust load Love numbers for frequency dependence
        dh[2, 1], dk, dl = pyTMD.earth.adjust_load_love_numbers(omega[i])
        # degree dependent factors for converting to crustal deformation
        # taking into account frequency dependence of load Love numbers
        dfactor = rad_e * (hl + dh)
        # summation over all spherical harmonic degrees
        p_cos = Plm.dot(dfactor * clm[c], dim="l")
        p_sin = Plm.dot(dfactor * slm[c], dim="l")
        # summation of cosine and sine harmonics
        tmp[c] = p_cos.dot(m_cos, dim="m") + p_sin.dot(m_sin, dim="m")

    # copy attributes from spherical harmonic dataset
    tmp.attrs.update(Ylms.attrs)
    # check if chunks were present in spherical harmonic dataset
    if hasattr(Ylms, "chunks") and bool(Ylms.chunks):
        tmp = tmp.chunk("auto")
    # return the spatial dataset of tidal constituents
    return tmp


def _greens_function(
    XYZ1: xr.Dataset,
    XYZ2: xr.Dataset,
    lmax: int | None = None,
    a_axis: float = _wgs84.a_axis,
    flat: float = _wgs84.flat,
    GM: float = _wgs84.GM,
    area: float | xr.DataArray | None = None,
    rho_w: float | xr.DataArray = 1025.0,
    lln: str = "han-wahr",
    reference: str = "CE",
    hinf: float | None = None,
    **kwargs,
):
    r"""
    Calculates the crustal deformation induced by ocean tide loading
    via Green's functions :cite:p:`Farrell:1972cm`

    Parameters
    ----------
    XYZ1: xr.Dataset
        Dataset with cartesian coordinates
    XYZ2: xr.Dataset
        Dataset containing tidal harmonic constants
    lmax: int or None, default None
        Upper bound of spherical harmonic degrees
    a_axis: float, default 6378136.3
        Semi-major axis of the Earth (meters)
    flat: float, default 1.0/298.257223563
        Ellipsoidal flattening
    GM: float, default 3.986004418e14
        Geocentric gravitational constant (m\ :sup:`3` s\ :sup:`-2`)
    rho_w: float or xarray.DataArray, default 1025.0
        Density of sea water (kg m\ :sup:`-3`)

        Can be spatially uniform or a grid matching the tidal dataset
    area: float, xr.DataArray or None, default None
        Grid cell areas for integrating tidal harmonics
    lln: str, default 'han-wahr'
        Name of the Load Love number dataset to use

            - ``'han-wahr'``: :cite:t:`Han:1995go`
            - ``'gegout'``: :cite:t:`Gegout:2010gc`
            - ``'wang-prem'``: :cite:t:`Wang:2012gc`
    reference: str, default 'CE'
        Reference frame of degree 1 load Love numbers

            - ``'CF'``: Center of Surface Figure
            - ``'CL'``: Center of Surface Lateral Figure
            - ``'CH'``: Center of Surface Height Figure
            - ``'CM'``: Center of Mass of Earth System
            - ``'CE'``: Center of Mass of Solid Earth
    hinf: float | None, default None
        Asymptotic value for the Love Number of Vertical Displacement
    kwargs: dict
        Keyword arguments for :py:func:`pyTMD.constituents.frequency`

    Returns
    -------
    ds: xr.Dataset
        Dataset containing tidal harmonic constants
    """
    # verify load Love number table name
    if lln not in _lln_table.keys():
        raise ValueError(f"Unknown load Love number dataset {lln}")
    # list of tidal constituents
    constituents = XYZ2.tmd.constituents
    # angular frequencies for constituents
    omega = pyTMD.constituents.frequency(constituents, **kwargs)

    # Earth parameters and physical constants
    # universal gravitational constant [N*m^2/kg^2]
    G = 6.67430e-11
    # average radius of the Earth with same volume as ellipsoid [m]
    rad_e = a_axis * np.power(1.0 - flat, 1.0 / 3.0)
    # average mass of the Earth
    m_e = GM / G
    # average density of the Earth [kg / m^3]
    rho_e = 0.75 * m_e / (np.pi * rad_e**3)
    # scale factor to convert a mass load to uplift
    scale = rad_e / m_e
    # verify scale factor values are close to expected
    if not np.allclose(scale, 6.371e6 / 5.972e24):
        raise ValueError(f"Unexpected scale factor: {scale}")

    # Earth's radius at the given latitude (meters)
    radius1 = pyTMD.math.radius(XYZ1.X, XYZ1.Y, XYZ1.Z)
    radius2 = pyTMD.math.radius(XYZ2.X, XYZ2.Y, XYZ2.Z)
    # angles between vectors
    scalar = pyTMD.math.scalar_product(
        XYZ1.X, XYZ1.Y, XYZ1.Z, XYZ2.X, XYZ2.Y, XYZ2.Z
    ) / (radius1 * radius2)
    alpha = np.arccos(scalar.clip(min=-1, max=1))
    # replace colocated points with spherical caps
    # with areas equal to those of the grid cells
    rad_cap = np.sqrt(area / np.pi) / rad_e
    alpha = alpha.clip(min=rad_cap)
    # calculate colatitude and longitude (radians)
    th1 = np.pi / 2.0 - np.arctan(XYZ1.Z / np.hypot(XYZ1.X, XYZ1.Y))
    th2 = np.pi / 2.0 - np.arctan(XYZ2.Z / np.hypot(XYZ2.X, XYZ2.Y))
    la1 = np.arctan2(XYZ1.Y, XYZ1.X)
    la2 = np.arctan2(XYZ2.Y, XYZ2.X)

    # calculate spherical harmonics for degree 2, order 1
    # to adjust Green's functions for frequency dependence
    l, m = 2, 1
    Ylms1, _ = pyTMD.math.sph_harm(l, th1, la1, m=m)
    Ylms2, _ = pyTMD.math.sph_harm(l, th2, la2, m=m)
    P21 = (8.0 * np.pi) * np.real(Ylms1 * np.conj(Ylms2)) / (2.0 * l + 1.0)

    # calculate Green's function for vertical displacement (unscaled)
    u = xr.zeros_like(alpha)
    u.values = pyTMD.earth.greens_function(
        alpha,
        lmax=lmax,
        lln=lln,
        reference=reference,
        hinf=hinf,
    )

    # create output dataset
    tmp = xr.Dataset(coords=XYZ1.coords)
    # for each constituent
    for i, c in enumerate(constituents):
        # get constituent and replace nans with 0
        data = XYZ2[c].fillna(0.0)
        # adjust load Love numbers for frequency dependence
        dh, dk, dl = pyTMD.earth.adjust_load_love_numbers(omega[i])
        # Green's function for converting to crustal deformation
        # taking into account frequency dependence of load Love numbers
        GF = scale * (u + dh * P21)
        # mass of the tide height change
        M = rho_w * data * area
        # summation over all tide points
        tmp[c] = xr.dot(GF, M, dim=XYZ2.dims, optimize=True)

    # copy attributes from input tidal dataset
    tmp.attrs.update(XYZ2.attrs)
    # add attributes for degree of truncation of Green's function
    if lmax is not None:
        tmp.attrs["max_degree"] = lmax
        tmp.attrs["max_order"] = lmax
    # add attributes for earth model and love numbers
    tmp.attrs["earth_love_numbers"] = lln
    tmp.attrs["reference_frame"] = reference
    # add attributes for earth and model parameters
    tmp.attrs["earth_radius"] = f"{rad_e:0.3f} m"
    tmp.attrs["earth_density"] = f"{rho_e:0.3f} kg/m^3"
    tmp.attrs["earth_inverse_flattening"] = f"{1.0 / flat:0.3f}"
    tmp.attrs["earth_gravity_constant"] = f"{GM:0.3f} m^3/s^2"
    # add attribute for seawater density (uniform or gridded)
    if np.ndim(rho_w) == 0:
        tmp.attrs["seawater_density"] = f"{float(rho_w):0.3f} kg/m^3"
    else:
        tmp.attrs["seawater_density"] = "gridded"
    # check if chunks were present in input tidal dataset
    if hasattr(XYZ2, "chunks") and bool(XYZ2.chunks):
        tmp = tmp.chunk("auto")
    # return the spatial dataset of tidal constituents
    return tmp


def _seawater_density(
    temperature: np.ndarray,
    salinity: np.ndarray,
):
    r"""
    Calculates the density of sea water using the EOS-80 model
    :cite:p:`UNESCO:1981um,Fofonoff:1983wi`

    Parameters
    ----------
    temperature: np.ndarray
        Sea water temperature (\ |degree| C)
    salinity: np.ndarray
        Sea water salinity (unitness)

    Returns
    -------
    rho_w: np.ndarray
        Sea water density (kg/m\ :sup:`3`)

    .. |degree|    unicode:: U+00B0 .. DEGREE SIGN
    """
    # fresh water density at atmospheric pressure as function of temperature
    # coefficients from equation 14 of Fofonoff (1983) derived from Bigg (1967)
    a = np.array(
        [999.842594, 6.793952e-2, -9.095290e-3, -1.120083e-6, 6.536332e-9]
    )
    # polynomial coefficients involving salinity
    b = np.array([0.824493, -4.0899e-3, 7.6438e-5, -8.2467e-7, 5.3875e-9])
    c = np.array([-5.72466e-3, 1.0227e-4, -1.6546e-6])
    d = np.array([4.8314e-4])
    # seawater density at atmospheric pressure
    # equation 13 of Fofonoff (1983)
    rho_w = (
        pyTMD.math.polynomial_sum(a, temperature)
        + np.power(salinity, 1.0) * pyTMD.math.polynomial_sum(b, temperature)
        + np.power(salinity, 1.5) * pyTMD.math.polynomial_sum(c, temperature)
        + np.power(salinity, 2.0) * pyTMD.math.polynomial_sum(d, temperature)
    )
    # return the sea water density
    return rho_w
