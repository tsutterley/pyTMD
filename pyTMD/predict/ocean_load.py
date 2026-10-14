#!/usr/bin/env python
"""
ocean_load.py
Written by Tyler Sutterley (06/2026)
Prediction routines for ocean, load and (long-period) equilibrium tides

REFERENCES:
    G. D. Egbert and S. Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).
    R. Ray and S. Erofeeva, "Long-period tidal variations in the length of day",
        Journal of Geophysical Research: Solid Earth, 119, (2014).

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    timescale: Python tools for time and astronomical calculations
        https://pypi.org/project/timescale/
    xarray: N-D labeled arrays and datasets in Python
        https://docs.xarray.dev/en/stable/

PROGRAM DEPENDENCIES:
    astro.py: computes the basic astronomical mean longitudes
    constituents.py: calculates constituent parameters and nodal arguments
    earth.py: calculates Earth parameters and Body Tide Love numbers
    interpolate.py: interpolation routines for spatial data
    math.py: Special functions of mathematical physics

UPDATE HISTORY:
    Updated 06/2026: moved function to find peaks in a tidal time series
        separated interpolation of admittances from inference functions
        added an admittance function for short-period minor constituents
    Updated 05/2026: verify unique frequencies in short period inference
        moved ellipsoid and love number parameters to earth module
        updated constituent parameters function can output arrays
    Updated 03/2026: simplify structure by splitting up IERS corrections
        and adding wrapper functions where appropriate
        set the maximum degree and order for the HW1995 catalog to 6
        clean up the ephemerides method of calculating solid earth tides
        clean up all of the ephemerides corrections for solid earth tides
        calculate tide-generating forces following Tamura (1982 and 1987)
        split up prediction functions into separate files
    Updated 02/2026: added attributes for constituents to output DataArrays
        do not infer minor constituents with frequencies equal to any major
        revert (again) load pole tides to a newer IERS convention definition
        but allow for both sign definitions based on the convention variable
        add function to calculate variations in Earth orientation parameters
    Updated 12/2025: added tidal LOD calculation from Ray and Erofeeva (2014)
    Updated 11/2025: update all prediction functions to use xarray Datasets
    Updated 09/2025: make permanent tide amplitude an input parameter
        can choose different tide potential catalogs for body tides
        generalize the calculation of body tides for degrees 3+
    Updated 08/2025: add simplified solid earth tide prediction function
        add correction of anelastic effects for long-period body tides
        use sign convention from IERS for complex body tide Love numbers
        include mantle anelastic effects when inferring long-period tides
        allow definition of nominal Love numbers for degree-2 constituents
        added option to include mantle anelastic effects for LPET predict
        switch time decimal in pole tides to nominal years of 365.25 days
        convert angles with numpy radians and degrees functions
        convert arcseconds to radians with asec2rad function in math.py
        return numpy arrays if cannot infer minor constituents
        use a vectorized linear interpolator for inferring from major tides
    Updated 07/2025: revert free-to-mean conversion to April 2023 version
        revert load pole tide to IERS 1996 convention definitions
        mask mean pole values prior to valid epoch of convention
    Updated 05/2025: pass keyword arguments to nodal corrections functions
    Updated 03/2025: changed argument for method calculating mean longitudes
    Updated 02/2025: verify dimensions of harmonic constants
    Updated 11/2024: use Love numbers for long-period tides when inferring
        move body tide Love/Shida numbers to arguments module
    Updated 10/2024: use PREM as the default Earth model for Love numbers
        more descriptive error message if cannot infer minor constituents
        updated calculation of long-period equilibrium tides
        added option to use Munk-Cartwright admittance interpolation for minor
    Updated 09/2024: verify order of minor constituents to infer
        fix to use case insensitive assertions of string argument values
        split infer minor function into short and long period calculations
        add two new functions to infer semi-diurnal and diurnal tides separately
    Updated 08/2024: minor nodal angle corrections in radians to match arguments
        include inference of eps2 and eta2 when predicting from GOT models
        add keyword argument to allow inferring specific minor constituents
        use nodal arguments for all non-OTIS model type cases
        add load pole tide function that exports in cartesian coordinates
        add ocean pole tide function that exports in cartesian coordinates
    Updated 07/2024: use normalize_angle from pyTMD astro module
        make number of days to convert tide time to MJD a variable
    Updated 02/2024: changed class name for ellipsoid parameters to datum
    Updated 01/2024: moved minor arguments calculation into new function
        moved constituent parameters function from predict to arguments
    Updated 12/2023: phase_angles function renamed to doodson_arguments
    Updated 09/2023: moved constituent parameters function within this module
    Updated 08/2023: changed ESR netCDF4 format to TMD3 format
    Updated 04/2023: using renamed astro mean_longitudes function
        using renamed arguments function for nodal corrections
        adding prediction routine for solid earth tides
        output solid earth tide corrections as combined XYZ components
    Updated 03/2023: add basic variable typing to function inputs
    Updated 12/2022: merged prediction functions into a single module
    Updated 05/2022: added ESR netCDF4 formats to list of model types
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 02/2021: replaced numpy bool to prevent deprecation warning
    Updated 09/2020: append output mask over each constituent
    Updated 08/2020: change time variable names to not overwrite functions
    Updated 07/2020: added function docstrings
    Updated 11/2019: output as numpy masked arrays instead of nan-filled arrays
    Updated 09/2019: added netcdf option to CORRECTIONS option
    Updated 08/2018: added correction option ATLAS for localized OTIS solutions
    Updated 07/2018: added option to use GSFC GOT nodal corrections
    Updated 09/2017: Rewritten in Python
"""

from __future__ import annotations

import logging
import numpy as np
import xarray as xr
import pyTMD.astro
import pyTMD.constituents
import pyTMD.interpolate
import pyTMD.math

__all__ = [
    "time_series",
    "infer_minor",
    "_infer_short_period",
    "_infer_semi_diurnal",
    "_infer_diurnal",
    "_infer_long_period",
    "minor_admittance",
    "_admittance_short_period",
    "_admittance_semi_diurnal",
    "_admittance_diurnal",
    "_admittance_long_period",
    "equilibrium_tide",
    "find_peaks",
]

# number of days between MJD and the tide epoch (1992-01-01T00:00:00)
_mjd_tide = 48622.0


def time_series(
    t: float | np.ndarray,
    ds: xr.Dataset,
    **kwargs,
):
    """
    Predict tides from ``Dataset`` at times

    Parameters
    ----------
    t: float or np.ndarray
        Days relative to 1992-01-01T00:00:00
    ds: xarray.Dataset
        Dataset containing tidal harmonic constants
    kwargs: dict
        Keyword arguments for :py:func:`pyTMD.constituents.arguments`

    Returns
    -------
    tpred: xarray.DataArray
        Predicted tidal time series
    """
    # set default keyword arguments
    kwargs.setdefault("corrections", "OTIS")
    # convert time to Modified Julian Days (MJD)
    MJD = t + _mjd_tide
    # list of constituents
    constituents = ds.tmd.constituents
    # load the nodal corrections
    pu, pf, G = pyTMD.constituents.arguments(MJD, constituents, **kwargs)
    # calculate constituent phase angles
    if kwargs["corrections"] in ("OTIS", "ATLAS", "TMD3"):
        # verify time is (at least) 1D
        t = np.atleast_1d(t)
        # load parameters for constituents
        _, p, o, _, _ = pyTMD.constituents._constituent_parameters(constituents)
        # broadcast parameters to time and constituent dimensions
        omega, phase0, t0 = np.broadcast_arrays(
            o[None, :], p[None, :], t[:, None], subok=True
        )
        # calculate phase angle from frequency and phase-0
        # convert angular frequency to radians per day
        theta = 86400.0 * omega * t0 + phase0 + pu
    else:
        # phase angle from arguments
        theta = np.radians(G) + pu
    # dataset of arguments
    arguments = xr.Dataset(
        data_vars=dict(
            u=(["time", "constituent"], pu),
            f=(["time", "constituent"], pf),
            G=(["time", "constituent"], G),
            theta=(["time", "constituent"], np.exp(1j * theta)),
        ),
        coords=dict(time=np.atleast_1d(MJD), constituent=constituents),
    )
    # convert Dataset to DataArray of complex tidal harmonics
    darr = ds.tmd.to_dataarray(constituents=constituents)
    # sum over tidal constituents
    tpred = (
        darr.real * arguments.f * arguments.theta.real
        - darr.imag * arguments.f * arguments.theta.imag
    ).sum(dim="constituent", skipna=False)
    # check if chunks are present
    if hasattr(tpred, "chunks") and bool(tpred.chunks):
        tpred = tpred.chunk(-1).compute()
    # copy units attribute
    tpred.attrs["units"] = ds[constituents[0]].attrs.get("units", None)
    tpred.attrs["constituents"] = constituents
    # return the predicted tides
    return tpred


# PURPOSE: infer the minor corrections from the major constituents
def infer_minor(
    t: float | np.ndarray,
    ds: xr.Dataset,
    **kwargs,
):
    """
    Infer the tidal values for minor constituents using their
    relation with major constituents
    :cite:p:`Doodson:1941td,Schureman:1958ty,Foreman:1989dt,Egbert:2002ge`

    Parameters
    ----------
    t: float or np.ndarray
        Days relative to 1992-01-01T00:00:00
    ds: xarray.Dataset
        Dataset containing major tidal harmonic constants
    deltat: float or np.ndarray, default 0.0
        Time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        Use nodal corrections from OTIS/ATLAS or GOT/FES models
    minor: list or None, default None
        Tidal constituent IDs of minor constituents for inference
    infer_long_period, bool, default True
        Try to infer long period tides from constituents
    raise_exception: bool, default False
        Raise a ``ValueError`` if major constituents are not found

    Returns
    -------
    tinfer: xr.DataArray
        Inferred tidal time series
    """
    # set default keyword arguments
    kwargs.setdefault("deltat", 0.0)
    kwargs.setdefault("corrections", "OTIS")
    kwargs.setdefault("infer_long_period", True)
    kwargs.setdefault("raise_exception", False)
    # list of minor constituents
    kwargs.setdefault("minor", None)
    # infer the minor tidal constituents
    tinfer = 0.0
    constituents = []
    species = []
    # infer short-period tides for minor constituents
    if kwargs["corrections"] in ("GOT",):
        species.extend(["semi_diurnal", "diurnal"])
    else:
        species.append("short_period")
    # infer long-period tides for minor constituents
    if kwargs["infer_long_period"]:
        species.append("long_period")
    # infer minor constituents for each species
    for s in species:
        result = _infer[s](t, ds, **kwargs)
        tinfer += result
        if hasattr(result, "constituents"):
            constituents.extend(result.constituents)
    # check if chunks are present
    if hasattr(tinfer, "chunks") and bool(tinfer.chunks):
        tinfer = tinfer.chunk(-1).compute()
    # update attributes for inferred constituents
    if hasattr(tinfer, "constituents"):
        tinfer.attrs["constituents"] = constituents
    # return the inferred values
    return tinfer


# PURPOSE: infer short-period minor constituents
def _infer_short_period(
    t: float | np.ndarray,
    ds: xr.Dataset,
    **kwargs,
):
    """
    Infer the tidal values for short-period minor constituents
    using their relation with major constituents
    :cite:p:`Egbert:2002ge,Ray:1999vm,Schureman:1958ty`

    For FES corrections, high precision spline coefficients
    are provided by ``pyfes`` :cite:p:`Lyard:2025tr`

    Parameters
    ----------
    t: float or np.ndarray
        Days relative to 1992-01-01T00:00:00
    ds: xarray.Dataset
        Dataset containing major tidal harmonic constants
    deltat: float or np.ndarray, default 0.0
        Time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        Use nodal corrections from OTIS/ATLAS or GOT/FES models
    minor: list or None, default None
        Tidal constituent IDs of minor constituents for inference
    raise_exception: bool, default False
        Raise a ``ValueError`` if major constituents are not found

    Returns
    -------
    tinfer: xr.DataArray
        Inferred tidal time series
    """
    # set default keyword arguments
    kwargs.setdefault("deltat", 0.0)
    kwargs.setdefault("corrections", "OTIS")
    kwargs.setdefault("raise_exception", False)
    # list of minor constituents
    kwargs.setdefault("minor", None)
    # major constituents used for inferring minor tides
    cindex = ["q1", "o1", "p1", "k1", "n2", "m2", "s2", "k2", "2n2"]
    # check that major constituents are in the dataset for inference
    nz = sum([(c in ds.tmd.constituents) for c in cindex])
    # angular frequencies for (available) major constituents
    omajor = pyTMD.constituents.frequency(ds.tmd.constituents, **kwargs)
    # raise exception or log error
    msg = "Not enough constituents to infer short-period tides"
    if (nz < 6) and kwargs["raise_exception"]:
        raise ValueError(msg)
    elif nz < 6:
        logging.debug(msg)
        return 0.0

    # complete list of minor constituents
    minor_constituents = [
        "2q1",
        "sigma1",
        "rho1",
        "m1b",
        "m1",
        "chi1",
        "pi1",
        "phi1",
        "theta1",
        "j1",
        "oo1",
        "2n2",
        "mu2",
        "nu2",
        "lambda2",
        "l2",
        "l2b",
        "t2",
    ]
    # FES models additionally infer eps2 and eta2
    if kwargs["corrections"] in ("FES",):
        minor_constituents.extend(["eps2", "eta2"])
    # possibly reduced list of minor constituents
    minor = kwargs["minor"] or minor_constituents
    # angular frequencies for inferred constituents
    omega = pyTMD.constituents.frequency(minor_constituents, **kwargs)
    # only add minor constituents that are not on the list of major values
    # and with frequencies not equal to any major constituent
    constituents = [
        m
        for i, m in enumerate(minor_constituents)
        if (m not in ds.tmd.constituents)
        and (m in minor)
        and (np.all(omega[i] != omajor))
    ]
    # if there are no constituents to infer
    msg = "No short-period tidal constituents to infer"
    if not any(constituents):
        logging.debug(msg)
        return 0.0

    # relationship between major and minor constituent complex amplitudes
    dmin = xr.Dataset()
    dmin["2q1"] = 0.263 * ds["q1"] - 0.0252 * ds["o1"]
    dmin["sigma1"] = 0.297 * ds["q1"] - 0.0264 * ds["o1"]
    dmin["rho1"] = 0.164 * ds["q1"] + 0.0048 * ds["o1"]
    dmin["m1b"] = 0.0140 * ds["o1"] + 0.0101 * ds["k1"]
    dmin["m1"] = 0.0389 * ds["o1"] + 0.0282 * ds["k1"]
    dmin["chi1"] = 0.0064 * ds["o1"] + 0.0060 * ds["k1"]
    dmin["pi1"] = 0.0030 * ds["o1"] + 0.0171 * ds["k1"]
    dmin["phi1"] = -0.0015 * ds["o1"] + 0.0152 * ds["k1"]
    dmin["theta1"] = -0.0065 * ds["o1"] + 0.0155 * ds["k1"]
    dmin["j1"] = -0.0389 * ds["o1"] + 0.0836 * ds["k1"]
    dmin["oo1"] = -0.0431 * ds["o1"] + 0.0613 * ds["k1"]
    dmin["2n2"] = 0.264 * ds["n2"] - 0.0253 * ds["m2"]
    dmin["mu2"] = 0.298 * ds["n2"] - 0.0264 * ds["m2"]
    dmin["nu2"] = 0.165 * ds["n2"] + 0.00487 * ds["m2"]
    dmin["lambda2"] = 0.0040 * ds["m2"] + 0.0074 * ds["s2"]
    dmin["l2"] = 0.0131 * ds["m2"] + 0.0326 * ds["s2"]
    dmin["l2b"] = 0.0033 * ds["m2"] + 0.0082 * ds["s2"]
    dmin["t2"] = 0.0585 * ds["s2"]
    # additional coefficients for FES models
    if kwargs["corrections"] in ("FES",):
        # spline coefficients for admittances
        mu2 = [0.069439968323, 0.351535557706, -0.046278307672]
        nu2 = [-0.006104695053, 0.156878802427, 0.006755704028]
        l2 = [0.077137765667, -0.051653455134, 0.027869916824]
        t2 = [0.180480173707, -0.020101177502, 0.008331518844]
        la2 = [0.016503557465, -0.013307812292, 0.007753383202]
        dmin["mu2"] = mu2[0] * ds["k2"] + mu2[1] * ds["n2"] + mu2[2] * ds["m2"]
        dmin["nu2"] = nu2[0] * ds["k2"] + nu2[1] * ds["n2"] + nu2[2] * ds["m2"]
        dmin["lambda2"] = (
            la2[0] * ds["k2"] + la2[1] * ds["n2"] + la2[2] * ds["m2"]
        )
        dmin["l2b"] = l2[0] * ds["k2"] + l2[1] * ds["n2"] + l2[2] * ds["m2"]
        dmin["t2"] = t2[0] * ds["k2"] + t2[1] * ds["n2"] + t2[2] * ds["m2"]
        dmin["eps2"] = 0.53285 * ds["2n2"] - 0.03304 * ds["n2"]
        dmin["eta2"] = -0.0034925 * ds["m2"] + 0.0831707 * ds["k2"]

    # convert time to Modified Julian Days (MJD)
    MJD = t + _mjd_tide
    # load the nodal corrections for minor constituents
    pu, pf, G = pyTMD.constituents.minor_arguments(
        MJD, deltat=kwargs["deltat"], corrections=kwargs["corrections"]
    )
    # phase angle from arguments
    theta = np.radians(G) + pu
    # dataset of minor arguments
    arguments = xr.Dataset(
        data_vars=dict(
            u=(["time", "constituent"], pu),
            f=(["time", "constituent"], pf),
            theta=(["time", "constituent"], np.exp(1j * theta)),
        ),
        coords=dict(time=np.atleast_1d(MJD), constituent=minor_constituents),
    )
    # convert Dataset to DataArray of complex tidal harmonics
    # (reduce list to only those constituents to infer)
    darr = dmin.tmd.to_dataarray(constituents=constituents)
    # select argument for constituents
    arg = arguments.sel(constituent=constituents)
    # sum over tidal constituents
    tinfer = (
        darr.real * arg.f * arg.theta.real - darr.imag * arg.f * arg.theta.imag
    ).sum(dim="constituent", skipna=False)
    # copy units attribute
    tinfer.attrs["units"] = ds["q1"].attrs.get("units", None)
    tinfer.attrs["constituents"] = constituents
    # return the inferred values
    return tinfer


# PURPOSE: infer semi-diurnal minor constituents
def _infer_semi_diurnal(
    t: float | np.ndarray,
    ds: xr.Dataset,
    **kwargs,
):
    """
    Infer the tidal values for semi-diurnal minor constituents
    using their relation with major constituents
    :cite:p:`Munk:1966go,Ray:1999vm,Cartwright:1971iz`

    Parameters
    ----------
    t: float or np.ndarray
        Days relative to 1992-01-01T00:00:00
    ds: xarray.Dataset
        Dataset containing major tidal harmonic constants
    deltat: float or np.ndarray, default 0.0
        Time correction for converting to Ephemeris Time (days)
    corrections: str, default 'GOT'
        Use nodal corrections from OTIS/ATLAS or GOT/FES models
    kwargs: dict
        Keyword arguments for
        :py:func:`pyTMD.predict._admittance_semi_diurnal`

    Returns
    -------
    tinfer: xr.DataArray
        Inferred tidal time series
    """
    # set default keyword arguments
    kwargs.setdefault("deltat", 0.0)
    kwargs.setdefault("corrections", "GOT")
    # convert time to Modified Julian Days (MJD)
    MJD = t + _mjd_tide
    # interpolate admittances for semi-diurnal minor constituents
    darr = _admittance_semi_diurnal(ds, **kwargs)
    # check if there are any constituents to infer
    if not hasattr(darr, "coords") or "constituent" not in darr.coords:
        return 0.0
    # list of constituents to infer
    constituents = darr.coords["constituent"].values
    # load the nodal corrections for minor constituents
    pu, pf, G = pyTMD.constituents.arguments(
        MJD,
        constituents,
        deltat=kwargs["deltat"],
        corrections=kwargs["corrections"],
    )
    # phase angle from arguments
    theta = np.radians(G) + pu
    # dataset of arguments
    arg = xr.Dataset(
        data_vars=dict(
            u=(["time", "constituent"], pu),
            f=(["time", "constituent"], pf),
            theta=(["time", "constituent"], np.exp(1j * theta)),
        ),
        coords=dict(time=np.atleast_1d(MJD), constituent=constituents),
    )
    # sum over tidal constituents
    tinfer = (
        darr.real * arg.f * arg.theta.real - darr.imag * arg.f * arg.theta.imag
    ).sum(dim="constituent", skipna=False)
    # copy units attribute
    tinfer.attrs["units"] = ds["n2"].attrs.get("units", None)
    tinfer.attrs["constituents"] = constituents
    # return the inferred values
    return tinfer


# PURPOSE: infer diurnal minor constituents
def _infer_diurnal(
    t: float | np.ndarray,
    ds: xr.Dataset,
    **kwargs,
):
    """
    Infer the tidal values for diurnal minor constituents
    using their relation with major constituents taking into
    account resonance due to free core nutation
    :cite:p:`Munk:1966go,Ray:2017jx,Wahr:1981if,Cartwright:1973em`

    Parameters
    ----------
    t: float or np.ndarray
        Days relative to 1992-01-01T00:00:00
    ds: xarray.Dataset
        Dataset containing major tidal harmonic constants
    deltat: float or np.ndarray, default 0.0
        Time correction for converting to Ephemeris Time (days)
    corrections: str, default 'GOT'
        Use nodal corrections from OTIS/ATLAS or GOT/FES models
    kwargs: dict
        Keyword arguments for
        :py:func:`pyTMD.predict._admittance_diurnal`

    Returns
    -------
    tinfer: xr.DataArray
        Inferred tidal time series
    """
    # set default keyword arguments
    kwargs.setdefault("deltat", 0.0)
    kwargs.setdefault("corrections", "GOT")
    # convert time to Modified Julian Days (MJD)
    MJD = t + _mjd_tide
    # interpolate admittances for diurnal minor constituents
    darr = _admittance_diurnal(ds, **kwargs)
    # check if there are any constituents to infer
    if not hasattr(darr, "coords") or "constituent" not in darr.coords:
        return 0.0
    # list of constituents to infer
    constituents = darr.coords["constituent"].values
    # load the nodal corrections for minor constituents
    pu, pf, G = pyTMD.constituents.arguments(
        MJD,
        constituents,
        deltat=kwargs["deltat"],
        corrections=kwargs["corrections"],
    )
    # phase angle from arguments
    theta = np.radians(G) + pu
    # dataset of arguments
    arg = xr.Dataset(
        data_vars=dict(
            u=(["time", "constituent"], pu),
            f=(["time", "constituent"], pf),
            theta=(["time", "constituent"], np.exp(1j * theta)),
        ),
        coords=dict(time=np.atleast_1d(MJD), constituent=constituents),
    )
    # sum over tidal constituents
    tinfer = (
        darr.real * arg.f * arg.theta.real - darr.imag * arg.f * arg.theta.imag
    ).sum(dim="constituent", skipna=False)
    # copy units attribute
    tinfer.attrs["units"] = ds["q1"].attrs.get("units", None)
    tinfer.attrs["constituents"] = constituents
    # return the inferred values
    return tinfer


# PURPOSE: infer long-period minor constituents
def _infer_long_period(
    t: float | np.ndarray,
    ds: xr.Dataset,
    **kwargs,
):
    """
    Infer the tidal values for long-period minor constituents
    using their relation with major constituents with option to
    take into account variations due to mantle anelasticity
    :cite:p:`Ray:1999vm,Ray:2014fu,Cartwright:1973em,Mathews:2002cr`

    Parameters
    ----------
    t: float or np.ndarray
        Days relative to 1992-01-01T00:00:00
    ds: xarray.Dataset
        Dataset containing major tidal harmonic constants
    deltat: float or np.ndarray, default 0.0
        Time correction for converting to Ephemeris Time (days)
    corrections: str, default 'GOT'
        Use nodal corrections from OTIS/ATLAS or GOT/FES models
    kwargs: dict
        Keyword arguments for
        :py:func:`pyTMD.predict._admittance_long_period`

    Returns
    -------
    tinfer: xr.DataArray
        Inferred tidal time series
    """
    # set default keyword arguments
    kwargs.setdefault("deltat", 0.0)
    kwargs.setdefault("corrections", "GOT")
    # convert time to Modified Julian Days (MJD)
    MJD = t + _mjd_tide
    # interpolate admittances for long-period minor constituents
    darr = _admittance_long_period(ds, **kwargs)
    # check if there are any constituents to infer
    if not hasattr(darr, "coords") or "constituent" not in darr.coords:
        return 0.0
    # list of constituents to infer
    constituents = darr.coords["constituent"].values
    # load the nodal corrections for minor constituents
    pu, pf, G = pyTMD.constituents.arguments(
        MJD,
        constituents,
        deltat=kwargs["deltat"],
        corrections=kwargs["corrections"],
    )
    # phase angle from arguments
    theta = np.radians(G) + pu
    # dataset of arguments
    arg = xr.Dataset(
        data_vars=dict(
            u=(["time", "constituent"], pu),
            f=(["time", "constituent"], pf),
            theta=(["time", "constituent"], np.exp(1j * theta)),
        ),
        coords=dict(time=np.atleast_1d(MJD), constituent=constituents),
    )
    # sum over tidal constituents
    tinfer = (
        darr.real * arg.f * arg.theta.real - darr.imag * arg.f * arg.theta.imag
    ).sum(dim="constituent", skipna=False)
    # copy units attribute
    tinfer.attrs["units"] = ds["node"].attrs.get("units", None)
    tinfer.attrs["constituents"] = constituents
    # return the inferred values
    return tinfer


# PURPOSE: interpolate admittances for minor constituents
def minor_admittance(
    ds: xr.Dataset,
    **kwargs,
):
    """
    Interpolate admittances from major constituents to a set of minor
    constituents :cite:p:`Doodson:1941td,Schureman:1958ty,Foreman:1989dt`

    Parameters
    ----------
    ds: xarray.Dataset
        Dataset containing major tidal harmonic constants
    corrections: str, default 'OTIS'
        Use nodal corrections from OTIS/ATLAS or GOT/FES models
    minor: list or None, default None
        Tidal constituent IDs of minor constituents for inference
    infer_long_period, bool, default True
        Try to interpolate long period tidal admittances
    raise_exception: bool, default False
        Raise a ``ValueError`` if major constituents are not found

    Returns
    -------
    other: xarray.Dataset
        Dataset containing interpolated minor tidal harmonic constants
    """
    # set default keyword arguments
    kwargs.setdefault("corrections", "OTIS")
    kwargs.setdefault("infer_long_period", True)
    kwargs.setdefault("raise_exception", False)
    # list of minor constituents
    kwargs.setdefault("minor", None)
    # interpolate the minor tidal constituents
    other = xr.Dataset(coords=ds.coords, attrs=ds.attrs)
    constituents = []
    species = []
    # interpolate short-period admittances
    if kwargs["corrections"] in ("GOT",):
        species.extend(["semi_diurnal", "diurnal"])
    else:
        species.append("short_period")
    # interpolate long-period admittances
    if kwargs["infer_long_period"]:
        species.append("long_period")
    # interpolate admittances for minor constituents of each species
    for s in species:
        darr = _admittance[s](ds, **kwargs)
        if hasattr(darr, "coords") and "constituent" in darr.coords:
            constituents.extend(darr.coords["constituent"].values)
            tmp = darr.to_dataset(dim="constituent")
            other = xr.merge([other, tmp], join="outer", compat="override")
    # add constituents attribute to the dataset
    other.attrs["constituents"] = constituents
    # return xarray dataset
    return other


# PURPOSE: interpolate admittances for short-period minor constituents
def _admittance_short_period(
    ds: xr.Dataset,
    **kwargs,
):
    """
    Interpolate admittances from short-period major constituents to a set of
    minor constituents :cite:p:`Egbert:2002ge,Ray:1999vm,Schureman:1958ty`

    For FES corrections, high precision spline coefficients
    are provided by ``pyfes`` :cite:p:`Lyard:2025tr`

    Parameters
    ----------
    ds: xarray.Dataset
        Dataset containing major tidal harmonic constants
    corrections: str, default 'OTIS'
        Use nodal corrections from OTIS/ATLAS or GOT/FES models
    minor: list or None, default None
        Tidal constituent IDs of minor constituents for inference
    raise_exception: bool, default False
        Raise a ``ValueError`` if major constituents are not found

    Returns
    -------
    darr: xr.DataArray
        Minor constituent harmonic constants
    """
    # set default keyword arguments
    kwargs.setdefault("corrections", "OTIS")
    kwargs.setdefault("raise_exception", False)
    # list of minor constituents
    kwargs.setdefault("minor", None)
    # major constituents used for inferring minor tides
    cindex = ["q1", "o1", "k1", "n2", "m2", "s2", "k2", "2n2"]
    # check that major constituents are in the dataset for inference
    nz = sum([(c in ds.tmd.constituents) for c in cindex])
    # angular frequencies for (available) major constituents
    omajor = pyTMD.constituents.frequency(ds.tmd.constituents, **kwargs)
    # raise exception or log error
    msg = "Not enough constituents to infer short-period tides"
    if (nz < 6) and kwargs["raise_exception"]:
        raise ValueError(msg)
    elif nz < 6:
        logging.debug(msg)
        return None
    # major constituents as a dataarray
    dmajor = ds.tmd.to_dataarray().rename(constituent="major")

    # complete list of minor constituents
    minor_constituents = [
        "2q1",
        "sigma1",
        "rho1",
        "m1b",
        "m1",
        "chi1",
        "pi1",
        "phi1",
        "theta1",
        "j1",
        "oo1",
        "2n2",
        "mu2",
        "nu2",
        "lambda2",
        "l2",
        "l2b",
        "t2",
    ]
    # FES models additionally infer eps2 and eta2
    if kwargs["corrections"] in ("FES",):
        minor_constituents.extend(["eps2", "eta2"])
    # possibly reduced list of minor constituents
    minor = kwargs["minor"] or minor_constituents
    # angular frequencies for minor constituents
    omega = pyTMD.constituents.frequency(minor_constituents, **kwargs)

    # relationship between major and minor constituent complex amplitudes
    S = {}
    # admittance interpolation coefficients
    # 2q1: 0.263 * q1 - 0.0252 * o1
    S["2q1"] = [0.263, -0.0252, 0, 0, 0, 0, 0, 0]
    # sigma1: 0.297 * q1 - 0.0264 * o1
    S["sigma1"] = [0.297, -0.0264, 0, 0, 0, 0, 0, 0]
    # rho1: 0.164 * q1 + 0.0048 * o1
    S["rho1"] = [0.164, 0.0048, 0, 0, 0, 0, 0, 0]
    # m1b: 0.0140 * o1 + 0.0101 * k1
    S["m1b"] = [0, 0.0140, 0.0101, 0, 0, 0, 0, 0]
    # m1: 0.0389 * o1 + 0.0282 * k1
    S["m1"] = [0, 0.0389, 0.0282, 0, 0, 0, 0, 0]
    # chi1: 0.0064 * o1 + 0.0060 * k1
    S["chi1"] = [0, 0.0064, 0.0060, 0, 0, 0, 0, 0]
    # pi1: 0.0030 * o1 + 0.0171 * k1
    S["pi1"] = [0, 0.0030, 0.0171, 0, 0, 0, 0, 0]
    # phi1: -0.0015 * o1 + 0.0152 * k1
    S["phi1"] = [0, -0.0015, 0.0152, 0, 0, 0, 0, 0]
    # theta1: -0.0065 * o1 + 0.0155 * k1
    S["theta1"] = [0, -0.0065, 0.0155, 0, 0, 0, 0, 0]
    # j1: -0.0389 * o1 + 0.0836 * k1
    S["j1"] = [0, -0.0389, 0.0836, 0, 0, 0, 0, 0]
    # oo1: -0.0431 * o1 + 0.0613 * k1
    S["oo1"] = [0, -0.0431, 0.0613, 0, 0, 0, 0, 0]
    # 2n2: 0.264 * n2 - 0.0253 * m2
    S["2n2"] = [0, 0, 0, 0.264, -0.0253, 0, 0, 0]
    # mu2: 0.298 * n2 - 0.0264 * m2
    S["mu2"] = [0, 0, 0, 0.298, -0.0264, 0, 0, 0]
    # nu2: 0.165 * n2 + 0.00487 * m2
    S["nu2"] = [0, 0, 0, 0.165, 0.00487, 0, 0, 0]
    # lambda2: 0.0040 * m2 + 0.0074 * s2
    S["lambda2"] = [0, 0, 0, 0, 0.0040, 0.0074, 0, 0]
    # l2: 0.0131 * m2 + 0.0326 * s2
    S["l2"] = [0, 0, 0, 0, 0.0131, 0.0326, 0, 0]
    # l2b: 0.0033 * m2 + 0.0082 * s2
    S["l2b"] = [0, 0, 0, 0, 0.0033, 0.0082, 0, 0]
    # t2: 0.0585 * s2
    S["t2"] = [0, 0, 0, 0, 0, 0.0585, 0, 0]
    # update coefficients for FES models
    if kwargs["corrections"] in ("FES",):
        # spline coefficients for admittances
        # semi_diurnal_minor: spl(1) * n2 + spl(2) * m2 + spl(3) * k2
        mu2 = [0.351535557706, -0.046278307672, 0.069439968323]
        nu2 = [0.156878802427, 0.006755704028, -0.006104695053]
        la2 = [-0.013307812292, 0.007753383202, 0.016503557465]
        l2 = [-0.051653455134, 0.027869916824, 0.077137765667]
        t2 = [-0.020101177502, 0.008331518844, 0.180480173707]
        S["mu2"] = [0, 0, 0, mu2[0], mu2[1], 0, mu2[2], 0]
        S["nu2"] = [0, 0, 0, nu2[0], nu2[1], 0, nu2[2], 0]
        S["lambda2"] = [0, 0, 0, la2[0], la2[1], 0, la2[2], 0]
        S["l2b"] = [0, 0, 0, l2[0], l2[1], 0, l2[2], 0]
        S["t2"] = [0, 0, 0, t2[0], t2[1], 0, t2[2], 0]
        # eps2: -0.03304 * n2 + 0.53285 * 2n2
        S["eps2"] = [0, 0, 0, -0.03304, 0, 0, 0, 0.53285]
        # eta2: -0.0034925 * m2 + 0.0831707 * k2
        S["eta2"] = [0, 0, 0, 0, -0.0034925, 0, 0.0831707, 0]

    # inferred constituents
    constituents = []
    # admittance interpolation coefficients
    coefficients = []
    # only add minor constituents that are not on the list of major values
    # and with frequencies not equal to any major constituent
    for i, c in enumerate(minor_constituents):
        if c in ds.tmd.constituents:
            continue
        elif c not in minor:
            continue
        elif np.any(omajor == omega[i]):
            continue
        # append to list
        constituents.append(c)
        # spline coefficients for constituent
        coefficients.append(S[c])

    # if there are no constituents to infer
    msg = "No short-period tidal constituents to infer"
    if not any(constituents):
        logging.debug(msg)
        return None

    # create dictionary of admittances
    admit = dict(dims=("constituent", "major"), coords={})
    admit["coords"]["constituent"] = dict(dims="constituent")
    admit["coords"]["constituent"]["data"] = constituents
    admit["coords"]["major"] = dict(dims="major")
    admit["coords"]["major"]["data"] = cindex
    admit["data"] = coefficients
    # convert admittances to xarray dataarray
    admit = xr.DataArray.from_dict(admit)
    # minor constituent complex amplitudes from major constituents
    darr = admit.dot(dmajor)
    return darr


# PURPOSE: interpolate admittances for semi-diurnal minor constituents
def _admittance_semi_diurnal(
    ds: xr.Dataset,
    **kwargs,
):
    """
    Interpolate admittances from semi-diurnal major constituents to a set of
    minor constituents :cite:p:`Munk:1966go,Ray:1999vm,Cartwright:1971iz`

    Parameters
    ----------
    ds: xarray.Dataset
        Dataset containing major tidal harmonic constants
    minor: list or None, default None
        Tidal constituent IDs of minor constituents for inference
    method: str, default 'linear'
        Method for interpolating between major constituents

            * ``'linear'``: linear interpolation
            * ``'admittance'``: Munk-Cartwright interpolation
    raise_exception: bool, default False
        Raise a ``ValueError`` if major constituents are not found

    Returns
    -------
    darr: xr.DataArray
        Minor constituent harmonic constants
    """
    # set default keyword arguments
    kwargs.setdefault("corrections", "GOT")
    kwargs.setdefault("method", "linear")
    kwargs.setdefault("raise_exception", False)
    # list of minor constituents
    kwargs.setdefault("minor", None)
    # validate interpolation method
    if kwargs["method"].lower() not in ("linear", "admittance"):
        raise ValueError("Invalid interpolation method")
    # major constituents used for inferring semi-diurnal minor tides
    # pivot waves listed in Table 6.7 of the 2010 IERS Conventions
    cindex = ["n2", "m2", "s2"]
    # check that major constituents are in the dataset for inference
    nz = sum([(c in ds.tmd.constituents) for c in cindex])
    # raise exception or log error
    msg = "Not enough constituents to infer semi-diurnal tides"
    if (nz < 3) and kwargs["raise_exception"]:
        raise ValueError(msg)
    elif nz < 3:
        logging.debug(msg)
        return None

    # angular frequencies for major constituents
    omajor = pyTMD.constituents.frequency(cindex, **kwargs)
    # Cartwright and Edden potential amplitudes for major constituents
    amajor = np.zeros(3)
    amajor[0] = 0.121006  # n2
    amajor[1] = 0.631931  # m2
    amajor[2] = 0.294019  # s2
    # "normalize" tide values
    dnorm = xr.Dataset()
    for i, c in enumerate(cindex):
        dnorm[c] = ds[c] / amajor[i]
    # major constituents as a dataarray
    z = dnorm.tmd.to_dataarray()

    # complete list of minor constituents
    minor_constituents = [
        "eps2",
        "2n2",
        "mu2",
        "nu2",
        "gamma2",
        "alpha2",
        "beta2",
        "delta2",
        "lambda2",
        "l2",
        "t2",
        "r2",
        "k2",
        "eta2",
    ]
    # possibly reduced list of minor constituents
    minor = kwargs["minor"] or minor_constituents
    # angular frequencies for minor constituents
    omin = pyTMD.constituents.frequency(minor_constituents, **kwargs)

    # Cartwright and Edden potential amplitudes for minor constituents
    amin = np.zeros(14)
    amin[0] = 0.004669  # eps2
    amin[1] = 0.016011  # 2n2
    amin[2] = 0.019316  # mu2
    amin[3] = 0.022983  # nu2
    amin[4] = 0.001902  # gamma2
    amin[5] = 0.002178  # alpha2
    amin[6] = 0.001921  # beta2
    amin[7] = 0.000714  # delta2
    amin[8] = 0.004662  # lambda2
    amin[9] = 0.017862  # l2
    amin[10] = 0.017180  # t2
    amin[11] = 0.002463  # r2
    amin[12] = 0.079924  # k2
    amin[13] = 0.004467  # eta

    # inferred constituents
    constituents = []
    # Cartwright and Edden potential amplitudes
    amplitude = []
    # angular frequencies
    omega = []
    # only add minor constituents that are not on the list of major values
    # and with frequencies not equal to any major constituent
    for i, c in enumerate(minor_constituents):
        if c in ds.tmd.constituents:
            continue
        elif c not in minor:
            continue
        elif np.any(omajor == omin[i]):
            continue
        # append to list
        constituents.append(c)
        # amplitude for constituent
        amplitude.append(amin[i])
        # angular frequencies for inferred constituents
        omega.append(omin[i])

    # if there are no constituents to infer
    msg = "No semi-diurnal tidal constituents to infer"
    if not any(constituents):
        logging.debug(msg)
        return None

    # convert minor constituent parameters to dataset
    arg = xr.Dataset(
        data_vars=dict(
            amplitude=(["constituent"], amplitude),
            omega=(["constituent"], omega),
        ),
        coords=dict(constituent=constituents),
    )

    # interpolate from major constituents
    if kwargs["method"].lower() == "linear":
        # linearly interpolate using constituent frequencies
        zmin = pyTMD.interpolate.interp1d(arg.omega.values, omajor, z)
        coords = z.coords.assign(dict(constituent=arg.constituent))
        zmin = xr.DataArray(zmin, dims=z.dims, coords=coords)
    elif kwargs["method"].lower() == "admittance":
        # admittance interpolation using Munk-Cartwright approach
        # coefficients for Munk-Cartwright admittance interpolation
        Ainv = xr.DataArray(
            [
                [3.3133, -4.2538, 1.9405],
                [-3.3133, 4.2538, -0.9405],
                [1.5018, -3.2579, 1.7561],
            ],
            coords=[np.arange(3), cindex],
            dims=["coefficient", "constituent"],
        )
        coef = Ainv.dot(z)
        # convert frequency to radians per 48 hours
        # following Munk and Cartwright (1966)
        f = np.exp(2.0 * 86400.0 * arg.omega * 1j)
        # calculate interpolated values for constituent
        zmin = coef[0] + coef[1] * f.real + coef[2] * f.imag
    # rescale tide values
    darr = arg.amplitude * zmin
    return darr


# PURPOSE: interpolate admittances for diurnal minor constituents
def _admittance_diurnal(
    ds: xr.Dataset,
    **kwargs,
):
    """
    Interpolate admittances from diurnal major constituents to a set of
    minor constituents taking into account resonance due to free core nutation
    :cite:p:`Munk:1966go,Ray:2017jx,Wahr:1981if,Cartwright:1973em`

    Parameters
    ----------
    ds: xarray.Dataset
        Dataset containing major tidal harmonic constants
    minor: list or None, default None
        Tidal constituent IDs of minor constituents for inference
    method: str, default 'linear'
        Method for interpolating between major constituents

            * ``'linear'``: linear interpolation
            * ``'admittance'``: Munk-Cartwright interpolation
    raise_exception: bool, default False
        Raise a ``ValueError`` if major constituents are not found

    Returns
    -------
    darr: xr.DataArray
        Minor constituent harmonic constants
    """
    # set default keyword arguments
    kwargs.setdefault("corrections", "GOT")
    kwargs.setdefault("method", "linear")
    kwargs.setdefault("raise_exception", False)
    # list of minor constituents
    kwargs.setdefault("minor", None)
    # validate interpolation method
    if kwargs["method"].lower() not in ("linear", "admittance"):
        raise ValueError("Invalid interpolation method")
    # major constituents used for inferring diurnal minor tides
    # pivot waves listed in Table 6.7 of the 2010 IERS Conventions
    cindex = ["q1", "o1", "k1"]
    # check that major constituents are in the dataset for inference
    nz = sum([(c in ds.tmd.constituents) for c in cindex])
    # raise exception or log error
    msg = "Not enough constituents to infer diurnal tides"
    if (nz < 3) and kwargs["raise_exception"]:
        raise ValueError(msg)
    elif nz < 3:
        logging.debug(msg)
        return None

    # angular frequencies for major constituents
    omajor = pyTMD.constituents.frequency(cindex, **kwargs)
    # Cartwright and Edden potential amplitudes for major constituents
    amajor = np.zeros(3)
    amajor[0] = 0.050184  # q1
    amajor[1] = 0.262163  # o1
    amajor[2] = 0.368731  # k1
    # "normalize" tide values
    dnorm = xr.Dataset()
    for i, c in enumerate(cindex):
        # Love numbers of degree 2 for constituent
        h2, k2, l2 = pyTMD.earth.love_numbers(omajor[i])
        # tilt factor: response with respect to the solid earth
        gamma_2 = 1.0 + k2 - h2
        dnorm[c] = ds[c] / (amajor[i] * gamma_2)
    # major constituents as a dataarray
    z = dnorm.tmd.to_dataarray()

    # complete list of minor constituents
    minor_constituents = [
        "2q1",
        "sigma1",
        "rho1",
        "tau1",
        "beta1",
        "m1a",
        "m1b",
        "chi1",
        "pi1",
        "p1",
        "psi1",
        "phi1",
        "theta1",
        "j1",
        "so1",
        "oo1",
        "ups1",
    ]
    # possibly reduced list of minor constituents
    minor = kwargs["minor"] or minor_constituents
    # angular frequencies for minor constituents
    omin = pyTMD.constituents.frequency(minor_constituents, **kwargs)

    # Cartwright and Edden potential amplitudes for minor constituents
    amin = np.zeros(17)
    amin[0] = 0.006638  # 2q1
    amin[1] = 0.008023  # sigma1
    amin[2] = 0.009540  # rho1
    amin[3] = 0.003430  # tau1
    amin[4] = 0.001941  # beta1
    amin[5] = 0.020604  # m1a
    amin[6] = 0.007420  # m1b
    amin[7] = 0.003925  # chi1
    amin[8] = 0.007125  # pi1
    amin[9] = 0.122008  # p1
    amin[10] = 0.002929  # psi1
    amin[11] = 0.005247  # phi1
    amin[12] = 0.003966  # theta1
    amin[13] = 0.020618  # j1
    amin[14] = 0.003417  # so1
    amin[15] = 0.011293  # oo1
    amin[16] = 0.002157  # ups1

    # inferred constituents
    constituents = []
    # Cartwright and Edden potential amplitudes
    amplitude = []
    # angular frequencies
    omega = []
    # tilt factors
    gamma_2 = []
    # only add minor constituents that are not on the list of major values
    # and with frequencies not equal to any major constituent
    for i, c in enumerate(minor_constituents):
        if c in ds.tmd.constituents:
            continue
        elif c not in minor:
            continue
        elif np.any(omajor == omin[i]):
            continue
        # append to list
        constituents.append(c)
        # amplitude for constituent
        amplitude.append(amin[i])
        # angular frequencies for inferred constituents
        omega.append(omin[i])
        # Love numbers of degree 2 for constituent
        h2, k2, l2 = pyTMD.earth.love_numbers(omin[i])
        # tilt factor: response with respect to the solid earth
        gamma_2.append(1.0 + k2 - h2)

    # if there are no constituents to infer
    msg = "No diurnal tidal constituents to infer"
    if not any(constituents):
        logging.debug(msg)
        return None

    # convert minor constituent parameters to dataset
    arg = xr.Dataset(
        data_vars=dict(
            amplitude=(["constituent"], amplitude),
            omega=(["constituent"], omega),
            gamma_2=(["constituent"], gamma_2),
        ),
        coords=dict(constituent=constituents),
    )

    # interpolate from major constituents
    if kwargs["method"].lower() == "linear":
        # linearly interpolate using constituent frequencies
        zmin = pyTMD.interpolate.interp1d(arg.omega.values, omajor, z)
        coords = z.coords.assign(dict(constituent=arg.constituent))
        zmin = xr.DataArray(zmin, dims=z.dims, coords=coords)
    elif kwargs["method"].lower() == "admittance":
        # admittance interpolation using Munk-Cartwright approach
        # coefficients for Munk-Cartwright admittance interpolation
        Ainv = xr.DataArray(
            [
                [3.1214, -3.8494, 1.728],
                [-3.1727, 3.9559, -0.7832],
                [1.438, -3.0297, 1.5917],
            ],
            coords=[np.arange(3), cindex],
            dims=["coefficient", "constituent"],
        )
        coef = Ainv.dot(z)
        # convert frequency to radians per 48 hours
        # following Munk and Cartwright (1966)
        f = np.exp(2.0 * 86400.0 * arg.omega * 1j)
        # calculate interpolated values for constituent
        zmin = coef[0] + coef[1] * f.real + coef[2] * f.imag
    # rescale tide values
    darr = arg.amplitude * arg.gamma_2 * zmin
    return darr


# PURPOSE: interpolate admittances for long-period minor constituents
def _admittance_long_period(
    ds: xr.Dataset,
    **kwargs,
):
    """
    Interpolate admittances from long-period major constituents to a set of
    minor constituents taking into account variations due to mantle anelasticity
    :cite:p:`Ray:1999vm,Ray:2014fu,Cartwright:1973em,Mathews:2002cr`

    Parameters
    ----------
    ds: xarray.Dataset
        Dataset containing major tidal harmonic constants
    minor: list or None, default None
        Tidal constituent IDs of minor constituents for inference
    include_anelasticity: bool, default False
        Compute Love numbers taking into account mantle anelasticity
    raise_exception: bool, default False
        Raise a ``ValueError`` if major constituents are not found

    Returns
    -------
    darr: xr.DataArray
        Minor constituent harmonic constants
    """
    kwargs.setdefault("corrections", "GOT")
    kwargs.setdefault("include_anelasticity", False)
    kwargs.setdefault("raise_exception", False)
    # list of minor constituents
    kwargs.setdefault("minor", None)
    # major constituents used for inferring long period minor tides
    # pivot waves listed in Table 6.7 of the 2010 IERS Conventions
    cindex = ["node", "mm", "mf"]
    # check that major constituents are in the dataset for inference
    nz = sum([(c in ds.tmd.constituents) for c in cindex])
    # raise exception or log error
    msg = "Not enough constituents to infer long-period tides"
    if (nz < 3) and kwargs["raise_exception"]:
        raise ValueError(msg)
    elif nz < 3:
        logging.debug(msg)
        return None

    # angular frequencies for major constituents
    omajor = pyTMD.constituents.frequency(cindex, **kwargs)
    # Cartwright and Edden potential amplitudes for major constituents
    amajor = np.zeros(3)
    amajor[0] = 0.027929  # node
    amajor[1] = 0.035184  # mm
    amajor[2] = 0.066607  # mf
    # "normalize" tide values
    dnorm = xr.Dataset()
    for i, c in enumerate(cindex):
        # complex Love numbers of degree 2 for long-period band
        if kwargs["include_anelasticity"]:
            # include variations largely due to mantle anelasticity
            h2, k2, l2 = pyTMD.earth.complex_love_numbers(omajor[i])
        else:
            # Love numbers for long-period tides (Wahr, 1981)
            h2, k2, l2 = pyTMD.earth.love_numbers(
                omajor[i], astype=np.complex128
            )
        # tilt factor: response with respect to the solid earth
        # use real components from Mathews et al. (2002)
        gamma_2 = 1.0 + k2.real - h2.real
        dnorm[c] = ds[c] / (amajor[i] * gamma_2)
    # major constituents as a dataarray
    z = dnorm.tmd.to_dataarray()

    # complete list of minor constituents
    minor_constituents = [
        "sa",
        "ssa",
        "sta",
        "msm",
        "msf",
        "mst",
        "mt",
        "msqm",
        "mq",
    ]
    # possibly reduced list of minor constituents
    minor = kwargs["minor"] or minor_constituents
    # angular frequencies for minor constituents
    omin = pyTMD.constituents.frequency(minor_constituents, **kwargs)

    # Cartwright and Edden potential amplitudes for minor constituents
    amin = np.zeros(9)
    amin[0] = 0.004922  # sa
    amin[1] = 0.030988  # ssa
    amin[2] = 0.001809  # sta
    amin[3] = 0.006728  # msm
    amin[4] = 0.005837  # msf
    amin[5] = 0.002422  # mst
    amin[6] = 0.012753  # mt
    amin[7] = 0.002037  # msqm
    amin[8] = 0.001687  # mq

    # inferred constituents
    constituents = []
    # Cartwright and Edden potential amplitudes
    amplitude = []
    # angular frequencies
    omega = []
    # tilt factors
    gamma_2 = []
    # only add minor constituents that are not on the list of major values
    # and with frequencies not equal to any major constituent
    for i, c in enumerate(minor_constituents):
        if c in ds.tmd.constituents:
            continue
        elif c not in minor:
            continue
        elif np.any(omajor == omin[i]):
            continue
        # append to list
        constituents.append(c)
        # amplitude for constituent
        amplitude.append(amin[i])
        # angular frequencies for inferred constituents
        omega.append(omin[i])
        # complex Love numbers of degree 2 for long-period band
        if kwargs["include_anelasticity"]:
            # include variations largely due to mantle anelasticity
            h2, k2, l2 = pyTMD.earth.complex_love_numbers(omin[i])
        else:
            # Love numbers for long-period tides (Wahr, 1981)
            h2, k2, l2 = pyTMD.earth.love_numbers(omin[i], astype=np.complex128)
        # tilt factor: response with respect to the solid earth
        # use real components from Mathews et al. (2002)
        gamma_2.append(1.0 + k2.real - h2.real)

    # if there are no constituents to infer
    msg = "No long-period tidal constituents to infer"
    if not any(constituents):
        logging.debug(msg)
        return None

    # convert minor constituent parameters to dataset
    arg = xr.Dataset(
        data_vars=dict(
            amplitude=(["constituent"], amplitude),
            omega=(["constituent"], omega),
            gamma_2=(["constituent"], gamma_2),
        ),
        coords=dict(constituent=constituents),
    )

    # linearly interpolate using constituent frequencies
    zmin = pyTMD.interpolate.interp1d(arg.omega.values, omajor, z)
    coords = z.coords.assign(dict(constituent=arg.constituent))
    zmin = xr.DataArray(zmin, dims=z.dims, coords=coords)
    # rescale tide values
    darr = arg.amplitude * arg.gamma_2 * zmin
    return darr


# dictionary of functions for inferring minor tidal constituents
_infer = {}
_infer["short_period"] = _infer_short_period
_infer["semi_diurnal"] = _infer_semi_diurnal
_infer["diurnal"] = _infer_diurnal
_infer["long_period"] = _infer_long_period
# dictionary of functions for interpolating admittances
_admittance = {}
_admittance["short_period"] = _admittance_short_period
_admittance["semi_diurnal"] = _admittance_semi_diurnal
_admittance["diurnal"] = _admittance_diurnal
_admittance["long_period"] = _admittance_long_period


# PURPOSE: estimate long-period equilibrium tides
def equilibrium_tide(
    t: np.ndarray,
    ds: xr.Dataset,
    **kwargs,
):
    """
    Compute the long-period equilibrium tides the summation of fifteen
    tidal spectral lines from Cartwright-Tayler-Edden tables
    :cite:p:`Cartwright:1971iz,Cartwright:1973em`

    Parameters
    ----------
    t: np.ndarray
        Days relative to 1992-01-01T00:00:00
    ds: xarray.Dataset
        Dataset with spatial coordinates
    deltat: float or np.ndarray, default 0.0
        Time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        Use nodal corrections from OTIS/ATLAS or GOT/FES models
    include_anelasticity: bool, default False
        Compute Love numbers taking into account mantle anelasticity
    constituents: list
        Long-period tidal constituent IDs

    Returns
    -------
    tpred: xr.DataArray
        Predicted tidal time series (meters)
    """
    # set default keyword arguments
    cindex = [
        "node",
        "sa",
        "ssa",
        "msm",
        "065.445",
        "mm",
        "065.465",
        "msf",
        "075.355",
        "mf",
        "mf+",
        "075.575",
        "mst",
        "mt",
        "085.465",
    ]
    kwargs.setdefault("constituents", cindex)
    kwargs.setdefault("deltat", 0.0)
    kwargs.setdefault("include_anelasticity", False)
    kwargs.setdefault("corrections", "OTIS")

    # number of constituents
    nc = len(cindex)
    # set function for astronomical longitudes
    # use ASTRO5 routines if not using an OTIS type model
    if kwargs["corrections"] in ("OTIS", "ATLAS", "TMD3", "netcdf"):
        method = "Cartwright"
    else:
        method = "ASTRO5"
    # convert time to Modified Julian Days (MJD)
    MJD = t + _mjd_tide
    # compute principal mean longitudes
    s, h, p, n, pp = pyTMD.astro.mean_longitudes(
        MJD + kwargs["deltat"], method=method
    )
    # initial time conversions
    hour = 24.0 * np.mod(MJD, 1)
    # convert from hours solar time into mean lunar time in degrees
    tau = 15.0 * hour - s + h
    # variable for multiples of 90 degrees (Ray technical note 2017)
    # full expansion of Equilibrium Tide includes some negative cosine
    # terms and some sine terms (Pugh and Woodworth, 2014)
    k = 90.0 + np.zeros_like(MJD)
    # convert to negative mean longitude of the ascending node (N')
    Np = pyTMD.math.normalize_angle(360.0 - n)
    # determine equilibrium arguments
    fargs = np.c_[tau, s, h, p, Np, pp, k]

    # Cartwright and Edden potential amplitudes (centimeters)
    # assemble long-period tide potential from 15 CTE terms greater than 1 mm
    amajor = xr.Dataset()
    # group 0,0
    # nodal term is included but not the constant term.
    amajor["node"] = 2.7929  # node
    amajor["sa"] = -0.4922  # sa
    amajor["ssa"] = -3.0988  # ssa
    # group 0,1
    amajor["msm"] = -0.6728  # msm
    amajor["065.445"] = 0.231
    amajor["mm"] = -3.5184  # mm
    amajor["065.465"] = 0.228
    # group 0,2
    amajor["msf"] = -0.5837  # msf
    amajor["075.355"] = -0.288
    amajor["mf"] = -6.6607  # mf
    amajor["mf+"] = -2.763  # mf+
    amajor["075.575"] = -0.258
    # group 0,3
    amajor["mst"] = -0.2422  # mst
    amajor["mt"] = -1.2753  # mt
    amajor["085.465"] = -0.528

    # set constituents to be iterable and lower case
    if isinstance(kwargs["constituents"], str):
        constituents = [kwargs["constituents"].lower()]
    else:
        constituents = [c.lower() for c in kwargs["constituents"]]

    # Doodson coefficients for 15 long-period terms
    coef = np.zeros((7, nc))
    # group 0,0
    coef[:, 0] = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]  # node
    coef[:, 1] = [0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0]  # sa
    coef[:, 2] = [0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0]  # ssa
    # group 0,1
    coef[:, 3] = [0.0, 1.0, -2.0, 1.0, 0.0, 0.0, 0.0]  # msm
    coef[:, 4] = [0.0, 1.0, 0.0, -1.0, -1.0, 0.0, 0.0]
    coef[:, 5] = [0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0]  # mm
    coef[:, 6] = [0.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0]
    # group 0,2
    coef[:, 7] = [0.0, 2.0, -2.0, 0.0, 0.0, 0.0, 0.0]  # msf
    coef[:, 8] = [0.0, 2.0, 0.0, -2.0, 0.0, 0.0, 0.0]
    coef[:, 9] = [0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # mf
    coef[:, 10] = [0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0]  # mf+
    coef[:, 11] = [0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0]
    # group 0,3
    coef[:, 12] = [0.0, 3.0, -2.0, 1.0, 0.0, 0.0, 0.0]  # mst
    coef[:, 13] = [0.0, 3.0, 0.0, -1.0, 0.0, 0.0, 0.0]  # mt
    coef[:, 14] = [0.0, 3.0, 0.0, -1.0, 1.0, 0.0, 0.0]

    # spherical harmonic degree and order
    l = 2
    m = 0
    # colatitude in radians
    theta = np.radians(90.0 - ds.y)
    # degree dependent normalization (4-pi)
    dfactor = np.sqrt((2.0 * l + 1.0) / (4.0 * np.pi))
    # 2nd degree Legendre polynomials
    Plm = pyTMD.math._assoc_legendre(l, m, np.cos(theta))
    P20 = dfactor * Plm.real

    # calculate tilt factors for each constituent
    gamma_2 = np.zeros(nc)
    for i, c in enumerate(cindex):
        # calculate angular frequencies of constituents
        omega = pyTMD.constituents._frequency(coef[:, i])
        # complex Love numbers of degree 2 for long-period band
        if kwargs["include_anelasticity"]:
            # include variations largely due to mantle anelasticity
            h2, k2, l2 = pyTMD.earth.complex_love_numbers(omega)
        else:
            # Love numbers for long-period tides (Wahr, 1981)
            h2, k2, l2 = pyTMD.earth.love_numbers(omega, astype=np.complex128)
        # tilt factor: response with respect to the solid earth
        # use real components from Mathews et al. (2002)
        gamma_2[i] = 1.0 + k2.real - h2.real

    # determine equilibrium arguments
    G = np.radians(np.dot(fargs, coef))
    # dataset of arguments
    arguments = xr.Dataset(
        data_vars=dict(
            G=(["time", "constituent"], G), gamma_2=(["constituent"], gamma_2)
        ),
        coords=dict(time=np.atleast_1d(MJD), constituent=cindex),
    )
    # reduce to selected constituents
    arg = arguments.sel(constituent=constituents)
    # convert dataset to dataarray of complex tidal elevations
    darr = amajor.tmd.to_dataarray(constituents=constituents)
    # sum equilibrium tide elevations
    tpred = (P20 * darr * arg.gamma_2 * np.cos(arg.G)).sum(
        dim="constituent", skipna=False
    )
    # add units attribute
    tpred.attrs["units"] = "centimeters"
    tpred.attrs["constituents"] = constituents
    # return the long-period equilibrium tides
    return tpred.tmd.to_units("meters")


def find_peaks(
    darr: xr.DataArray,
    dim: str = "time",
    **kwargs,
):
    """
    Find peaks in an xarray ``DataArray`` using a first order
    differentiation method

    Parameters
    ----------
    darr: xarray.DataArray
        Input ``DataArray`` containing a signal with peaks
    dim: str, default 'time'
        Dimension along which to find peaks
    kwargs: dict
        Keyword arguments for ``xarray.DataArray.differentiate``

    Returns
    -------
    high_peaks: xarray.DataArray
        Boolean array indicating locations of high peaks
    low_peaks: xarray.DataArray
        Boolean array indicating locations of low peaks
    """
    # differentiate to calculate high and low peaks
    diff = darr.differentiate(dim, **kwargs)
    # look for zero crossings in the derivative to find peaks
    # compare the sign of the derivative with the next step
    sign = np.sign(diff)
    next_sign = sign.shift({dim: -1})
    # get the zero crossings to find the high and low peaks
    # checking the gradient of the sign change gives the peak type
    high_peaks = (sign >= 0) & (next_sign < 0)
    low_peaks = (sign <= 0) & (next_sign > 0)
    # return the peaks
    return (high_peaks, low_peaks)
