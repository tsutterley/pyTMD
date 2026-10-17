#!/usr/bin/env python
"""
solid_earth.py
Written by Tyler Sutterley (06/2026)
Prediction routines for solid Earth (body) tides

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
    math.py: Special functions of mathematical physics
    spatial.py: utilities for working with geospatial data

UPDATE HISTORY:
    Updated 06/2026: standardize use of lambda (lmda) to denote longitudes
    Updated 05/2026: use numpy hypot function to calculate magnitudes
        moved ellipsoid and love number parameters to earth module
        refactored body tide function to improve computational times
    Updated 04/2026: use xarray dot product for calculating constituent phases
    Updated 03/2026: use table of body tide love numbers for degrees 4+
    Written 03/2026: split up prediction functions into separate files
"""

from __future__ import annotations

import numpy as np
import xarray as xr
import pyTMD.astro
import pyTMD.constituents
import pyTMD.math
import pyTMD.spatial

# get tide potential tables and ellipsoid parameters
from pyTMD.constituents import _tide_potential_table
from pyTMD.earth import _iers

__all__ = [
    "body_tide",
    "solid_earth_tide",
    "_tide_potential_table",
    "_out_of_phase",
    "_out_of_phase_diurnal",
    "_out_of_phase_semidiurnal",
    "_latitude_dependence",
    "_latitude_dependence_diurnal",
    "_latitude_dependence_semidiurnal",
    "_frequency_dependence",
    "_frequency_dependence_diurnal",
    "_frequency_dependence_long_period",
    "_free_to_mean",
]

# number of days between MJD and the tide epoch (1992-01-01T00:00:00)
_mjd_tide = 48622.0


# PURPOSE: estimate solid Earth tides due to gravitational attraction
# using a simplified approach based on Cartwright and Tayler (1971)
def body_tide(
    t: np.ndarray,
    ds: xr.Dataset,
    deltat: float | np.ndarray = 0.0,
    method: str = "ASTRO5",
    tide_system: str = "tide_free",
    catalog: str = "CTE1973",
    **kwargs,
):
    """
    Compute the solid Earth tides due to the gravitational
    attraction of the Moon and Sun using the approach of
    :cite:t:`Cartwright:1971iz` adjusting the degree-2 Love numbers
    for a near-diurnal frequency dependence :cite:p:`Mathews:1995go`

    Parameters
    ----------
    t: np.ndarray
        Days relative to 1992-01-01T00:00:00
    ds: xarray.Dataset
        Dataset with spatial coordinates
    deltat: float or np.ndarray, default 0.0
        Time correction for converting to Ephemeris Time (days)
    method: str, default 'ASTRO5'
        Method for computing the mean longitudes

            - ``'Cartwright'``
            - ``'Meeus'``
            - ``'ASTRO5'``
            - ``'IERS'``
    tide_system: str, default 'tide_free'
        Output permanent tide system

            - ``'tide_free'``: no permanent direct and indirect tidal potentials
            - ``'mean_tide'``: permanent tidal potentials (direct and indirect)
    catalog: str, default 'CTE1973'
        Name of the tide potential catalog

            - ``'CTE1973'``: :cite:t:`Cartwright:1973em`
            - ``'HW1995'``: :cite:t:`Hartmann:1995jp`
            - ``'T1987'``: :cite:t:`Tamura:1987tp`
            - ``'W1990'``: Woodworth updates to ``'CTE1973'``
    lmax: int, default 6
        Maximum degree of spherical harmonic expansion

        Will be based on the maximum degree available in the catalog
    include_planets: bool, default False
        Include tide potentials from planetary bodies
    h2: float or None, default None
        Degree-2 Love number of vertical displacement
    l2: float or None, default None
        Degree-2 Love (Shida) number of horizontal displacement
    h3: float, default 0.291
        Degree-3 Love number of vertical displacement
    l3: float, default 0.015
        Degree-3 Love (Shida) number of horizontal displacement

    Returns
    -------
    zeta: xr.Dataset
        Solid Earth tide (meters)
    """
    # set default keyword arguments
    # maximum degree of spherical harmonic expansion
    kwargs.setdefault("lmax", 6)
    # include contributions from planets
    kwargs.setdefault("include_planets", False)
    # nominal Love and Shida numbers for degrees 2 and 3
    kwargs.setdefault("h2", None)
    kwargs.setdefault("l2", None)
    kwargs.setdefault("h3", 0.291)
    kwargs.setdefault("l3", 0.015)
    # check if user has provided degree-2 Love numbers
    user_degree_2 = (kwargs["h2"] is not None) and (kwargs["l2"] is not None)
    # validate method and output tide system
    if method.lower() not in ("cartwright", "meeus", "astro5", "iers"):
        raise ValueError("Invalid method for computing mean longitudes")
    if tide_system.lower() not in ("mean_tide", "tide_free"):
        raise ValueError("Invalid permanent tide system")
    if catalog not in _tide_potential_table:
        raise ValueError("Unknown tide potential catalog")

    # convert dates to Modified Julian Days
    MJD = t + _mjd_tide

    # check if tide catalog includes planetary contributions
    if catalog == "HW1995":
        # catalog includes planetary contributions
        # and harmonics up to degree and order 6
        include_planets = True
        lmax = np.minimum(6, kwargs["lmax"])
    elif catalog == "T1987":
        # catalog includes planetary contributions
        # and harmonics up to degree and order 4
        include_planets = True
        lmax = np.minimum(4, kwargs["lmax"])
    else:
        # older catalogs without planetary contributions
        # and harmonics up to degree and order 3
        include_planets = False
        lmax = np.minimum(3, kwargs["lmax"])

    # parse tide potential table for constituents
    table = _tide_potential_table[catalog]
    tide_pot = pyTMD.constituents._parse_tide_potential_table(
        table,
        skiprows=1,
        columns=1,
        include_degree=True,
        include_planets=include_planets,
    )
    # number of constituents in the catalog
    nrows = len(tide_pot)

    # compute principal mean longitudes
    # convert dates into Ephemeris Time
    s, h, p, n, pp = pyTMD.astro.mean_longitudes(MJD + deltat, method=method)
    # initial time conversions
    hour = 24.0 * np.mod(MJD, 1)
    # convert from hours solar time into mean lunar time in degrees
    tau = 15.0 * hour - s + h
    # variable for multiples of 90 degrees (Ray technical note 2017)
    # full expansion of Equilibrium Tide includes some negative cosine
    # terms and some sine terms (Pugh and Woodworth, 2014)
    k = 90.0 + np.zeros_like(MJD)

    # astronomical and planetary mean longitudes
    args = ["tau", "s", "h", "p", "n", "pp", "k"]
    # verify that the catalog includes planetary contributions
    if kwargs["include_planets"] and include_planets:
        # calculate planetary mean longitudes
        # me: Mercury, ve: Venus, ma: Mars, ju: Jupiter, sa: Saturn
        me, ve, ma, ju, sa = pyTMD.astro.planetary_longitudes(MJD + deltat)
        arguments = np.c_[tau, s, h, p, n, pp, k, me, ve, ma, ju, sa]
        args.extend(["me", "ve", "ma", "ju", "sa"])
    else:
        arguments = np.c_[tau, s, h, p, n, pp, k]
    # number of arguments
    nargs = len(args)
    # allocate array for Doodson coefficients
    coef = np.zeros((nargs, nrows))

    # longitudes and colatitudes in radians
    lmda = np.radians(ds.x)
    th = np.radians(90.0 - ds.y)
    # latitude dependence of Love/Shida numbers for degree 2
    dep2 = 1.0 - 1.5 * np.sin(th) ** 2

    # create dictionary of tide potential parameters for each constituent
    CTE = dict(dims=("time", "constituent"), coords={}, data_vars={})
    # time and constituent coordinates
    for coord in ["time", "constituent"]:
        CTE["coords"][coord] = dict(dims=coord)
    CTE["coords"]["time"]["data"] = np.atleast_1d(MJD)
    CTE["coords"]["constituent"]["data"] = np.arange(nrows)
    # constituent parameters
    for field in ["degree", "order", "amplitude", "hl", "ll"]:
        CTE["data_vars"][field] = dict(dims=("constituent"))
    # constituent equilibrium phases
    for field in ["phase"]:
        CTE["data_vars"][field] = dict(dims=("time", "constituent"))
    # spherical harmonic degree and order
    CTE["data_vars"]["degree"]["data"] = np.zeros(nrows)
    CTE["data_vars"]["order"]["data"] = np.zeros(nrows)
    # constituent potential amplitudes (meters)
    CTE["data_vars"]["amplitude"]["data"] = np.zeros(nrows)
    # complex love numbers
    CTE["data_vars"]["hl"]["data"] = np.zeros(nrows, dtype=np.complex128)
    CTE["data_vars"]["ll"]["data"] = np.zeros(nrows, dtype=np.complex128)

    # for each line in the table
    for i, line in enumerate(tide_pot):
        # spherical harmonic degree
        l = line["l"]
        # skip if degree is above the specified expansion limit
        if l > lmax:
            continue
        # spherical harmonic dependence (order)
        TAU = line["tau"]
        # update Doodson coefficients for constituent
        coef[0, i] = TAU
        coef[1, i] = line["s"]
        coef[2, i] = line["h"]
        coef[3, i] = line["p"]
        # convert N for ascending lunar node (from N')
        coef[4, i] = -1.0 * line["n"]
        coef[5, i] = line["pp"]
        # use cosines for (l + tau) even
        # and sines for (l + tau) odd
        coef[6, i] = -1.0 * np.mod(l + TAU, 2)
        # include planetary contributions
        if kwargs["include_planets"] and include_planets:
            # coefficients including planetary terms
            coef[7, i] = line["lme"]
            coef[8, i] = line["lve"]
            coef[9, i] = line["lma"]
            coef[10, i] = line["lju"]
            coef[11, i] = line["lsa"]
        # calculate angular frequency of constituents
        # if not including planets the coefficients will be zero
        # (and not contribute to the angular frequency)
        omega = pyTMD.constituents._frequency(
            coef[:, i], method=method, include_planets=kwargs["include_planets"]
        )
        # skip the permanent tide if using a mean-tide system
        if (omega == 0) and (tide_system.lower() == "mean_tide"):
            continue
        # save degree and order for constituent
        CTE["data_vars"]["degree"]["data"][i] = l
        CTE["data_vars"]["order"]["data"][i] = TAU
        # potential amplitude for constituent
        CTE["data_vars"]["amplitude"]["data"][i] = line["Hs1"]
        # calculate Love numbers for degree
        if (l == 2) and user_degree_2:
            # user-defined Love numbers for all constituents
            hl = np.complex128(kwargs["h2"])
            ll = np.complex128(kwargs["l2"])
        elif (l == 2) and (method.lower() == "iers"):
            # IERS: including both in-phase and out-of-phase components
            # 1) using resonance formula for tides in the diurnal band
            # 2) adjusting some long-period tides for anelastic effects
            hl, _, ll = pyTMD.earth.complex_love_numbers(omega, method=method)
        elif l == 2:
            # use resonance formula for tides in the diurnal band
            hl, _, ll = pyTMD.earth.love_numbers(
                omega, method=method, astype=np.complex128
            )
        else:
            # extract the body tide love numbers for degree
            hb, _, lb = pyTMD.earth.degree_love_numbers(l)
            # use nominal Love numbers for all other degrees
            hl = np.complex128(kwargs.get(f"h{l}", hb))
            ll = np.complex128(kwargs.get(f"l{l}", lb))
        # save Love numbers for constituent
        CTE["data_vars"]["hl"]["data"][i] = hl.copy()
        CTE["data_vars"]["ll"]["data"][i] = ll.copy()
    # calculate constituent phase using equilibrium arguments
    # this calculates over all constituents and time steps
    G = pyTMD.math.normalize_angle(np.dot(arguments, coef))
    CTE["data_vars"]["phase"]["data"] = np.exp(1j * np.radians(G))
    # convert to xarray Dataset from the constituent data dictionary
    CTE = xr.Dataset.from_dict(CTE)

    # allocate for output body tide estimates (meters)
    # latitudinal, longitudinal and radial components
    zeta = xr.Dataset()
    # initialize output body tides
    for key in ["R", "N", "E"]:
        zeta[key] = xr.zeros_like(th * CTE["time"])

    # order of operations for each spherical harmonic degree and order:
    # 1) extract parameters for constituents with degree l and order m
    # 2) extract or calculate the complex Love numbers for degree and order
    # 3) calculate spherical harmonic functions and their derivatives
    # 4) rotate the spherical harmonics by the constituent phase and
    # scale by the Cartwright-Tayler-Edden (CTE) amplitude
    # 5) convert the potentials to displacements and sum over constituents
    # 6) add to totals (latitudinal, longitudinal and radial components)
    # note: looping over constituents is slow due to repeated rotations and
    #       dimensional broadcasting needed for the spherical harmonics
    # here we split the calculations and do a summation over constituents
    for l in range(2, lmax + 1):
        for m in range(l + 1):
            # skip if degree or order are not included in the catalog
            if l not in CTE.degree or m not in CTE.order:
                continue
            # extract parameters for constituents with degree l and order m
            params = CTE.where((CTE.degree == l) & (CTE.order == m), drop=True)
            # calculate or extract complex Love numbers for degree and order
            # (includes frequency dependence for degree 2)
            if (l == 2) and (method.lower() == "iers"):
                # include (complex) latitudinal dependence for degree 2
                hl = params["hl"] - (0.615e-3 + 0.122e-4j) * dep2
                ll = params["ll"] + (0.19334e-3 - 0.3819e-5j) * dep2
            elif l == 2 and not user_degree_2:
                # include (nominal) latitudinal dependence for degree 2
                hl = params["hl"] - 0.0006 * dep2
                ll = params["ll"] + 0.0002 * dep2
            else:
                # nominal Love numbers
                hl, ll = params["hl"], params["ll"]
            # spherical harmonic functions and derivatives
            # will need to be rotated by constituent phase
            Ylm, dYlm = pyTMD.math.sph_harm(l, th, lmda, m=m)
            # rotate spherical harmonics and scale by amplitude
            S = params["amplitude"] * Ylm * params["phase"]
            dS = params["amplitude"] * dYlm * params["phase"]
            # convert potentials for constituent and sum over constituents
            # add to totals (latitudinal, longitudinal and radial components)
            zeta["N"] += (ll.real * dS.real - ll.imag * dS.imag).sum(
                dim="constituent"
            )
            zeta["E"] -= m * (ll.real * S.imag - ll.imag * S.real).sum(
                dim="constituent"
            )
            zeta["R"] += (hl.real * S.real - hl.imag * S.imag).sum(
                dim="constituent"
            )

    # add units attributes to output dataset
    for var in zeta.data_vars:
        zeta[var].attrs["units"] = "meters"
    # return the body tides
    return zeta


# PURPOSE: estimate solid Earth tides due to gravitational attraction
def solid_earth_tide(
    t: np.ndarray,
    XYZ: xr.Dataset,
    SXYZ: xr.Dataset,
    LXYZ: xr.Dataset,
    deltat: float = 0.0,
    a_axis: float = _iers.a_axis,
    tide_system: str = "tide_free",
    **kwargs,
):
    """
    Compute the solid Earth tides in Cartesian coordinates
    due to the gravitational attraction of the Moon and Sun
    :cite:p:`Mathews:1991kv,Mathews:1997js,Ries:1992ip,Wahr:1981ea`

    Parameters
    ----------
    t: np.ndarray
        Days relative to 1992-01-01T00:00:00
    XYZ: xr.Dataset
        Dataset with cartesian coordinates
    SXYZ: xr.Dataset
        Dataset with Earth-centered Earth-fixed coordinates of the Sun
    LXYZ: xr.Dataset
        Dataset with Earth-centered Earth-fixed coordinates of the Moon
    deltat: float or np.ndarray, default 0.0
        Time correction for converting to Ephemeris Time (days)
    a_axis: float, default 6378136.3
        Semi-major axis of the Earth (meters)
    tide_system: str, default 'tide_free'
        Permanent tide system for the output solid Earth tide

            - ``'tide_free'``: no permanent direct and indirect tidal potentials
            - ``'mean_tide'``: permanent tidal potentials (direct and indirect)
    lmax: int, default 3
        Maximum degree of spherical harmonic expansion
    h2: float, default 0.6078
        Degree-2 Love number of vertical displacement
    l2: float, default 0.0847
        Degree-2 Love (Shida) number of horizontal displacement
    h3: float, default 0.292
        Degree-3 Love number of vertical displacement
    l3: float, default 0.015
        Degree-3 Love (Shida) number of horizontal displacement
    mass_ratio_solar: float, default 332946.0482
        Mass ratio between the Earth and the Sun
    mass_ratio_lunar: float, default 0.0123000371
        Mass ratio between the Earth and the Moon

    Returns
    -------
    dxt: xr.Dataset
        Solid Earth tide displacements (meters)
    """
    # set default keyword arguments
    # maximum degree of spherical harmonic expansion
    kwargs.setdefault("lmax", 3)
    # nominal Love and Shida numbers for degrees 2 and 3
    kwargs.setdefault("h2", 0.6078)
    kwargs.setdefault("l2", 0.0847)
    kwargs.setdefault("h3", 0.292)
    kwargs.setdefault("l3", 0.015)
    # mass ratios between earth and Sun/Moon
    kwargs.setdefault("mass_ratio_solar", 332946.0482)
    kwargs.setdefault("mass_ratio_lunar", 0.0123000371)
    # validate output tide system
    if tide_system.lower() not in ("mean_tide", "tide_free"):
        raise ValueError("Invalid permanent tide system")
    # convert time to Modified Julian Days (MJD)
    MJD = t + _mjd_tide
    # radius of the point on the Earth's surface
    radius = pyTMD.math.radius(XYZ["X"], XYZ["Y"], XYZ["Z"])
    # sine of geocentric latitude
    sinphi = XYZ["Z"] / radius
    # distance between the Earth and the Sun/Moon
    solar_radius = pyTMD.math.radius(SXYZ["X"], SXYZ["Y"], SXYZ["Z"])
    lunar_radius = pyTMD.math.radius(LXYZ["X"], LXYZ["Y"], LXYZ["Z"])
    # cosine of angles between vectors of the point and the Sun/Moon
    solar_scalar = pyTMD.math.scalar_product(
        XYZ["X"], XYZ["Y"], XYZ["Z"], SXYZ["X"], SXYZ["Y"], SXYZ["Z"]
    ) / (radius * solar_radius)
    lunar_scalar = pyTMD.math.scalar_product(
        XYZ["X"], XYZ["Y"], XYZ["Z"], LXYZ["X"], LXYZ["Y"], LXYZ["Z"]
    ) / (radius * lunar_radius)
    # unit vectors for dimensions
    unit_vector = XYZ / radius
    solar_unit_vector = SXYZ / solar_radius
    lunar_unit_vector = LXYZ / lunar_radius
    # factors for Sun and Moon using IAU estimates of mass ratios
    K_solar = kwargs["mass_ratio_solar"] * a_axis**3 / solar_radius**2
    K_lunar = kwargs["mass_ratio_lunar"] * a_axis**3 / lunar_radius**2
    # factors for degree 2
    F2_solar = K_solar * (a_axis / solar_radius)
    F2_lunar = K_lunar * (a_axis / lunar_radius)
    # allocate for output displacements
    dx_solar = xr.Dataset()
    dx_lunar = xr.Dataset()
    for var in ["X", "Y", "Z"]:
        dx_solar[var] = xr.zeros_like(solar_scalar)
        dx_lunar[var] = xr.zeros_like(lunar_scalar)
    # compute total displacement (Mathews et al. 1997)
    # from the tide-generating potentials
    # for each spherical harmonic degree
    for l in range(2, kwargs["lmax"] + 1):
        # extract the body tide love numbers for degree
        hb, kb, lb = pyTMD.earth.degree_love_numbers(l)
        # get the degree-dependent Love numbers
        hl = kwargs.get(f"h{l}", hb)
        ll = kwargs.get(f"l{l}", lb)
        # include latitudinal dependence for degree 2
        # from equations 5 and 6 of Mathews et al., (1997)
        if l == 2:
            hl -= 0.0006 * (1.5 * sinphi**2 - 0.5)
            ll += 0.0002 * (1.5 * sinphi**2 - 0.5)
        # update gravitational parameters for degree
        K_solar *= a_axis / solar_radius
        K_lunar *= a_axis / lunar_radius
        # legendre polynomial for degree
        Pl_solar, dPl_solar = pyTMD.math.legendre(l, solar_scalar)
        Pl_lunar, dPl_lunar = pyTMD.math.legendre(l, lunar_scalar)
        # divide differential by u
        # ignore divide by zero and invalid value warnings
        with np.errstate(divide="ignore", invalid="ignore"):
            dPl_solar /= np.sqrt(1 - solar_scalar**2)
            dPl_lunar /= np.sqrt(1 - lunar_scalar**2)
        # add solar and lunar terms for degree
        dx_solar += K_solar * (
            hl * Pl_solar * unit_vector
            + ll * dPl_solar * solar_scalar * unit_vector
            - ll * dPl_solar * solar_unit_vector
        )
        dx_lunar += K_lunar * (
            hl * Pl_lunar * unit_vector
            + ll * dPl_lunar * lunar_scalar * unit_vector
            - ll * dPl_lunar * lunar_unit_vector
        )
    # sum solar and lunar components
    dxt = dx_solar + dx_lunar
    # corrections for out-of-phase portions of the Love and Shida numbers
    dxt += _out_of_phase(XYZ, SXYZ, LXYZ, F2_solar, F2_lunar)
    # corrections for the latitudinal dependence (diurnal and semi-diurnal)
    dxt += _latitude_dependence(XYZ, SXYZ, LXYZ, F2_solar, F2_lunar)
    # corrections for the frequency dependence (diurnal and long-period)
    dxt += _frequency_dependence(XYZ, MJD, deltat=deltat)
    # convert the permanent tide system if specified
    if tide_system.lower() == "mean_tide":
        # compute new h2 and l2 (Mathews et al., 1997)
        h2 = kwargs["h2"] - 0.0006 * (1.5 * sinphi**2 - 0.5)
        l2 = kwargs["l2"] + 0.0002 * (1.5 * sinphi**2 - 0.5)
        dxt += _free_to_mean(XYZ, h2, l2)
    # add units attributes to output dataset
    for var in dxt.data_vars:
        dxt[var].attrs["units"] = "meters"
    # return the solid earth tide
    return dxt


def _out_of_phase(
    XYZ: xr.Dataset,
    SXYZ: xr.Dataset,
    LXYZ: xr.Dataset,
    F2_solar: np.ndarray,
    F2_lunar: np.ndarray,
):
    """
    Wrapper function to compute the out-of-phase corrections induced
    by mantle anelasticity :cite:p:`Petit:2010tp`

    Parameters
    ----------
    XYZ: xr.Dataset
        Dataset with cartesian coordinates
    SXYZ: xr.Dataset
        Dataset with Earth-centered Earth-fixed coordinates of the Sun
    LXYZ: xr.Dataset
        Dataset with Earth-centered Earth-fixed coordinates of the Moon
    F2_solar: np.ndarray
        Factors for the Sun
    F2_lunar: np.ndarray
        Factors for the Moon

    Returns
    -------
    D: xr.Dataset
        Solid Earth tide corrections
    """
    # compute diurnal and semi-diurnal corrections separately
    # for both the Sun and Moon
    D = _out_of_phase_diurnal(XYZ, SXYZ, F2_solar)
    D += _out_of_phase_diurnal(XYZ, LXYZ, F2_lunar)
    D += _out_of_phase_semidiurnal(XYZ, SXYZ, F2_solar)
    D += _out_of_phase_semidiurnal(XYZ, LXYZ, F2_lunar)
    # return the out-of-phase corrections
    return D


def _out_of_phase_diurnal(
    XYZ: xr.Dataset,
    LSXYZ: xr.Dataset,
    F2: np.ndarray,
    dh2: float = -0.0025,
    dl2: float = -0.0007,
):
    """
    Computes the out-of-phase corrections induced by mantle
    anelasticity in the diurnal band :cite:p:`Petit:2010tp`

    Parameters
    ----------
    XYZ: xr.Dataset
        Dataset with cartesian coordinates
    LSXYZ: xr.Dataset
        Dataset with Earth-centered Earth-fixed coordinates of the Sun or Moon
    F2: np.ndarray
        Factors for the Sun or Moon
    dh2: float, default -0.0025
        Love number correction for the diurnal band
    dl2: float, default -0.0007
        Shida number correction for the diurnal band

    Returns
    -------
    D: xr.Dataset
        Solid Earth tide corrections
    """
    # compute the normalized position vector of coordinates
    radius = pyTMD.math.radius(XYZ["X"], XYZ["Y"], XYZ["Z"])
    # sine and cosine of (geocentric) latitude
    sinphi = XYZ["Z"] / radius
    cosphi = np.hypot(XYZ["X"], XYZ["Y"]) / radius
    # double angle formulas of cosine/sine latitude
    sin2phi = 2.0 * sinphi * cosphi
    cos2phi = cosphi**2 - sinphi**2
    # sine and cosine of longitude
    sinla = XYZ["Y"] / cosphi / radius
    cosla = XYZ["X"] / cosphi / radius
    # compute the normalized position vector of the Sun/Moon
    lunisolar_radius = pyTMD.math.radius(LSXYZ["X"], LSXYZ["Y"], LSXYZ["Z"])
    # sine and cosine of Solar/Lunar declinations
    lunisolar_sinphi = LSXYZ["Z"] / lunisolar_radius
    lunisolar_cosphi = np.hypot(LSXYZ["X"], LSXYZ["Y"]) / lunisolar_radius
    # double angle formulas of sine Solar/Lunar declinations
    lunisolar_sin2phi = 2.0 * lunisolar_cosphi * lunisolar_sinphi
    # sine and cosine of Solar/Lunar hour angles
    lunisolar_sinla = LSXYZ["Y"] / lunisolar_cosphi / lunisolar_radius
    lunisolar_cosla = LSXYZ["X"] / lunisolar_cosphi / lunisolar_radius
    # calculate offsets
    # equation 19 from Mathews et al. (1997)
    DR = (-0.75 * dh2 * F2 * sin2phi * lunisolar_sin2phi) * (
        sinla * lunisolar_cosla - cosla * lunisolar_sinla
    )
    # equation 20 from Mathews et al. (1997)
    DN = (-1.5 * dl2 * F2 * cos2phi * lunisolar_sin2phi) * (
        sinla * lunisolar_cosla - cosla * lunisolar_sinla
    )
    DE = (-1.5 * dl2 * F2 * sinphi * lunisolar_sin2phi) * (
        cosla * lunisolar_cosla + sinla * lunisolar_sinla
    )
    # output corrections
    D = xr.Dataset()
    # compute corrections in cartesian coordinates
    D["X"] = DR * cosla * cosphi - DE * sinla - DN * cosla * sinphi
    D["Y"] = DR * sinla * cosphi + DE * cosla - DN * sinla * sinphi
    D["Z"] = DR * sinphi + DN * cosphi
    # return the corrections
    return D


def _out_of_phase_semidiurnal(
    XYZ: xr.Dataset,
    LSXYZ: xr.Dataset,
    F2: np.ndarray,
    dh2: float = -0.0022,
    dl2: float = -0.0007,
):
    """
    Computes the out-of-phase corrections induced by mantle
    anelasticity in the semi-diurnal band :cite:p:`Petit:2010tp`

    Parameters
    ----------
    XYZ: xr.Dataset
        Dataset with cartesian coordinates
    LSXYZ: xr.Dataset
        Dataset with Earth-centered Earth-fixed coordinates of the Sun or Moon
    F2: np.ndarray
        Factors for the Sun or Moon
    dh2: float, default -0.0022
        Love number correction for the semi-diurnal band
    dl2: float, default -0.0007
        Shida number correction for the semi-diurnal band

    Returns
    -------
    D: xr.Dataset
        Solid Earth tide corrections
    """
    # compute the normalized position vector of coordinates
    radius = pyTMD.math.radius(XYZ["X"], XYZ["Y"], XYZ["Z"])
    # sine and cosine of (geocentric) latitude
    sinphi = XYZ["Z"] / radius
    cosphi = np.hypot(XYZ["X"], XYZ["Y"]) / radius
    # sine and cosine of longitude
    sinla = XYZ["Y"] / cosphi / radius
    cosla = XYZ["X"] / cosphi / radius
    # double angle formulas of cosine/sine longitude
    cos2la = cosla**2 - sinla**2
    sin2la = 2.0 * cosla * sinla
    # compute the normalized position vector of the Sun/Moon
    lunisolar_radius = pyTMD.math.radius(LSXYZ["X"], LSXYZ["Y"], LSXYZ["Z"])
    # cosine of Solar/Lunar declinations
    lunisolar_cosphi = np.hypot(LSXYZ["X"], LSXYZ["Y"]) / lunisolar_radius
    # sine and cosine of Solar/Lunar hour angles
    lunisolar_sinla = LSXYZ["Y"] / lunisolar_cosphi / lunisolar_radius
    lunisolar_cosla = LSXYZ["X"] / lunisolar_cosphi / lunisolar_radius
    # double angle formulas of cosine/sine Solar/Lunar hour angles
    lunisolar_cos2la = lunisolar_cosla**2 - lunisolar_sinla**2
    lunisolar_sin2la = 2.0 * lunisolar_cosla * lunisolar_sinla
    # calculate offsets
    # equation 21 from Mathews et al. (1997)
    DR = (-0.75 * dh2 * F2 * cosphi**2 * lunisolar_cosphi**2) * (
        sin2la * lunisolar_cos2la - cos2la * lunisolar_sin2la
    )
    # equation 22 from Mathews et al. (1997)
    DN = (1.5 * dl2 * F2 * cosphi * sinphi * lunisolar_cosphi**2) * (
        sin2la * lunisolar_cos2la - cos2la * lunisolar_sin2la
    )
    DE = (-1.5 * dl2 * F2 * cosphi * lunisolar_cosphi**2) * (
        cos2la * lunisolar_cos2la + sin2la * lunisolar_sin2la
    )
    # output corrections
    D = xr.Dataset()
    # compute corrections in cartesian coordinates
    D["X"] = DR * cosla * cosphi - DE * sinla - DN * cosla * sinphi
    D["Y"] = DR * sinla * cosphi + DE * cosla - DN * sinla * sinphi
    D["Z"] = DR * sinphi + DN * cosphi
    # return the corrections
    return D


def _latitude_dependence(
    XYZ: xr.Dataset,
    SXYZ: xr.Dataset,
    LXYZ: xr.Dataset,
    F2_solar: np.ndarray,
    F2_lunar: np.ndarray,
):
    r"""
    Wrapper function to compute the latitudinal dependent corrections
    given by L\ :sup:`1` for both the diurnal and semi-diurnal bands
    :cite:p:`Petit:2010tp`

    Parameters
    ----------
    XYZ: xr.Dataset
        Dataset with cartesian coordinates
    SXYZ: xr.Dataset
        Dataset with Earth-centered Earth-fixed coordinates of the Sun
    LXYZ: xr.Dataset
        Dataset with Earth-centered Earth-fixed coordinates of the Moon
    F2_solar: np.ndarray
        Factors for the Sun
    F2_lunar: np.ndarray
        Factors for the Moon

    Returns
    -------
    D: xr.Dataset
        Solid Earth tide corrections
    """
    # compute diurnal and semi-diurnal corrections separately
    # for both the Sun and Moon
    D = _latitude_dependence_diurnal(XYZ, SXYZ, F2_solar)
    D += _latitude_dependence_diurnal(XYZ, LXYZ, F2_lunar)
    D += _latitude_dependence_semidiurnal(XYZ, SXYZ, F2_solar)
    D += _latitude_dependence_semidiurnal(XYZ, LXYZ, F2_lunar)
    # return the latitudinal dependent corrections
    return D


def _latitude_dependence_diurnal(
    XYZ: xr.Dataset,
    LSXYZ: xr.Dataset,
    F2: np.ndarray,
    L1: float = 0.0012,
):
    r"""
    Computes the corrections induced by the latitudinal
    dependence of the diurnal band :cite:p:`Petit:2010tp`

    Parameters
    ----------
    XYZ: xr.Dataset
        Dataset with cartesian coordinates
    LSXYZ: xr.Dataset
        Dataset with Earth-centered Earth-fixed coordinates of the Sun or Moon
    F2: np.ndarray
        Factors for the Sun or Moon
    L1: float, default 0.0012
        Love/Shida number correction for the diurnal band

    Returns
    -------
    D: xr.Dataset
        Solid Earth tide corrections
    """
    # compute the normalized position vector of coordinates
    radius = pyTMD.math.radius(XYZ["X"], XYZ["Y"], XYZ["Z"])
    # sine and cosine of (geocentric) latitude
    sinphi = XYZ["Z"] / radius
    cosphi = np.hypot(XYZ["X"], XYZ["Y"]) / radius
    # double angle formulas of cosine latitude
    cos2phi = cosphi**2 - sinphi**2
    # sine and cosine of longitude
    sinla = XYZ["Y"] / cosphi / radius
    cosla = XYZ["X"] / cosphi / radius
    # compute the normalized position vector of the Sun/Moon
    lunisolar_radius = pyTMD.math.radius(LSXYZ["X"], LSXYZ["Y"], LSXYZ["Z"])
    # sine and cosine of Solar/Lunar declinations
    lunisolar_sinphi = LSXYZ["Z"] / lunisolar_radius
    lunisolar_cosphi = np.hypot(LSXYZ["X"], LSXYZ["Y"]) / lunisolar_radius
    # double angle formulas of sin Solar/Lunar declinations
    lunisolar_sin2phi = 2.0 * lunisolar_sinphi * lunisolar_cosphi
    # sine and cosine of Solar/Lunar hour angles
    lunisolar_sinla = LSXYZ["Y"] / lunisolar_cosphi / lunisolar_radius
    lunisolar_cosla = LSXYZ["X"] / lunisolar_cosphi / lunisolar_radius
    # calculate offsets for the diurnal band
    # equation 25 from Mathews et al. (1997)
    DN = (-1.5 * L1 * F2 * sinphi**2 * lunisolar_sin2phi) * (
        cosla * lunisolar_cosla + sinla * lunisolar_sinla
    )
    DE = (1.5 * L1 * F2 * sinphi * cos2phi * lunisolar_sin2phi) * (
        sinla * lunisolar_cosla - cosla * lunisolar_sinla
    )
    # output corrections
    D = xr.Dataset()
    # compute corrections in cartesian coordinates
    D["X"] = -DE * sinla - DN * cosla * sinphi
    D["Y"] = DE * cosla - DN * sinla * sinphi
    D["Z"] = DN * cosphi
    # return the corrections
    return D


def _latitude_dependence_semidiurnal(
    XYZ: xr.Dataset,
    LSXYZ: xr.Dataset,
    F2: np.ndarray,
    L1: float = 0.0024,
):
    r"""
    Computes the corrections induced by the latitudinal
    dependence of the semi-diurnal band :cite:p:`Petit:2010tp`

    Parameters
    ----------
    XYZ: xr.Dataset
        Dataset with cartesian coordinates
    LSXYZ: xr.Dataset
        Dataset with Earth-centered Earth-fixed coordinates of the Sun or Moon
    F2: np.ndarray
        Factors for the Sun or Moon
    L1: float, default 0.0024
        Love/Shida number correction for the semi-diurnal band

    Returns
    -------
    D: xr.Dataset
        Solid Earth tide corrections
    """
    # compute the normalized position vector of coordinates
    radius = pyTMD.math.radius(XYZ["X"], XYZ["Y"], XYZ["Z"])
    # sine and cosine of (geocentric) latitude
    sinphi = XYZ["Z"] / radius
    cosphi = np.hypot(XYZ["X"], XYZ["Y"]) / radius
    # sine and cosine of longitude
    sinla = XYZ["Y"] / cosphi / radius
    cosla = XYZ["X"] / cosphi / radius
    # double angle formulas of cos/sin longitude
    cos2la = cosla**2 - sinla**2
    sin2la = 2.0 * cosla * sinla
    # compute the normalized position vector of the Sun/Moon
    lunisolar_radius = pyTMD.math.radius(LSXYZ["X"], LSXYZ["Y"], LSXYZ["Z"])
    # cosine of Solar/Lunar declinations
    lunisolar_cosphi = np.hypot(LSXYZ["X"], LSXYZ["Y"]) / lunisolar_radius
    # sine and cosine of Solar/Lunar hour angles
    lunisolar_sinla = LSXYZ["Y"] / lunisolar_cosphi / lunisolar_radius
    lunisolar_cosla = LSXYZ["X"] / lunisolar_cosphi / lunisolar_radius
    # double angle formulas of cosine/sine Solar/Lunar hour angles
    lunisolar_cos2la = lunisolar_cosla**2 - lunisolar_sinla**2
    lunisolar_sin2la = 2.0 * lunisolar_cosla * lunisolar_sinla
    # calculate offsets for the semi-diurnal band
    # equation 26 from Mathews et al. (1997)
    DN = (-1.5 * L1 * F2 * sinphi * cosphi * lunisolar_cosphi**2) * (
        cos2la * lunisolar_cos2la + sin2la * lunisolar_sin2la
    )
    DE = (-1.5 * L1 * F2 * sinphi**2 * cosphi * lunisolar_cosphi**2) * (
        sin2la * lunisolar_cos2la - cos2la * lunisolar_sin2la
    )
    # output corrections
    D = xr.Dataset()
    # compute corrections in cartesian coordinates
    D["X"] = -DE * sinla - DN * cosla * sinphi
    D["Y"] = DE * cosla - DN * sinla * sinphi
    D["Z"] = DN * cosphi
    # return the corrections
    return D


def _frequency_dependence(
    XYZ: xr.Dataset,
    MJD: np.ndarray,
    deltat: float | np.ndarray = 0.0,
):
    """
    Wrapper function to compute the frequency dependent in-phase and
    out-of-phase corrections :cite:p:`Petit:2010tp`

    Parameters
    ----------
    XYZ: xr.Dataset
        Dataset with cartesian coordinates
    MJD: np.ndarray
        Modified Julian Day (MJD)
    deltat: float or np.ndarray, default 0.0
        Time correction for converting to Ephemeris Time (days)

    Returns
    -------
    D: xr.Dataset
        Solid Earth tide corrections
    """
    # Love/Shida number corrections (diurnal and long-period)
    # compute diurnal and long-period corrections separately
    D = _frequency_dependence_diurnal(XYZ, MJD, deltat=deltat)
    D += _frequency_dependence_long_period(XYZ, MJD, deltat=deltat)
    # return the frequency dependent corrections
    return D


def _frequency_dependence_diurnal(
    XYZ: xr.Dataset,
    MJD: np.ndarray,
    deltat: float | np.ndarray = 0.0,
):
    """
    Computes the frequency dependent in-phase and out-of-phase corrections
    of the diurnal band :cite:p:`Petit:2010tp`

    Parameters
    ----------
    XYZ: xr.Dataset
        Dataset with cartesian coordinates
    MJD: np.ndarray
        Modified Julian Day (MJD)
    deltat: float or np.ndarray, default 0.0
        Time correction for converting to Ephemeris Time (days)

    Returns
    -------
    D: xr.Dataset
        Solid Earth tide corrections
    """
    # Corrections for Frequency Dependence of Diurnal Tides
    # reduced version of table 7.3a from IERS conventions
    columns = [
        "tau",
        "s",
        "h",
        "p",
        "np",
        "ps",
        "dR_ip",
        "dR_op",
        "dT_ip",
        "dT_op",
    ]
    table = xr.DataArray(
        np.array(
            [
                [1, -3, 0, 2, 0, 0, -0.01, 0.0, 0.0, 0.0],
                [1, -3, 2, 0, 0, 0, -0.01, 0.0, 0.0, 0.0],
                [1, -2, 0, 1, -1, 0, -0.02, 0.0, 0.0, 0.0],
                [1, -2, 0, 1, 0, 0, -0.08, 0.0, -0.01, 0.01],
                [1, -2, 2, -1, 0, 0, -0.02, 0.0, 0.0, 0.0],
                [1, -1, 0, 0, -1, 0, -0.10, 0.0, 0.0, 0.0],
                [1, -1, 0, 0, 0, 0, -0.51, 0.0, -0.02, 0.03],
                [1, -1, 2, 0, 0, 0, 0.01, 0.0, 0.0, 0.0],
                [1, 0, -2, 1, 0, 0, 0.01, 0.0, 0.0, 0.0],
                [1, 0, 0, -1, 0, 0, 0.02, 0.0, 0.0, 0.0],
                [1, 0, 0, 1, 0, 0, 0.06, 0.0, 0.0, 0.0],
                [1, 0, 0, 1, 1, 0, 0.01, 0.0, 0.0, 0.0],
                [1, 0, 2, -1, 0, 0, 0.01, 0.0, 0.0, 0.0],
                [1, 1, -3, 0, 0, 1, -0.06, 0.0, 0.0, 0.0],
                [1, 1, -2, 0, -1, 0, 0.01, 0.0, 0.0, 0.0],
                [1, 1, -2, 0, 0, 0, -1.23, -0.07, 0.06, 0.01],
                [1, 1, -1, 0, 0, -1, 0.02, 0.0, 0.0, 0.0],
                [1, 1, -1, 0, 0, 1, 0.04, 0.0, 0.0, 0.0],
                [1, 1, 0, 0, -1, 0, -0.22, 0.01, 0.01, 0.0],
                [1, 1, 0, 0, 0, 0, 12.00, -0.80, -0.67, -0.03],
                [1, 1, 0, 0, 1, 0, 1.73, -0.12, -0.10, 0.0],
                [1, 1, 0, 0, 2, 0, -0.04, 0.0, 0.0, 0.0],
                [1, 1, 1, 0, 0, -1, -0.50, -0.01, 0.03, 0.0],
                [1, 1, 1, 0, 0, 1, 0.01, 0.0, 0.0, 0.0],
                [1, 0, 1, 0, 1, -1, -0.01, 0.0, 0.0, 0.0],
                [1, 1, 2, -2, 0, 0, -0.01, 0.0, 0.0, 0.0],
                [1, 1, 2, 0, 0, 0, -0.11, 0.01, 0.01, 0.0],
                [1, 2, -2, 1, 0, 0, -0.01, 0.0, 0.0, 0.0],
                [1, 2, 0, -1, 0, 0, -0.02, 0.0, 0.0, 0.0],
                [1, 3, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0],
                [1, 3, 0, 0, 1, 0, 0.0, 0.0, 0.0, 0.0],
            ]
        ),
        dims=["constituent", "argument"],
        coords=dict(argument=columns),
    )
    coef = table.to_dataset(dim="argument")
    # get phase angles (Doodson arguments)
    TAU, S, H, P, ZNS, PS = pyTMD.astro.doodson_arguments(MJD + deltat)
    # dataset of arguments
    arguments = xr.Dataset(
        data_vars=dict(
            tau=(["time"], TAU),
            s=(["time"], S),
            h=(["time"], H),
            p=(["time"], P),
            np=(["time"], ZNS),
            ps=(["time"], PS),
        ),
        coords=dict(time=np.atleast_1d(MJD)),
    )
    # compute the normalized position vector of coordinates
    radius = pyTMD.math.radius(XYZ["X"], XYZ["Y"], XYZ["Z"])
    # sine and cosine of (geocentric) latitude
    sinphi = XYZ["Z"] / radius
    cosphi = np.hypot(XYZ["X"], XYZ["Y"]) / radius
    sin2phi = 2.0 * sinphi * cosphi
    cos2phi = cosphi**2 - sinphi**2
    # sine and cosine of longitude
    sinla = XYZ["Y"] / cosphi / radius
    cosla = XYZ["X"] / cosphi / radius
    # compute phase angle of tide potential (Greenwich)
    thetaf = (
        arguments.tau * coef["tau"]
        + arguments.s * coef["s"]
        + arguments.h * coef["h"]
        + arguments.p * coef["p"]
        + arguments.np * coef["np"]
        + arguments.ps * coef["ps"]
    )
    # calculate sine and cosine of phase (local hour angle)
    sphase = np.sin(thetaf) * cosla + np.cos(thetaf) * sinla
    cphase = np.cos(thetaf) * cosla - np.sin(thetaf) * sinla
    # calculate offsets in local coordinates
    # equations 27 and 28 from Mathews et al. (1997)
    DR = sin2phi * (coef["dR_ip"] * sphase + coef["dR_op"] * cphase)
    DN = cos2phi * (coef["dT_ip"] * sphase + coef["dT_op"] * cphase)
    DE = sinphi * (coef["dT_ip"] * cphase - coef["dT_op"] * sphase)
    # compute corrections (Mathews et al. 1997)
    # rotate to cartesian coordinates
    DX = (DR * cosla * cosphi - DE * sinla - DN * cosla * sinphi).sum(
        dim="constituent", skipna=False
    )
    DY = (DR * sinla * cosphi + DE * cosla - DN * sinla * sinphi).sum(
        dim="constituent", skipna=False
    )
    DZ = (DR * sinphi + DN * cosphi).sum(dim="constituent", skipna=False)
    # convert from millimeters to meters
    D = xr.Dataset()
    D["X"] = 1e-3 * DX
    D["Y"] = 1e-3 * DY
    D["Z"] = 1e-3 * DZ
    # return the corrections
    return D


def _frequency_dependence_long_period(
    XYZ: xr.Dataset,
    MJD: np.ndarray,
    deltat: float | np.ndarray = 0.0,
):
    """
    Computes the frequency dependent in-phase and out-of-phase corrections
    induced by mantle anelasticity in the long-period band
    :cite:p:`Petit:2010tp`

    Parameters
    ----------
    XYZ: xr.Dataset
        Dataset with cartesian coordinates
    MJD: np.ndarray
        Modified Julian Day (MJD)
    deltat: float or np.ndarray, default 0.0
        Time correction for converting to Ephemeris Time (days)

    Returns
    -------
    D: xr.Dataset
        Solid Earth tide corrections
    """
    # Corrections for Frequency Dependence of Long-Period Tides
    # reduced version of table 7.3b from IERS conventions
    columns = [
        "tau",
        "s",
        "h",
        "p",
        "np",
        "ps",
        "dR_ip",
        "dR_op",
        "dT_ip",
        "dT_op",
    ]
    table = xr.DataArray(
        np.array(
            [
                [0, 0, 0, 0, 1, 0, 0.47, 0.23, 0.16, 0.07],
                [0, 0, 2, 0, 0, 0, -0.20, -0.12, -0.11, -0.05],
                [0, 1, 0, -1, 0, 0, -0.11, -0.08, -0.09, -0.04],
                [0, 2, 0, 0, 0, 0, -0.13, -0.11, -0.15, -0.07],
                [0, 2, 0, 0, 1, 0, -0.05, -0.05, -0.06, -0.03],
            ]
        ),
        dims=["constituent", "argument"],
        coords=dict(argument=columns),
    )
    coef = table.to_dataset(dim="argument")
    # get phase angles (Doodson arguments)
    TAU, S, H, P, ZNS, PS = pyTMD.astro.doodson_arguments(MJD + deltat)
    # dataset of arguments
    arguments = xr.Dataset(
        data_vars=dict(
            tau=(["time"], TAU),
            s=(["time"], S),
            h=(["time"], H),
            p=(["time"], P),
            np=(["time"], ZNS),
            ps=(["time"], PS),
        ),
        coords=dict(time=np.atleast_1d(MJD)),
    )
    # compute the normalized position vector of coordinates
    radius = pyTMD.math.radius(XYZ["X"], XYZ["Y"], XYZ["Z"])
    # sine and cosine of (geocentric) latitude
    sinphi = XYZ["Z"] / radius
    cosphi = np.hypot(XYZ["X"], XYZ["Y"]) / radius
    sin2phi = 2.0 * cosphi * sinphi
    # sine and cosine of longitude
    sinla = XYZ["Y"] / cosphi / radius
    cosla = XYZ["X"] / cosphi / radius
    # compute phase angle of tide potential (Greenwich)
    thetaf = (
        arguments.tau * coef["tau"]
        + arguments.s * coef["s"]
        + arguments.h * coef["h"]
        + arguments.p * coef["p"]
        + arguments.np * coef["np"]
        + arguments.ps * coef["ps"]
    )
    # calculate sine and cosine of phase (local hour angle)
    # note that zonal harmonics have no longitudinal dependence
    sphase = np.sin(thetaf)
    cphase = np.cos(thetaf)
    # calculate offsets in local coordinates
    # equations 31 and 32 from Mathews et al. (1997)
    DR = (1.5 * sinphi**2 - 0.5) * (
        coef["dT_ip"] * sphase + coef["dR_ip"] * cphase
    )
    DN = sin2phi * (coef["dT_op"] * sphase + coef["dR_op"] * cphase)
    # compute corrections (Mathews et al. 1997)
    # rotate to cartesian coordinates
    DX = (DR * cosla * cosphi - DN * cosla * sinphi).sum(
        dim="constituent", skipna=False
    )
    DY = (DR * sinla * cosphi - DN * sinla * sinphi).sum(
        dim="constituent", skipna=False
    )
    DZ = (DR * sinphi + DN * cosphi).sum(dim="constituent", skipna=False)
    # convert from millimeters to meters
    D = xr.Dataset()
    D["X"] = 1e-3 * DX
    D["Y"] = 1e-3 * DY
    D["Z"] = 1e-3 * DZ
    # return the corrections
    return D


def _free_to_mean(
    XYZ: xr.Dataset,
    h2: float | np.ndarray,
    l2: float | np.ndarray,
    H0: float = -0.31460,
):
    """
    Calculate offsets for converting the permanent tide from
    a tide-free to a mean-tide state :cite:p:`Mathews:1997js`

    Parameters
    ----------
    XYZ: xr.Dataset
        Dataset with cartesian coordinates
    h2: float or np.ndarray
        Degree-2 Love number of vertical displacement
    l2: float or np.ndarray
        Degree-2 Love (Shida) number of horizontal displacement
    H0: float, default -0.31460
        Mean amplitude of the permanent tide (meters)

    Returns
    -------
    D: xr.Dataset
        free-to-mean tide offset
    """
    # compute the normalized position vector of coordinates
    radius = pyTMD.math.radius(XYZ["X"], XYZ["Y"], XYZ["Z"])
    # sine and cosine of (geocentric) latitude
    sinphi = XYZ["Z"] / radius
    cosphi = np.hypot(XYZ["X"], XYZ["Y"]) / radius
    # sine and cosine of longitude
    sinla = XYZ["Y"] / cosphi / radius
    cosla = XYZ["X"] / cosphi / radius
    # spherical harmonic normalization of degree 2 (order 0)
    dfactor = np.sqrt(5.0 / (4.0 * np.pi))
    # in Mathews et al. (1997):
    # dR0=-0.1196 m with h2=0.6026
    # dN0=-0.0165 m with l2=0.0831
    dR0 = dfactor * h2 * H0
    dN0 = dfactor * l2 * H0
    # calculate offsets in local coordinates
    DR = 0.5 * dR0 * (3.0 * sinphi**2 - 1.0)
    DN = 3.0 * dN0 * cosphi * sinphi
    # compute as an additive correction (Mathews et al. 1997)
    D = xr.Dataset()
    D["X"] = -DR * cosla * cosphi + DN * cosla * sinphi
    D["Y"] = -DR * sinla * cosphi + DN * sinla * sinphi
    D["Z"] = -DR * sinphi - DN * cosphi
    # return the corrections
    return D
