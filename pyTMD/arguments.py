#!/usr/bin/env python
u"""
arguments.py
Written by Tyler Sutterley (12/2024)
Calculates the nodal corrections for tidal constituents
Modification of ARGUMENTS fortran subroutine by Richard Ray 03/1999

CALLING SEQUENCE:
    pu, pf, G = arguments(MJD, constituents)

INPUTS:
    MJD: Modified Julian Day of input date
    constituents: tidal constituent IDs

OUTPUTS:
    pu, pf: nodal corrections for the constituents
    G: phase correction in degrees

OPTIONS:
    deltat: time correction for converting to Ephemeris Time (days)
    corrections: use nodal corrections from OTIS, FES or GOT models

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

PROGRAM DEPENDENCIES:
    astro.py: computes the basic astronomical mean longitudes
    math.py: Special functions of mathematical physics

REFERENCES:
    A. T. Doodson and H. Warburg, "Admiralty Manual of Tides", HMSO, (1941).
    P. Schureman, "Manual of Harmonic Analysis and Prediction of Tides"
        US Coast and Geodetic Survey, Special Publication, 98, (1958).
    M. G. G. Foreman and R. F. Henry, "The harmonic analysis of tidal model
        time series", Advances in Water Resources, 12, (1989).
    G. D. Egbert and S. Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).

UPDATE HISTORY:
    Updated 12/2024: added function to calculate tidal aliasing periods
    Updated 11/2024: allow variable case for Doodson number formalisms
        fix species in constituent parameters for complex tides
        move body tide Love/Shida numbers from predict module
    Updated 10/2024: can convert Doodson numbers formatted as strings
        update Doodson number conversions to follow Cartwright X=10 convention
        add function to parse Cartwright/Tayler/Edden tables
        add functions to calculate UKHO Extended Doodson numbers for constituents
    Updated 09/2024: add function to calculate tidal angular frequencies
    Updated 08/2024: add support for constituents in PERTH5 tables
        add back nodal arguments from PERTH3 for backwards compatibility
    Updated 01/2024: add function to create arguments coefficients table
        add function to calculate the arguments for minor constituents
        include multiples of 90 degrees as variable following Ray 2017
        add function to calculate the Doodson numbers for constituents
        add option to return None and not raise error for Doodson numbers
        moved constituent parameters function from predict to arguments
        added more constituent parameters for OTIS/ATLAS predictions
    Updated 12/2023: made keyword argument for selecting M1 coefficients
    Updated 08/2023: changed ESR netCDF4 format to TMD3 format
    Updated 04/2023: using renamed astro mean_longitudes function
        function renamed from original load_nodal_corrections
    Updated 03/2023: add basic variable typing to function inputs
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 05/2022: added ESR netCDF4 formats to list of model types
        changed keyword arguments to camel case
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 12/2020: fix k1 for FES models
    Updated 08/2020: change time variable names to not overwrite functions
        update nodal corrections for FES models
    Updated 07/2020: added function docstrings.  add shallow water constituents
    Updated 09/2019: added netcdf option to corrections option
    Updated 08/2018: added correction option ATLAS for localized OTIS solutions
    Updated 07/2018: added option to use GSFC GOT nodal corrections
    Updated 09/2017: Rewritten in Python
    Rewritten in Matlab by Lana Erofeeva 01/2003
    Written by Richard Ray 03/1999
"""
from __future__ import annotations

import re
import json
import pathlib
import numpy as np
import pyTMD.astro
from pyTMD.utilities import get_data_path

__all__ = [
    "arguments",
    "minor_arguments",
    "coefficients_table",
    "doodson_number",
    "nodal",
    "frequency",
    "aliasing_period",
    "_arguments_table",
    "_minor_table",
    "_constituent_parameters",
    "_love_numbers",
    "_parse_tide_potential_table",
    "_to_doodson_number",
    "_to_extended_doodson",
    "_from_doodson_number",
    "_from_extended_doodson"
]

def arguments(
        MJD: np.ndarray,
        constituents: list | np.ndarray,
        **kwargs
    ):
    """
    Calculates the nodal corrections for tidal constituents
    :cite:p:`Doodson:1941td` :cite:p:`Schureman:1958ty` :cite:p:`Foreman:1989dt` :cite:p:`Egbert:2002ge`

    Parameters
    ----------
    MJD: np.ndarray
        modified Julian day of input date
    constituents: list
        tidal constituent IDs
    deltat: float or np.ndarray, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS, FES or GOT models
    M1: str, default 'perth5'
        coefficients to use for M1 tides

                - ``'Doodson'``
                - ``'Ray'``
                - ``'perth5'``

    Returns
    -------
    pu: np.ndarray
        nodal angle correction
    pf: np.ndarray
        nodal factor correction
    G: np.ndarray
        phase correction in degrees
    """
    # set default keyword arguments
    kwargs.setdefault('deltat', 0.0)
    kwargs.setdefault('corrections', 'OTIS')
    kwargs.setdefault('M1', 'perth5')

    # set function for astronomical longitudes
    # use ASTRO5 routines if not using an OTIS type model
    ASTRO5 = kwargs['corrections'] not in ('OTIS','ATLAS','TMD3','netcdf')
    # convert from Modified Julian Dates into Ephemeris Time
    s, h, p, n, pp = pyTMD.astro.mean_longitudes(MJD + kwargs['deltat'],
        ASTRO5=ASTRO5)

    # number of temporal values
    nt = len(np.atleast_1d(MJD))
    # initial time conversions
    hour = 24.0*np.mod(MJD, 1)
    # convert from hours solar time into mean lunar time in degrees
    tau = 15.0*hour - s + h
    # variable for multiples of 90 degrees (Ray technical note 2017)
    k = 90.0 + np.zeros((nt))

    # determine equilibrium arguments
    fargs = np.c_[tau, s, h, p, n, pp, k]
    G = np.dot(fargs, coefficients_table(constituents, **kwargs))

    # set nodal corrections
    # determine nodal corrections f and u for each model type
    pu, pf = nodal(n, p, constituents, **kwargs)

    # return values as tuple
    return (pu, pf, G)

def minor_arguments(
        MJD: np.ndarray,
        **kwargs
    ):
    """
    Calculates the nodal corrections for minor tidal constituents
    in order to infer their values :cite:p:`Doodson:1941td` :cite:p:`Schureman:1958ty`
    :cite:p:`Foreman:1989dt` :cite:p:`Egbert:2002ge`


    Parameters
    ----------
    MJD: np.ndarray
        modified Julian day of input date
    deltat: float or np.ndarray, default 0.0
        time correction for converting to Ephemeris Time (days)
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS, FES or GOT models

    Returns
    -------
    pu: np.ndarray
        nodal angle correction
    pf: np.ndarray
        nodal factor correction
    G: np.ndarray
        phase correction in degrees
    """
    # set default keyword arguments
    kwargs.setdefault('deltat', 0.0)
    kwargs.setdefault('corrections', 'OTIS')

    # degrees to radians
    dtr = np.pi/180.0
    # set function for astronomical longitudes
    # use ASTRO5 routines if not using an OTIS type model
    ASTRO5 = kwargs['corrections'] not in ('OTIS','ATLAS','TMD3','netcdf')
    # convert from Modified Julian Dates into Ephemeris Time
    s, h, p, n, pp = pyTMD.astro.mean_longitudes(MJD + kwargs['deltat'],
        ASTRO5=ASTRO5)

    # number of temporal values
    nt = len(np.atleast_1d(MJD))
    # initial time conversions
    hour = 24.0*np.mod(MJD, 1)
    # convert from hours solar time into mean lunar time in degrees
    tau = 15.0*hour - s + h
    # variable for multiples of 90 degrees (Ray technical note 2017)
    k = 90.0 + np.zeros((nt))

    # determine equilibrium arguments
    fargs = np.c_[tau, s, h, p, n, pp, k]
    arg = np.dot(fargs, _minor_table())

    # determine nodal corrections f and u
    sinn = np.sin(n*dtr)
    cosn = np.cos(n*dtr)
    sin2n = np.sin(2.0*n*dtr)
    cos2n = np.cos(2.0*n*dtr)

    # nodal factor corrections for minor constituents
    f = np.ones((nt, 20))
    f[:,0] = np.sqrt((1.0 + 0.189*cosn - 0.0058*cos2n)**2 +
        (0.189*sinn - 0.0058*sin2n)**2)# 2Q1
    f[:,1] = f[:,0]# sigma1
    f[:,2] = f[:,0]# rho1
    f[:,3] = np.sqrt((1.0 + 0.185*cosn)**2 + (0.185*sinn)**2)# M12
    f[:,4] = np.sqrt((1.0 + 0.201*cosn)**2 + (0.201*sinn)**2)# M11
    f[:,5] = np.sqrt((1.0 + 0.221*cosn)**2 + (0.221*sinn)**2)# chi1
    f[:,9] = np.sqrt((1.0 + 0.198*cosn)**2 + (0.198*sinn)**2)# J1
    f[:,10] = np.sqrt((1.0 + 0.640*cosn + 0.134*cos2n)**2 +
        (0.640*sinn + 0.134*sin2n)**2)# OO1
    f[:,11] = np.sqrt((1.0 - 0.0373*cosn)**2 + (0.0373*sinn)**2)# 2N2
    f[:,12] = f[:,11]# mu2
    f[:,13] = f[:,11]# nu2
    f[:,15] = f[:,11]# L2
    f[:,16] = np.sqrt((1.0 + 0.441*cosn)**2 + (0.441*sinn)**2)# L2

    # nodal angle corrections for minor constituents
    u = np.zeros((nt, 20))
    u[:,0] = np.arctan2(0.189*sinn - 0.0058*sin2n,
        1.0 + 0.189*cosn - 0.0058*sin2n)# 2Q1
    u[:,1] = u[:,0]# sigma1
    u[:,2] = u[:,0]# rho1
    u[:,3] = np.arctan2( 0.185*sinn, 1.0 + 0.185*cosn)# M12
    u[:,4] = np.arctan2(-0.201*sinn, 1.0 + 0.201*cosn)# M11
    u[:,5] = np.arctan2(-0.221*sinn, 1.0 + 0.221*cosn)# chi1
    u[:,9] = np.arctan2(-0.198*sinn, 1.0 + 0.198*cosn)# J1
    u[:,10] = np.arctan2(-0.640*sinn - 0.134*sin2n,
        1.0 + 0.640*cosn + 0.134*cos2n)# OO1
    u[:,11] = np.arctan2(-0.0373*sinn, 1.0 - 0.0373*cosn)# 2N2
    u[:,12] = u[:,11]# mu2
    u[:,13] = u[:,11]# nu2
    u[:,15] = u[:,11]# L2
    u[:,16] = np.arctan2(-0.441*sinn, 1.0 + 0.441*cosn)# L2

    if kwargs['corrections'] in ('FES',):
        # additional astronomical terms for FES models
        II = np.arccos(0.913694997 - 0.035692561*np.cos(n*dtr))
        at1 = np.arctan(1.01883*np.tan(n*dtr/2.0))
        at2 = np.arctan(0.64412*np.tan(n*dtr/2.0))
        xi = -at1 - at2 + n*dtr
        xi = np.arctan2(np.sin(xi), np.cos(xi))
        nu = at1 - at2
        I2 = np.tan(II/2.0)
        Ra1 = np.sqrt(1.0 - 12.0*(I2**2)*np.cos(2.0*(p - xi)) + 36.0*(I2**4))
        P2 = np.sin(2.0*(p - xi))
        Q2 = 1.0/(6.0*(I2**2)) - np.cos(2.0*(p - xi))
        R = np.arctan(P2/Q2)

        # nodal factor corrections for minor constituents
        f[:,0] = np.sin(II)*(np.cos(II/2.0)**2)/0.38 # 2Q1
        f[:,1] = f[:,0] # sigma1
        f[:,2] = f[:,0] # rho1
        f[:,3] = f[:,0] # M12
        f[:,4] = np.sin(2.0*II)/0.7214 # M11
        f[:,5] = f[:,4] # chi1
        f[:,9] = f[:,4] # J1
        f[:,10] = np.sin(II)*np.power(np.sin(II/2.0),2.0)/0.01640 # OO1
        f[:,11] = np.power(np.cos(II/2.0),4.0)/0.9154 # 2N2
        f[:,12] = f[:,11] # mu2
        f[:,13] = f[:,11] # nu2
        f[:,14] = f[:,11] # lambda2
        f[:,15] = f[:,11]*Ra1 # L2
        f[:,18] = f[:,11] # eps2
        f[:,19] = np.power(np.sin(II),2.0)/0.1565 # eta2

        # nodal angle corrections for minor constituents
        u[:,0] = (2.0*xi - nu) # 2Q1
        u[:,1] = u[:,0] # sigma1
        u[:,2] = u[:,0] # rho1
        u[:,3] = u[:,0] # M12
        u[:,4] = -nu # M11
        u[:,5] = u[:,4] # chi1
        u[:,9] = u[:,4] # J1
        u[:,10] = (-2.0*xi - nu) # OO1
        u[:,11] = (2.0*xi - 2.0*nu) # 2N2
        u[:,12] = u[:,11] # mu2
        u[:,13] = u[:,11] # nu2
        u[:,14] = (2.0*xi - 2.0*nu) # lambda2
        u[:,15] = (2.0*xi - 2.0*nu - R)# L2
        u[:,18] = u[:,12] # eps2
        u[:,19] = -2.0*nu # eta2
    elif kwargs['corrections'] in ('GOT',):
        f[:,18] = f[:,11] # eps2
        f[:,19] = np.sqrt((1.0 + 0.436*cosn)**2 + (0.436*sinn)**2) # eta2
        u[:,18] = u[:,11] # eps2
        u[:,19] = np.arctan(-0.436*sinn/(1.0 + 0.436*cosn)) # eta2

    # return values as tuple
    return (u, f, arg)

# JSON file of Doodson coefficients
_coefficients_table = get_data_path(['data','doodson.json'])

def coefficients_table(
        constituents: list | tuple | np.ndarray | str,
        **kwargs
    ):
    """
    Doodson table coefficients for tidal constituents
    :cite:p:`Doodson:1921kt` :cite:p:`Doodson:1941td`

    Parameters
    ----------
    constituents: list, tuple, np.ndarray or str
        tidal constituent IDs
    corrections: str, default 'OTIS'
        use coefficients from OTIS, FES or GOT models
    file: str or pathlib.Path, default `coefficients.json`
        JSON file of Doodson coefficients

    Returns
    -------
    coef: np.ndarray
        Doodson coefficients (Cartwright numbers) for each constituent
    """
    # set default keyword arguments
    kwargs.setdefault('corrections', 'OTIS')
    kwargs.setdefault('file', _coefficients_table)

    # verify coefficients table path
    table = pathlib.Path(kwargs['file']).expanduser().absolute()
    # modified Doodson coefficients for constituents
    # using 7 index variables: tau, s, h, p, n, pp, k
    # tau: mean lunar time
    # s: mean longitude of moon
    # h: mean longitude of sun
    # p: mean longitude of lunar perigee
    # n: mean longitude of ascending lunar node
    # pp: mean longitude of solar perigee
    # k: 90-degree phase
    with table.open(mode='r', encoding='utf8') as fid:
        coefficients = json.load(fid)

    # # Without p'
    # coefficients['sa'] = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    # coefficients['sta'] = [0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0]
    # set s1 coefficients
    if kwargs['corrections'] in ('OTIS','ATLAS','TMD3','netcdf'):
        coefficients['s1'] = [1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 1.0]

    # set constituents to be iterable
    if isinstance(constituents, str):
        constituents = [constituents]
    # allocate for output coefficients
    nc = len(constituents)
    coef = np.zeros((7, nc))
    # for each constituent of interest
    for i, c in enumerate(constituents):
        try:
            coef[:,i] = coefficients[c]
        except KeyError:
            raise ValueError(f'Unsupported constituent: {c}')

    # return Doodson coefficients for constituents
    return coef

def doodson_number(
        constituents: str | list | np.ndarray,
        **kwargs
    ):
    """
    Calculates the Doodson or Cartwright number for
    tidal constituents :cite:p:`Doodson:1921kt`

    Parameters
    ----------
    constituents: str, list or np.ndarray
        tidal constituent ID(s)
    corrections: str, default 'OTIS'
        use arguments from OTIS, FES or GOT models
    formalism: str, default 'Doodson'
        constituent identifier formalism

            - ``'Cartwright'``
            - ``'Doodson'``
            - ``'Extended'``
    raise_error: bool, default True
        Raise exception if constituent is unsupported

    Returns
    -------
    numbers: float, np.ndarray or dict
        Doodson or Cartwright number for each constituent
    """
    # set default keyword arguments
    kwargs.setdefault('corrections', 'OTIS')
    kwargs.setdefault('formalism', 'Doodson')
    kwargs.setdefault('raise_error', True)
    # validate inputs
    assert kwargs['formalism'].title() in ('Cartwright', 'Doodson','Extended'), \
        f'Unknown formalism {kwargs["formalism"]}'
    # get the coefficients of coefficients
    if isinstance(constituents, str):
        # try to get the Doodson coefficients for constituent
        try:
            coefficients = coefficients_table(constituents.lower(), **kwargs)
        except ValueError as exc:
            if kwargs['raise_error']:
                raise ValueError(f'Unsupported constituent {constituents}')
            else:
                return None
        # extract identifier in formalism
        if (kwargs['formalism'].title() == 'Cartwright'):
            # extract Cartwright number
            numbers = np.array(coefficients[:6,0])
        elif (kwargs['formalism'].title() == 'Doodson'):
            # convert from coefficients to Doodson number
            numbers = _to_doodson_number(coefficients[:,0], **kwargs)
        elif (kwargs['formalism'].title() == 'Extended'):
            # convert to extended Doodson number in UKHO format
            numbers = _to_extended_doodson(coefficients[:,0], **kwargs)
    else:
        # output dictionary with Doodson numbers
        numbers = {}
        # for each input constituent
        for i,c in enumerate(constituents):
            # try to get the Doodson coefficients for constituent
            try:
                coefficients = coefficients_table(c.lower(), **kwargs)
            except ValueError as exc:
                if kwargs['raise_error']:
                    raise ValueError(f'Unsupported constituent {c}')
                else:
                    numbers[c] = None
                    continue
            # convert from coefficients to Doodson number
            if (kwargs['formalism'].title() == 'Cartwright'):
                # extract Cartwright number
                numbers[c] = np.array(coefficients[:6,0])
            elif (kwargs['formalism'].title() == 'Doodson'):
                # convert from coefficients to Doodson number
                numbers[c] = _to_doodson_number(coefficients[:,0], **kwargs)
            elif (kwargs['formalism'].title() == 'Extended'):
                # convert to extended Doodson number in UKHO format
                numbers[c] = _to_extended_doodson(coefficients[:,0], **kwargs)
    # return the Doodson or Cartwright number
    return numbers

# PURPOSE: compute the nodal corrections
def nodal(
        n: np.ndarray,
        p: np.ndarray,
        constituents: list | tuple | np.ndarray | str,
        **kwargs
    ):
    """
    Calculates the nodal corrections for tidal constituents
    :cite:p:`Doodson:1941td` :cite:p:`Schureman:1958ty` :cite:p:`Foreman:1989dt` :cite:p:`Ray:1999vm`

    Calculates factors for compound tides using recursion

    Parameters
    ----------
    n: np.ndarray
        mean longitude of ascending lunar node (degrees)
    p: np.ndarray
        mean longitude of lunar perigee (degrees)
    constituents: list, tuple, np.ndarray or str
        tidal constituent IDs
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS, FES or GOT models
    M1: str, default 'perth5'
        coefficients to use for M1 tides

                - ``'Doodson'``
                - ``'Ray'``
                - ``'perth5'``

    Returns
    -------
    f: np.ndarray
        nodal factor correction
    u: np.ndarray
        nodal angle correction
    """
    # set default keyword arguments
    kwargs.setdefault('corrections', 'OTIS')
    kwargs.setdefault('M1', 'perth5')
    # set correction type
    OTIS_TYPE = kwargs['corrections'] in ('OTIS','ATLAS','TMD3','netcdf')
    FES_TYPE = kwargs['corrections'] in ('FES',)
    PERTH3_TYPE = kwargs['corrections'] in ('perth3',)

    # degrees to radians
    dtr = np.pi/180.0
    # trigonometric factors for nodal corrections
    sinn = np.sin(n*dtr)
    cosn = np.cos(n*dtr)
    sin2n = np.sin(2.0*n*dtr)
    cos2n = np.cos(2.0*n*dtr)
    sin3n = np.sin(3.0*n*dtr)
    sinp  = np.sin(p*dtr)
    cosp  = np.cos(p*dtr)
    sin2p = np.sin(2.0*p*dtr)
    cos2p = np.cos(2.0*p*dtr)

    # set constituents to be iterable
    if isinstance(constituents, str):
        constituents = [constituents]

    # set nodal corrections
    nt = len(np.atleast_1d(n))
    nc = len(constituents)
    # nodal factor correction
    f = np.zeros((nt, nc))
    # nodal angle correction
    u = np.zeros((nt, nc))

    # additional astronomical terms for FES models
    II = np.arccos(0.913694997 - 0.035692561*np.cos(n*dtr))
    at1 = np.arctan(1.01883*np.tan(n*dtr/2.0))
    at2 = np.arctan(0.64412*np.tan(n*dtr/2.0))
    xi = -at1 - at2 + n*dtr
    xi = np.arctan2(np.sin(xi), np.cos(xi))
    nu = at1 - at2
    I2 = np.tan(II/2.0)
    Ra1 = np.sqrt(1.0 - 12.0*(I2**2)*np.cos(2.0*(p - xi)) + 36.0*(I2**4))
    P2 = np.sin(2.0*(p - xi))
    Q2 = 1.0/(6.0*(I2**2)) - np.cos(2.0*(p - xi))
    R = np.arctan(P2/Q2)
    P_prime = np.sin(2.0*II)*np.sin(nu)
    Q_prime = np.sin(2.0*II)*np.cos(nu) + 0.3347
    nu_prime = np.arctan(P_prime/Q_prime)
    P_sec = (np.sin(II)**2)*np.sin(2.0*nu)
    Q_sec = (np.sin(II)**2)*np.cos(2.0*nu) + 0.0727
    nu_sec = 0.5*np.arctan(P_sec/Q_sec)

    # compute standard nodal corrections f and u
    for i, c in enumerate(constituents):
        if c in ('msf','tau1','p1','theta1','lambda2','s2') and OTIS_TYPE:
            term1 = 0.0
            term2 = 1.0
        elif c in ('p1','s2') and (FES_TYPE or PERTH3_TYPE):
            term1 = 0.0
            term2 = 1.0
        elif c in ('mm','msm') and OTIS_TYPE:
            term1 = 0.0
            term2 = 1.0 - 0.130*cosn
        elif c in ('mm','msm') and FES_TYPE:
            term1 = 0.0
            term2 = (2.0/3.0 - np.power(np.sin(II),2.0))/0.5021
        elif c in ('mm','msm'):
            term1 = -0.0534*sin2p - 0.0219*np.sin((2.0*p-n)*dtr)
            term2 = 1.0 - 0.1308*cosn - 0.0534*cos2p - 0.0219*np.cos((2.0*p-n)*dtr)
        elif c in ('mf','msqm','msp','mq','mtm') and OTIS_TYPE:
            f[:,i] = 1.043 + 0.414*cosn
            u[:,i] = dtr*(-23.7*sinn + 2.7*sin2n - 0.4*sin3n)
            continue
        elif c in ('mf','msqm','msp','mq','mt','mtm') and FES_TYPE:
            f[:,i] = np.power(np.sin(II),2.0)/0.1578
            u[:,i] = -2.0*xi
            continue
        elif c in ('mf','msqm','msp','mq'):
            term1 = -0.04324*sin2p - 0.41465*sinn - 0.03873*sin2n
            term2 = 1.0 + 0.04324*cos2p + 0.41465*cosn + 0.03873*cos2n
        elif c in ('mt',) and OTIS_TYPE:
            term1 = -0.203*sinn - 0.040*sin2n
            term2 = 1.0 + 0.203*cosn + 0.040*cos2n
        elif c in ('mt','mtm',):
            term1 = -0.018*sin2p - 0.4145*sinn - 0.040*sin2n
            term2 = 1.0 + 0.018*cos2p + 0.4145*cosn + 0.040*cos2n
        elif c in ('msf',) and FES_TYPE:
            f[:,i] = 1.0
            u[:,i] = (2.0*xi - 2.0*nu)
            continue
        elif c in ('msf',):
            # linear tide and not compound
            term1 = 0.137*sinn
            term2 = 1.0
        elif c in ('mst',):
            term1 = -0.380*sin2p - 0.413*sinn - 0.037*sin2n
            term2 = 1.0 + 0.380*cos2p + 0.413*cosn + 0.037*cos2n
        elif c in ('o1','so3','op2') and OTIS_TYPE:
            term1 = 0.189*sinn - 0.0058*sin2n
            term2 = 1.0 + 0.189*cosn - 0.0058*cos2n
            f[:,i] = np.sqrt(term1**2 + term2**2) # O1
            u[:,i] = dtr*(10.8*sinn - 1.3*sin2n + 0.2*sin3n)
            continue
        elif c in ('o1','so3','op2','2q1','q1','rho1','sigma1') and FES_TYPE:
            f[:,i] = np.sin(II)*(np.cos(II/2.0)**2)/0.38
            u[:,i] = (2.0*xi - nu)
            continue
        elif c in ('q1','o1') and PERTH3_TYPE:
            f[:,i] = 1.009 + 0.187*cosn - 0.015*cos2n
            u[:,i] = dtr*(10.8*sinn - 1.3*sin2n)
            continue
        elif c in ('o1','so3','op2'):
            term1 = 0.1886*sinn - 0.0058*sin2n - 0.0065*sin2p
            term2 = 1.0 + 0.1886*cosn - 0.0058*cos2n - 0.0065*cos2p
        elif c in ('2q1','q1','rho1','sigma1') and OTIS_TYPE:
            f[:,i] = np.sqrt((1.0 + 0.188*cosn)**2 + (0.188*sinn)**2)
            u[:,i] = np.arctan(0.189*sinn/(1.0 + 0.189*cosn))
            continue
        elif c in ('2q1','q1','rho1','sigma1'):
            term1 = 0.1886*sinn
            term2 = 1.0 + 0.1886*cosn
        elif c in ('tau1',):
            term1 = 0.219*sinn
            term2 = 1.0 - 0.219*cosn
        elif c in ('beta1',):
            term1 = 0.226*sinn
            term2 = 1.0 + 0.226*cosn
        elif c in ('m1',) and (kwargs['M1'] == 'Doodson'):
            # A. T. Doodson's coefficients for M1 tides
            term1 = sinp + 0.2*np.sin((p-n)*dtr)
            term2 = 2.0*cosp + 0.4*np.cos((p-n)*dtr)
        elif c in ('m1',) and (kwargs['M1'] == 'Ray'):
            # R. Ray's coefficients for M1 tides (perth3)
            term1 = 0.64*sinp + 0.135*np.sin((p-n)*dtr)
            term2 = 1.36*cosp + 0.267*np.cos((p-n)*dtr)
        elif c in ('m1',) and (kwargs['M1'] == 'perth5'):
            # assumes M1 argument includes p
            term1 = -0.2294*sinn - 0.3594*sin2p - 0.0664*np.sin((2.0*p-n)*dtr)
            term2 = 1.0 + 0.1722*cosn + 0.3594*cos2p + 0.0664*np.cos((2.0*p-n)*dtr)
        elif c in ('chi1',) and OTIS_TYPE:
            term1 = -0.221*sinn
            term2 = 1.0 + 0.221*cosn
        elif c in ('chi1','theta1','j1') and FES_TYPE:
            f[:,i] = np.sin(2.0*II) / 0.7214
            u[:,i] = -nu
            continue
        elif c in ('chi1',):
            term1 = -0.250*sinn
            term2 = 1.0 + 0.193*cosn
        elif c in ('p1',):
            term1 = -0.0112*sinn
            term2 = 1.0 - 0.0112*cosn
        elif c in ('k1','sk3','2sk5') and OTIS_TYPE:
            term1 = -0.1554*sinn + 0.0029*sin2n
            term2 = 1.0 + 0.1158*cosn - 0.0029*cos2n
        elif c in ('k1','sk3','2sk5') and FES_TYPE:
            temp1 = 0.8965*np.power(np.sin(2.0*II),2.0)
            temp2 = 0.6001*np.sin(2.0*II)*np.cos(nu)
            f[:,i] = np.sqrt(temp1 + temp2 + 0.1006)
            u[:,i] = -nu_prime
            continue
        elif c in ('k1',) and PERTH3_TYPE:
            f[:,i] = 1.006 + 0.115*cosn - 0.009*cos2n
            u[:,i] = dtr*(-8.9*sinn + 0.7*sin2n)
            continue
        elif c in ('k1','sk3','2sk5'):
            term1 = -0.1554*sinn + 0.0031*sin2n
            term2 = 1.0 + 0.1158*cosn - 0.0028*cos2n
        elif c in ('j1','theta1'):
            term1 = -0.227*sinn
            term2 = 1.0 + 0.169*cosn
        elif c in ('oo1','ups1') and OTIS_TYPE:
            term1 = -0.640*sinn - 0.134*sin2n
            term2 = 1.0 + 0.640*cosn + 0.134*cos2n
        elif c in ('oo1','ups1') and FES_TYPE:
            f[:,i] = np.sin(II)*np.power(np.sin(II/2.0),2.0)/0.01640
            u[:,i] = -2.0*xi - nu
            continue
        elif c in ('oo1','ups1'):
            term1 = -0.640*sinn - 0.134*sin2n - 0.150*sin2p
            term2 = 1.0 + 0.640*cosn + 0.134*cos2n + 0.150*cos2p
        elif c in ('m2','2n2','mu2','n2','nu2','lambda2','ms4','eps2','2sm6',
                '2sn6','mp1','mp3','sn4') and FES_TYPE:
            f[:,i] = np.power(np.cos(II/2.0),4.0)/0.9154
            u[:,i] = 2.0*xi - 2.0*nu
            continue
        elif c in ('m2','n2') and PERTH3_TYPE:
            f[:,i] = 1.000 - 0.037*cosn
            u[:,i] = dtr*(-2.1*sinn)
            continue
        elif c in ('m2','2n2','mu2','n2','nu2','lambda2','ms4','eps2','2sm6',
                '2sn6','mp1','mp3','sn4'):
            term1 = -0.03731*sinn + 0.00052*sin2n
            term2 = 1.0 - 0.03731*cosn + 0.00052*cos2n
        elif c in ('l2','sl4') and OTIS_TYPE:
            term1 = -0.25*sin2p - 0.11*np.sin((2.0*p-n)*dtr) - 0.04*sinn
            term2 = 1.0 - 0.25*cos2p - 0.11*np.cos((2.0*p - n)*dtr) - 0.04*cosn
        elif c in ('l2','sl4') and FES_TYPE:
            f[:,i] = Ra1*np.power(np.cos(II/2.0),4.0)/0.9154
            u[:,i] = 2.0*xi - 2.0*nu - R
            continue
        elif c in ('l2','sl4'):
            term1 = -0.25*sin2p - 0.11*np.sin((2.0*p-n)*dtr) - 0.037*sinn
            term2 = 1.0 - 0.25*cos2p - 0.11*np.cos((2.0*p-n)*dtr) - 0.037*cosn
        elif c in ('l2b',):
            # for when l2 is split into two constituents
            term1 = 0.441*sinn
            term2 = 1.0 + 0.441*cosn
        elif c in ('k2','sk4','2sk6','kp1') and OTIS_TYPE:
            term1 = -0.3108*sinn - 0.0324*sin2n
            term2 = 1.0 + 0.2852*cosn + 0.0324*cos2n
        elif c in ('k2','sk4','2sk6','kp1') and FES_TYPE:
            term1 = 19.0444 * np.power(np.sin(II),4.0)
            term2 = 2.7702 * np.power(np.sin(II),2.0) * np.cos(2.0*nu)
            f[:,i] = np.sqrt(term1 + term2 + 0.0981)
            u[:,i] = -2.0*nu_sec
            continue
        elif c in ('k2',) and PERTH3_TYPE:
            f[:,i] = 1.024 + 0.286*cosn + 0.008*cos2n
            u[:,i] = dtr*(-17.7*sinn + 0.7*sin2n)
            continue
        elif c in ('k2','sk4','2sk6','kp1'):
            term1 = -0.3108*sinn - 0.0324*sin2n
            term2 = 1.0 + 0.2853*cosn + 0.0324*cos2n
        elif c in ('gamma2',):
            term1 = 0.147*np.sin(2.0*(n-p)*dtr)
            term2 = 1.0 + 0.147*np.cos(2.0*(n-p)*dtr)
        elif c in ('delta2',):
            term1 = 0.505*sin2p + 0.505*sinn - 0.165*sin2n
            term2 = 1.0 - 0.505*cos2p - 0.505*cosn + 0.165*cos2n
        elif c in ('eta2','zeta2') and FES_TYPE:
            f[:,i] = np.power(np.sin(II),2.0)/0.1565
            u[:,i] = -2.0*nu
            continue
        elif c in ('eta2','zeta2'):
            term1 = -0.436*sinn
            term2 = 1.0 + 0.436*cosn
        elif c in ('s2',):
            term1 = 0.00225*sinn
            term2 = 1.0 + 0.00225*cosn
        elif c in ("m1'",):
            # Linear 3rd degree terms
            term1 = -0.01815*sinn
            term2 = 1.0 - 0.27837*cosn
        elif c in ("q1'",):
            # Linear 3rd degree terms
            term1 = 0.3915*sinn + 0.033*sin2n + 0.061*sin2p
            term2 = 1.0 + 0.3915*cosn + 0.033*cos2n + 0.06*cos2p
        elif c in ("j1'",):
            # Linear 3rd degree terms
            term1 = -0.438*sinn - 0.033*sin2n
            term2 = 1.0 + 0.372*cosn + 0.033*cos2n
        elif c in ("2n2'",):
            # Linear 3rd degree terms
            term1 = 0.166*sinn
            term2 = 1.0 + 0.166*cosn
        elif c in ("n2'",):
            # Linear 3rd degree terms
            term1 = 0.1705*sinn - 0.0035*sin2n - 0.0176*sin2p
            term2 = 1.0 + 0.1705*cosn - 0.0035*cos2n - 0.0176*cos2p
        elif c in ("l2'"):
            # Linear 3rd degree terms
            term1 = -0.2495*sinn
            term2 = 1.0 + 0.1315*cosn
        elif c in ('m3',) and FES_TYPE:
            f[:,i] = np.power(np.cos(II/2.0), 6.0) / 0.8758
            u[:,i] = (3.0*xi - 3.0*nu)
            continue
        elif c in ('m3','e3'):
            # Linear 3rd degree terms
            term1 = -0.05644*sinn
            term2 = 1.0 - 0.05644*cosn
        elif c in ('j3','f3'):
            term1 = -0.464*sinn - 0.052*sin2n
            term2 = 1.0 + 0.387*cosn + 0.052*cos2n
        elif c in ('l3',):
            term1 = -0.373*sin2p - 0.164*np.sin((2.0*p-n)*dtr)
            term2 = 1.0 - 0.373*cos2p - 0.164*np.cos((2.0*p-n)*dtr)
        elif c in ('mfdw',):
            # special test of Doodson-Warburg formula
            f[:,i] = 1.043 + 0.414*cosn
            u[:,i] = dtr*(-23.7*sinn + 2.7*sin2n - 0.4*sin3n)
            continue
        elif c in ('so1','2so3','2po1'):
            # compound tides calculated using recursion
            parents = ['o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]
            u[:,i] = -utmp[:,0]
            continue
        elif c in ('o3',):
            # compound tides calculated using recursion
            parents = ['o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3
            u[:,i] = 3.0*utmp[:,0]
            continue
        elif c in ('2k2'):
            # compound tides calculated using recursion
            parents = ['k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2
            u[:,i] = 2.0*utmp[:,0]
            continue
        elif c in ('tk1'):
            # compound tides calculated using recursion
            parents = ['k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]
            u[:,i] = -utmp[:,0]
            continue
        elif c in ('2oop1'):
            # compound tides calculated using recursion
            parents = ['oo1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2
            u[:,i] = 2.0*utmp[:,0]
            continue
        elif c in ('oq2'):
            # compound tides calculated using recursion
            parents = ['o1','q1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('2oq1'):
            # compound tides calculated using recursion
            parents = ['o1','q1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('ko2'):
            # compound tides calculated using recursion
            parents = ['o1','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('opk1',):
            # compound tides calculated using recursion
            parents = ['o1','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] - utmp[:,1]
            continue
        elif c in ('2ook1',):
            # compound tides calculated using recursion
            parents = ['oo1','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('kj2',):
            # compound tides calculated using recursion
            parents = ['k1','j1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('kjq1'):
            # compound tides calculated using recursion
            parents = ['k1','j1','q1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1] * ftmp[:,2]
            u[:,i] = utmp[:,0] + utmp[:,1] - utmp[:,2]
            continue
        elif c in ('k3',):
            # compound tides calculated using recursion
            parents = ['k1','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in('m4','mn4','mns2','2ms2','mnus2','mmus2','2ns2','n4','mnu4',
                'mmu4','2mt6','2ms6','msn6','mns6','2mr6','msmu6','2mp3','2ms3',
                '2mp5','2msp7','2(ms)8','2ms8'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2
            u[:,i] = 2.0*utmp[:,0]
            continue
        elif c in ('msn2','snm2','nsm2'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2
            u[:,i] = 0.0
            continue
        elif c in ('mmun2','2mn2'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3
            u[:,i] = utmp[:,0]
            continue
        elif c in ('2sm2',):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]
            u[:,i] = -utmp[:,0]
            continue
        elif c in ('m6','2mn6','2mnu6','2mmu6','2nm6','mnnu6','mnmu6','3ms8',
                '3mp7','2msn8','3ms5','3mp5','3ms4','3m2s2','3m2s10','2mn2s2'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3
            u[:,i] = 3.0*utmp[:,0]
            continue
        elif c in ('m8','ma8','3mn8','3mnu8','3mmu8','2mn8','2(mn):8','3msn10',
                '4ms10','2(mn)S10','4m2s12'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**4
            u[:,i] = 4.0*utmp[:,0]
            continue
        elif c in ('m10','4mn10','5ms12','4msn12','4mns12'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**5
            u[:,i] = 5.0*utmp[:,0]
            continue
        elif c in ('m12','5mn12','6ms14','5msn14'):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**6
            u[:,i] = 6.0*utmp[:,0]
            continue
        elif c in ('m14',):
            # compound tides calculated using recursion
            parents = ['m2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**7
            u[:,i] = 7.0*utmp[:,0]
            continue
        elif c in ('mo3','no3','mso5'):
            # compound tides calculated using recursion
            parents = ['m2','o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('no1','nso3'):
            # compound tides calculated using recursion
            parents = ['m2','o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] - utmp[:,1]
            continue
        elif c in ('mq3','nq3'):
            # compound tides calculated using recursion
            parents = ['m2','q1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('2mq3',):
            # compound tides calculated using recursion
            parents = ['m2','q1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('2no3',):
            # compound tides calculated using recursion
            parents = ['m2','o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('2mo5','2no5','mno5','2mso7','2(ms):o9'):
            # compound tides calculated using recursion
            parents = ['m2','o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('2mno7','3mo7'):
            # compound tides calculated using recursion
            parents = ['m2','o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3 * ftmp[:,1]
            u[:,i] = 3.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('mk3','nk3','msk5','nsk5'):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('mnk5','2mk5','2nk5','2msk7'):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('2mk3',):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('3mk7','2mnk7','2nmk7','3nk7','3msk9'):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3 * ftmp[:,1]
            u[:,i] = 3.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('3msk7',):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3 * ftmp[:,1]
            u[:,i] = 3.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('4mk9','3mnk9','2m2nk9','2(mn):k9','3nmk9','4msk11'):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**4 * ftmp[:,1]
            u[:,i] = 4.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('3km5',):
            # compound tides calculated using recursion
            parents = ['m2','k1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]**3
            u[:,i] = utmp[:,0] + 3.0*utmp[:,1]
            continue
        elif c in ('mk4','nk4','mks2'):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('msk2','2smk4','msk6','snk6'):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] - utmp[:,1]
            continue
        elif c in ('mnk6','2mk6','2msk8','msnk8'):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('mnk2','2mk2'):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('mkn2','nkm2'):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = utmp[:,1]
            continue
        elif c in ('skm2',):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = -utmp[:,0] + utmp[:,1]
            continue
        elif c in ('3mk8','2mnk8'):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3 * ftmp[:,1]
            u[:,i] = 3.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('m2(ks):2',):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]**2
            u[:,i] = utmp[:,0] + 2.0*utmp[:,1]
            continue
        elif c in ('2ms2k2',):
            # compound tides calculated using recursion
            parents = ['m2','k2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]**2
            u[:,i] = 2.0*utmp[:,0] - 2.0*utmp[:,1]
            continue
        elif c in ('mko5','msko7'):
            # compound tides calculated using recursion
            parents = ['m2','k2','o1']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1] * ftmp[:,2]
            u[:,i] = utmp[:,0] + utmp[:,1] + utmp[:,2]
            continue
        elif c in ('ml4','msl6'):
            # compound tides calculated using recursion
            parents = ['m2','l2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0] * ftmp[:,1]
            u[:,i] = utmp[:,0] + utmp[:,1]
            continue
        elif c in ('2ml2',):
            # compound tides calculated using recursion
            parents = ['m2','l2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] - utmp[:,1]
            continue
        elif c in ('2ml6','2ml2s2','2mls4','2msl8'):
            # compound tides calculated using recursion
            parents = ['m2','l2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**2 * ftmp[:,1]
            u[:,i] = 2.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('2nmls6','3mls6','2mnls6','3ml8','2mnl8','3msl10'):
            # compound tides calculated using recursion
            parents = ['m2','l2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**3 * ftmp[:,1]
            u[:,i] = 3.0*utmp[:,0] + utmp[:,1]
            continue
        elif c in ('4msl12',):
            # compound tides calculated using recursion
            parents = ['m2','l2']
            utmp, ftmp = nodal(n, p, parents, **kwargs)
            f[:,i] = ftmp[:,0]**4 * ftmp[:,1]
            u[:,i] = 4.0*utmp[:,0] + utmp[:,1]
            continue
        else:
            # default for linear tides
            term1 = 0.0
            term2 = 1.0

        # calculate factors for linear tides
        # and parent waves in compound tides
        f[:,i] = np.sqrt(term1**2 + term2**2)
        u[:,i] = np.arctan2(term1, term2)

    # return corrections for constituents
    return (u, f)

def frequency(
        constituents: list | np.ndarray,
        **kwargs
    ):
    """
    Calculates the angular frequency for tidal constituents :cite:p:`Ray:1999vm`

    Parameters
    ----------
    constituents: list
        tidal constituent IDs
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS, FES or GOT models
    M1: str, default 'perth5'
        coefficients to use for M1 tides

                - ``'Doodson'``
                - ``'Ray'``
                - ``'perth5'``

    Returns
    -------
    omega: np.ndarray
        angular frequency in radians per second
    """
    # set default keyword arguments
    kwargs.setdefault('corrections', 'OTIS')
    kwargs.setdefault('M1', 'perth5')
    # set function for astronomical longitudes
    # use ASTRO5 routines if not using an OTIS type model
    ASTRO5 = kwargs['corrections'] not in ('OTIS','ATLAS','TMD3','netcdf')
    # Modified Julian Dates at J2000
    MJD = np.array([51544.5, 51544.55])
    # time interval in seconds
    deltat = 86400.0*(MJD[1] - MJD[0])
    # calculate the mean longitudes of the sun and moon
    s, h, p, n, pp = pyTMD.astro.mean_longitudes(MJD, ASTRO5=ASTRO5)

    # number of temporal values
    nt = len(np.atleast_1d(MJD))
    # initial time conversions
    hour = 24.0*np.mod(MJD, 1)
    # convert from hours solar time into mean lunar time in degrees
    tau = 15.0*hour - s + h
    # variable for multiples of 90 degrees (Ray technical note 2017)
    k = 90.0 + np.zeros((nt))

    # determine equilibrium arguments
    fargs = np.c_[tau, s, h, p, n, pp, k]
    rates = (fargs[1,:] - fargs[0,:])/deltat
    fd = np.dot(rates, coefficients_table(constituents, **kwargs))
    # convert to radians per second
    omega = 2.0*np.pi*fd/360.0
    return omega

def aliasing_period(
        constituents: list | np.ndarray,
        sampling: float | np.ndarray,
        **kwargs
    ):
    """
    Calculates the tidal aliasing for a repeat period

    Parameters
    ----------
    constituents: list
        tidal constituent IDs
    sampling: float
        sampling repeat period in seconds
    corrections: str, default 'OTIS'
        use nodal corrections from OTIS, FES or GOT models
    M1: str, default 'perth5'
        coefficients to use for M1 tides

                - ``'Doodson'``
                - ``'Ray'``
                - ``'perth5'``

    Returns
    -------
    period: np.ndarray
        tidal aliasing period in seconds
    """
    # get the angular frequency for tidal constituents
    omega = frequency(constituents, **kwargs)
    # convert to cycles per second
    f = omega/(2.0*np.pi)
    # calculate the sampling frequency
    fs = 1.0/sampling
    # calculate the aliasing period
    period = 1.0/pyTMD.math.aliasing(f, fs)
    # reutrn the aliasing period
    return period

def _arguments_table(**kwargs):
    """
    Arguments table for tidal constituents :cite:p:`Doodson:1921kt` :cite:p:`Doodson:1941td`

    Parameters
    ----------
    corrections: str, default 'OTIS'
        use arguments from OTIS, FES or GOT models

    Returns
    -------
    coef: np.ndarray
        Doodson coefficients (Cartwright numbers) for each constituent
    """
    # set default keyword arguments
    kwargs.setdefault('corrections', 'OTIS')

    # constituents array (not all are included in tidal program)
    cindex = ['sa', 'ssa', 'mm', 'msf', 'mf', 'mt', 'alpha1', '2q1', 'sigma1',
        'q1', 'rho1', 'o1', 'tau1', 'm1', 'chi1', 'pi1', 'p1', 's1', 'k1',
        'psi1', 'phi1', 'theta1', 'j1', 'oo1', '2n2', 'mu2', 'n2', 'nu2', 'm2a',
        'm2', 'm2b', 'lambda2', 'l2', 't2', 's2', 'r2', 'k2', 'eta2', 'mns2',
        '2sm2', 'm3', 'mk3', 's3', 'mn4', 'm4', 'ms4', 'mk4', 's4', 's5', 'm6',
        's6', 's7', 's8', 'm8', 'mks2', 'msqm', 'mtm', 'n4', 'eps2', 'z0']
    # modified Doodson coefficients for constituents
    # using 7 index variables: tau, s, h, p, n, pp, k
    # tau: mean lunar time
    # s: mean longitude of moon
    # h: mean longitude of sun
    # p: mean longitude of lunar perigee
    # n: mean longitude of ascending lunar node
    # pp: mean longitude of solar perigee
    # k: 90-degree phase
    coef = coefficients_table(cindex, **kwargs)
    # return the coefficient table
    return coef

def _minor_table(**kwargs):
    """
    Arguments table for minor tidal constituents :cite:p:`Doodson:1921kt` :cite:p:`Doodson:1941td`

    Returns
    -------
    coef: np.ndarray
        Doodson coefficients (Cartwright numbers) for each constituent
    """
    # modified Doodson coefficients for constituents
    # using 7 index variables: tau, s, h, p, n, pp, k
    # tau: mean lunar time
    # s: mean longitude of moon
    # h: mean longitude of sun
    # p: mean longitude of lunar perigee
    # n: mean longitude of ascending lunar node
    # pp: mean longitude of solar perigee
    # k: 90-degree phase
    minor = ['2q1', 'sigma1', 'rho1', 'm1b', 'm1a', 'chi1', 'pi1',
        'phi1', 'theta1', 'j1', 'oo1', '2n2', 'mu2', 'nu2', 'lambda2',
        'l2a', 'l2b', 't2', 'eps2', 'eta2']
    coef = coefficients_table(minor, **kwargs)
    # return the coefficient table
    return coef

def _constituent_parameters(c: str, **kwargs):
    """
    Loads parameters for a given tidal constituent :cite:p:`Egbert:2002ge`

    Parameters
    ----------
    c: str
        tidal constituent ID
    raise_error: bool, default False
        Raise exception if constituent is unsupported

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
    # default keyword arguments
    kwargs.setdefault('raise_error', False)
    # constituents array that are included in tidal program
    cindex = ['m2', 's2', 'k1', 'o1', 'n2', 'p1', 'k2', 'q1', '2n2', 'mu2',
        'nu2', 'l2', 't2', 'j1', 'm1', 'oo1', 'rho1', 'mf', 'mm', 'ssa',
        'm4', 'ms4', 'mn4', 'm6', 'm8', 'mk3', 's6', '2sm2', '2mk3',
        'msf', 'sa', 'mt', '2q1']
    # species type (spherical harmonic dependence of quadrupole potential)
    _species = np.array([2, 2, 1, 1, 2, 1, 2, 1, 2, 2, 2, 2, 2, 1, 1,
        1, 1, 0, 0, 0, 4, 4, 4, 6, 8, 3, 6, 2, 3, 0, 0, 0, 1])
    # Load Love numbers
    # alpha = correction factor for first order load tides
    _alpha = np.array([0.693, 0.693, 0.736, 0.695, 0.693, 0.706, 0.693,
        0.695, 0.693, 0.693, 0.693, 0.693, 0.693, 0.695, 0.695, 0.695, 0.695,
        0.693, 0.693, 0.693, 0.693, 0.693, 0.693, 0.693, 0.693, 0.693, 0.693,
        0.693, 0.693, 0.693, 0.693, 0.693, 0.693])
    # omega: angular frequency of constituent, in radians
    _omega = np.array([1.405189e-04, 1.454441e-04, 7.292117e-05, 6.759774e-05,
        1.378797e-04, 7.252295e-05, 1.458423e-04, 6.495854e-05, 1.352405e-04,
        1.355937e-04, 1.382329e-04, 1.431581e-04, 1.452450e-04, 7.556036e-05,
        7.025945e-05, 7.824458e-05, 6.531174e-05, 0.053234e-04, 0.026392e-04,
        0.003982e-04, 2.810377e-04, 2.859630e-04, 2.783984e-04, 4.215566e-04,
        5.620755e-04, 2.134402e-04, 4.363323e-04, 1.503693e-04, 2.081166e-04,
        4.925200e-06, 1.990970e-07, 7.962619e-06, 6.231934e-05])
    # Astronomical arguments (relative to t0 = 1 Jan 0:00 1992)
    # phases for each constituent are referred to the time when the phase of
    # the forcing for that constituent is zero on the Greenwich meridian
    _phase = np.array([1.731557546, 0.000000000, 0.173003674, 1.558553872,
        6.050721243, 6.110181633, 3.487600001, 5.877717569, 4.086699633,
        3.463115091, 5.427136701, 0.553986502, 0.050398470, 2.137025284,
        2.436575000, 1.929046130, 5.254133027, 1.756042456, 1.964021610,
        3.487600001, 3.463115091, 1.731557546, 1.499093481, 5.194672637,
        6.926230184, 1.904561220, 0.000000000, 4.551627762, 3.290111417,
        4.551627762, 6.232786837, 3.720064066, 3.91369596])
    # amplitudes of equilibrium tide in meters
    # _amplitude = np.array([0.242334,0.112743,0.141565,0.100661,0.046397,
    _amplitude = np.array([0.2441, 0.112743, 0.141565, 0.100661, 0.046397,
        0.046848, 0.030684, 0.019273, 0.006141, 0.007408, 0.008811, 0.006931,
        0.006608, 0.007915, 0.007915, 0.004338, 0.003661, 0.042041, 0.022191,
        0.019567, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.003681, 0.003104,
        0.008044, 0.002565])

    # map between input constituent and cindex
    j = [j for j,val in enumerate(cindex) if (val == c.lower())]
    # set the values for the constituent
    if j:
        amplitude, = _amplitude[j]
        phase, = _phase[j]
        omega, = _omega[j]
        alpha, = _alpha[j]
        species, = _species[j]
    elif kwargs['raise_error']:
        raise ValueError(f'Unsupported constituent {c}')
    else:
        amplitude = 0.0
        phase = 0.0
        omega = 0.0
        alpha = 0.0
        species = 0
    # return the values for the constituent
    return (amplitude, phase, omega, alpha, species)

def _love_numbers(
        omega: np.ndarray,
        model: str = 'PREM'
    ):
    """
    Compute the body tide Love/Shida numbers for a given
    frequency :cite:p:`Wahr:1981ea` :cite:p:`Wahr:1981if` :cite:p:`Mathews:1995go`

    Parameters
    ----------
    omega: np.ndarray
        angular frequency (radians per second)
    model: str, default 'PREM'
        Earth model to use for Love numbers

            - '1066A'
            - 'PEM-C'
            - 'C2'
            - 'PREM'

    Returns
    -------
    h2: float
        Degree-2 Love number of vertical displacement
    k2: float
        Degree-2 Love number of gravitational potential
    l2: float
        Degree-2 Love (Shida) number of horizontal displacement
    """
    # free core nutation frequencies (cycles per sidereal day) and
    # Love number parameters from Wahr (1981) table 6
    # and Mathews et al. (1995) table 3
    if (model == '1066A'):
        fcn = 1.0021714
        h0, h1 = np.array([6.03e-1, -2.46e-3])
        k0, k1 = np.array([2.98e-1, -1.23e-3])
        l0, l1 = np.array([8.42e-2, 7.81e-5])
    elif (model == 'PEM-C'):
        fcn = 1.0021771
        h0, h1 = np.array([6.02e-1, -2.46e-3])
        k0, k1 = np.array([2.98e-1, -1.24e-3])
        l0, l1 = np.array([8.39e-2, 7.69e-5])
    elif (model == 'C2'):
        fcn = 1.0021844
        h0, h1 = np.array([6.02e-1, -2.45e-3])
        k0, k1 = np.array([2.98e-1, -1.23e-3])
        l0, l1 = np.array([8.46e-2, 7.58e-5])
    elif (model == 'PREM'):
        fcn = 1.0023214
        h0, h1 = np.array([5.994e-1, -2.532e-3])
        k0, k1 = np.array([2.962e-1, -1.271e-3])
        l0, l1 = np.array([8.378e-2, 7.932e-5])
    else:
        raise ValueError(f'Unknown Earth model: {model}')
    # Love numbers for different frequency bands
    if (omega > 1e-4):
        # tides in the semi-diurnal band
        h2 = 0.609
        k2 = 0.302
        l2 = 0.0852
    elif (omega < 2e-5):
        # tides in the long period band
        h2 = 0.606
        k2 = 0.299
        l2 = 0.0840
    else:
        # use resonance formula for tides in the diurnal band
        # frequency of the o1 tides (radians/second)
        omega_o1, = frequency('o1')
        # frequency of free core nutation (radians/second)
        omega_fcn = fcn*7292115e-11
        # Love numbers for frequency using equation 4.18 of Wahr (1981)
        # (simplification to use only the free core nutation term)
        ratio = (omega - omega_o1)/(omega_fcn - omega)
        h2 = h0 + h1*ratio
        k2 = k0 + k1*ratio
        l2 = l0 + l1*ratio
    # return the Love numbers for frequency
    return (h2, k2, l2)

# Cartwright and Tayler (1971) table with 3rd-degree values
_ct1971_table_5 = get_data_path(['data','ct1971_tab5.txt'])
# Cartwright and Edden (1973) table with updated values
_ce1973_table_1 = get_data_path(['data','ce1973_tab1.txt'])

def _parse_tide_potential_table(table: str | pathlib.Path):
    """Parse tables of tide-generating potential from
    :cite:p:`Cartwright:1971iz` and :cite:p:`Cartwright:1973em`

    Parameters
    ----------
    table: str or pathlib.Path
        table of tide-generating potentials

    Returns
    -------
    CTE: float
        Cartwright-Tayler-Edden table values
    """
    # verify table path
    table = pathlib.Path(table).expanduser().absolute()
    with table.open(mode='r', encoding='utf8') as f:
        file_contents = f.readlines()
    # number of lines in the file
    file_lines = len(file_contents)
    # tau: coefficient for mean lunar time
    # s: coefficient for mean longitude of moon
    # h: coefficient for mean longitude of sun
    # p: coefficient for mean longitude of lunar perigee
    # n: coefficient for mean longitude of ascending lunar node
    # pp: coefficient for mean longitude of solar perigee
    # Hs1: amplitude for epoch span 1 (1861-09-21 to 1879-09-22)
    # Hs2: amplitude for epoch span 2 (1915-05-16 to 1933-05-22)
    # Hs3: amplitude for epoch span 2 (1951-05-23 to 1969-05-22)
    # DO: Doodson number for coefficient
    # Hs0: Doodson scaled amplitude for 1900
    names = ('tau','s','h','p','n','pp','Hs1','Hs2','Hs3','DO','Hs0')
    formats = ('i','i','i','i','i','i','f','f','f','U7','f')
    dtype = np.dtype({'names':names, 'formats':formats})
    CTE = np.zeros((file_lines), dtype=dtype)
    for i,line in enumerate(file_contents):
        # drop last column with values from Doodson (1921)
        CTE[i] = np.array(tuple(line.split()[:11]), dtype=dtype)
    # return the table values
    return CTE

def _to_doodson_number(coef: list | np.ndarray, **kwargs):
    """
    Converts Cartwright numbers into a Doodson number

    Parameters
    ----------
    coef: list or np.ndarray
        Doodson coefficients (Cartwright numbers) for constituent
    raise_error: bool, default True
        Raise exception if constituent is unsupported

    Returns
    -------
    DO: float or string
        Doodson number for constituent
    """
    # default keyword arguments
    kwargs.setdefault('raise_error', True)
    # assert length and verify array
    coef = np.array(coef[:6]).astype(int)
    # add 5 to values following Doodson convention (prevent negatives)
    coef[1:] += 5
    # check for unsupported constituents
    if (np.any(coef < 0) or np.any(coef > 10)) and kwargs['raise_error']:
        raise ValueError('Unsupported constituent')
    elif (np.any(coef < 0) or np.any(coef > 10)):
        return None
    elif np.any(coef == 10):
        # convert to string and replace 10 with X (Cartwright convention)
        DO = [str(v).replace('10','X') for v in coef]
        # convert to Doodson number
        return np.str_('{0}{1}{2}.{3}{4}{5}'.format(*DO))
    else:
        # convert to single number and round off floating point errors
        DO = np.sum([v*10**(2-o) for o,v in enumerate(coef)])
        return np.round(DO, decimals=3)

def _to_extended_doodson(coef: list | np.ndarray, **kwargs):
    """
    Converts Cartwright numbers into an UKHO Extended Doodson number

    Parameters
    ----------
    coef: list or np.ndarray
        Doodson coefficients (Cartwright numbers) for constituent

    Returns
    -------
    XDO: string
        Extended Doodson number for constituent
    """
    # assert length and verify array
    coef = np.array(coef).astype(int)
    # digits for UKHO Extended Doodson number
    # Z = 0
    # A - P = 1 to 15
    # R - Y = -8 to -1
    digits = 'RSTUVWXYZABCDEFGHIJKLMNOP'
    XDO = ''.join([digits[v+8] for v in coef])
    return np.str_(XDO)

def _from_doodson_number(DO: str | float | np.ndarray, **kwargs):
    """
    Converts Doodson numbers into Cartwright numbers

    Parameters
    ----------
    DO: str, float or np.ndarray
        Doodson number for constituent

    Returns
    -------
    coef: np.ndarray
        Doodson coefficients (Cartwright numbers) for constituent
    """
    # convert from Doodson number to Cartwright numbers
    coef = [c.replace('X', '10') for c in re.findall(r'\w', str(DO).zfill(7))]
    coef = np.array(coef, dtype=int)
    # remove 5 from values following Doodson convention
    coef[1:] -= 5
    return coef

def _from_extended_doodson(XDO: str | np.str_, **kwargs):
    """
    Converts UKHO Extended Doodson number into Cartwright numbers

    Parameters
    ----------
    XDO: string
        Extended Doodson number for constituent

    Returns
    -------
    coef: np.ndarray
        Doodson coefficients (Cartwright numbers) for constituent
    """
    # digits for UKHO Extended Doodson number
    # Z = 0
    # A - P = 1 to 15
    # R - Y = -8 to -1
    digits = 'RSTUVWXYZABCDEFGHIJKLMNOP'
    # convert from extended Doodson number to Cartwright numbers
    coef = np.array([(digits.index(c)-8) for c in str(XDO)], dtype=int)
    return coef
