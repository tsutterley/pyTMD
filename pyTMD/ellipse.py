#!/usr/bin/env python
u"""
ellipse.py
Written by Tyler Sutterley (01/2024)
Expresses the amplitudes and phases for the u and v components in terms of
    four ellipse parameters using Foreman's formula

CALLING SEQUENCE:
    umajor,uminor,uincl,uphase = pyTMD.ellipse.ellipse(u,v)

INPUTS:
    u: zonal current (EW)
    v: meridional current (NS)

OUTPUTS:
    umajor: amplitude of the semimajor semi-axis
    uminor: amplitude of the semiminor semi-axis
    uincl: angle of inclination of the northern semimajor semi-axis
    uphase: phase lag of the maximum current behind the maximum tidal potential
        of the individual constituent

REFERENCE:
    M. G. G. Foreman and R. F. Henry, "The harmonic analysis of tidal model time
        series", Advances in Water Resources, 12(3), 109-120, (1989).
        https://doi.org/10.1016/0309-1708(89)90017-1

UPDATE HISTORY:
    Updated 01/2024: added inverse function to get currents from parameters
        use complex algebra to calculate tidal ellipse parameters
    Updated 09/2023: renamed to ellipse.py (from tidal_ellipse.py)
    Updated 03/2023: add basic variable typing to function inputs
    Updated 04/2022: updated docstrings to numpy documentation format
    Written 07/2020
"""
from __future__ import annotations

import numpy as np

def ellipse(u: np.ndarray, v: np.ndarray):
    """
    Expresses the amplitudes and phases for the u and v components in terms of
    four ellipse parameters using Foreman's formula [1]_

    Parameters
    ----------
    u: np.ndarray
        zonal current (EW)
    v: np.ndarray
        meridional current (NS)

    Returns
    -------
    umajor: np.ndarray
        amplitude of the semimajor semi-axis
    uminor: np.ndarray
        amplitude of the semiminor semi-axis
    uincl: np.ndarray
        angle of inclination of the northern semimajor semi-axis
    uphase: np.ndarray
        phase lag of the maximum current behind the maximum tidal potential
        of the individual constituent

    References
    ----------
    .. [1] M. G. G. Foreman and R. F. Henry, "The harmonic analysis of tidal
        model time series," *Advances in Water Resources*, 12(3), 109--120,
        (1989). `doi: 10.1016/0309-1708(89)90017-1
        <https://doi.org/10.1016/0309-1708(89)90017-1>`_
    """
    # validate inputs
    u = np.atleast_1d(u)
    v = np.atleast_1d(v)
    # wp, wm: complex radius of positively and negatively rotating vectors
    wp = (u + 1j*v)/2.0
    wm = np.conj(u - 1j*v)/2.0
    # ap, am: amplitudes of positively and negatively rotating vectors
    ap = np.abs(wp)
    am = np.abs(wm)
    # ep, em: phases of positively and negatively rotating vectors
    ep = np.angle(wp, deg=True)
    em = np.angle(wm, deg=True)
    # determine the amplitudes of the semimajor and semiminor axes
    # using Foreman's formula
    umajor = (ap + am)
    uminor = (ap - am)
    # determine the inclination and phase using Foreman's formula
    uincl = (em + ep)/2.0
    uphase = (em - ep)/2.0
    # adjust orientation of ellipse
    k = (uincl//180.0)
    uincl -= 180.0*k
    uphase += 180.0*k
    uphase = np.mod(uphase, 360.0)
    # return values
    return (umajor, uminor, uincl, uphase)

def inverse(
        umajor: np.ndarray,
        uminor: np.ndarray,
        uincl: np.ndarray,
        uphase: np.ndarray
    ):
    """
    Calculates currents u, v using the four tidal ellipse
    parameters from Foreman's formula [1]_

    Parameters
    ----------
    umajor: np.ndarray
        amplitude of the semimajor semi-axis
    uminor: np.ndarray
        amplitude of the semiminor semi-axis
    uincl: np.ndarray
        angle of inclination of the northern semimajor semi-axis
    uphase: np.ndarray
        phase lag of the maximum current behind the maximum tidal potential
        of the individual constituent

    Returns
    -------
    u: np.ndarray
        zonal current (EW)
    v: np.ndarray
        meridional current (NS)

    References
    ----------
    .. [1] M. G. G. Foreman and R. F. Henry, "The harmonic analysis of tidal
        model time series," *Advances in Water Resources*, 12(3), 109--120,
        (1989). `doi: 10.1016/0309-1708(89)90017-1
        <https://doi.org/10.1016/0309-1708(89)90017-1>`_
    """
    # validate inputs
    umajor = np.atleast_1d(umajor)
    uminor = np.atleast_1d(uminor)
    # convert inclination and phase to radians
    uincl = np.atleast_1d(uincl)*np.pi/180.0
    uphase = np.atleast_1d(uphase)*np.pi/180.0
    # ep, em: phases of positively and negatively rotating vectors
    ep = (uincl - uphase)
    em = (uincl + uphase)
    # ap, am: amplitudes of positively and negatively rotating vectors
    ap = (umajor + uminor)/2.0
    am = (umajor - uminor)/2.0
    # wp, wm: complex radius of positively and negatively rotating vectors
    wp = ap * np.exp(1j*ep)
    wm = am * np.exp(1j*em)
    # calculate complex currents
    u = wp + np.conj(wm)
    v = -1j*(wp - np.conj(wm))
    # return values
    return (u, v)
