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
    # change to polar coordinates
    t1p = u.real - v.imag
    t2p = v.real + u.imag
    t1m = u.real + v.imag
    t2m = v.real - u.imag
    # ap, am: amplitudes of positively and negatively rotating vectors
    ap = np.sqrt(t1p**2 + t2p**2)/2.0
    am = np.sqrt(t1m**2 + t2m**2)/2.0
    # ep, em: phases of positively and negatively rotating vectors
    ep = 180.0*np.arctan2(t2p, t1p)/np.pi
    ep[ep < 0.0] += 360.0
    em = 180.0*np.arctan2(t2m, t1m)/np.pi
    em[em < 0.0] += 360.0
    # determine the amplitudes of the semimajor and semiminor axes
    # using Foreman's formula
    umajor = (ap + am)
    uminor = (ap - am)
    # determine the inclination and phase using Foreman's formula
    uincl = 0.5 * (em + ep)
    uincl[uincl > 180.0] -= 180.0
    uphase = -0.5*(ep - em)
    uphase[uphase < 0.0] += 360.0
    uphase[uphase >= 360.0] -= 360.0
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
    ep[ep < 0] += 2*np.pi
    ep[ep > 2*np.pi] -= 2*np.pi
    em = (uincl + uphase)
    em[em < 0] += 2*np.pi
    em[em > 2*np.pi] -= 2*np.pi
    # ap, am: amplitudes of positively and negatively rotating vectors
    ap = (umajor + uminor)/2.0
    am = (umajor - uminor)/2.0
    # currents in polar coordinates
    rp = np.tan(ep)
    rm = np.tan(em)
    t1p = 2.0*ap/np.sqrt(1.0 + rp**2)
    t2p = rp*t1p
    t1m = 2.0*am/np.sqrt(1.0 + rm**2)
    t2m = rm*t1m
    # adjust quadrants
    i1, = np.nonzero(np.arctan(rm) < 0)
    t1m[i1] *= -1.0
    t2m[i1] *= -1.0
    i2, = np.nonzero(np.isclose(np.arctan(rp) + np.pi, ep))
    t1p[i2] *= -1.0
    t2p[i2] *= -1.0
    # calculate currents
    u = (t1m + t1p)/2.0 + 1j*(t2p - t2m)/2.0
    v = (t2m + t2p)/2.0 + 1j*(t1m - t1p)/2.0
    # return values
    return (u, v)
