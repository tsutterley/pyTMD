#!/usr/bin/env python
u"""
tidal_ellipse.py (03/2023)
Expresses the amplitudes and phases for the u and v components in terms of
    four ellipse parameters using Foreman's formula

CALLING SEQUENCE:
    umajor,uminor,uincl,uphase = tidal_ellipse(u,v)

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
    Updated 03/2023: add basic variable typing to function inputs
    Updated 04/2022: updated docstrings to numpy documentation format
    Written 07/2020
"""
from __future__ import annotations

import numpy as np

def tidal_ellipse(u: np.ndarray, v: np.ndarray):
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
    umajor: float
        amplitude of the semimajor semi-axis
    uminor: float
        amplitude of the semiminor semi-axis
    uincl: float
        angle of inclination of the northern semimajor semi-axis
    uphase: float
        phase lag of the maximum current behind the maximum tidal potential
        of the individual constituent

    References
    ----------
    .. [1] M. G. G. Foreman and R. F. Henry, "The harmonic analysis of tidal
        model time series," *Advances in Water Resources*, 12(3), 109--120,
        (1989). `doi: 10.1016/0309-1708(89)90017-1
        <https://doi.org/10.1016/0309-1708(89)90017-1>`_
    """
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
