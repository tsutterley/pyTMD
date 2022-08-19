#!/usr/bin/env python
u"""
compute_equilibrium_tide.py (04/2022)
Calculates the long-period equilibrium ocean tides
Uses the summation of fifteen tidal spectral lines from the
    Cartwright-Tayler-Edden tables to compute the long-period tides

INPUTS:
    t: days relative to Jan 1, 1992 (48622 MJD)
    lat: latitudes in degrees

OUTPUTS:
    lpet: long-period equilibrium tide in meters

REFERENCES:
    Cartwright & Tayler, Geophys. J. R.A.S., 23, 45, 1971.
    Cartwright & Edden, Geophys. J. R.A.S., 33, 253, 1973.

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 11/2020: check sizes of inputs to check if using grid or drift data
        drift: drift buoys or satellite/airborne altimetry (time per data point)
        grid: spatial grids or images (single time for all data points)
    Updated 08/2020: use Load love number value of 0.693 for gamma_2
    Written 08/2020
"""
import numpy as np

def compute_equilibrium_tide(t, lat):
    """
    Parameters
    ----------
    t: float
        time (days relative to January 1, 1992)
    lat: float
        latitudes (degrees)

    Returns
    -------
    lpet: float
        long-period equilibrium tide in meters

    References
    ----------
    .. [1] Cartwright & Tayler, Geophys. J. R.A.S., 23, 45, 1971.
    .. [2] Cartwright & Edden, Geophys. J. R.A.S., 33, 253, 1973.
    """
    # longitude of moon
    # longitude of sun
    # longitude of lunar perigee
    # longitude of ascending lunar node
    PHC = np.array([290.21,280.12,274.35,343.51])
    DPD = np.array([13.1763965,0.9856473,0.1114041,0.0529539])

    # number of input points
    nt = len(np.atleast_1d(t))
    nlat = len(np.atleast_1d(lat))
    # compute 4 principal mean longitudes in radians at delta time (SHPN)
    SHPN = np.zeros((4,nt))
    for N in range(4):
        # convert time from days relative to 1992-01-01 to 1987-01-01
        ANGLE = PHC[N] + (t + 1826.0)*DPD[N]
        SHPN[N,:] = np.pi*np.mod(ANGLE, 360.0)/180.0

    # assemble long-period tide potential from 15 CTE terms greater than 1 mm
    # nodal term is included but not the constant term.
    PH = np.zeros((nt))
    Z = np.zeros((nt))
    Z += 2.79*np.cos(SHPN[3,:]) - 0.49*np.cos(SHPN[1,:] - \
        283.0*np.pi/180.0) - 3.10*np.cos(2.0*SHPN[1,:])
    PH += SHPN[0,:]
    Z += -0.67*np.cos(PH - 2.0*SHPN[1,:] + SHPN[2,:]) - \
        (3.52 - 0.46*np.cos(SHPN[3,:]))*np.cos(PH - SHPN[2,:])
    PH += SHPN[0,:]
    Z += - 6.66*np.cos(PH) - 2.76*np.cos(PH + SHPN[3,:]) - \
        0.26 * np.cos(PH + 2.*SHPN[3,:]) - 0.58 * np.cos(PH - 2.*SHPN[1,:]) - \
        0.29 * np.cos(PH - 2.*SHPN[2,:])
    PH += SHPN[0,:]
    Z += - 1.27*np.cos(PH - SHPN[2,:]) - \
        0.53*np.cos(PH - SHPN[2,:] + SHPN[3,:]) - \
        0.24*np.cos(PH - 2.0*SHPN[1,:] + SHPN[2,:])

    # Multiply by gamma_2 * normalization * P20(lat)
    gamma_2 = 0.693
    P20 = 0.5*(3.0*np.sin(lat*np.pi/180.0)**2 - 1.0)
    # calculate long-period equilibrium tide and convert to meters
    if (nlat != nt):
        lpet = gamma_2*np.sqrt((4.0+1.0)/(4.0*np.pi))*np.outer(P20,Z/100.0)
    else:
        lpet = gamma_2*np.sqrt((4.0+1.0)/(4.0*np.pi))*P20*(Z/100.0)
    # return the long-period equilibrium tides
    return lpet
