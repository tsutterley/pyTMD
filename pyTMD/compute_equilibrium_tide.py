"""
compute_equilibrium_tide.py (08/2020)
Calculates the long-period equilibrium ocean tides
Fifteen tidal spectral lines from the Cartwright-Tayler-Edden
    tables are summed over to compute the long-period tides

INPUTS:
    t: days relative to Jan 1, 1992 (48622 MJD)
    lat: latitudes in degrees

OUTPUTS:
    lpet: long-period equilibrium tide in meters

REFERENCES:
    Cartwright & Tayler, Geophys. J. R.A.S., 23, 45, 1971.
    Cartwright & Edden, Geophys. J. R.A.S., 33, 253, 1973.
"""
import numpy as np

def compute_equilibrium_tide(t, lat):
    #-- longitude of moon
    #-- longitude of sun
    #-- longitude of lunar perigee
    #-- longitude of ascending lunar node
    PHC = np.array([290.21,280.12,274.35,343.51])
    DPD = np.array([13.1763965,0.9856473,0.1114041,0.0529539])

    #-- convert time from days relative to 1992-01-01 to 1987-01-01
    #-- Compute 4 principal mean longitudes in radians at delta time
    npts = len(t)
    SHPN = np.zeros((4,npts))
    for N in range(4):
        ANGLE = PHC[N] + (t + 1826.0)*DPD[N]
        SHPN[N,:] = np.pi*np.mod(ANGLE, 360.0)/180.0

    #-- Assemble long-period tide potential from 15 CTE terms > 1 mm.
    #-- Nodal term is included but not the constant term.
    PH = np.zeros((npts))
    Z = np.zeros((npts))
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

    #-- Multiply by gamma_2 * sqrt(5/4 pi) * P20(lat)
    lpet = 0.437*Z*(1.5*np.sin(lat*np.pi/180.0)**2 - 0.5)/100.0
    return lpet
