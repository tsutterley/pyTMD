#!/usr/bin/env python
u"""
map_ll_tides.py
Original IDL program mapll.pro written by Eric Rignot
Adapted by Tyler Sutterley (09/2017)

Converts from Polar Stereographic (X,Y) coordinates to geodetic latitude and
    longitude for the polar regions.

CALLING SEQUENCE:
    xy = map_ll_tides(alon, alat, HEM='S', FORMAT='dict')
    x,y = map_ll_tides(alon, alat, HEM='S', FORMAT='tuple')
    xy = map_ll_tides(alon, alat, HEM='S', FORMAT='zip')

INPUTS:
    alon: longitude
    alat: latitude
    slat: latitude_of_origin (N: 70, S: -71)
    sn: sign of latitude (N: +1, S: -1)
    xlam: longitude_of_center (N: 45, S: 0)

OUTPUTS:
    x: horizontal coordinate polar stereographic projection (easting) in km
    y: vertical coordinate polar stereographic projection (northing) in km

OPTIONS:
    FORMAT: format of output coordinates
        'dict': python dictionary with keys 'x' and 'y'
        'tuple': python tuple with first variable x and second variable y
        'zip': aggregated coordinate pairs

NOTES:
    Snyder, J P (1982) Map Projections used by the U.S. Geological Survey
        Forward formulas for the ellipsoid.  Geological Survey Bulletin 1532,
        U.S. Government Printing Office.
        See JPL Technical Memorandum 3349-85-101 for further details.

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python (http://www.numpy.org)

UPDATE HISTORY:
    Forked 09/2017: using ecc2, slat, sn and xlam for tidal models
    Updated 08/2017: can enter slat, sn and xlam directly into HEM for custom
    Updated 05/2017: added comments and updated header text
        added option FORMAT to change the format of the output variables
    Updated 04/2015: added option for Bamber 5km DEM
    Updated 06/2014: updated comments
    Updated 02/2014: minor update to if statements
    Updated 05/2013: converted to python
"""

def map_ll_tides(alon, alat, slat, sn, xlam, FORMAT='tuple'):
    import numpy as np

    #-- semimajor axis of the ellipsoid
    rad_e = 6378.137#-- WGS84 [km]
    #-- ecc^2 value from TMD tidal conversion algorithm (ecc = 0.08181919)
    ecc2 = 6.694379852e-3
    # ecc2 = 0.00669437999015#-- WGS84 (ecc = 0.0818191908426215)
    ecc = np.sqrt(ecc2)

    #-- convert signs and rotate
    rlat = sn*alat*np.pi/180.0
    rlon = -(alon-xlam)*np.pi/180.0

    t = np.tan(np.pi/4.0 - rlat/2.0) / ((1.0 - ecc*np.sin(rlat)) / \
        (1.0 + ecc*np.sin(rlat)))**(ecc/2.0)

    if ((sn*slat) == 90.0):
        rho = 2.0*rad_e/((1.0+ecc)**(1.0+ecc)*(1.0-ecc)**(1.0-ecc))**(ecc/2.0)
    else:
        SL = sn*slat*np.pi/180.0
        tc = np.tan(np.pi/4.0 - SL/2.0) / ((1.0 - ecc*np.sin(SL)) /
            (1.0 + ecc*np.sin(SL)))**(ecc/2.0)
        #-- m at standard latitude
        cm = np.cos(SL) / np.sqrt(1.0-ecc2*(np.sin(SL)**2.0))
        rho = rad_e*cm/tc

    #-- polar stereographic x and y
    x = -rho*t*np.sin(rlon)
    y = rho*t*np.cos(rlon)

    #-- return coordinates in output format (default python dictionary)
    if (FORMAT == 'dict'):
        return dict(x=x,y=y)
    elif (FORMAT == 'tuple'):
        return (x,y)
    elif (FORMAT == 'zip'):
        return zip(x,y)
