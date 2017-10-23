#!/usr/bin/env python
u"""
map_xy_tides.py
Original IDL program mapxy.pro written by Eric Rignot
Adapted by Tyler Sutterley (09/2017)

Converts from Polar Stereographic (X,Y) coordinates to geodetic latitude and
	longitude for the polar regions

CALLING SEQUENCE:
	ll = map_xy_tides(x, y, slat, sn, xlam,  FORMAT='dict')
	ln,lt = map_xy_tides(x, y, slat, sn, xlam,  FORMAT='tuple')
	ll = map_xy_tides(x, y, slat, sn, xlam, FORMAT='zip')

INPUTS:
	x: horizontal coordinate polar stereographic projection (easting) in km
	y: vertical coordinate polar stereographic projection (northing) in km
	slat: latitude_of_origin (N: 70, S: -71)
	sn: sign of latitude (N: +1, S: -1)
	xlam: longitude_of_center (N: 45, S: 0)

OUTPUTS:
	alon: longitude
	alat: latitude

OPTIONS:
	FORMAT: format of output coordinates
		'dict': python dictionary with keys 'lon' and 'lat'
		'tuple': python tuple with first variable alon and second variable alat
		'zip': aggregated coordinate pairs

NOTES:
	Snyder, J P (1982) Map Projections used by the U.S. Geological Survey
		Inverse formulas for the ellipsoid.  Geological Survey Bulletin 1532,
		U.S. Government Printing Office.
		See JPL Technical Memorandum 3349-85-101 for further details.
	Adams (1921) Latitude developments connected with geodesy and
		cartography with tables

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python (http://www.numpy.org)

UPDATE HISTORY:
	Forked 09/2017: using ecc2, slat, sn and xlam for tidal models
	Updated 08/2017: can enter slat, sn and xlam directly into HEM for custom
	Updated 05/2017: added comments and updated header text
		added option FORMAT to change the format of the output variables
	Updated 04/2015: added option for Bamber 5km DEM
	Updated 06/2014: updated comments for references and expanded
		Adams series solution to order 4
	Updated 02/2014: minor update to if statements
	Updated 06/2013: converted to python
"""

def map_xy_tides(x, y, slat, sn, xlam, FORMAT='tuple'):
	import numpy as np

	#-- semimajor axis of the ellipsoid
	rad_e = 6378.137#-- WGS84 [km]
	#-- ecc^2 value from TMD tidal conversion algorithm (ecc = 0.08181919)
	ecc2 = 6.694379852e-3
	# ecc2 = 0.00669437999015#-- WGS84 (ecc = 0.0818191908426215)
	ecc = np.sqrt(ecc2)

	#-- size of the input array
	sz = np.ndim(x)
	x = x*1.0
	y = y*1.0

	#-- distance rho
	rho = np.sqrt(x**2 + y**2)

	#-- if points are exactly on the pole:
	#-- make y slightly off pole to prevent division error
	if (sz == 0):
		y = 1e-20 if (rho == 0.0) else y
	else:
		#-- making y slightly off pole to prevent division error
		if np.count_nonzero(rho == 0.0):
			ind, = np.nonzero(rho == 0.0)
			y[ind] = 1e-20

	#-- m at standard latitude
	SL = sn*slat*np.pi/180.0
	cm = np.cos(SL) / np.sqrt(1.0 - ecc2*(np.sin(SL)**2))
	t = np.tan(np.pi/4.0 - SL/2.0) / ((1.0 - ecc*np.sin(SL)) / \
		(1.0 + ecc*np.sin(SL)))**(ecc/2.0)

	if ((sn*slat) == 90.0):
	    t = rho*np.sqrt((1.0+ecc)**(1.0+ecc)*(1.0-ecc)**(1.0-ecc))/(2.0*rad_e)
	else:
		t = rho*t/(rad_e*cm)

	#-- conformal latitude
	chi = (np.pi/2.0)-2.0*np.arctan(t)
	#-- inverse formula for calculating lat in terms in chi
	#-- Snyder equation 3-5 and Adams (1921) Page 85 phi-chi
	#-- series solution to avoid iterating to convergence
	a1 = (ecc2/2.0) + (5.0*ecc2**2/24.0) + (ecc2**3/12.0) + (13.0*ecc2**4/360.0)
	a2 = (7.0*ecc2**2/48.0) + (29.0*ecc2**3/240.0) + (811.0*ecc2**4/11520.0)
	a3 = (7.0*ecc2**3/120.0) + (81.0*ecc2**4/1120.0)
	a4 = (4279.0*ecc2**4/161280.0)
	alat = chi + a1*np.sin(2.0*chi) + a2*np.sin(4.0*chi) + a3*np.sin(6.0*chi) + \
		a4*np.sin(8.0*chi)
	alat = (sn*alat)*(180.0/np.pi)

	#-- arctan2 computes the arctangent of xpr/ypr with range (-pi,pi)
	#-- this portion was modified from the original in the tidal mapxy.m code
	alon = np.arctan2(-x,y)*(180.0/np.pi) + xlam

	#-- Fixing longitudes to be 0:360
	#-- input data is points
	if (sz == 0):
		if (alon < 0):
			alon = alon+360.
		elif (alon > 360.0):
			alon = alon-360.
	else: #-- input data is arrays
		if np.count_nonzero(alon < 0.0):
			ind, = np.nonzero(alon < 0.0)
			alon[ind]=alon[ind]+360.0
		if np.count_nonzero(alon > 360.0):
			ind, = np.nonzero(alon > 360.0)
			alon[ind]=alon[ind]-360.0

	#-- return coordinates in output format (default python dictionary)
	if (FORMAT == 'dict'):
		return dict(lon=alon,lat=alat)
	elif (FORMAT == 'tuple'):
		return (alon,alat)
	elif (FORMAT == 'zip'):
		return zip(alon,alat)
