#!/usr/bin/env python
u"""
convert_xy_ll.py (09/2017)
Wrapper function to convert lat/lon points to and from projected coordinates

CALLING SEQUENCE:
	x,y = convert_xy_ll(lon,lat,PROJ,'F')
	lon,lat = convert_xy_ll(x,y,PROJ,'B')

INPUTS:
	i1: longitude ('F') or projection easting x ('B')
	i2: latitude ('F') or projection northing y ('B')
	PROJ: spatial reference system code for coordinate transformations
		https://www.epsg-registry.org/
	BF: backwards ('B') or forward ('F') translations

OUTPUTS:
	o1: projection easting x ('F') or longitude ('B')
	o2: projection northing y ('F') or latitude ('B')

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python
		http://www.numpy.org
		http://www.scipy.org/NumPy_for_Matlab_Users

PROGRAM DEPENDENCIES:
	map_ll_tides.py: converts from latitude/longitude to polar stereographic
	map_xy_tides.py: converts from polar stereographic to latitude/longitude

UPDATE HISTORY:
	Written 09/2017
"""
import numpy as np
from pyTMD.map_ll_tides import map_ll_tides
from pyTMD.map_xy_tides import map_xy_tides

def convert_xy_ll(i1,i2,PROJ,BF):
	#-- python dictionary with conversion functions
	conversion_functions = {}
	conversion_functions['3031'] = xy_ll_EPSG3031
	conversion_functions['CATS2008'] = xy_ll_CATS2008
	conversion_functions['3976'] = xy_ll_EPSG3976
	conversion_functions['3996'] = xy_ll_EPSG3996
	conversion_functions['4326'] = pass_values
	#-- check that PROJ for conversion was entered correctly
	if PROJ not in conversion_functions.keys():
		raise Exception('PROJ:{0} conversion function not found'.format(PROJ))
	#-- run conversion program and return values
	o1,o2 = conversion_functions[PROJ](i1,i2,BF)
	return (o1,o2)

#-- wrapper function for models in EPSG 3031 (Antarctic Polar Stereographic)
def xy_ll_EPSG3031(i1,i2,BF):
	#-- convert lat/lon to Polar-Stereographic x/y
	if (BF.upper() == 'F'):
		o1,o2 = map_ll_tides(i1,i2,-71.0,-1.0,0.0)
	#-- convert Polar-Stereographic x/y to lat/lon
	elif (BF.upper() == 'B'):
		o1,o2 = map_xy_tides(i1,i2,-71.0,-1.0,0.0)
	#-- return the output variables
	return (o1,o2)

#-- wrapper function for CATS2008 tide models
def xy_ll_CATS2008(i1,i2,BF):
	#-- convert lat/lon to Polar-Stereographic x/y
	if (BF.upper() == 'F'):
		o1,o2 = map_ll_tides(i1,i2,-71.0,-1.0,-70.0)
	#-- convert Polar-Stereographic x/y to lat/lon
	elif (BF.upper() == 'B'):
		o1,o2 = map_xy_tides(i1,i2,-71.0,-1.0,-70.0)
	#-- return the output variables
	return (o1,o2)

#-- wrapper function for models in EPSG 3976 (NSIDC Sea Ice Stereographic South)
def xy_ll_EPSG3976(i1,i2,BF):
	#-- convert lat/lon to Polar-Stereographic x/y
	if (BF.upper() == 'F'):
		o1,o2 = map_ll_tides(i1,i2,-70.0,-1.0,0.0)
	#-- convert Polar-Stereographic x/y to lat/lon
	elif (BF.upper() == 'B'):
		o1,o2 = map_xy_tides(i1,i2,-70.0,-1.0,0.0)
	#-- return the output variables
	return (o1,o2)

#-- wrapper function for models in EPSG 3996 (IBCAO Polar Stereographic North)
def xy_ll_EPSG3996(i1,i2,BF):
	#-- convert lat/lon to Polar-Stereographic x/y
	if (BF.upper() == 'F'):
		# o1,o2 = map_ll_tides(i1,i2,75.0,1.0,0.0)
		o1 = (90.0-i2)*111.7*np.cos(i1/180.0*np.pi)
		o2 = (90.0-i2)*111.7*np.sin(i1/180.0*np.pi)
	#-- convert Polar-Stereographic x/y to lat/lon
	elif (BF.upper() == 'B'):
		# o1,o2 = map_xy_tides(i1,i2,75.0,1.0,0.0)
		o1=90.0-np.sqrt(i1**2+i2**2)/111.7
		o2=np.arctan2(i2,i1)*180/np.pi
		ii,=np.nonzero(o1<0)
		o1[ii] += 360.0
	#-- return the output variables
	return (o1,o2)

#-- wrapper function to pass lat/lon values
def pass_values(i1,i2,BF):
	return (i1,i2)
