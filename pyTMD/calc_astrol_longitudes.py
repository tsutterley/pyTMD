#!/usr/bin/env python
u"""
calc_astrol_longitudes.py (09/2017)
Modification of ASTROL fortran subroutine by Richard Ray 03/1999

Computes the basic astronomical mean longitudes: s, h, p, N
Formulae are for the period 1990--2010 and were derived by David Cartwright
Note N is not N', i.e. N is decreasing with time.

CALLING SEQUENCE:
	s,h,p,N = calc_astrol_longitudes(time)

INPUTS:
	time: modified julian day of input date

OUTPUTS:
	s: mean longitude of moon (degrees)
	h: mean longitude of sun (degrees)
	p: mean longitude of lunar perigee (degrees)
	N: mean longitude of ascending lunar node (degrees)

OPTIONS:
	MEEUS: use additional coefficients from Meeus Astronomical Algorithms

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python
		http://www.numpy.org
		http://www.scipy.org/NumPy_for_Matlab_Users

UPDATE HISTORY:
	Updated 09/2017: added option MEEUS to use additional coefficients
		from Meeus Astronomical Algorithms to calculate mean longitudes
	Updated 09/2017: Rewritten in Python
	Rewritten in Matlab by Lana Erofeeva 2003
	Written by Richard Ray 12/1990
"""
import numpy as np

def polynomial_sum(coefficients, time):
	return sum([c * (time ** i) for i,c in enumerate(coefficients)])

def calc_astrol_longitudes(time, MEEUS=False):
	circle = 360.0
	if MEEUS:
		#-- convert from MJD to days relative to 2000-01-01
		T = time - 51544.5
		#-- mean longitude of moon
		lunar_longitude = [218.3164591, 13.17639647754579, -9.9454632e-13,
			3.8086292e-20, -8.6184958e-27]
		s = polynomial_sum(lunar_longitude,T)
		#-- mean longitude of sun
		solar_longitude = [280.46645, 0.985647360164271, 2.2727347e-13]
		h = polynomial_sum(solar_longitude,T)
		#-- mean longitude of lunar perigee
		lunar_perigee = [83.3532430, 0.11140352391786447, -7.7385418e-12,
			-2.5636086e-19, 2.95738836e-26]
		p = polynomial_sum(lunar_perigee,T)
		#-- mean longitude of ascending lunar node
		lunar_node = [125.0445550, -0.052953762762491446, 1.55628359e-12,
			4.390675353e-20, -9.26940435e-27]
		N = polynomial_sum(lunar_node,T)
	else:
		#-- convert from MJD to days relative to 2000-01-01
		T = time - 51544.4993
		#-- mean longitude of moon
		s = 218.3164 + 13.17639648 * T
		#-- mean longitude of sun
		h = 280.4661 + 0.98564736 * T
		#-- mean longitude of lunar perigee
		p =  83.3535 + 0.11140353 * T
		#-- mean longitude of ascending lunar node
		N = 125.0445 - 0.05295377 * T

	#-- take the modulus of each
	s = s % circle
	h = h % circle
	p = p % circle
	N = N % circle

	#-- return as tuple
	return (s,h,p,N)
