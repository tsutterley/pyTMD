#!/usr/bin/env python
u"""
infer_minor_corrections.py (09/2017)
Return correction for 16 minor constituents based on Richard Ray's perth2 code

CALLING SEQUENCE:
	dh = infer_minor_corrections(c)

INPUTS:
	constituents: tidal constituent IDs
	zmajor: Complex HC for GIVEN constituents/points
	time: days relative to Jan 1, 1992 (48622mjd)

OUTPUT:
	dh: height from minor constituents

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python
		http://www.numpy.org
		http://www.scipy.org/NumPy_for_Matlab_Users

PROGRAM DEPENDENCIES:
	calc_astrol_longitudes.py: computes the basic astronomical mean longitudes

UPDATE HISTORY:
	Updated 09/2017: Rewritten in Python
"""
import numpy as np
from calc_astrol_longitudes import calc_astrol_longitudes

def infer_minor_corrections(time, zmajor, constituents):
	dtr = np.pi/180.0
	PP = 282.8
	#-- number of constituents
	npts,nc = np.shape(zmajor)
	#-- allocate for output elevation correction
	dh = np.zeros((npts))
	#-- convert time from days relative to Jan 1, 1992 to modified julian days
	time_mjd = 48622.0 + time
	cindex = ['q1','o1','p1','k1','n2','m2','s2','k2']
	#-- re-order zmaj to correspond to cindex
	z8 = np.zeros((npts,8),dtype=np.complex64)
	ni = 0
	for i,c in enumerate(cindex):
		j = [j for j,val in enumerate(constituents) if val == c]
		if j:
			j1, = j
			z8[:,i] = zmajor[:,j1]
			ni += 1

	if (ni < 6):
		raise Exception('Not enough constituents for inference')

	zmin = np.zeros((npts,18),dtype=np.complex64)
	zmin[:,0] = 0.263*z8[:,0] - 0.0252*z8[:,1]#-- 2Q1
	zmin[:,1] = 0.297*z8[:,0] - 0.0264*z8[:,1]#-- sigma1
	zmin[:,2] = 0.164*z8[:,0] + 0.0048*z8[:,1]#-- rho1 +
	zmin[:,3] = 0.0140*z8[:,1] + 0.0101*z8[:,3]#-- M1
	zmin[:,4] = 0.0389*z8[:,1] + 0.0282*z8[:,3]#-- M1
	zmin[:,5] = 0.0064*z8[:,1] + 0.0060*z8[:,3]#-- chi1
	zmin[:,6] = 0.0030*z8[:,1] + 0.0171*z8[:,3]#-- pi1
	zmin[:,7] = -0.0015*z8[:,1] + 0.0152*z8[:,3]#-- phi1
	zmin[:,8] = -0.0065*z8[:,1] + 0.0155*z8[:,3]#-- theta1
	zmin[:,9] = -0.0389*z8[:,1] + 0.0836*z8[:,3]#-- J1 +
	zmin[:,10] = -0.0431*z8[:,1] + 0.0613*z8[:,3]#-- OO1 +
	zmin[:,11] = 0.264*z8[:,4] - 0.0253*z8[:,5]#-- 2N2 +
	zmin[:,12] = 0.298*z8[:,4] - 0.0264*z8[:,5]#-- mu2 +
	zmin[:,13] = 0.165*z8[:,4] + 0.00487*z8[:,5]#-- nu2 +
	zmin[:,14] = 0.0040*z8[:,5] + 0.0074*z8[:,6]#-- lambda2
	zmin[:,15] = 0.0131*z8[:,5] + 0.0326*z8[:,6]#-- L2 +
	zmin[:,16] = 0.0033*z8[:,5] + 0.0082*z8[:,6]#-- L2 +
	zmin[:,17] = 0.0585*z8[:,6]#-- t2 +

	hour = (time % 1)*24.0
	t1 = 15.0*hour
	t2 = 30.0*hour
	S,H,P,omega = calc_astrol_longitudes(time_mjd)

	arg = np.zeros((npts,18))
	arg[:,0] = t1 - 4.0*S + H + 2.0*P - 90.#-- 2Q1
	arg[:,1] = t1 - 4.0*S + 3.0*H - 90.#-- sigma1
	arg[:,2] = t1 - 3.0*S + 3.0*H - P - 90.#-- rho1
	arg[:,3] = t1 - S + H - P + 90.#-- M1
	arg[:,4] = t1 - S + H + P + 90.#-- M1
	arg[:,5] = t1 - S + 3.0*H - P + 90.#-- chi1
	arg[:,6] = t1 - 2.0*H + PP - 90.#-- pi1
	arg[:,7] = t1 + 3.0*H + 90.#-- phi1
	arg[:,8] = t1 + S - H + P + 90.#-- theta1
	arg[:,9] = t1 + S + H - P + 90.#-- J1
	arg[:,10] = t1 + 2.0*S + H + 90.#-- OO1
	arg[:,11] = t2 - 4.0*S + 2.0*H + 2.0*P#-- 2N2
	arg[:,12] = t2 - 4.0*S + 4.0*H#-- mu2
	arg[:,13] = t2 - 3.0*S + 4.0*H - P#-- nu2
	arg[:,14] = t2 - S + P + 180.0#-- lambda2
	arg[:,15] = t2 - S + 2.0*H - P + 180.0#-- L2
	arg[:,16] = t2 - S + 2.0*H + P#-- L2
	arg[:,17] = t2 - H + PP#-- t2

	#-- determine nodal corrections f and u
	sinn = np.sin(omega*dtr)
	cosn = np.cos(omega*dtr)
	sin2n = np.sin(2.0*omega*dtr)
	cos2n = np.cos(2.0*omega*dtr)

	f = np.ones((npts,18))
	f[:,0] = np.sqrt((1.0 + 0.189*cosn - 0.0058*cos2n)**2 +
		(0.189*sinn - 0.0058*sin2n)**2)
	f[:,1] = f[:,0]
	f[:,2] = f[:,0]
	f[:,3] = np.sqrt((1.0 + 0.185*cosn)**2 + (0.185*sinn)**2)
	f[:,4] = np.sqrt((1.0 + 0.201*cosn)**2 + (0.201*sinn)**2)
	f[:,5] = np.sqrt((1.0 + 0.221*cosn)**2 + (0.221*sinn)**2)
	f[:,9] = np.sqrt((1.0 + 0.198*cosn)**2 + (0.198*sinn)**2)
	f[:,10] = np.sqrt((1.0 + 0.640*cosn + 0.134*cos2n)**2 +
		(0.640*sinn + 0.134*sin2n)**2)
	f[:,11] = np.sqrt((1.0 - 0.0373*cosn)**2 + (0.0373*sinn)**2)
	f[:,12] = f[:,11]
	f[:,13] = f[:,11]
	f[:,15] = f[:,11]
	f[:,16] = np.sqrt((1.0 + 0.441*cosn)**2 + (0.441*sinn)**2)

	u = np.zeros((npts,18))
	u[:,0] = np.arctan2(0.189*sinn - 0.0058*sin2n,
		1.0 + 0.189*cosn - 0.0058*sin2n)/dtr
	u[:,1] = u[:,0]
	u[:,2] = u[:,0]
	u[:,3] = np.arctan2( 0.185*sinn, 1.0 + 0.185*cosn)/dtr
	u[:,4] = np.arctan2(-0.201*sinn, 1.0 + 0.201*cosn)/dtr
	u[:,5] = np.arctan2(-0.221*sinn, 1.0 + 0.221*cosn)/dtr
	u[:,9] = np.arctan2(-0.198*sinn, 1.0 + 0.198*cosn)/dtr
	u[:,10] = np.arctan2(-0.640*sinn - 0.134*sin2n,
		1.0 + 0.640*cosn + 0.134*cos2n)/dtr
	u[:,11] = np.arctan2(-0.0373*sinn, 1.0 - 0.0373*cosn)/dtr
	u[:,12] = u[:,11]
	u[:,13] = u[:,11]
	u[:,15] = u[:,11]
	u[:,16] = np.arctan2(-0.441*sinn, 1.0 + 0.441*cosn)/dtr

	#-- sum over all tides
	for k in range(18):
		th = (arg[:,k] + u[:,k])*dtr
		dh += zmin.real[:,k]*f[:,k]*np.cos(th)-zmin.imag[:,k]*f[:,k]*np.sin(th)
	#-- return the inferred elevation
	return dh
