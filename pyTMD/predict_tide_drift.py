#!/usr/bin/env python
u"""
predict_tide_drift.py (09/2017)
Predict tidal elevation at multiple times and locations using harmonic constants

CALLING SEQUENCE:
	ht = predict_tide_drift(time,hc,con)

INPUTS:
	time: days relative to Jan 1, 1992 (48622mjd)
	hc: harmonic constant vector (complex)
	constituents: tidal constituent IDs

OUTPUT:
	ht: time series reconstructed using the nodal corrections

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python
		http://www.numpy.org
		http://www.scipy.org/NumPy_for_Matlab_Users

PROGRAM DEPENDENCIES:
	load_constituent.py: loads parameters for a given tidal constitent
	load_nodal_corrections.py: loads nodal corrections for tidal constituents

UPDATE HISTORY:
	Updated 09/2017: Rewritten in Python
"""
import numpy as np
from load_constituent import load_constituent
from load_nodal_corrections import load_nodal_corrections

def predict_tide_drift(time,hc,constituents):
	nt = len(time)
	#-- load the nodal corrections
	pu,pf,G = load_nodal_corrections(time + 48622.0, constituents)
	#-- allocate for output time series
	ht = np.zeros((nt))
	#-- for each constituent
	for k,c in enumerate(constituents):
		#-- load parameters for each constituent
		amp,ph,omega,alpha,species = load_constituent(c)
		#-- add component for constituent to output tidal elevation
		th = omega*time*86400.0 + ph + pu[:,k]
		ht += pf[:,k]*hc.real[:,k]*np.cos(th) - pf[:,k]*hc.imag[:,k]*np.sin(th)
	#-- return the tidal elevation
	return ht
