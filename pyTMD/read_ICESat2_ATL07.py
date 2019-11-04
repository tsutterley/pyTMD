#!/usr/bin/env python
u"""
read_ICESat2_ATL07.py (11/2019)
Read ICESat-2 ATL07 (Sea Ice Height) data files

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python
		http://www.numpy.org
		http://www.scipy.org/NumPy_for_Matlab_Users
	h5py: Python interface for Hierarchal Data Format 5 (HDF5)
		http://h5py.org

UPDATE HISTORY:
	Written 11/2019
"""
from __future__ import print_function

import os
import re
import h5py
import numpy as np

#-- PURPOSE: read ICESat-2 ATL07 HDF5 data files
def read_HDF5_ATL07(FILENAME, ATTRIBUTES=False, VERBOSE=False):
	#-- Open the HDF5 file for reading
	fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

	#-- Output HDF5 file information
	if VERBOSE:
		print(fileID.filename)
		print(list(fileID.keys()))

	#-- allocate python dictionaries for ICESat-2 ATL07 variables and attributes
	IS2_atl07_mds = {}
	IS2_atl07_attrs = {} if ATTRIBUTES else None

	#-- read each input beam within the file
	IS2_atl07_beams = []
	for gtx in [k for k in fileID.keys() if bool(re.match(r'gt\d[lr]',k))]:
		#-- check if subsetted beam contains sea ice data
		try:
			fileID[gtx]['sea_ice_segments']['height_segment_id']
		except KeyError:
			pass
		else:
			IS2_atl07_beams.append(gtx)

	#-- read each input beam within the file
	for gtx in IS2_atl07_beams:
		IS2_atl07_mds[gtx] = {}
		IS2_atl07_mds[gtx]['sea_ice_segments'] = {}
		IS2_atl07_mds[gtx]['sea_ice_segments']['geolocation'] = {}
		IS2_atl07_mds[gtx]['sea_ice_segments']['geophysical'] = {}
		IS2_atl07_mds[gtx]['sea_ice_segments']['heights'] = {}
		IS2_atl07_mds[gtx]['sea_ice_segments']['stats'] = {}
		#-- get each HDF5 variable
		#-- ICESat-2 sea_ice_segments Group
		for key,val in fileID[gtx]['sea_ice_segments'].items():
			if isinstance(val, h5py.Dataset):
				IS2_atl07_mds[gtx]['sea_ice_segments'][key] = val[:]
			elif isinstance(val, h5py.Group):
				for k,v in val.items():
					IS2_atl07_mds[gtx]['sea_ice_segments'][key][k] = v[:]

		#-- Getting attributes of included variables
		if ATTRIBUTES:
			#-- Getting attributes of ICESat-2 ATL07 beam variables
			IS2_atl07_attrs[gtx] = {}
			IS2_atl07_attrs[gtx]['sea_ice_segments'] = {}
			IS2_atl07_attrs[gtx]['sea_ice_segments']['geolocation'] = {}
			IS2_atl07_attrs[gtx]['sea_ice_segments']['geophysical'] = {}
			IS2_atl07_attrs[gtx]['sea_ice_segments']['heights'] = {}
			IS2_atl07_attrs[gtx]['sea_ice_segments']['stats'] = {}
			#-- Global Group Attributes for ATL07 beam
			for att_name,att_val in fileID[gtx].attrs.items():
				IS2_atl07_attrs[gtx][att_name] = att_val
			for key,val in fileID[gtx]['sea_ice_segments'].items():
				IS2_atl07_attrs[gtx]['sea_ice_segments'][key] = {}
				for att_name,att_val in val.attrs.items():
					IS2_atl07_attrs[gtx]['sea_ice_segments'][key][att_name] = att_val
				if isinstance(val, h5py.Group):
					for k,v in val.items():
						IS2_atl07_attrs[gtx]['sea_ice_segments'][key][k] = {}
						for att_name,att_val in v.attrs.items():
							IS2_atl07_attrs[gtx]['sea_ice_segments'][key][k][att_name] = att_val

	#-- ICESat-2 orbit_info Group
	IS2_atl07_mds['orbit_info'] = {}
	for key,val in fileID['orbit_info'].items():
		IS2_atl07_mds['orbit_info'][key] = val[:]
	#-- ICESat-2 quality_assessment Group
	IS2_atl07_mds['quality_assessment'] = {}
	for key,val in fileID['quality_assessment'].items():
		if isinstance(val, h5py.Dataset):
			IS2_atl07_mds['quality_assessment'][key] = val[:]
		elif isinstance(val, h5py.Group):
			IS2_atl07_mds['quality_assessment'][key] = {}
			for k,v in val.items():
				IS2_atl07_mds['quality_assessment'][key][k] = v[:]

	#-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
	#-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01:T00:00:00Z UTC)
	#-- Add this value to delta time parameters to compute full gps_seconds
	#-- could alternatively use the Julian day of the ATLAS SDP epoch: 2458119.5
	#-- and add leap seconds since 2018-01-01:T00:00:00Z UTC (ATLAS SDP epoch)
	IS2_atl07_mds['ancillary_data'] = {}
	IS2_atl07_attrs['ancillary_data'] = {} if ATTRIBUTES else None
	for key in ['atlas_sdp_gps_epoch']:
		#-- get each HDF5 variable
		IS2_atl07_mds['ancillary_data'][key] = fileID['ancillary_data'][key][:]
		#-- Getting attributes of group and included variables
		if ATTRIBUTES:
			#-- Variable Attributes
			IS2_atl07_attrs['ancillary_data'][key] = {}
			for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
				IS2_atl07_attrs['ancillary_data'][key][att_name] = att_val

	#-- sea ice ancillary information (processing flags and parameters)
	for cal in ('fine_surface_finding','sea_ice','surface_classification'):
		IS2_atl07_mds['ancillary_data'][cal] = {}
		IS2_atl07_attrs['ancillary_data'][cal] = {} if ATTRIBUTES else None
		for key,val in fileID['ancillary_data'][cal].items():
			#-- get each HDF5 variable
			IS2_atl07_mds['ancillary_data'][cal][key] = val[:]
			#-- Getting attributes of group and included variables
			if ATTRIBUTES:
				#-- Variable Attributes
				IS2_atl07_attrs['ancillary_data'][cal][key] = {}
				for att_name,att_val in val.attrs.items():
					IS2_atl07_attrs['ancillary_data'][cal][key][att_name] = att_val

	#-- get each global attribute and the attributes for orbit and quality
	if ATTRIBUTES:
		#-- ICESat-2 HDF5 global attributes
		for att_name,att_val in fileID.attrs.items():
			IS2_atl07_attrs[att_name] = att_name
		#-- ICESat-2 orbit_info Group
		IS2_atl07_attrs['orbit_info'] = {}
		for key,val in fileID['orbit_info'].items():
			IS2_atl07_attrs['orbit_info'][key] = {}
			for att_name,att_val in val.attrs.items():
				IS2_atl07_attrs['orbit_info'][key][att_name]= att_val
		#-- ICESat-2 quality_assessment Group
		IS2_atl07_attrs['quality_assessment'] = {}
		for key,val in fileID['quality_assessment'].items():
			IS2_atl07_attrs['quality_assessment'][key] = {}
			for att_name,att_val in val.attrs.items():
				IS2_atl07_attrs['quality_assessment'][key][att_name]= att_val
			if isinstance(val, h5py.Group):
				for k,v in val.items():
					IS2_atl07_attrs['quality_assessment'][key][k] = {}
					for att_name,att_val in v.attrs.items():
						IS2_atl07_attrs['quality_assessment'][key][k][att_name]= att_val

	#-- Closing the HDF5 file
	fileID.close()
	#-- Return the datasets and variables
	return (IS2_atl07_mds,IS2_atl07_attrs,IS2_atl07_beams)
