#!/usr/bin/env python
u"""
read_ATM1b_QFIT_binary.py
Written by Tyler Sutterley (02/2020)

Reads Level-1b Airborne Topographic Mapper (ATM) QFIT binary data products
	http://nsidc.org/data/docs/daac/icebridge/ilatm1b/docs/ReadMe.qfit.txt

Can be the following ATM QFIT file types:
	ILATM1B: Airborne Topographic Mapper QFIT Elevation
	BLATM1B: Pre-Icebridge Airborne Topographic Mapper QFIT Elevation
	ILNSA1B: Narrow Swath Airborne Topographic Mapper QFIT Elevation

Based on the QFIT C reader provided on NSIDC
	ftp://sidads.colorado.edu/pub/tools/icebridge/qfit/c/

INPUTS:
	full_filename: full path to ATM QFIT .qi file (can have tilde-prefix)

OUTPUTS:
	Data variables for the given input .qi file format listed below
		outputs are scaled from the inputs listed in the ReadMe.qfit.txt file

	10-word format (used prior to 2006):
		time:			Relative Time (seconds from start of data file)
		latitude:		Laser Spot Latitude (degrees)
		longitude:		Laser Spot Longitude (degrees)
		elevation:		Elevation above WGS84 ellipsoid (meters)
		xmt_sigstr:		Start Pulse Signal Strength (relative)
		rcv_sigstr:		Reflected Laser Signal Strength (relative)
		azimuth:		Scan Azimuth (degrees)
		pitch:			Pitch (degrees)
		roll:			Roll (degrees)
		time_hhmmss:	GPS Time packed (example: 153320.1000 = 15h 33m 20.1s)
		time_J2000:		Time converted to seconds since 2000-01-01 12:00:00 UTC

	12-word format (in use since 2006):
		time:			Relative Time (seconds from start of data file)
		latitude:		Laser Spot Latitude (degrees)
		longitude:		Laser Spot Longitude (degrees)
		elevation:		Elevation above WGS84 ellipsoid (meters)
		xmt_sigstr:		Start Pulse Signal Strength (relative)
		rcv_sigstr:		Reflected Laser Signal Strength (relative)
		azimuth:		Scan Azimuth (degrees)
		pitch:			Pitch (degrees)
		roll:			Roll (degrees)
		gps_pdop:		GPS PDOP (dilution of precision)
		pulse_width:	Laser received pulse width (digitizer samples)
		time_hhmmss:	GPS Time packed (example: 153320.1000 = 15h 33m 20.1s)
		time_J2000:		Time converted to seconds since 2000-01-01 12:00:00 UTC

	14-word format (used in some surveys between 1997 and 2004):
		time:			Relative Time (seconds from start of data file)
		latitude:		Laser Spot Latitude (degrees)
		longitude:		Laser Spot Longitude (degrees)
		elevation:		Elevation above WGS84 ellipsoid (meters)
		xmt_sigstr:		Start Pulse Signal Strength (relative)
		rcv_sigstr:		Reflected Laser Signal Strength (relative)
		azimuth:		Scan Azimuth (degrees)
		pitch:			Pitch (degrees)
		roll:			Roll (degrees)
		passive_sig:	Passive Signal (relative)
		pass_foot_lat:	Passive Footprint Latitude (degrees)
		pass_foot_long:	Passive Footprint Longitude (degrees)
		pass_foot_synth_elev:	Passive Footprint Synthesized Elevation (meters)
		time_hhmmss:	GPS Time packed (example: 153320.1000 = 15h 33m 20.1s)
		time_J2000:		Time converted to seconds since 2000-01-01 12:00:00 UTC

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python
		http://www.numpy.org
		http://www.scipy.org/NumPy_for_Matlab_Users

UPDATE HISTORY:
	Updated 02/2020: using python3 division for calculating record counts
	Updated 10/2018: updated GPS time calculation for calculating leap seconds
	Updated 01/2018: simplified regex for extracting YYMMSS from filenames
	Updated 10/2017: value as integer if big-endian (was outputting as list)
		calculate and output time as J2000 in addition to packed hhmmss
	Updated 06/2017: read and output ATM QFIT file headers
	Written 05/2017
"""
from __future__ import print_function, division

import os
import re
import numpy as np
from pyTMD.count_leap_seconds import count_leap_seconds

#-- PURPOSE: get the record length and endianness of the input QFIT file
def get_record_length(fid):
	#-- assume initially big endian (all input data 32-bit integers)
	dtype = np.dtype('>i4')
	value, = np.fromfile(fid, dtype=dtype, count=1)
	fid.seek(0)
	#-- swap to little endian and reread first line
	if (value > 100):
		dtype = np.dtype('<i4')
		value, = np.fromfile(fid, dtype=dtype, count=1)
		fid.seek(0)
	#-- get the number of variables
	n_blocks = value//dtype.itemsize
	#-- read past first record
	np.fromfile(fid, dtype=dtype, count=n_blocks)
	#-- return the number of variables and the endianness
	return (n_blocks, dtype)

#-- PURPOSE: get length and text of ATM1b file headers
def read_ATM1b_QFIT_header(fid, n_blocks, dtype):
	header_count = 0
	header_text = ''
	value = np.full((n_blocks), -1, dtype=np.int32)
	while (value[0] < 0):
		#-- read past first record
		line = fid.read(n_blocks*dtype.itemsize)
		value = np.fromstring(line, dtype=dtype, count=n_blocks)
		header_text += line[dtype.itemsize:]
		header_count += dtype.itemsize*n_blocks
	#-- rewind file to previous record and remove last record from header text
	fid.seek(header_count)
	header_text = header_text[:-dtype.itemsize*n_blocks]
	return header_count, header_text.replace('\x00','').rstrip()

#-- PURPOSE: read ATM L1b variables from a QFIT binary file
def read_ATM1b_QFIT_records(fid,n_blocks,n_records,dtype,date,SUBSETTER=None):
	#-- 10 word format = 0
	#-- 12 word format = 1
	#-- 14 word format = 2
	w = (n_blocks-10)//2
	#-- scaling factors for each variable for the 3 word formats (14 max)
	scaling_table = [
		[1e3, 1e6, 1e6, 1e3, 1, 1, 1e3, 1e3, 1e3, 1e3],
		[1e3, 1e6, 1e6, 1e3, 1, 1, 1e3, 1e3, 1e3, 1.0e1, 1, 1e3],
		[1e3, 1e6, 1e6, 1e3, 1, 1, 1e3, 1e3, 1e3, 1, 1e6, 1e6, 1e3, 1e3]]
	#-- input variable names for the 3 word formats (14 max)
	variable_table = []
	#-- 10 word format
	variable_table.append(['rel_time','latitude','longitude','elevation',
		'xmt_sigstr','rcv_sigstr','azimuth','pitch','roll','time_hhmmss'])
	#-- 12 word format
	variable_table.append(['rel_time','latitude','longitude','elevation',
		'xmt_sigstr','rcv_sigstr','azimuth','pitch','roll',
		'gps_pdop','pulse_width','time_hhmmss'])
	#-- 14 word format
	variable_table.append(['rel_time','latitude','longitude','elevation',
		'xmt_sigstr','rcv_sigstr','azimuth','pitch','roll','passive_sig',
		'pass_foot_lat','pass_foot_long','pass_foot_synth_elev','time_hhmmss'])
	#-- input variable data types for the 3 word formats (14 max)
	dtype_table = []
	#-- 10 word format
	dtype_table.append(['f','f','f','f','i','i','f','f','f','f'])
	#-- 12 word format
	dtype_table.append(['f','f','f','f','i','i','f','f','f','f','i','f'])
	#-- 14 word format
	dtype_table.append(['f','f','f','f','i','i','f','f','f','i','f','f','f','f'])
	#-- dictionary with output variables
	ATM_L1b_input = {}
	for n,d in zip(variable_table[w],dtype_table[w]):
		ATM_L1b_input[n] = np.zeros((n_records), dtype=np.dtype(d))
	#-- hour, minute and second from time_hhmmss
	hour = np.zeros((n_records),dtype=np.float)
	minute = np.zeros((n_records),dtype=np.float)
	second = np.zeros((n_records),dtype=np.float)
	#-- for each record in the ATM Level-1b file
	for r in range(n_records):
		#-- set binary to point if using input subsetter
		if SUBSETTER is not None:
			fid.seek(SUBSETTER[r])
		#-- input data record r
		i = np.fromfile(fid,dtype=dtype,count=n_blocks)
		#-- read variable and scale to output format
		for v,n,d,s in zip(i,variable_table[w],dtype_table[w],scaling_table[w]):
			ATM_L1b_input[n][r] = v.astype(d)/s
		#-- unpack GPS time
		time_hhmmss = '{0:010.3f}'.format(ATM_L1b_input['time_hhmmss'][r])
		hour[r] = np.float(time_hhmmss[:2])
		minute[r] = np.float(time_hhmmss[2:4])
		second[r] = np.float(time_hhmmss[4:])
	#-- leap seconds for converting from GPS time to UTC
	S = calc_GPS_to_UTC(date[0],date[1],date[2],hour,minute,second)
	#-- calculation of Julian day
	JD = calc_julian_day(date[0],date[1],date[2],hour,minute,second-S)
	#-- converting to J2000 seconds
	ATM_L1b_input['time_J2000'] = (JD - 2451545.0)*86400.0
	#-- return the input data dictionary
	return ATM_L1b_input

#-- PURPOSE: calculate the Julian day from calendar date
#-- http://scienceworld.wolfram.com/astronomy/JulianDate.html
def calc_julian_day(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND):
	JD = 367.*YEAR - np.floor(7.*(YEAR + np.floor((MONTH+9.)/12.))/4.) - \
		np.floor(3.*(np.floor((YEAR + (MONTH - 9.)/7.)/100.) + 1.)/4.) + \
		np.floor(275.*MONTH/9.) + DAY + 1721028.5 + HOUR/24. + MINUTE/1440. + \
		SECOND/86400.
	return np.array(JD,dtype=np.float)

#-- PURPOSE: calculate the number of leap seconds between GPS time (seconds
#-- since Jan 6, 1980 00:00:00) and UTC
def calc_GPS_to_UTC(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND):
	GPS = 367.*YEAR - np.floor(7.*(YEAR + np.floor((MONTH+9.)/12.))/4.) - \
		np.floor(3.*(np.floor((YEAR + (MONTH - 9.)/7.)/100.) + 1.)/4.) + \
		np.floor(275.*MONTH/9.) + DAY + 1721028.5 - 2444244.5
	GPS_Time = GPS*86400.0 + HOUR*3600.0 + MINUTE*60.0 + SECOND
	return count_leap_seconds(GPS_Time)

#-- PURPOSE: get shape of ATM Level-1b binary file without reading data
def ATM1b_QFIT_shape(full_filename):
	#-- read the input file to get file information
	fd = os.open(os.path.expanduser(full_filename),os.O_RDONLY)
	file_info = os.fstat(fd)
	#-- open the filename in binary read mode
	fid = os.fdopen(fd, 'rb')
	#-- get the number of variables and the endianness of the file
	n_blocks,dtype = get_record_length(fid)
	MAXARG = 14
	#-- check that the number of blocks per record is less than MAXARG
	if (n_blocks > MAXARG):
		raise Exception('ERROR: Unexpected number of variables')
	#-- read over header text
	header_count,header_text = read_ATM1b_QFIT_header(fid, n_blocks, dtype)
	#-- number of records within file
	n_records = (file_info.st_size-header_count)//n_blocks//dtype.itemsize
	#-- close the input file
	fid.close()
	#-- return the data shape
	return n_records

#-- PURPOSE: read ATM Level-1b QFIT binary file
def read_ATM1b_QFIT_binary(full_filename, SUBSETTER=None):
	#-- read the input file to get file information
	fd = os.open(os.path.expanduser(full_filename),os.O_RDONLY)
	file_info = os.fstat(fd)
	#-- open the filename in binary read mode
	fid = os.fdopen(fd, 'rb')

	#-- regular expression pattern for extracting parameters
	rx = re.compile('(BLATM1B|ILATM1B)_(\d+)(.*?)\.qi$', re.VERBOSE)
	#-- extract mission and other parameters from filename
	MISSION,YYMMDD,AUX=rx.findall(os.path.basename(full_filename)).pop()
	#-- early date strings omitted century and millenia (e.g. 93 for 1993)
	if (len(YYMMDD) == 6):
		ypre,month,day = np.array([YYMMDD[:2],YYMMDD[2:4],YYMMDD[4:]],dtype='f')
		year = (ypre + 1900.0) if (ypre >= 90) else (ypre + 2000.0)
	elif (len(YYMMDD) == 8):
		year,month,day = np.array([YYMMDD[:4],YYMMDD[4:6],YYMMDD[6:]],dtype='f')

	#-- get the number of variables and the endianness of the file
	n_blocks,dtype = get_record_length(fid)
	MAXARG = 14
	#-- check that the number of blocks per record is less than MAXARG
	if (n_blocks > MAXARG):
		raise Exception('ERROR: Unexpected number of variables')
	#-- read over header text
	header_count,header_text = read_ATM1b_QFIT_header(fid, n_blocks, dtype)

	#-- number of records to read with and without input subsetter
	if SUBSETTER is None:
		#-- number of records within file (file size - header size)
		n_records = (file_info.st_size-header_count)//n_blocks//dtype.itemsize
	else:
		#-- number of records in subsetter
		n_records = len(SUBSETTER)
		#-- convert from data point indices into binary variable indices
		SUBSETTER = header_count + dtype.itemsize*(np.array(SUBSETTER)*n_blocks)

	#-- read input data
	ATM_L1b_input = read_ATM1b_QFIT_records(fid, n_blocks, n_records, dtype,
		[year, month, day], SUBSETTER=SUBSETTER)

	#-- close the input file
	fid.close()
	#-- return the data and header text
	return ATM_L1b_input, header_text
