#!/usr/bin/env python
u"""
plot_tide_forecasts.py
Written by Tyler Sutterley (11/2019)
Plots the daily tidal displacements for a given location

Uses OTIS format tidal solutions provided by Ohio State University and ESR
	http://volkov.oce.orst.edu/tides/region.html
	https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
	ftp://ftp.esr.org/pub/datasets/tmd/
or Global Tide Model (GOT) solutions provided by Richard Ray at GSFC

COMMAND LINE OPTIONS:
	-D X, --directory=X: Working data directory with tide models
	-C X, --coordinates=X: latitude and longitude of point to forecast
	--date=X: date to forecast in ISO format (YYYY-MM-DD)
	-T X, --tide=X: Tide model to use
		CATS0201
		CATS2008
		CATS2008_load
		TPXO9-atlas
		TPXO9.1
		TPXO8-atlas
		TPXO7.2
		TPXO7.2_load
		AODTM-5
		AOTIM-5
		AOTIM-5-2018
		GOT4.7
		GOT4.7_load
		GOT4.8
		GOT4.8_load
		GOT4.10
		GOT4.10_load

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python
		http://www.numpy.org
		http://www.scipy.org/NumPy_for_Matlab_Users
	scipy: Scientific Tools for Python
		http://www.scipy.org/
	pyproj: Python interface to PROJ library
		https://pypi.org/project/pyproj/
	netCDF4: Python interface to the netCDF C library
	 	https://unidata.github.io/netcdf4-python/netCDF4/index.html
	matplotlib: Python 2D plotting library
		http://matplotlib.org/
		https://github.com/matplotlib/matplotlib

PROGRAM DEPENDENCIES:
	calc_astrol_longitudes.py: computes the basic astronomical mean longitudes
	calc_delta_time.py: calculates difference between universal and dynamic time
	convert_xy_ll.py: convert lat/lon points to and from projected coordinates
	load_constituent.py: loads parameters for a given tidal constituent
	load_nodal_corrections.py: load the nodal corrections for tidal constituents
	infer_minor_corrections.py: return corrections for 16 minor constituents
	read_tide_model.py: extract tidal harmonic constants from OTIS tide models
	read_netcdf_model.py: extract tidal harmonic constants from netcdf models
	read_GOT_model.py: extract tidal harmonic constants from GSFC GOT models
	predict_tidal_ts.py: predict tidal elevations using harmonic constants

UPDATE HISTORY:
	Updated 11/2019: added AOTIM-5-2018 tide model (2018 update to 2004 model)
	Updated 09/2019: added TPXO9_atlas reading from netcdf4 tide files
	Updated 08/2018: added correction option ATLAS for localized OTIS solutions
	Written 07/2018 for public release
"""
from __future__ import print_function

import sys
import os
import time
import getopt
import numpy as np
import matplotlib.pyplot as plt
from pyTMD.calc_delta_time import calc_delta_time
from pyTMD.infer_minor_corrections import infer_minor_corrections
from pyTMD.predict_tidal_ts import predict_tidal_ts
from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.read_netcdf_model import extract_netcdf_constants
from pyTMD.read_GOT_model import extract_GOT_constants

#-- PURPOSE: calculate the Modified Julian Day (MJD) from calendar date
#-- http://scienceworld.wolfram.com/astronomy/JulianDate.html
def calc_modified_julian_day(YEAR, MONTH, DAY):
	MJD = 367.*YEAR - np.floor(7.*(YEAR + np.floor((MONTH+9.)/12.))/4.) - \
		np.floor(3.*(np.floor((YEAR + (MONTH - 9.)/7.)/100.) + 1.)/4.) + \
		np.floor(275.*MONTH/9.) + DAY + 1721028.5 - 2400000.5
	return np.array(MJD,dtype=np.float)

#-- PURPOSE: read HDF5 data from merge_HDF5_triangle_files.py
#-- compute tides at points and times using tidal model driver algorithms
def plot_tide_forecasts(tide_dir, LON, LAT, DATE, TIDE_MODEL=''):
	#-- select between tide models
	if (TIDE_MODEL == 'CATS0201'):
		grid_file = os.path.join(tide_dir,'cats0201_tmd','grid_CATS')
		model_file = os.path.join(tide_dir,'cats0201_tmd','h0_CATS02_01')
		reference = 'https://mail.esr.org/polar_tide_models/Model_CATS0201.html'
		model_format = 'OTIS'
		EPSG = '4326'
		type = 'z'
	elif (TIDE_MODEL == 'CATS2008'):
		grid_file = os.path.join(tide_dir,'CATS2008','grid_CATS2008a_opt')
		model_file = os.path.join(tide_dir,'CATS2008','hf.CATS2008.out')
		reference = ('https://www.esr.org/research/polar-tide-models/'
			'list-of-polar-tide-models/cats2008/')
		model_format = 'OTIS'
		EPSG = 'CATS2008'
		type = 'z'
	elif (TIDE_MODEL == 'CATS2008_load'):
		grid_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','grid_CATS2008a_opt')
		model_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','h_CATS2008a_SPOTL_load')
		reference = ('https://www.esr.org/research/polar-tide-models/'
			'list-of-polar-tide-models/cats2008/')
		model_format = 'OTIS'
		EPSG = 'CATS2008'
		type = 'z'
	elif (TIDE_MODEL == 'TPXO9-atlas'):
		model_directory = os.path.join(tide_dir,'TPXO9_atlas')
		grid_file = 'grid_tpxo9_atlas.nc.gz'
		model_files = ['h_q1_tpxo9_atlas_30.nc.gz','h_o1_tpxo9_atlas_30.nc.gz',
			'h_p1_tpxo9_atlas_30.nc.gz','h_k1_tpxo9_atlas_30.nc.gz',
			'h_n2_tpxo9_atlas_30.nc.gz','h_m2_tpxo9_atlas_30.nc.gz',
			'h_s2_tpxo9_atlas_30.nc.gz','h_k2_tpxo9_atlas_30.nc.gz',
			'h_m4_tpxo9_atlas_30.nc.gz','h_ms4_tpxo9_atlas_30.nc.gz',
			'h_mn4_tpxo9_atlas_30.nc.gz','h_2n2_tpxo9_atlas_30.nc.gz']
		reference = 'http://volkov.oce.orst.edu/tides/tpxo9_atlas.html'
		model_format = 'netcdf'
		type = 'z'
		SCALE = 1.0/1000.0
	elif (TIDE_MODEL == 'TPXO9.1'):
		grid_file = os.path.join(tide_dir,'TPXO9.1','DATA','grid_tpxo9')
		model_file = os.path.join(tide_dir,'TPXO9.1','DATA','h_tpxo9.v1')
		reference = 'http://volkov.oce.orst.edu/tides/global.html'
		model_format = 'OTIS'
		EPSG = '4326'
		type = 'z'
	elif (TIDE_MODEL == 'TPXO8-atlas'):
		grid_file = os.path.join(tide_dir,'tpxo8_atlas','grid_tpxo8atlas_30_v1')
		model_file = os.path.join(tide_dir,'tpxo8_atlas','hf.tpxo8_atlas_30_v1')
		reference = 'http://volkov.oce.orst.edu/tides/tpxo8_atlas.html'
		model_format = 'ATLAS'
		EPSG = '4326'
		type = 'z'
	elif (TIDE_MODEL == 'TPXO7.2'):
		grid_file = os.path.join(tide_dir,'TPXO7.2_tmd','grid_tpxo7.2')
		model_file = os.path.join(tide_dir,'TPXO7.2_tmd','h_tpxo7.2')
		reference = 'http://volkov.oce.orst.edu/tides/global.html'
		model_format = 'OTIS'
		EPSG = '4326'
		type = 'z'
	elif (TIDE_MODEL == 'TPXO7.2_load'):
		grid_file = os.path.join(tide_dir,'TPXO7.2_load','grid_tpxo6.2')
		model_file = os.path.join(tide_dir,'TPXO7.2_load','h_tpxo7.2_load')
		reference = 'http://volkov.oce.orst.edu/tides/global.html'
		model_format = 'OTIS'
		EPSG = '4326'
		type = 'z'
	elif (TIDE_MODEL == 'AODTM-5'):
		grid_file = os.path.join(tide_dir,'aodtm5_tmd','grid_Arc5km')
		model_file = os.path.join(tide_dir,'aodtm5_tmd','h0_Arc5km.oce')
		reference = ('https://www.esr.org/research/polar-tide-models/'
			'list-of-polar-tide-models/aodtm-5/')
		model_format = 'OTIS'
		EPSG = 'PSNorth'
		type = 'z'
	elif (TIDE_MODEL == 'AOTIM-5'):
		grid_file = os.path.join(tide_dir,'aotim5_tmd','grid_Arc5km')
		model_file = os.path.join(tide_dir,'aotim5_tmd','h_Arc5km.oce')
		reference = ('https://www.esr.org/research/polar-tide-models/'
			'list-of-polar-tide-models/aotim-5/')
		model_format = 'OTIS'
		EPSG = 'PSNorth'
		type = 'z'
	elif (TIDE_MODEL == 'AOTIM-5-2018'):
		grid_file = os.path.join(tide_dir,'Arc5km2018','grid_Arc5km2018')
		model_file = os.path.join(tide_dir,'Arc5km2018','h_Arc5km2018')
		reference = ('https://www.esr.org/research/polar-tide-models/'
			'list-of-polar-tide-models/aotim-5/')
		model_format = 'OTIS'
		EPSG = 'PSNorth'
		type = 'z'
	elif (TIDE_MODEL == 'GOT4.7'):
		model_directory = os.path.join(tide_dir,'GOT4.7','grids_oceantide')
		model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
			'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
		c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
		reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
			'MiscPubs/19990089548_1999150788.pdf')
		model_format = 'GOT'
		SCALE = 1.0/100.0
	elif (TIDE_MODEL == 'GOT4.7_load'):
		model_directory = os.path.join(tide_dir,'GOT4.7','grids_loadtide')
		model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
			'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
			's1load.d.gz','m4load.d.gz']
		c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
		reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
			'MiscPubs/19990089548_1999150788.pdf')
		model_format = 'GOT'
		SCALE = 1.0/1000.0
	elif (TIDE_MODEL == 'GOT4.8'):
		model_directory = os.path.join(tide_dir,'got4.8','grids_oceantide')
		model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
			'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
		c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
		reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
			'MiscPubs/19990089548_1999150788.pdf')
		model_format = 'GOT'
		SCALE = 1.0/100.0
	elif (TIDE_MODEL == 'GOT4.8_load'):
		model_directory = os.path.join(tide_dir,'got4.8','grids_loadtide')
		model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
			'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
			's1load.d.gz','m4load.d.gz']
		c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
		reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
			'MiscPubs/19990089548_1999150788.pdf')
		model_format = 'GOT'
		SCALE = 1.0/1000.0
	elif (TIDE_MODEL == 'GOT4.10'):
		model_directory = os.path.join(tide_dir,'GOT4.10c','grids_oceantide')
		model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
			'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
		c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
		reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
			'MiscPubs/19990089548_1999150788.pdf')
		model_format = 'GOT'
		SCALE = 1.0/100.0
	elif (TIDE_MODEL == 'GOT4.10_load'):
		model_directory = os.path.join(tide_dir,'GOT4.10c','grids_loadtide')
		model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
			'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
			's1load.d.gz','m4load.d.gz']
		c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
		reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
			'MiscPubs/19990089548_1999150788.pdf')
		model_format = 'GOT'
		SCALE = 1.0/1000.0

	#-- calculate the modified Julian day from the calendar date
	MJD = calc_modified_julian_day(DATE.tm_year, DATE.tm_mon, DATE.tm_mday)
	TIME = np.arange(7*1440)/1440.0

	#-- read tidal constants and interpolate to grid points
	if model_format in ('OTIS','ATLAS'):
		amp,ph,D,c = extract_tidal_constants(np.array([LON]), np.array([LAT]),
			grid_file,model_file,EPSG,type,METHOD='spline',GRID=model_format)
		deltat = np.zeros_like(MJD)
	elif (model_format == 'netcdf'):
		amp,ph,D,c = extract_netcdf_constants(np.array([LON]), np.array([LAT]),
			model_directory, grid_file, model_files, type, METHOD='spline',
			SCALE=SCALE)
		deltat = np.zeros_like(MJD)
	elif (model_format == 'GOT'):
		amp,ph = extract_GOT_constants(np.array([LON]), np.array([LAT]),
			model_directory, model_files, METHOD='spline', SCALE=SCALE)
		delta_file = os.path.join(tide_dir,'deltat.data')
		deltat = calc_delta_time(delta_file, MJD + TIME)

	#-- calculate complex phase in radians for Euler's
	cph = -1j*ph*np.pi/180.0
	#-- calculate constituent oscillation
	hc = amp*np.exp(cph)

	#-- convert time from MJD to days relative to Jan 1, 1992 (48622 MJD)
	#-- predict tidal elevations at time 1 and infer minor corrections
	TIDE = predict_tidal_ts(MJD + TIME - 48622.0, hc, c,
		DELTAT=deltat, CORRECTIONS=model_format)
	MINOR = infer_minor_corrections(MJD + TIME - 48622.0, hc, c,
		DELTAT=deltat, CORRECTIONS=model_format)
	TIDE.data[:] += MINOR.data[:]
	#-- convert to centimeters
	TIDE.data[:] *= 100.0

	#-- differentiate to calculate high and low tides
	diff = np.zeros_like(TIME, dtype=np.float)
	#-- forward differentiation for starting point
	diff[0] = TIDE.data[1] - TIDE.data[0]
	#-- backward differentiation for end point
	diff[-1] = TIDE.data[-1] - TIDE.data[-2]
	#-- centered differentiation for all others
	diff[1:-1] = (TIDE.data[2:] - TIDE.data[0:-2])/2.0
	#-- indices of high and low tides
	htindex, = np.nonzero((np.sign(diff[0:-1]) >= 0) & (np.sign(diff[1:]) < 0))
	ltindex, = np.nonzero((np.sign(diff[0:-1]) <= 0) & (np.sign(diff[1:]) > 0))

	#-- create plot with tidal displacements, high and low tides and dates
	fig,ax1 = plt.subplots(num=1)
	ax1.plot(TIME*24.0,TIDE.data,'k')
	ax1.plot(TIME[htindex]*24.0,TIDE.data[htindex],'r*')
	ax1.plot(TIME[ltindex]*24.0,TIDE.data[ltindex],'b*')
	for h in range(24,192,24):
		ax1.axvline(h,color='gray',lw=0.5,ls='dashed',dashes=(11,5))
	ax1.set_xlim(0,7*24)
	ax1.set_ylabel('{0} Tidal Displacement [cm]'.format(TIDE_MODEL))
	args = (DATE.tm_year,DATE.tm_mon,DATE.tm_mday)
	ax1.set_xlabel('Time from {0:4d}-{1:02d}-{2:02d} UTC [Hours]'.format(*args))
	ax1.set_title(u'{0:0.6f}\u00b0N {1:0.6f}\u00b0W'.format(LAT,LON))
	fig.subplots_adjust(left=0.10,right=0.98,bottom=0.10,top=0.95)
	plt.show()

#-- PURPOSE: help module to describe the optional input parameters
def usage():
	print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
	print(' -D X, --directory=X\tWorking data directory')
	print(' -C X, --coordinates=X\tLatitude and longitude of point')
	print(' --date=X\t\tdate to forecast in ISO format (YYYY-MM-DD)')
	print(' -T X, --tide=X\t\tTide model to use in correction\n')

#-- Main program that calls plot_tide_forecasts()
def main():
	#-- Read the system arguments listed after the program
	long_options = ['help','directory=','coordinates=','date=','tide=']
	optlist,arglist = getopt.getopt(sys.argv[1:], 'hD:C:T:', long_options)

	#-- set data directory
	tide_dir = os.getcwd()
	#-- tide model to use
	TIDE_MODEL = 'GOT4.10'
	#-- coordinates to use
	LAT,LON = (32.93301304,242.7294513)
	#-- use current date
	DATE = time.gmtime()
	for opt, arg in optlist:
		if opt in ('-h','--help'):
			usage()
			sys.exit()
		elif opt in ("-D","--directory"):
			tide_dir = os.path.expanduser(arg)
		elif opt in ("-C","--coordinates"):
			LAT,LON = np.array(arg.split(','),dtype=np.float)
		elif opt in ("--date"):
			DATE = time.strptime(arg,'%Y-%m-%d')
		elif opt in ("-T","--tide"):
			TIDE_MODEL = arg

	#-- run tidal elevation program
	plot_tide_forecasts(tide_dir, LON, LAT, DATE, TIDE_MODEL=TIDE_MODEL)

#-- run main program
if __name__ == '__main__':
	main()
