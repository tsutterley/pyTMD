#!/usr/bin/env python
u"""
compute_tidal_elevations.py
Written by Tyler Sutterley (07/2018)
Calculates tidal elevations for an input csv filea

Uses OTIS format tidal solutions provided by Ohio State University and ESR
	http://volkov.oce.orst.edu/tides/region.html
	https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
	ftp://ftp.esr.org/pub/datasets/tmd/
or Global Tide Model (GOT) solutions provided by Richard Ray at GSFC

INPUTS:
	csv file with columns:
		Modified Julian Day (days since 1858-11-17 at 00:00:00)
		latitude: degrees
		longitude: degrees
		elevation (height above or below WGS84 ellipsoid)

COMMAND LINE OPTIONS:
	-D X, --directory=X: Working data directory
	-T X, --tide=X: Tide model to use in correction
	-M X, --mode=X: Permission mode of output file

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python
		http://www.numpy.org
		http://www.scipy.org/NumPy_for_Matlab_Users
	scipy: Scientific Tools for Python
		http://www.scipy.org/

PROGRAM DEPENDENCIES:
	calc_astrol_longitudes.py: computes the basic astronomical mean longitudes
	calc_delta_time.py: calculates difference between universal and dynamic time
	convert_xy_ll.py: convert lat/lon points to and from projected coordinates
	map_ll_tides.py: converts from latitude/longitude to polar stereographic
	map_xy_tides.py: converts from polar stereographic to latitude/longitude
	load_constituent.py: loads parameters for a given tidal constituent
	load_nodal_corrections.py: load the nodal corrections for tidal constituents
	infer_minor_corrections.py: return corrections for 16 minor constituents
	read_tide_model.py: extract tidal harmonic constants out of a tidal model
	read_GOT_model.py: extract tidal harmonic constants from GSFC GOT models
	predict_tide_drift.py: predict tidal elevations using harmonic constants

UPDATE HISTORY:
	Updated 07/2018: added GSFC Global Ocean Tides (GOT) models
	Written 10/2017 for public release
"""
from __future__ import print_function

import sys
import os
import getopt
import numpy as np
import pyTMD.calc_delta_time
import pyTMD.infer_minor_corrections
import pyTMD.predict_tide_drift
from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.read_GOT_model import extract_GOT_constants

#-- PURPOSE: read HDF5 data from merge_HDF5_triangle_files.py
#-- compute tides at points and times using tidal model driver algorithms
def compute_tidal_elevations(tide_dir, input_file, output_file,
	TIDE_MODEL='', MODE=0775):

	#-- read input *.csv file to extract MJD, latitude, longitude and elevation
	dtype = dict(names=('MJD','lat','lon','h'),formats=('f8','f8','f8','f8'))
	dinput = np.loadtxt(input_file, delimiter=',', dtype=dtype)

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
		EPSG = '3032'
		type = 'z'
	elif (TIDE_MODEL == 'CATS2008_load'):
		grid_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','grid_CATS2008a_opt')
		model_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','h_CATS2008a_SPOTL_load')
		reference = ('https://www.esr.org/research/polar-tide-models/'
			'list-of-polar-tide-models/cats2008/')
		model_format = 'OTIS'
		EPSG = '3032'
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
		EPSG = '3996'
		type = 'z'
	elif (TIDE_MODEL == 'AOTIM-5'):
		grid_file = os.path.join(tide_dir,'aotim5_tmd','grid_Arc5km')
		model_file = os.path.join(tide_dir,'aotim5_tmd','h_Arc5km.oce')
		reference = ('https://www.esr.org/research/polar-tide-models/'
			'list-of-polar-tide-models/aotim-5/')
		model_format = 'OTIS'
		EPSG = '3996'
		type = 'z'
	elif (TIDE_MODEL == 'GOT4.7'):
		model_directory = os.path.join(tide_dir,'GOT4.7','grids_oceantide')
		model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
			'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
		c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
		reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
			'MiscPubs/19990089548_1999150788.pdf')
		model_format = 'GOT'
	elif (TIDE_MODEL == 'GOT4.7_load'):
		model_directory = os.path.join(tide_dir,'GOT4.7','grids_loadtide')
		model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
			'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
			's1load.d.gz','m4load.d.gz']
		c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
		reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
			'MiscPubs/19990089548_1999150788.pdf')
		model_format = 'GOT'
	elif (TIDE_MODEL == 'GOT4.8'):
		model_directory = os.path.join(tide_dir,'got4.8','grids_oceantide')
		model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
			'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
		c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
		reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
			'MiscPubs/19990089548_1999150788.pdf')
		model_format = 'GOT'
	elif (TIDE_MODEL == 'GOT4.8_load'):
		model_directory = os.path.join(tide_dir,'got4.8','grids_loadtide')
		model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
			'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
			's1load.d.gz','m4load.d.gz']
		c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
		reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
			'MiscPubs/19990089548_1999150788.pdf')
		model_format = 'GOT'
	elif (TIDE_MODEL == 'GOT4.10'):
		model_directory = os.path.join(tide_dir,'GOT4.10c','grids_oceantide')
		model_files = ['q1.d.gz','o1.d.gz','p1.d.gz','k1.d.gz','n2.d.gz',
			'm2.d.gz','s2.d.gz','k2.d.gz','s1.d.gz','m4.d.gz']
		c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
		reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
			'MiscPubs/19990089548_1999150788.pdf')
		model_format = 'GOT'
	elif (TIDE_MODEL == 'GOT4.10_load'):
		model_directory = os.path.join(tide_dir,'GOT4.10c','grids_loadtide')
		model_files = ['q1load.d.gz','o1load.d.gz','p1load.d.gz','k1load.d.gz',
			'n2load.d.gz','m2load.d.gz','s2load.d.gz','k2load.d.gz',
			's1load.d.gz','m4load.d.gz']
		c = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
		reference = ('https://denali.gsfc.nasa.gov/personal_pages/ray/'
			'MiscPubs/19990089548_1999150788.pdf')
		model_format = 'GOT'

	#-- read tidal constants and interpolate to grid points
	if (model_format == 'OTIS'):
		amp,ph,D,c = extract_tidal_constants(dinput['lon'], dinput['lat'],
			grid_file, model_file, EPSG, type, METHOD='linear')
		deltat = 0.0
	elif (model_format == 'GOT'):
		amp,ph = extract_GOT_constants(dinput['lon'], dinput['lat'],
			model_directory, model_files, METHOD='linear')
		delta_file = os.path.join(tide_dir,'deltat.data')
		deltat = pyTMD.calc_delta_time(delta_file,dinput['MJD'])

	#-- convert phase to radians
	cph = -1j*ph*np.pi/180.0
	hc = amp*np.exp(cph)

	#-- convert time from MJD to days relative to Jan 1, 1992 (48622 MJD)
	#-- predict tidal elevations at time 1 and infer minor corrections
	tide = pyTMD.predict_tide_drift(dinput['MJD'] - 48622.0, hc, c,
		DELTAT=deltat, CORRECTIONS=model_format)
	tide += pyTMD.infer_minor_corrections(dinput['MJD'] - 48622.0, hc, c,
		DELTAT=deltat, CORRECTIONS=model_format)
	tide *= 100.0
	#-- replace nan values with fill value
	fill_value = -9999.0
	inan, = np.nonzero(np.isnan(tide))
	tide[inan] = fill_value

	#-- output to file
	with open(output_file,'w') as f:
		for d,lt,ln,td in zip(dinput['MJD'],dinput['lat'],dinput['lon'],tide):
			f.write('{0:0.6f},{1:0.2f},{2:0.2f},{3:0.2f}\n'.format(d,lt,ln,td))
	#-- change the permissions level to MODE
	os.chmod(output_file, MODE)

#-- PURPOSE: help module to describe the optional input parameters
def usage():
	print('\nHelp: {}'.format(os.path.basename(sys.argv[0])))
	print(' -D X, --directory=X\tWorking data directory')
	print(' -T X, --tide=X\t\tTide model to use in correction')
	print(' -M X, --mode=X\t\tPermission mode of output file\n')

#-- Main program that calls compute_tidal_elevations()
def main():
	#-- Read the system arguments listed after the program
	long_options = ['help','directory=','tide=','mode=']
	optlist,arglist = getopt.getopt(sys.argv[1:], 'hD:T:M:', long_options)

	#-- set data directory containing the tidal data
	tide_dir = None
	#-- tide model to use
	TIDE_MODEL = 'CATS2008'
	#-- permissions mode of the local files (number in octal)
	MODE = 0775
	for opt, arg in optlist:
		if opt in ('-h','--help'):
			usage()
			sys.exit()
		elif opt in ("-D","--directory"):
			tide_dir = os.path.expanduser(arg)
		elif opt in ("-T","--tide"):
			TIDE_MODEL = arg
		elif opt in ("-M","--mode"):
			MODE = int(arg, 8)

	#-- enter input and output files as system argument
	if not arglist:
		raise Exception('No System Arguments Listed')

	#-- tilde-expand input and output csv files
	input_file = os.path.expanduser(arglist[0])
	output_file = os.path.expanduser(arglist[1])
	#-- set base directory from the input file if not current set
	tide_dir = os.path.dirname(input_file) if (tide_dir is None) else tide_dir

	#-- run tidal elevation program for input *.csv file
	compute_tidal_elevations(tide_dir, input_file, output_file,
		TIDE_MODEL=TIDE_MODEL, MODE=MODE)

#-- run main program
if __name__ == '__main__':
	main()
