#!/usr/bin/env python
u"""
read_iers_EOP.py
Written by Tyler Sutterley (07/2020)
Provides the daily earth orientation parameters (EOP) from IERS
    http://www.usno.navy.mil/USNO/earth-orientation/eo-products/weekly
Data format: http://maia.usno.navy.mil/ser7/readme.finals

INPUTS:
    input_file: full path to IERS EOP "finals" file

OUTPUTS:
    MJD: modified julian date of EOP measurements
    x: Angular coordinate x [arcsec]
    y: Angular coordinate y [arcsec]
    flag: IERS (I) or Prediction (P) flag for polar motion values

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

REFERENCE:
    Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
        IERS Technical Note No. 36, BKG (2010)

UPDATE HISTORY:
    Updated 07/2020: added function docstrings
    Written 09/2017
"""
import re
import numpy as np

#-- calculates the daily earth orientation parameters (EOP) from IERS
def read_iers_EOP(input_file):
    """
    Calculates the daily earth orientation parameters (EOP) from IERS

    Arguments
    ---------
    input_file: full path to IERS EOP "finals" file

    Returns
    -------
    MJD: modified Julian date of EOP measurements
    x: Angular coordinate x [arcsec]
    y: Angular coordinate y [arcsec]
    flag: IERS (I) or Prediction (P) flag for polar motion values
    """
    #-- read data file splitting at line breaks
    with open(input_file,'r') as f:
        file_contents = f.read().splitlines()
    #-- number of data lines
    n_lines = len(file_contents)
    dinput = {}
    dinput['MJD'] = np.zeros((n_lines))
    dinput['x'] = np.zeros((n_lines))
    dinput['y'] = np.zeros((n_lines))
    #-- for each line in the file
    flag = 'I'
    counter = 0
    while (flag == 'I'):
        line = file_contents[counter]
        i = 2+2+2+1; j = i+8
        dinput['MJD'][counter] = np.float(line[i:j])
        i = j+1
        flag = line[i]
        i += 2; j = i+9
        dinput['x'][counter] = np.float(line[i:j])
        i = j+10; j = i+9
        dinput['y'][counter] = np.float(line[i:j])
        counter += 1
    #-- reduce to data values
    dinput['MJD'] = dinput['MJD'][:counter]
    dinput['x'] = dinput['x'][:counter]
    dinput['y'] = dinput['y'][:counter]
    #-- return the date, flag and polar motion values
    return dinput
