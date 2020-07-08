#!/usr/bin/env python
u"""
calc_iers_mean_pole.py
Written by Tyler Sutterley (07/2020)
Calculates the mean pole coordinates x and y are obtained by Gaussian-weighted
    average of the IERS C01 pole coordinates following the provided readme
    ftp://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.readme
Should provide updated values similar to
    ftp://hpiers.obspm.fr/iers/eop/eopc01/mean-pole.tab
Coordinates (particularly at end of the time series) will change with more dates

INPUTS:
    input_file: full path to eopc01.1900-now.dat file provided by IERS

OUTPUTS:
    T: date [decimal-years]
    xm: mean pole coordinate x [arcsec]
    ym: mean pole coordinate y [arcsec]

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 07/2020: added function docstrings
    Written 10/2017
"""
from __future__ import print_function

import sys
import os
import time
import numpy as np

#-- PURPOSE: read table of IERS pole coordinates and calculate Gaussian-average
def calc_iers_mean_pole(input_file):
    """
    Calculates the mean pole coordinates x and y are obtained by a
    Gaussian-weighted average of the IERS pole coordinates

    Arguments
    ---------
    input_file: pole coordinates file
    """
    #-- read the pole coordinates file (e.g. eopc01.1900-now.dat)
    dinput = np.loadtxt(os.path.expanduser(input_file))
    DIRECTORY = os.path.dirname(os.path.expanduser(input_file))
    ndat,ncol = np.shape(dinput)
    Ti = dinput[:,0]
    Xi = dinput[:,1]
    Yi = dinput[:,3]
    xm = np.zeros((ndat))
    ym = np.zeros((ndat))
    #-- output file with mean pole coordinates
    today = time.strftime('%Y-%m-%d',time.localtime())
    fid = open(os.path.join(DIRECTORY,'mean_pole_{0}.tab'.format(today)),'w')
    for i,T in enumerate(Ti):
        #-- mean pole is Gaussian Weight of all dates with a = 3.40 years.
        Wi = np.exp(-0.5*((Ti-T)/3.4)**2)
        xm[i] = np.sum(Wi*Xi)/np.sum(Wi)
        ym[i] = np.sum(Wi*Yi)/np.sum(Wi)
        print('{0:6.2f} {1:11.7f} {2:11.7f}'.format(T,xm[i],ym[i]),file=fid)
    #-- close the output file
    fid.close()

#-- run program with input file from system argument
if __name__ == '__main__':
    calc_iers_mean_pole(os.path.expanduser(sys.argv[1]))
