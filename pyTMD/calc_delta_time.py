#!/usr/bin/env python
u"""
calc_delta_time.py
Written by Tyler Sutterley (11/2019)
Calculates the difference between universal time and dynamical time (TT - UT1)
    following Richard Ray's PERTH3 algorithms

INPUTS:
    delta_file from
        http://maia.usno.navy.mil/ser7/deltat.data
        ftp://cddis.nasa.gov/products/iers/deltat.data
    iMJD: Modified Julian Day of times to interpolate

UPDATE HISTORY:
    Updated 11/2019: pad input time dimension if entering a single value
    Updated 07/2018: linearly extrapolate if using dates beyond the table
    Written 07/2018
"""
import os
import numpy as np

#-- PURPOSE: calculate the Modified Julian Day (MJD) from calendar date
#-- http://scienceworld.wolfram.com/astronomy/JulianDate.html
def calc_modified_julian_day(YEAR, MONTH, DAY):
    MJD = 367.*YEAR - np.floor(7.*(YEAR + np.floor((MONTH+9.)/12.))/4.) - \
        np.floor(3.*(np.floor((YEAR + (MONTH - 9.)/7.)/100.) + 1.)/4.) + \
        np.floor(275.*MONTH/9.) + DAY + 1721028.5 - 2400000.5
    return np.array(MJD,dtype=np.float)

#-- interpolate delta time
def calc_delta_time(delta_file,iMJD):
    #-- change dimensions if entering a single value
    if (np.ndim(iMJD) == 0):
        iMJD = np.array([iMJD])
    #-- number of output points
    npts = len(iMJD)
    #-- allocate for delta time before and after the measurement
    delta1 = np.zeros((npts))
    delta2 = np.zeros((npts))
    deltat = np.zeros((npts))
    time1 = np.zeros((npts))
    time2 = np.zeros((npts))
    #-- read delta time file
    dinput = np.loadtxt(os.path.expanduser(delta_file))
    #-- calculate julian days and convert to MJD
    MJD = calc_modified_julian_day(dinput[:,0],dinput[:,1],dinput[:,2])
    ii, = np.nonzero(iMJD < MJD[-1])
    jj, = np.nonzero(iMJD >= MJD[-1])
    for i in ii:
        indice, = np.nonzero((MJD[:-1] < iMJD[i]) & (MJD[1:] >= iMJD[i]))
        delta1[i] = dinput[indice,3]
        time1[i] = MJD[indice]
        delta2[i] = dinput[indice+1,3]
        time2[i] = MJD[indice+1]
    #-- linearly interpolate to date and convert to days
    dt = (iMJD[ii] - time1[ii])/(time2[ii] - time1[ii])
    deltat[ii] = ((1.0-dt)*delta1[ii] + dt*delta2[ii])/86400.0
    #-- for dates beyond the maximum date in the table: extrapolate to date
    dt = (dinput[-1,3]-dinput[-2,3])/(MJD[-1] - MJD[-2])
    deltat[jj] = (iMJD[jj] - MJD[-1])*dt/86400.0
    #-- return the delta time in days after removing singleton dimensions
    return np.squeeze(deltat)
