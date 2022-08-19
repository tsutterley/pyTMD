#!/usr/bin/env python
u"""
calc_delta_time.py
Written by Tyler Sutterley (04/2022)
Calculates the difference between dynamic time and universal time (TT - UT1)
    following Richard Ray's PERTH3 algorithms

INPUTS:
    delta_file from
        http://maia.usno.navy.mil/ser7/deltat.data
        ftp://cddis.nasa.gov/products/iers/deltat.data
    idays: input times to interpolate (days since 1992-01-01T00:00:00)

OUTPUTS:
    deltat: delta time estimates at the output times in days

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/

REFERENCES:
    Jean Meeus, Astronomical Algorithms, 2nd edition, 1998.

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 08/2020: using builtin time operations, interpolate with tide time
        convert output units to be in days
    Updated 07/2020: added function docstrings. scipy interpolating splines
    Updated 11/2019: pad input time dimension if entering a single value
    Updated 07/2018: linearly extrapolate if using dates beyond the table
    Written 07/2018
"""
import os
import numpy as np
import scipy.interpolate
import pyTMD.time

# PURPOSE: calculate the difference between universal time and dynamical time
# by interpolating a delta time file to a given date
def calc_delta_time(delta_file, idays):
    """
    Calculates the difference between universal time (UT) and
    dynamical time (TT) [Meeus1998]_

    Parameters
    ----------
    delta_file: str
        file containing the delta times
    idays: float
        input times to interpolate (days since 1992-01-01T00:00:00)

    Returns
    -------
    deltat: float
        delta time at the input time

    References
    ----------
    .. [Meeus1998] J. Meeus, *Astronomical Algorithms*, 2nd edition, 477 pp., (1998).
    """
    # read delta time file
    dinput = np.loadtxt(os.path.expanduser(delta_file))
    # calculate Julian days and then convert to days since 1992-01-01T00:00:00
    days = pyTMD.time.convert_calendar_dates(dinput[:,0],dinput[:,1],dinput[:,2],
        epoch=(1992,1,1,0,0,0))
    # use scipy interpolating splines to interpolate delta times
    spl = scipy.interpolate.UnivariateSpline(days,dinput[:,3],k=1,s=0,ext=0)
    # return the delta time for the input date converted to days
    return spl(idays)/86400.0
