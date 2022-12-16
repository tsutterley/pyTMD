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
    Updated 12/2022: refactor delta time as part of time module
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 08/2020: using builtin time operations, interpolate with tide time
        convert output units to be in days
    Updated 07/2020: added function docstrings. scipy interpolating splines
    Updated 11/2019: pad input time dimension if entering a single value
    Updated 07/2018: linearly extrapolate if using dates beyond the table
    Written 07/2018
"""
import warnings
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
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.time instead",DeprecationWarning)
    # call renamed version to not break workflows
    return pyTMD.time.interpolate_delta_time(delta_file, idays)
