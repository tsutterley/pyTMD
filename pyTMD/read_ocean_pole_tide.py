#!/usr/bin/env python
u"""
read_ocean_pole_tide.py
Written by Tyler Sutterley (12/2022)

Reads ocean pole load tide coefficients provided by IERS
http://maia.usno.navy.mil/conventions/2010/2010_official/chapter7/tn36_c7.pdf
http://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/icc7.pdf

IERS 0.5x0.5 map of ocean pole tide coefficients:
ftp://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/additional_info/
    opoleloadcoefcmcor.txt.gz

OUTPUTS:
    ur: radial ocean pole tide coefficients
    un: north ocean pole tide coefficients
    ue: east ocean pole tide coefficients
    glon: ocean grid longitude
    glat: ocean grid latitude

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

REFERENCES:
    S Desai, "Observing the pole tide with satellite altimetry", Journal of
        Geophysical Research: Oceans, 107(C11), 2002. doi: 10.1029/2001JC001224
    S Desai, J Wahr and B Beckley "Revisiting the pole tide for and from
        satellite altimetry", Journal of Geodesy, 89(12), p1233-1243, 2015.
        doi: 10.1007/s00190-015-0848-7

UPDATE HISTORY:
    Updated 12/2022: refactored ocean pole tide read program under io
    Updated 04/2022: updated docstrings to numpy documentation format
        use longcomplex data format to be windows compliant
    Updated 07/2021: added check that ocean pole tide file is accessible
    Updated 03/2021: replaced numpy bool/int to prevent deprecation warnings
    Updated 08/2020: output north load and east load deformation components
    Updated 07/2020: added function docstrings
    Updated 12/2018: Compatibility updates for Python3
    Written 09/2017
"""
import warnings
import pyTMD.io

# PURPOSE: read real and imaginary ocean pole tide coefficients
def read_ocean_pole_tide(input_file):
    """
    Read real and imaginary ocean pole tide coefficients

    Parameters
    ----------
    input_file: str
        IERS 0.5x0.5 map of ocean pole tide coefficients

    Returns
    -------
    ur: float
        radial ocean pole tide coefficients
    un: float
        north ocean pole tide coefficients
    ue: float
        east ocean pole tide coefficients
    glon: float
        ocean grid longitude
    glat: float
        ocean grid latitude

    References
    ----------
    .. [1] S Desai, "Observing the pole tide with satellite altimetry", *Journal of
        Geophysical Research: Oceans*, 107(C11), (2002).
        `doi: 10.1029/2001JC001224 <https://doi.org/10.1029/2001JC001224>`_
    .. [2] S Desai, J Wahr and B Beckley "Revisiting the pole tide for and from
        satellite altimetry", *Journal of Geodesy*, 89(12), p1233-1243, (2015).
        `doi: 10.1007/s00190-015-0848-7 <https://doi.org/10.1007/s00190-015-0848-7>`_
    """
    # raise warnings for deprecation of module
    warnings.filterwarnings("always")
    warnings.warn("Deprecated. Please use pyTMD.io instead",DeprecationWarning)
    return pyTMD.io.ocean_pole_tide(input_file)
