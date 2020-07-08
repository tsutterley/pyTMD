#!/usr/bin/env python
u"""
read_ocean_pole_tide.py (07/2020)
Reads ocean pole load tide coefficients provided by IERS
http://maia.usno.navy.mil/conventions/2010/2010_official/chapter7/tn36_c7.pdf
http://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/icc7.pdf

IERS 0.5x0.5 map of ocean pole tide coefficients:
ftp://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/additional_info/
    opoleloadcoefcmcor.txt.gz

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
    Updated 07/2020: added function docstrings
    Updated 12/2018: Compatibility updates for Python3
    Written 09/2017
"""
import re
import gzip
import numpy as np

#-- PURPOSE: read real and imaginary ocean pole tide coefficients
def read_ocean_pole_tide(input_file):
    """
    Read real and imaginary ocean pole tide coefficients

    Arguments
    ---------
    input_file: IERS 0.5x0.5 map of ocean pole tide coefficients

    Returns
    -------
    ur: ocean pole tide coefficients
    glon: ocean grid longitude
    glat: ocean grid latitude
    """
    #-- read GZIP ocean pole tide file
    with gzip.open(input_file,'rb') as f:
        file_contents = f.read().splitlines()

    #-- counts the number of lines in the header
    count = 0
    #-- Reading over header text
    HEADER = True
    while HEADER:
        #-- file line at count
        line = file_contents[count]
        #-- find --------- within line to set HEADER flag to False when found
        HEADER = not bool(re.match(b'---------',line))
        #-- add 1 to counter
        count += 1

    #-- grid parameters and dimensions
    dlon,dlat = (0.50,0.50)
    glon = np.arange(dlon/2.0,360+dlon/2.0,dlon)
    glat = np.arange(90.0-dlat/2.0,-90.0-dlat/2.0,-dlat)
    nlon = len(glon)
    nlat = len(glat)
    #-- allocate for output grid maps
    ur = np.zeros((nlon,nlat),dtype=np.complex128)
    #-- read lines of file and add to output variables
    for i,line in enumerate(file_contents[count:]):
        ln,lt,urr,uri,unr,uni,uer,uei = np.array(line.split(), dtype='f8')
        ilon = np.int(ln/dlon)
        ilat = np.int((90.0-lt)/dlat)
        ur[ilon,ilat] = urr + 1j*uri

    #-- extend matrix for bilinear interpolation
    glon = extend_array(glon,dlon)
    ur = extend_matrix(ur)
    #-- return values
    return (ur,glon,glat)

#-- PURPOSE: wrapper function to extend an array
def extend_array(input_array,step_size):
    """
    Wrapper function to extend an array

    Arguments
    ---------
    input_array: array to extend
    step_size: step size between elements of array

    Returns
    -------
    temp: extended array
    """
    n = len(input_array)
    temp = np.zeros((n+2),dtype=input_array.dtype)
    #-- extended array [x-1,x0,...,xN,xN+1]
    temp[0] = input_array[0] - step_size
    temp[1:-1] = input_array[:]
    temp[-1] = input_array[-1] + step_size
    return temp

#-- PURPOSE: wrapper function to extend a matrix
def extend_matrix(input_matrix):
    """
    Wrapper function to extend a matrix

    Arguments
    ---------
    input_matrix: matrix to extend

    Returns
    -------
    temp: extended matrix
    """
    nx,ny = np.shape(input_matrix)
    temp = np.zeros((nx+2,ny),dtype=input_matrix.dtype)
    temp[0,:] = input_matrix[-1,:]
    temp[1:-1,:] = input_matrix[:,:]
    temp[-1,:] = input_matrix[0,:]
    return temp
