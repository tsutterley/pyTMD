#!/usr/bin/env python
u"""
output_otis_tides.py
Written by Tyler Sutterley (07/2020)
Writes OTIS-format tide files for use in the tide program
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 07/2020: added function docstrings
    Written 08/2018
"""
import os
import struct
import numpy as np

#-- PURPOSE: output grid file in OTIS format
def output_otis_grid(grid_file,xlim,ylim,hz,mz,iob,dt):
    """
    Writes OTIS-format grid files

    Arguments
    ---------
    grid_file: output OTIS grid file name
    xlim: x-coordinate grid-cell edges of output grid
    ylim: y-coordinate grid-cell edges of output grid
    hz: bathymety
    mz: land/water mask
    iob: open boundary index
    dt: time step
    """
    #-- open this way for files
    fid = open(grid_file,'wb')
    nob = len(iob)
    ny,nx = np.shape(hz)
    reclen = 32
    fid.write(struct.pack('>i4',reclen))
    fid.write(struct.pack('>i4',nx))
    fid.write(struct.pack('>i4',ny))
    ylim.tofile(fid,format='>f4')
    xlim.tofile(fid,format='>f4')
    fid.write(struct.pack('>f4',dt))
    fid.write(struct.pack('>i4',nob))
    fid.write(struct.pack('>i4',reclen))
    if (nob == 0):
        fid.write(struct.pack('>i4',4))
        fid.write(struct.pack('>i4',0))
        fid.write(struct.pack('>i4',4))
    else:
        reclen = 8*nob
        fid.write(struct.pack('>i4',reclen))
        iob.tofile(fid,format='>i4')
        fid.write(struct.pack('>i4',reclen))
    reclen = 4*nx*ny
    #-- write depth and mask data to file
    fid.write(struct.pack('>i4',reclen))
    hz.tofile(fid,format='>f4')
    for m in range(ny):
        hz[m,:].tofile(fid,format='>f4')
    fid.write(struct.pack('>i4',reclen))
    fid.write(struct.pack('>i4',reclen))
    for m in range(ny):
        mz[m,:].tofile(fid,format='>i4')
    fid.write(struct.pack('>i4',reclen))
    #-- close the output OTIS file
    fid.close()

#-- PURPOSE: output elevation file in OTIS format
def output_otis_elevation(elevation_file,h,xlim,ylim,constituents):
    """
    Writes OTIS-format elevation files

    Arguments
    ---------
    elevation_file: output OTIS elevation file name
    h: Eulerian form of tidal height oscillation
    xlim: x-coordinate grid-cell edges of output grid
    ylim: y-coordinate grid-cell edges of output grid
    constituents: tidal constituent IDs
    """
    fid = open(elevation_file,'wb')
    ny,nx,nc = np.shape(h)
    #-- length of header: allow for 4 character >i4 c_id strings
    header_length = 4*(7 + nc)
    fid.write(struct.pack('>i4',header_length))
    fid.write(struct.pack('>i4',nx))
    fid.write(struct.pack('>i4',ny))
    fid.write(struct.pack('>i4',nc))
    ylim.tofile(fid,format='>f4')
    xlim.tofile(fid,format='>f4')
    for c in constituents:
        fid.write('{0:4}'.format(c))
    fid.write(struct.pack('>i4',header_length))
    #-- write each constituent to file
    constituent_header = 8*nx*ny
    for ic in range(nc):
        fid.write(struct.pack('>i4',constituent_header))
        for m in range(ny):
            temp = np.zeros((2*nx),dtype='>f4')
            temp[0:2*nx-1:2] = h.real[m,:,ic]
            temp[1:2*nx:2] = h.imag[m,:,ic]
            temp.tofile(fid,format='>f4')
        fid.write(struct.pack('>i4',constituent_header))
    #-- close the output OTIS file
    fid.close()

#-- PURPOSE: output transport file in OTIS format
def output_otis_transport(transport_file,u,v,xlim,ylim,constituents):
    """
    Writes OTIS-format elevation files

    Arguments
    ---------
    transport_file: output OTIS transport file name
    u: Eulerian form of tidal zonal transport oscillation
    v: Eulerian form of tidal meridional transport oscillation
    xlim: x-coordinate grid-cell edges of output grid
    ylim: y-coordinate grid-cell edges of output grid
    constituents: tidal constituent IDs
    """
    fid = open(transport_file,'wb')
    ny,nx,nc = np.shape(u)
    #-- length of header: allow for 4 character >i4 c_id strings
    header_length = 4*(7 + nc)
    fid.write(struct.pack('>i4',header_length))
    fid.write(struct.pack('>i4',nx))
    fid.write(struct.pack('>i4',ny))
    fid.write(struct.pack('>i4',nc))
    ylim.tofile(fid,format='>f4')
    xlim.tofile(fid,format='>f4')
    for c in constituents:
        fid.write('{0:4}'.format(c))
    fid.write(struct.pack('>i4',header_length))
    #-- write each constituent to file
    constituent_header = 2*8*nx*ny
    for ic in range(nc):
        fid.write(struct.pack('>i4',constituent_header))
        for m in range(ny):
            temp = np.zeros((4*nx),dtype='>f4')
            temp[0:4*nx-3:4] = u.real[m,:,ic]
            temp[1:4*nx-2:4] = u.imag[m,:,ic]
            temp[2:4*nx-1:4] = v.real[m,:,ic]
            temp[3:4*nx:4] = v.imag[m,:,ic]
            temp.tofile(fid,format='>f4')
        fid.write(struct.pack('>i4',constituent_header))
    #-- close the output OTIS file
    fid.close()
