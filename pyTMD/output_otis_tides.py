#!/usr/bin/env python
u"""
output_otis_tides.py
Written by Tyler Sutterley (04/2022)
Writes OTIS-format tide files for use in the tide program
    http://volkov.oce.orst.edu/tides/region.html
    https://www.esr.org/research/polar-tide-models/list-of-polar-tide-models/
    ftp://ftp.esr.org/pub/datasets/tmd/

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html

UPDATE HISTORY:
    Updated 04/2022: updated docstrings to numpy documentation format
    Updated 09/2020: python3 compatibility updates for struct and utf8 encoding
    Updated 07/2020: added function docstrings
    Written 08/2018
"""
import os
import struct
import numpy as np

# PURPOSE: output grid file in OTIS format
def output_otis_grid(FILE, xlim, ylim, hz, mz, iob, dt):
    """
    Writes OTIS-format grid files

    Parameters
    ----------
    FILE: str
        output OTIS grid file name
    xlim: float
        x-coordinate grid-cell edges of output grid
    ylim: float
        y-coordinate grid-cell edges of output grid
    hz:float
        bathymetry
    mz: int
        land/water mask
    iob: int
        open boundary index
    dt: float
        time step
    """
    # open this way for files
    fid = open(os.path.expanduser(FILE), 'wb')
    nob = len(iob)
    ny,nx = np.shape(hz)
    reclen = 32
    fid.write(struct.pack('>i',reclen))
    fid.write(struct.pack('>i',nx))
    fid.write(struct.pack('>i',ny))
    ylim.tofile(fid,format='>f4')
    xlim.tofile(fid,format='>f4')
    fid.write(struct.pack('>f',dt))
    fid.write(struct.pack('>i',nob))
    fid.write(struct.pack('>i',reclen))
    if (nob == 0):
        fid.write(struct.pack('>i',4))
        fid.write(struct.pack('>i',0))
        fid.write(struct.pack('>i',4))
    else:
        reclen = 8*nob
        fid.write(struct.pack('>i',reclen))
        iob.tofile(fid,format='>i4')
        fid.write(struct.pack('>i',reclen))
    reclen = 4*nx*ny
    # write depth and mask data to file
    fid.write(struct.pack('>i',reclen))
    hz.tofile(fid,format='>f4')
    for m in range(ny):
        hz[m,:].tofile(fid,format='>f4')
    fid.write(struct.pack('>i',reclen))
    fid.write(struct.pack('>i',reclen))
    for m in range(ny):
        mz[m,:].tofile(fid,format='>i4')
    fid.write(struct.pack('>i',reclen))
    # close the output OTIS file
    fid.close()

# PURPOSE: output elevation file in OTIS format
def output_otis_elevation(FILE, h, xlim, ylim, constituents):
    """
    Writes OTIS-format elevation files

    Parameters
    ----------
    FILE: str
        output OTIS elevation file name
    h: complex
        Eulerian form of tidal height oscillation
    xlim: float
        x-coordinate grid-cell edges of output grid
    ylim: float
        y-coordinate grid-cell edges of output grid
    constituents: list
        tidal constituent IDs
    """
    fid = open(os.path.expanduser(FILE), 'wb')
    ny,nx,nc = np.shape(h)
    # length of header: allow for 4 character >i c_id strings
    header_length = 4*(7 + nc)
    fid.write(struct.pack('>i',header_length))
    fid.write(struct.pack('>i',nx))
    fid.write(struct.pack('>i',ny))
    fid.write(struct.pack('>i',nc))
    ylim.tofile(fid,format='>f4')
    xlim.tofile(fid,format='>f4')
    for c in constituents:
        fid.write(c.ljust(4).encode('utf8'))
    fid.write(struct.pack('>i',header_length))
    # write each constituent to file
    constituent_header = 8*nx*ny
    for ic in range(nc):
        fid.write(struct.pack('>i',constituent_header))
        for m in range(ny):
            temp = np.zeros((2*nx),dtype='>f')
            temp[0:2*nx-1:2] = h.real[m,:,ic]
            temp[1:2*nx:2] = h.imag[m,:,ic]
            temp.tofile(fid,format='>f4')
        fid.write(struct.pack('>i',constituent_header))
    # close the output OTIS file
    fid.close()

# PURPOSE: output transport file in OTIS format
def output_otis_transport(FILE, u, v, xlim, ylim, constituents):
    """
    Writes OTIS-format transport files

    Parameters
    ----------
    FILE: str
        output OTIS transport file name
    u: complex
        Eulerian form of tidal zonal transport oscillation
    v: complex
        Eulerian form of tidal meridional transport oscillation
    xlim: float
        x-coordinate grid-cell edges of output grid
    ylim: float
        y-coordinate grid-cell edges of output grid
    constituents: list
        tidal constituent IDs
    """
    fid = open(os.path.expanduser(FILE), 'wb')
    ny,nx,nc = np.shape(u)
    # length of header: allow for 4 character >i c_id strings
    header_length = 4*(7 + nc)
    fid.write(struct.pack('>i',header_length))
    fid.write(struct.pack('>i',nx))
    fid.write(struct.pack('>i',ny))
    fid.write(struct.pack('>i',nc))
    ylim.tofile(fid,format='>f4')
    xlim.tofile(fid,format='>f4')
    for c in constituents:
        fid.write(c.ljust(4).encode('utf8'))
    fid.write(struct.pack('>i',header_length))
    # write each constituent to file
    constituent_header = 2*8*nx*ny
    for ic in range(nc):
        fid.write(struct.pack('>i',constituent_header))
        for m in range(ny):
            temp = np.zeros((4*nx),dtype='>f')
            temp[0:4*nx-3:4] = u.real[m,:,ic]
            temp[1:4*nx-2:4] = u.imag[m,:,ic]
            temp[2:4*nx-1:4] = v.real[m,:,ic]
            temp[3:4*nx:4] = v.imag[m,:,ic]
            temp.tofile(fid,format='>f4')
        fid.write(struct.pack('>i',constituent_header))
    # close the output OTIS file
    fid.close()
