#!/usr/bin/env python
u"""
reduce_OTIS_files.py
Written by Tyler Sutterley (07/2020)
Read OTIS-format tidal files and reduce to a regional subset

COMMAND LINE OPTIONS:
    -D X, --directory=X: working data directory
    -T X, --tide=X: Tide model to use
    -B X, --bounds=X: Grid Bounds (xmin,xmax,ymin,ymax)
    -M X, --mode=X: permissions mode of the output files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    read_tide_model.py: extract tidal harmonic constants out of a tidal model
    convert_ll_xy.py: converts lat/lon points to and from projected coordinates
    output_otis_tides.py: writes OTIS-format tide files

UPDATE HISTORY:
    Updated 07/2020: renamed coordinate conversion program
    Updated 02/2020: changed CATS2008 grid to match version on U.S. Antarctic
        Program Data Center http://www.usap-dc.org/view/dataset/601235
    Updated 11/2019: added AOTIM-5-2018 tide model (2018 update to 2004 model)
    Written 08/2018
"""
from __future__ import print_function

import sys
import os
import getopt
import numpy as np
from pyTMD.convert_ll_xy import convert_ll_xy
from pyTMD.read_tide_model import *
from pyTMD.output_otis_tides import *

#-- PURPOSE: read 1 degree land sea mask and create masks for specific regions
def make_regional_OTIS_files(tide_dir, TIDE_MODEL, BOUNDS, MODE=0o775):
    #-- select between tide models
    if (TIDE_MODEL == 'CATS0201'):
        grid_file = os.path.join(tide_dir,'cats0201_tmd','grid_CATS')
        z_file = os.path.join(tide_dir,'cats0201_tmd','h0_CATS02_01')
        uv_file = os.path.join(tide_dir,'cats0201_tmd','UV0_CATS02_01')
        reference = 'https://mail.esr.org/polar_tide_models/Model_CATS0201.html'
        model_format = 'OTIS'
        EPSG = '4326'
    elif (TIDE_MODEL == 'CATS2008'):
        grid_file = os.path.join(tide_dir,'CATS2008','grid_CATS2008')
        z_file = os.path.join(tide_dir,'CATS2008','hf.CATS2008.out')
        uv_file = os.path.join(tide_dir,'CATS2008','uv.CATS2008.out')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/cats2008/')
        model_format = 'OTIS'
        EPSG = 'CATS2008'
    elif (TIDE_MODEL == 'CATS2008_load'):
        grid_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','grid_CATS2008a_opt')
        z_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','h_CATS2008a_SPOTL_load')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/cats2008/')
        model_format = 'OTIS'
        EPSG = 'CATS2008'
    elif (MODEL == 'TPXO9_atlas'):
        grid_file = os.path.join(tide_dir,'tpxo9_atlas','grid_tpxo9atlas_30_v1')
        z_file = os.path.join(tide_dir,'tpxo9_atlas','hf.tpxo9_atlas_30_v1')
        uv_file = os.path.join(tide_dir,'tpxo9_atlas','uv.tpxo9_atlas_30_v1')
        reference = 'http://volkov.oce.orst.edu/tides/tpxo9_atlas.html'
        model_format = 'ATLAS'
        EPSG = '4326'
    elif (TIDE_MODEL == 'TPXO9.1'):
        grid_file = os.path.join(tide_dir,'TPX09.1','DATA','grid_tpxo9')
        z_file = os.path.join(tide_dir,'TPX09.1','DATA','h_tpxo9.v1')
        uv = os.path.join(tide_dir,'TPX09.1','DATA','h_tpxo9.v1')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        model_format = 'OTIS'
        EPSG = '4326'
    elif (TIDE_MODEL == 'TPXO8-atlas'):
        grid_file = os.path.join(tide_dir,'tpxo8_atlas','grid_tpxo8atlas_30_v1')
        z_file = os.path.join(tide_dir,'tpxo8_atlas','hf.tpxo8_atlas_30_v1')
        uv_file = os.path.join(tide_dir,'tpxo8_atlas','uv.tpxo8_atlas_30_v1')
        reference = 'http://volkov.oce.orst.edu/tides/tpxo8_atlas.html'
        model_format = 'ATLAS'
        EPSG = '4326'
    elif (TIDE_MODEL == 'TPXO7.2'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_tmd','grid_tpxo7.2')
        z_file = os.path.join(tide_dir,'TPXO7.2_tmd','h_tpxo7.2')
        uv_file = os.path.join(tide_dir,'TPXO7.2_tmd','u_tpxo7.2')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        model_format = 'OTIS'
        EPSG = '4326'
    elif (TIDE_MODEL == 'TPXO7.2_load'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_load','grid_tpxo6.2')
        z_file = os.path.join(tide_dir,'TPXO7.2_load','h_tpxo7.2_load')
        uv_file = os.path.join(tide_dir,'TPXO7.2_load','u_tpxo7.2_load')
        reference = 'http://volkov.oce.orst.edu/tides/global.html'
        model_format = 'OTIS'
        EPSG = '4326'
    elif (TIDE_MODEL == 'AODTM-5'):
        grid_file = os.path.join(tide_dir,'aodtm5_tmd','grid_Arc5km')
        z_file = os.path.join(tide_dir,'aodtm5_tmd','h0_Arc5km.oce')
        uv_file = os.path.join(tide_dir,'aodtm5_tmd','u0_Arc5km.oce')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aodtm-5/')
        model_format = 'OTIS'
        EPSG = 'PSNorth'
    elif (TIDE_MODEL == 'AOTIM-5'):
        grid_file = os.path.join(tide_dir,'aotim5_tmd','grid_Arc5km')
        z_file = os.path.join(tide_dir,'aotim5_tmd','h_Arc5km.oce')
        uv_file = os.path.join(tide_dir,'aotim5_tmd','u_Arc5km.oce')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aotim-5/')
        model_format = 'OTIS'
        EPSG = 'PSNorth'
    elif (TIDE_MODEL == 'AOTIM-5-2018'):
        grid_file = os.path.join(tide_dir,'aotim5_tmd','grid_Arc5km')
        z_file = os.path.join(tide_dir,'aotim5_tmd','h_Arc5km2018')
        uv_file = os.path.join(tide_dir,'aotim5_tmd','u_Arc5km2018_T')
        reference = ('https://www.esr.org/research/polar-tide-models/'
            'list-of-polar-tide-models/aotim-5/')
        model_format = 'OTIS'
        EPSG = 'PSNorth'

    #-- read the OTIS-format tide grid file
    if (model_format == 'ATLAS'):
        #-- if reading a global solution with localized solutions
        x0,y0,hz0,mz0,iob,dt,pmask,local = read_atlas_grid(grid_file)
        xi,yi,hz = combine_atlas_model(x0,y0,hz0,pmask,local,VARIABLE='depth')
        mz = create_atlas_mask(x0,y0,mz0,local,VARIABLE='depth')
    else:
        #-- if reading a pure global solution
        xi,yi,hz,mz,iob,dt = read_tide_grid(grid_file)

    #-- convert bounds from latitude/longitude to model coordinates
    x,y = convert_ll_xy(BOUNDS[:2],BOUNDS[2:],EPSG,'F')
    #-- find indices to reduce to xmin,xmax,ymin,ymax
    gridx,gridy = np.meshgrid(xi,yi)
    indy,indx = np.nonzero((gridx >= x[0]) & (gridx <= x[1]) &
        (gridy >= y[0]) & (gridy <= y[1]))
    nx = np.count_nonzero((xi >= x[0]) & (xi <= x[1]))
    ny = np.count_nonzero((yi >= y[0]) & (yi <= y[1]))
    #-- calculate new grid limits and convert back to grid-cell edges
    dx = np.abs(xi[1]-xi[0])
    dy = np.abs(yi[1]-yi[0])
    xlim = np.array([xi[indx[0]]-dx/2.0,xi[indx[-1]]+dx/2.0],dtype='>f4')
    ylim = np.array([yi[indy[0]]-dy/2.0,yi[indy[-1]]+dy/2.0],dtype='>f4')
    #-- reduce grid and mask to new bounds
    hz1 = np.zeros((ny,nx),dtype='>f4')
    mz1 = np.zeros((ny,nx),dtype='>i4')
    hz1[:,:] = hz[indy,indx].reshape(ny,nx)
    mz1[:,:] = mz[indy,indx].reshape(ny,nx)
    #-- output reduced grid to file
    new_grid_file = create_unique_filename('{0}.reduced'.format(grid_file))
    output_otis_grid(new_grid_file,xlim,ylim,hz1,mz1,iob,dt)
    #-- change the permissions level to MODE
    os.chmod(new_grid_file, MODE)

    #-- read each constituent
    constituents,nc = read_constituents(z_file)
    z1 = np.zeros((ny,nx,nc),dtype=np.complex64)
    u1 = np.zeros((ny,nx,nc),dtype=np.complex64)
    v1 = np.zeros((ny,nx,nc),dtype=np.complex64)
    for i,c in enumerate(constituents):
        #-- read constituent from elevation file
        if (model_format == 'ATLAS'):
            z0,zlocal = read_atlas_elevation(z_file,i,c)
            xi,yi,z = combine_atlas_model(x0,y0,z0,pmask,zlocal,VARIABLE='z')
        else:
            z = read_elevation_file(z_file,i)
        #-- reduce elevation to new bounds
        z1[:,:,i] = z[indy,indx].reshape(ny,nx)
        #-- read constituent from transport file
        if (model_format == 'ATLAS'):
            u0,v0,uvlocal = read_atlas_transport(uv_file,i,c)
            xi,yi,u = combine_atlas_model(x0,y0,u0,pmask,uvlocal,VARIABLE='u')
            xi,yi,v = combine_atlas_model(x0,y0,v0,pmask,uvlocal,VARIABLE='v')
        else:
            u,v = read_transport_file(uv_file,i)
        #-- reduce transport components to new bounds
        u1[:,:,i] = u[indy,indx].reshape(ny,nx)
        v1[:,:,i] = v[indy,indx].reshape(ny,nx)
    #-- output reduced elevation and transport components
    new_z_file = create_unique_filename('{0}.reduced'.format(z_file))
    new_uv_file = create_unique_filename('{0}.reduced'.format(uv_file))
    output_otis_elevation(new_z_file,z1,xlim,ylim,constituents)
    output_otis_transport(new_uv_file,u1,v1,xlim,ylim,constituents)
    #-- change the permissions level to MODE
    os.chmod(new_z_file, MODE)
    os.chmod(new_uv_file, MODE)

#-- PURPOSE: create a unique filename adding a numerical instance if existing
def create_unique_filename(filename):
    #-- split filename into fileBasename and fileExtension
    fileBasename, fileExtension = os.path.splitext(filename)
    #-- create counter to add to the end of the filename if existing
    counter = 1
    while counter:
        try:
            #-- open file descriptor only if the file doesn't exist
            fd = os.open(filename, os.O_CREAT | os.O_EXCL | os.O_RDWR)
        except OSError:
            pass
        else:
            #-- close the file descriptor and return the filename
            os.close(fd)
            return filename
        #-- new filename adds counter the between fileBasename and fileExtension
        filename = '{0}{1}_{2:d}'.format(fileBasename, fileExtension, counter)
        counter += 1

#-- PURPOSE: help module to describe the optional input parameters
def usage():
    print('\nHelp: {0}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tWorking data directory')
    print(' -T X, --tide=X\t\tTide model to use in correction')
    print(' -B X, --bounds=X\tGrid Bounds (xmin,xmax,ymin,ymax)')
    print(' -M X, --mode=X\t\tPermission mode of directories and files\n')

#-- This is the main part of the program that calls the individual modules
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','tide=','bounds=','mode=']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hD:T:B:M:',long_options)

    #-- command line parameters
    tide_dir = os.getcwd()
    BOUNDS = None
    #-- tide model to use
    TIDE_MODEL = 'TPXO9.1'
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-D","--directory"):
            tide_dir = os.path.expanduser(arg)
        elif opt in ("-T","--tide"):
            TIDE_MODEL = arg
        elif opt in ("-B","--bounds"):
            BOUNDS = np.array(arg.split(','),dtype=np.float)
        elif opt in ("-M","--mode"):
            MODE = int(arg,8)

    #-- verify model before running program
    model_list = ['CATS0201','CATS2008','CATS2008_load','TPXO9-atlas','TPXO9.1',
        'TPXO8-atlas','TPXO7.2','TPXO7.2_load','AODTM-5','AOTIM-5',
        'AOTIM-5-2018']
    assert TIDE_MODEL in model_list, 'Unlisted tide model'

    #-- run program
    make_regional_OTIS_files(tide_dir, TIDE_MODEL, BOUNDS, MODE=MODE)

#-- run main program
if __name__ == '__main__':
    main()
