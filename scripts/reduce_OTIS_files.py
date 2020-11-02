#!/usr/bin/env python
u"""
reduce_OTIS_files.py
Written by Tyler Sutterley (10/2020)
Read OTIS-format tidal files and reduce to a regional subset

COMMAND LINE OPTIONS:
    -D X, --directory X: working data directory
    -T X, --tide X: Tide model to use
    -B X, --bounds X: Grid Bounds (xmin,xmax,ymin,ymax)
    --projection X: spatial projection of bounds as EPSG code or PROJ4 string
        4326: latitude and longitude coordinates on WGS84 reference ellipsoid
    -M X, --mode X: permissions mode of the output files

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
    Updated 10/2020: using argparse to set command line parameters
    Updated 09/2020: can use projected coordinates for output model bounds
        compatibility updates for python3
    Updated 07/2020: renamed coordinate conversion program
    Updated 02/2020: changed CATS2008 grid to match version on U.S. Antarctic
        Program Data Center http://www.usap-dc.org/view/dataset/601235
    Updated 11/2019: added AOTIM-5-2018 tide model (2018 update to 2004 model)
    Written 08/2018
"""
from __future__ import print_function

import sys
import os
import pyproj
import argparse
import numpy as np
from pyTMD.convert_ll_xy import convert_ll_xy
from pyTMD.read_tide_model import *
from pyTMD.output_otis_tides import *

#-- PURPOSE: read 1 degree land sea mask and create masks for specific regions
def make_regional_OTIS_files(tide_dir, TIDE_MODEL, TYPE=['z','uv'],
    BOUNDS=4*[None], PROJECTION='4326', MODE=0o775):
    model_file = {}
    new_model_file = {}
    #-- select between tide models
    if (TIDE_MODEL == 'CATS0201'):
        grid_file = os.path.join(tide_dir,'cats0201_tmd','grid_CATS')
        model_file['z'] = os.path.join(tide_dir,'cats0201_tmd','h0_CATS02_01')
        model_file['uv'] = os.path.join(tide_dir,'cats0201_tmd','UV0_CATS02_01')
        model_format = 'OTIS'
        EPSG = '4326'
    elif (TIDE_MODEL == 'CATS2008'):
        grid_file = os.path.join(tide_dir,'CATS2008','grid_CATS2008')
        model_file['z'] = os.path.join(tide_dir,'CATS2008','hf.CATS2008.out')
        model_file['uv'] = os.path.join(tide_dir,'CATS2008','uv.CATS2008.out')
        model_format = 'OTIS'
        EPSG = 'CATS2008'
    elif (TIDE_MODEL == 'CATS2008_load'):
        grid_file = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','grid_CATS2008a_opt')
        model_file['z'] = os.path.join(tide_dir,'CATS2008a_SPOTL_Load','h_CATS2008a_SPOTL_load')
        model_format = 'OTIS'
        EPSG = 'CATS2008'
    elif (MODEL == 'TPXO9_atlas'):
        grid_file = os.path.join(tide_dir,'tpxo9_atlas','grid_tpxo9atlas_30_v1')
        model_file['z'] = os.path.join(tide_dir,'tpxo9_atlas','hf.tpxo9_atlas_30_v1')
        model_file['uv'] = os.path.join(tide_dir,'tpxo9_atlas','uv.tpxo9_atlas_30_v1')
        model_format = 'ATLAS'
        EPSG = '4326'
    elif (TIDE_MODEL == 'TPXO9.1'):
        grid_file = os.path.join(tide_dir,'TPXO9.1','DATA','grid_tpxo9')
        model_file['z'] = os.path.join(tide_dir,'TPXO9.1','DATA','h_tpxo9.v1')
        model_file['uv'] = os.path.join(tide_dir,'TPXO9.1','DATA','u_tpxo9.v1')
        model_format = 'OTIS'
        EPSG = '4326'
    elif (TIDE_MODEL == 'TPXO8-atlas'):
        grid_file = os.path.join(tide_dir,'tpxo8_atlas','grid_tpxo8atlas_30_v1')
        model_file['z'] = os.path.join(tide_dir,'tpxo8_atlas','hf.tpxo8_atlas_30_v1')
        model_file['uv'] = os.path.join(tide_dir,'tpxo8_atlas','uv.tpxo8_atlas_30_v1')
        model_format = 'ATLAS'
        EPSG = '4326'
    elif (TIDE_MODEL == 'TPXO7.2'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_tmd','grid_tpxo7.2')
        model_file['z'] = os.path.join(tide_dir,'TPXO7.2_tmd','h_tpxo7.2')
        model_file['uv'] = os.path.join(tide_dir,'TPXO7.2_tmd','u_tpxo7.2')
        model_format = 'OTIS'
        EPSG = '4326'
    elif (TIDE_MODEL == 'TPXO7.2_load'):
        grid_file = os.path.join(tide_dir,'TPXO7.2_load','grid_tpxo6.2')
        model_file['z'] = os.path.join(tide_dir,'TPXO7.2_load','h_tpxo7.2_load')
        model_format = 'OTIS'
        EPSG = '4326'
    elif (TIDE_MODEL == 'AODTM-5'):
        grid_file = os.path.join(tide_dir,'aodtm5_tmd','grid_Arc5km')
        model_file['z'] = os.path.join(tide_dir,'aodtm5_tmd','h0_Arc5km.oce')
        model_file['uv'] = os.path.join(tide_dir,'aodtm5_tmd','u0_Arc5km.oce')
        model_format = 'OTIS'
        EPSG = 'PSNorth'
    elif (TIDE_MODEL == 'AOTIM-5'):
        grid_file = os.path.join(tide_dir,'aotim5_tmd','grid_Arc5km')
        model_file['z'] = os.path.join(tide_dir,'aotim5_tmd','h_Arc5km.oce')
        model_file['uv'] = os.path.join(tide_dir,'aotim5_tmd','u_Arc5km.oce')
        model_format = 'OTIS'
        EPSG = 'PSNorth'
    elif (TIDE_MODEL == 'AOTIM-5-2018'):
        grid_file = os.path.join(tide_dir,'aotim5_tmd','grid_Arc5km')
        model_file['z'] = os.path.join(tide_dir,'aotim5_tmd','h_Arc5km2018')
        model_file['uv'] = os.path.join(tide_dir,'aotim5_tmd','u_Arc5km2018_T')
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

    #-- converting bounds x,y from projection to latitude/longitude
    try:
        crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(int(PROJECTION)))
    except (ValueError,pyproj.exceptions.CRSError):
        crs1 = pyproj.CRS.from_string(PROJECTION)
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    xbox = np.array([BOUNDS[0],BOUNDS[1],BOUNDS[1],BOUNDS[0],BOUNDS[0]])
    ybox = np.array([BOUNDS[2],BOUNDS[2],BOUNDS[3],BOUNDS[3],BOUNDS[2]])
    lon,lat = transformer.transform(xbox,ybox)

    #-- convert bounds from latitude/longitude to model coordinates
    x,y = convert_ll_xy(lon,lat,EPSG,'F')

    #-- find indices to reduce to xmin,xmax,ymin,ymax
    gridx,gridy = np.meshgrid(xi,yi)
    indy,indx = np.nonzero((gridx >= x.min()) & (gridx <= x.max()) &
        (gridy >= y.min()) & (gridy <= y.max()))
    nx = np.count_nonzero((xi >= x.min()) & (xi <= x.max()))
    ny = np.count_nonzero((yi >= y.min()) & (yi <= y.max()))
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
    new_grid_file = create_unique_filename(grid_file)
    output_otis_grid(new_grid_file,xlim,ylim,hz1,mz1,iob,dt)
    #-- change the permissions level to MODE
    os.chmod(new_grid_file, MODE)

    #-- combine ATLAS sub-grids into single output grid
    #-- reduce elevation files to bounds
    if 'z' in model_file.keys():
        #-- read each constituent
        constituents,nc = read_constituents(model_file['z'])
        z1 = np.zeros((ny,nx,nc),dtype=np.complex64)
        for i,c in enumerate(constituents):
            #-- read constituent from elevation file
            if (model_format == 'ATLAS'):
                z0,zlocal=read_atlas_elevation(model_file['z'],i,c)
                xi,yi,z=combine_atlas_model(x0,y0,z0,pmask,zlocal,VARIABLE='z')
            else:
                z=read_elevation_file(model_file['z'],i)
            #-- reduce elevation to new bounds
            z1[:,:,i] = z[indy,indx].reshape(ny,nx)
        #-- output reduced elevation components
        new_model_file['z'] = create_unique_filename(model_file['z'])
        output_otis_elevation(new_model_file['z'],z1,xlim,ylim,constituents)
        #-- change the permissions level to MODE
        os.chmod(new_model_file['z'], MODE)

    #-- combine ATLAS sub-grids into single output grid
    #-- reduce transport files to bounds
    if 'uv' in model_file.keys():
        #-- read each constituent
        constituents,nc = read_constituents(model_file['uv'])
        u1 = np.zeros((ny,nx,nc),dtype=np.complex64)
        v1 = np.zeros((ny,nx,nc),dtype=np.complex64)
        for i,c in enumerate(constituents):
            #-- read constituent from transport file
            if (model_format == 'ATLAS'):
                u0,v0,uvlocal=read_atlas_transport(model_file['uv'],i,c)
                xi,yi,u=combine_atlas_model(x0,y0,u0,pmask,uvlocal,VARIABLE='u')
                xi,yi,v=combine_atlas_model(x0,y0,v0,pmask,uvlocal,VARIABLE='v')
            else:
                u,v=read_transport_file(model_file['uv'],i)
            #-- reduce transport components to new bounds
            u1[:,:,i] = u[indy,indx].reshape(ny,nx)
            v1[:,:,i] = v[indy,indx].reshape(ny,nx)
        #-- output reduced transport components
        new_model_file['uv'] = create_unique_filename(model_file['uv'])
        output_otis_transport(new_model_file['uv'],u1,v1,xlim,ylim,constituents)
        #-- change the permissions level to MODE
        os.chmod(new_model_file['uv'], MODE)

#-- PURPOSE: create a unique filename adding a numerical instance if existing
def create_unique_filename(filename):
    #-- split filename into fileBasename and fileExtension
    fileBasename, fileExtension = os.path.splitext(filename)
    fileExtension = '' if (fileExtension in ('.out','.oce')) else fileExtension
    #-- replace extension with reduced flag
    filename = '{0}{1}{2}'.format(fileBasename, fileExtension, '.reduced')
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
        #-- new filename adds counter
        args = (fileBasename, fileExtension, '.reduced', counter)
        filename = '{0}{1}{2}_{3:d}'.format(*args)
        counter += 1

#-- This is the main part of the program that calls the individual modules
def main():
    #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Read OTIS-format tidal files and reduce to a regional
            subset
            """
    )
    #-- command line options
    #-- set data directory containing the tidal data
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- tide model to use
    model_choices = ('CATS0201','CATS2008','CATS2008_load','TPXO9-atlas',
        'TPXO9.1','TPXO8-atlas','TPXO7.2','TPXO7.2_load','AODTM-5','AOTIM-5',
        'AOTIM-5-2018')
    parser.add_argument('--tide','-T',
        metavar='TIDE', type=str, default='TPXO9.1',
        choices=model_choices,
        help='Tide model to use')
    #-- spatial projection (EPSG code or PROJ4 string)
    parser.add_argument('--projection','-P',
        type=str, default='4326',
        help='Spatial projection as EPSG code or PROJ4 string')
    #-- bounds for reducing model (xmin,xmax,ymin,ymax)
    parser.add_argument('--bounds','-B',
        metavar=('xmin','xmax','ymin','ymax'), type=float, nargs=4,
        help='Grid bounds for reducing model')
    #-- permissions mode of output reduced files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of the output files')
    args = parser.parse_args()

    #-- run program
    make_regional_OTIS_files(args.directory, args.tide, BOUNDS=args.bounds,
        PROJECTION=args.projection, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
