"""
test_model.py (04/2023)
Tests the reading of model definition files
"""
from __future__ import annotations

import io
import pytest
import shutil
import inspect
import pathlib
import pyTMD.io

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent

def test_definition_CATS2008():
    """Tests the reading of the CATS2008 model definition file
    """
    m = pyTMD.io.model().from_file(filepath.joinpath('model_CATS2008.def'))
    # test read variables
    assert m.format == 'OTIS'
    assert m.name == 'CATS2008'
    assert m.model_file == pathlib.Path('CATS2008/hf.CATS2008.out')
    assert m.grid_file == pathlib.Path('CATS2008/grid_CATS2008')
    assert m.projection == 'CATS2008'
    assert m.type == 'z'
    assert m.variable == 'tide_ocean'
    # test properties
    assert m.atl03 == 'tide_ocean'
    assert m.atl06 == 'tide_ocean'
    assert m.atl07 == 'height_segment_ocean'
    assert m.atl10 == 'height_segment_ocean'
    assert m.atl11 == 'tide_ocean'
    assert m.atl12 == 'tide_ocean_seg'
    assert m.gla12 == 'd_ocElv'
    assert m.long_name == 'ocean_tide_elevation'

def test_definition_FES():
    """Tests the reading of the FES2014 model definition file
    """
    m = pyTMD.io.model().from_file(filepath.joinpath('model_FES2014.def'))
    # test read variables
    assert m.format == 'FES'
    assert m.name == 'FES2014'
    model_files = ['fes2014/ocean_tide/2n2.nc.gz',
        'fes2014/ocean_tide/eps2.nc.gz', 'fes2014/ocean_tide/j1.nc.gz',
        'fes2014/ocean_tide/k1.nc.gz', 'fes2014/ocean_tide/k2.nc.gz',
        'fes2014/ocean_tide/l2.nc.gz', 'fes2014/ocean_tide/la2.nc.gz',
        'fes2014/ocean_tide/m2.nc.gz', 'fes2014/ocean_tide/m3.nc.gz',
        'fes2014/ocean_tide/m4.nc.gz', 'fes2014/ocean_tide/m6.nc.gz',
        'fes2014/ocean_tide/m8.nc.gz', 'fes2014/ocean_tide/mf.nc.gz',
        'fes2014/ocean_tide/mks2.nc.gz', 'fes2014/ocean_tide/mm.nc.gz',
        'fes2014/ocean_tide/mn4.nc.gz', 'fes2014/ocean_tide/ms4.nc.gz',
        'fes2014/ocean_tide/msf.nc.gz', 'fes2014/ocean_tide/msqm.nc.gz',
        'fes2014/ocean_tide/mtm.nc.gz', 'fes2014/ocean_tide/mu2.nc.gz',
        'fes2014/ocean_tide/n2.nc.gz', 'fes2014/ocean_tide/n4.nc.gz',
        'fes2014/ocean_tide/nu2.nc.gz', 'fes2014/ocean_tide/o1.nc.gz',
        'fes2014/ocean_tide/p1.nc.gz', 'fes2014/ocean_tide/q1.nc.gz',
        'fes2014/ocean_tide/r2.nc.gz', 'fes2014/ocean_tide/s1.nc.gz',
        'fes2014/ocean_tide/s2.nc.gz', 'fes2014/ocean_tide/s4.nc.gz',
        'fes2014/ocean_tide/sa.nc.gz', 'fes2014/ocean_tide/ssa.nc.gz',
        'fes2014/ocean_tide/t2.nc.gz']
    # assert that all model files are in the model definition
    for f in model_files:
        assert pathlib.Path(f) in m.model_file
    # assert that all constituents are in the model definition
    constituents = ['2n2','eps2','j1','k1','k2','l2',
                'lambda2','m2','m3','m4','m6','m8','mf','mks2','mm',
                'mn4','ms4','msf','msqm','mtm','mu2','n2','n4','nu2',
                'o1','p1','q1','r2','s1','s2','s4','sa','ssa','t2']
    assert m.constituents == constituents
    assert m.type == 'z'
    assert m.scale == 1.0/100.0
    assert m.variable == 'tide_ocean'
    assert m.compressed is True
    # check validity of parsed constituents
    parsed_constituents = [pyTMD.io.model.parse_file(f) for f in model_files]
    assert parsed_constituents == constituents
    # test derived properties
    assert m.atl03 == 'tide_ocean'
    assert m.atl06 == 'tide_ocean'
    assert m.atl07 == 'height_segment_ocean'
    assert m.atl10 == 'height_segment_ocean'
    assert m.atl11 == 'tide_ocean'
    assert m.atl12 == 'tide_ocean_seg'
    assert m.gla12 == 'd_ocElv'
    assert m.long_name == 'ocean_tide_elevation'

# PURPOSE: test glob file functionality
def test_definition_FES_glob():
    """Tests the reading of the FES2014 model definition file
    with glob file searching
    """
    # get model parameters
    m = pyTMD.io.model().from_file(filepath.joinpath('model_FES2014.def'))
    model_files = ['fes2014/ocean_tide/2n2.nc.gz',
        'fes2014/ocean_tide/eps2.nc.gz', 'fes2014/ocean_tide/j1.nc.gz',
        'fes2014/ocean_tide/k1.nc.gz', 'fes2014/ocean_tide/k2.nc.gz',
        'fes2014/ocean_tide/l2.nc.gz', 'fes2014/ocean_tide/la2.nc.gz',
        'fes2014/ocean_tide/m2.nc.gz', 'fes2014/ocean_tide/m3.nc.gz',
        'fes2014/ocean_tide/m4.nc.gz', 'fes2014/ocean_tide/m6.nc.gz',
        'fes2014/ocean_tide/m8.nc.gz', 'fes2014/ocean_tide/mf.nc.gz',
        'fes2014/ocean_tide/mks2.nc.gz', 'fes2014/ocean_tide/mm.nc.gz',
        'fes2014/ocean_tide/mn4.nc.gz', 'fes2014/ocean_tide/ms4.nc.gz',
        'fes2014/ocean_tide/msf.nc.gz', 'fes2014/ocean_tide/msqm.nc.gz',
        'fes2014/ocean_tide/mtm.nc.gz', 'fes2014/ocean_tide/mu2.nc.gz',
        'fes2014/ocean_tide/n2.nc.gz', 'fes2014/ocean_tide/n4.nc.gz',
        'fes2014/ocean_tide/nu2.nc.gz', 'fes2014/ocean_tide/o1.nc.gz',
        'fes2014/ocean_tide/p1.nc.gz', 'fes2014/ocean_tide/q1.nc.gz',
        'fes2014/ocean_tide/r2.nc.gz', 'fes2014/ocean_tide/s1.nc.gz',
        'fes2014/ocean_tide/s2.nc.gz', 'fes2014/ocean_tide/s4.nc.gz',
        'fes2014/ocean_tide/sa.nc.gz', 'fes2014/ocean_tide/ssa.nc.gz',
        'fes2014/ocean_tide/t2.nc.gz']
    # create temporary files for testing glob functionality
    for model_file in model_files:
        local = filepath.joinpath(model_file)
        local.parent.mkdir(parents=True, exist_ok=True)
        local.touch(exist_ok=True)
    # create model definition file
    fid = io.StringIO()
    attrs = ['name','format','compressed','type','scale','version']
    for attr in attrs:
        val = getattr(m,attr)
        if isinstance(val,list):
            fid.write('{0}\t{1}\n'.format(attr,','.join(val)))
        else:
            fid.write('{0}\t{1}\n'.format(attr,val))
    # append glob strings for model file
    glob_string = r'fes2014/ocean_tide/*.nc.gz'
    fid.write('{0}\t{1}\n'.format('model_file',glob_string))
    fid.seek(0)
    # use model definition file as input
    model = pyTMD.io.model(directory=filepath).from_file(fid)
    for attr in attrs:
        assert getattr(model,attr) == getattr(m,attr)
    # verify that the model files and constituents match
    assert (len(model.model_file) == len(model_files))
    for f in model_files:
        assert pathlib.Path(filepath).joinpath(f) in model.model_file
    for c in m.constituents:
        assert c in model.constituents
    # close the glob definition file
    fid.close()
    # clean up model
    shutil.rmtree(filepath.joinpath('fes2014'))

def test_definition_FES_currents():
    """Tests the reading of the FES2014 model definition file for currents
    """
    m = pyTMD.io.model().from_file(filepath.joinpath('model_FES2014_currents.def'))
    # test read variables
    assert m.format == 'FES'
    assert m.name == 'FES2014'
    model_files_u = ['fes2014/eastward_velocity/2n2.nc.gz',
        'fes2014/eastward_velocity/eps2.nc.gz', 'fes2014/eastward_velocity/j1.nc.gz',
        'fes2014/eastward_velocity/k1.nc.gz', 'fes2014/eastward_velocity/k2.nc.gz',
        'fes2014/eastward_velocity/l2.nc.gz', 'fes2014/eastward_velocity/la2.nc.gz',
        'fes2014/eastward_velocity/m2.nc.gz', 'fes2014/eastward_velocity/m3.nc.gz',
        'fes2014/eastward_velocity/m4.nc.gz', 'fes2014/eastward_velocity/m6.nc.gz',
        'fes2014/eastward_velocity/m8.nc.gz', 'fes2014/eastward_velocity/mf.nc.gz',
        'fes2014/eastward_velocity/mks2.nc.gz', 'fes2014/eastward_velocity/mm.nc.gz',
        'fes2014/eastward_velocity/mn4.nc.gz', 'fes2014/eastward_velocity/ms4.nc.gz',
        'fes2014/eastward_velocity/msf.nc.gz', 'fes2014/eastward_velocity/msqm.nc.gz',
        'fes2014/eastward_velocity/mtm.nc.gz', 'fes2014/eastward_velocity/mu2.nc.gz',
        'fes2014/eastward_velocity/n2.nc.gz', 'fes2014/eastward_velocity/n4.nc.gz',
        'fes2014/eastward_velocity/nu2.nc.gz', 'fes2014/eastward_velocity/o1.nc.gz',
        'fes2014/eastward_velocity/p1.nc.gz', 'fes2014/eastward_velocity/q1.nc.gz',
        'fes2014/eastward_velocity/r2.nc.gz', 'fes2014/eastward_velocity/s1.nc.gz',
        'fes2014/eastward_velocity/s2.nc.gz', 'fes2014/eastward_velocity/s4.nc.gz',
        'fes2014/eastward_velocity/sa.nc.gz', 'fes2014/eastward_velocity/ssa.nc.gz',
        'fes2014/eastward_velocity/t2.nc.gz']
    model_files_v = ['fes2014/northward_velocity/2n2.nc.gz',
        'fes2014/northward_velocity/eps2.nc.gz', 'fes2014/northward_velocity/j1.nc.gz',
        'fes2014/northward_velocity/k1.nc.gz', 'fes2014/northward_velocity/k2.nc.gz',
        'fes2014/northward_velocity/l2.nc.gz', 'fes2014/northward_velocity/la2.nc.gz',
        'fes2014/northward_velocity/m2.nc.gz', 'fes2014/northward_velocity/m3.nc.gz',
        'fes2014/northward_velocity/m4.nc.gz', 'fes2014/northward_velocity/m6.nc.gz',
        'fes2014/northward_velocity/m8.nc.gz', 'fes2014/northward_velocity/mf.nc.gz',
        'fes2014/northward_velocity/mks2.nc.gz', 'fes2014/northward_velocity/mm.nc.gz',
        'fes2014/northward_velocity/mn4.nc.gz', 'fes2014/northward_velocity/ms4.nc.gz',
        'fes2014/northward_velocity/msf.nc.gz', 'fes2014/northward_velocity/msqm.nc.gz',
        'fes2014/northward_velocity/mtm.nc.gz', 'fes2014/northward_velocity/mu2.nc.gz',
        'fes2014/northward_velocity/n2.nc.gz', 'fes2014/northward_velocity/n4.nc.gz',
        'fes2014/northward_velocity/nu2.nc.gz', 'fes2014/northward_velocity/o1.nc.gz',
        'fes2014/northward_velocity/p1.nc.gz', 'fes2014/northward_velocity/q1.nc.gz',
        'fes2014/northward_velocity/r2.nc.gz', 'fes2014/northward_velocity/s1.nc.gz',
        'fes2014/northward_velocity/s2.nc.gz', 'fes2014/northward_velocity/s4.nc.gz',
        'fes2014/northward_velocity/sa.nc.gz', 'fes2014/northward_velocity/ssa.nc.gz',
        'fes2014/northward_velocity/t2.nc.gz']
    # assert that all model files are in the model definition
    for f in model_files_u:
        assert pathlib.Path(f) in m.model_file['u']
    for f in model_files_v:
        assert pathlib.Path(f) in m.model_file['v']
    # assert that all constituents are in the model definition
    constituents = ['2n2','eps2','j1','k1','k2','l2',
                'lambda2','m2','m3','m4','m6','m8','mf','mks2','mm',
                'mn4','ms4','msf','msqm','mtm','mu2','n2','n4','nu2',
                'o1','p1','q1','r2','s1','s2','s4','sa','ssa','t2']
    assert m.constituents == constituents
    assert m.type == ['u','v']
    assert m.scale == 1.0
    assert m.compressed is True
    # check validity of parsed constituents
    parsed_constituents = [pyTMD.io.model.parse_file(f) for f in model_files_u]
    assert parsed_constituents == constituents
    # test derived properties
    assert m.long_name['u'] == 'zonal_tidal_current'
    assert m.long_name['v'] == 'meridional_tidal_current'

# PURPOSE: test glob file functionality
def test_definition_FES_currents_glob():
    """Tests the reading of the FES2014 model definition file
    with glob file searching for currents
    """
    # get model parameters
    m = pyTMD.io.model().from_file(filepath.joinpath('model_FES2014_currents.def'))
    model_files_u = ['fes2014/eastward_velocity/2n2.nc.gz',
        'fes2014/eastward_velocity/eps2.nc.gz', 'fes2014/eastward_velocity/j1.nc.gz',
        'fes2014/eastward_velocity/k1.nc.gz', 'fes2014/eastward_velocity/k2.nc.gz',
        'fes2014/eastward_velocity/l2.nc.gz', 'fes2014/eastward_velocity/la2.nc.gz',
        'fes2014/eastward_velocity/m2.nc.gz', 'fes2014/eastward_velocity/m3.nc.gz',
        'fes2014/eastward_velocity/m4.nc.gz', 'fes2014/eastward_velocity/m6.nc.gz',
        'fes2014/eastward_velocity/m8.nc.gz', 'fes2014/eastward_velocity/mf.nc.gz',
        'fes2014/eastward_velocity/mks2.nc.gz', 'fes2014/eastward_velocity/mm.nc.gz',
        'fes2014/eastward_velocity/mn4.nc.gz', 'fes2014/eastward_velocity/ms4.nc.gz',
        'fes2014/eastward_velocity/msf.nc.gz', 'fes2014/eastward_velocity/msqm.nc.gz',
        'fes2014/eastward_velocity/mtm.nc.gz', 'fes2014/eastward_velocity/mu2.nc.gz',
        'fes2014/eastward_velocity/n2.nc.gz', 'fes2014/eastward_velocity/n4.nc.gz',
        'fes2014/eastward_velocity/nu2.nc.gz', 'fes2014/eastward_velocity/o1.nc.gz',
        'fes2014/eastward_velocity/p1.nc.gz', 'fes2014/eastward_velocity/q1.nc.gz',
        'fes2014/eastward_velocity/r2.nc.gz', 'fes2014/eastward_velocity/s1.nc.gz',
        'fes2014/eastward_velocity/s2.nc.gz', 'fes2014/eastward_velocity/s4.nc.gz',
        'fes2014/eastward_velocity/sa.nc.gz', 'fes2014/eastward_velocity/ssa.nc.gz',
        'fes2014/eastward_velocity/t2.nc.gz']
    model_files_v = ['fes2014/northward_velocity/2n2.nc.gz',
        'fes2014/northward_velocity/eps2.nc.gz', 'fes2014/northward_velocity/j1.nc.gz',
        'fes2014/northward_velocity/k1.nc.gz', 'fes2014/northward_velocity/k2.nc.gz',
        'fes2014/northward_velocity/l2.nc.gz', 'fes2014/northward_velocity/la2.nc.gz',
        'fes2014/northward_velocity/m2.nc.gz', 'fes2014/northward_velocity/m3.nc.gz',
        'fes2014/northward_velocity/m4.nc.gz', 'fes2014/northward_velocity/m6.nc.gz',
        'fes2014/northward_velocity/m8.nc.gz', 'fes2014/northward_velocity/mf.nc.gz',
        'fes2014/northward_velocity/mks2.nc.gz', 'fes2014/northward_velocity/mm.nc.gz',
        'fes2014/northward_velocity/mn4.nc.gz', 'fes2014/northward_velocity/ms4.nc.gz',
        'fes2014/northward_velocity/msf.nc.gz', 'fes2014/northward_velocity/msqm.nc.gz',
        'fes2014/northward_velocity/mtm.nc.gz', 'fes2014/northward_velocity/mu2.nc.gz',
        'fes2014/northward_velocity/n2.nc.gz', 'fes2014/northward_velocity/n4.nc.gz',
        'fes2014/northward_velocity/nu2.nc.gz', 'fes2014/northward_velocity/o1.nc.gz',
        'fes2014/northward_velocity/p1.nc.gz', 'fes2014/northward_velocity/q1.nc.gz',
        'fes2014/northward_velocity/r2.nc.gz', 'fes2014/northward_velocity/s1.nc.gz',
        'fes2014/northward_velocity/s2.nc.gz', 'fes2014/northward_velocity/s4.nc.gz',
        'fes2014/northward_velocity/sa.nc.gz', 'fes2014/northward_velocity/ssa.nc.gz',
        'fes2014/northward_velocity/t2.nc.gz']
    # create temporary files for testing glob functionality
    for model_file in model_files_u:
        local = filepath.joinpath(model_file)
        local.parent.mkdir(parents=True, exist_ok=True)
        local.touch(exist_ok=True)
    for model_file in model_files_v:
        local = filepath.joinpath(model_file)
        local.parent.mkdir(parents=True, exist_ok=True)
        local.touch(exist_ok=True)
    # create model definition file
    fid = io.StringIO()
    attrs = ['name','format','compressed','type','scale','version']
    for attr in attrs:
        val = getattr(m,attr)
        if isinstance(val,list):
            fid.write('{0}\t{1}\n'.format(attr,','.join(val)))
        else:
            fid.write('{0}\t{1}\n'.format(attr,val))
    # append glob strings for model file
    eastward = r'fes2014/eastward_velocity/*.nc.gz'
    northward = r'fes2014/northward_velocity/*.nc.gz'
    fid.write('{0}\t{1};{2}\n'.format('model_file',eastward,northward))
    fid.seek(0)
    # use model definition file as input
    model = pyTMD.io.model(directory=filepath).from_file(fid)
    for attr in attrs:
        assert getattr(model,attr) == getattr(m,attr)
    # verify that the model files and constituents match
    assert (len(model.model_file['u']) == len(model_files_u))
    assert (len(model.model_file['v']) == len(model_files_v))
    for f in model_files_u:
        assert pathlib.Path(filepath).joinpath(f) in model.model_file['u']
    for f in model_files_v:
        assert pathlib.Path(filepath).joinpath(f) in model.model_file['v']
    for c in m.constituents:
        assert c in model.constituents
    # close the glob definition file
    fid.close()
    # clean up model
    shutil.rmtree(filepath.joinpath('fes2014'))

def test_definition_GOT():
    """Tests the reading of the GOT4.10 model definition file
    """
    m = pyTMD.io.model().from_file(filepath.joinpath('model_GOT4.10.def'))
    # test read variables
    assert m.format == 'GOT'
    assert m.name == 'GOT4.10'
    model_files = ['GOT4.10c/grids_loadtide/k1load.d.gz',
        'GOT4.10c/grids_loadtide/k2load.d.gz',
        'GOT4.10c/grids_loadtide/m2load.d.gz',
        'GOT4.10c/grids_loadtide/m4load.d.gz',
        'GOT4.10c/grids_loadtide/n2load.d.gz',
        'GOT4.10c/grids_loadtide/o1load.d.gz',
        'GOT4.10c/grids_loadtide/p1load.d.gz',
        'GOT4.10c/grids_loadtide/q1load.d.gz',
        'GOT4.10c/grids_loadtide/s1load.d.gz',
        'GOT4.10c/grids_loadtide/s2load.d.gz']
    # assert that all model files are in the model definition
    for f in model_files:
        assert pathlib.Path(f) in m.model_file
    assert m.type == 'z'
    assert m.scale == 1.0/1000.0
    assert m.variable == 'tide_load'
    assert m.compressed is True
    # test derived properties
    assert m.atl03 == 'tide_load'
    assert m.atl06 == 'tide_load'
    assert m.atl07 == 'height_segment_load'
    assert m.atl10 == 'height_segment_load'
    assert m.atl11 == 'tide_load'
    assert m.atl12 == 'tide_load_seg'
    assert m.gla12 == 'd_ldElv'
    assert m.long_name == 'load_tide_elevation'

# PURPOSE: test glob file functionality
def test_definition_GOT_glob():
    """Tests the reading of the GOT4.10 model definition file
    with glob file searching
    """
    # get model parameters
    m = pyTMD.io.model().from_file(filepath.joinpath('model_GOT4.10.def'))
    model_files = ['GOT4.10c/grids_loadtide/k1load.d.gz',
        'GOT4.10c/grids_loadtide/k2load.d.gz',
        'GOT4.10c/grids_loadtide/m2load.d.gz',
        'GOT4.10c/grids_loadtide/m4load.d.gz',
        'GOT4.10c/grids_loadtide/n2load.d.gz',
        'GOT4.10c/grids_loadtide/o1load.d.gz',
        'GOT4.10c/grids_loadtide/p1load.d.gz',
        'GOT4.10c/grids_loadtide/q1load.d.gz',
        'GOT4.10c/grids_loadtide/s1load.d.gz',
        'GOT4.10c/grids_loadtide/s2load.d.gz']
    # create temporary files for testing glob functionality
    for model_file in model_files:
        local = filepath.joinpath(model_file)
        local.parent.mkdir(parents=True, exist_ok=True)
        local.touch(exist_ok=True)
    # create model definition file
    fid = io.StringIO()
    attrs = ['name','format','compressed','type','scale']
    for attr in attrs:
        val = getattr(m,attr)
        if isinstance(val,list):
            fid.write('{0}\t{1}\n'.format(attr,','.join(val)))
        else:
            fid.write('{0}\t{1}\n'.format(attr,val))
    # append glob strings for model file
    glob_string = r'GOT4.10c/grids_loadtide/*.d.gz'
    fid.write('{0}\t{1}\n'.format('model_file',glob_string))
    fid.seek(0)
    # use model definition file as input
    model = pyTMD.io.model(directory=filepath).from_file(fid)
    for attr in attrs:
        assert getattr(model,attr) == getattr(m,attr)
    # verify that the model files match
    assert (len(model.model_file) == len(model_files))
    for f in model_files:
        assert pathlib.Path(filepath).joinpath(f) in model.model_file
    # close the glob definition file
    fid.close()
    # clean up model
    shutil.rmtree(filepath.joinpath('GOT4.10c'))

def test_definition_TPXO9():
    """Tests the reading of the TPXO9-atlas-v5 model definition file
    """
    m = pyTMD.io.model().from_file(filepath.joinpath('model_TPXO9-atlas-v5.def'))
    # test read variables
    assert m.format == 'netcdf'
    assert m.name == 'TPXO9-atlas-v5'
    model_files = ['TPXO9_atlas_v5/h_2n2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_k1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_k2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_m2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_m4_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_mf_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_mm_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_mn4_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_ms4_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_n2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_o1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_p1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_q1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_s1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_s2_tpxo9_atlas_30_v5.nc']
    assert m.grid_file == pathlib.Path('TPXO9_atlas_v5/grid_tpxo9_atlas_30_v5.nc')
    # assert that all model files are in the model definition
    for f in model_files:
        assert pathlib.Path(f) in m.model_file
    assert m.type == 'z'
    assert m.scale == 1.0/100.0
    assert m.variable == 'tide_ocean'
    assert m.compressed is False
    # test derived properties
    assert m.atl03 == 'tide_ocean'
    assert m.atl06 == 'tide_ocean'
    assert m.atl07 == 'height_segment_ocean'
    assert m.atl10 == 'height_segment_ocean'
    assert m.atl11 == 'tide_ocean'
    assert m.atl12 == 'tide_ocean_seg'
    assert m.gla12 == 'd_ocElv'
    assert m.long_name == 'ocean_tide_elevation'

# PURPOSE: test glob file functionality
def test_definition_TPXO9_glob():
    """Tests the reading of the TPXO9-atlas-v5 model definition file
    with glob file searching
    """
    m = pyTMD.io.model().from_file(filepath.joinpath('model_TPXO9-atlas-v5.def'))
    # test read variables
    assert m.format == 'netcdf'
    assert m.name == 'TPXO9-atlas-v5'
    model_files = ['TPXO9_atlas_v5/h_2n2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_k1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_k2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_m2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_m4_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_mf_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_mm_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_mn4_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_ms4_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_n2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_o1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_p1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_q1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_s1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/h_s2_tpxo9_atlas_30_v5.nc']
    # create temporary files for testing glob functionality
    for model_file in model_files:
        local = filepath.joinpath(model_file)
        local.parent.mkdir(parents=True, exist_ok=True)
        local.touch(exist_ok=True)
    # create temporary grid file
    local = filepath.joinpath(m.grid_file)
    local.touch(exist_ok=True)
    # create model definition file
    fid = io.StringIO()
    attrs = ['name','format','compressed','type','scale']
    for attr in attrs:
        val = getattr(m,attr)
        if isinstance(val,list):
            fid.write('{0}\t{1}\n'.format(attr,','.join(val)))
        else:
            fid.write('{0}\t{1}\n'.format(attr,val))
    # append glob strings for model file
    glob_string = r'TPXO9_atlas_v5/h*.nc'
    fid.write('{0}\t{1}\n'.format('model_file',glob_string))
    fid.write('{0}\t{1}\n'.format('grid_file',m.grid_file))
    fid.seek(0)
    # use model definition file as input
    model = pyTMD.io.model(directory=filepath).from_file(fid)
    for attr in attrs:
        assert getattr(model,attr) == getattr(m,attr)
    # verify that the model files match
    assert (len(model.model_file) == len(model_files))
    for f in model_files:
        assert pathlib.Path(filepath).joinpath(f) in model.model_file
    # close the glob definition file
    fid.close()
    # clean up model
    shutil.rmtree(filepath.joinpath('TPXO9_atlas_v5'))

def test_definition_TPXO9_currents():
    """Tests the reading of the TPXO9-atlas-v5 model definition file for currents
    """
    m = pyTMD.io.model().from_file(filepath.joinpath('model_TPXO9-atlas-v5_currents.def'))
    # test read variables
    assert m.format == 'netcdf'
    assert m.name == 'TPXO9-atlas-v5'
    model_files = ['TPXO9_atlas_v5/u_2n2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_k1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_k2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_m2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_m4_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_mf_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_mm_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_mn4_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_ms4_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_n2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_o1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_p1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_q1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_s1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_s2_tpxo9_atlas_30_v5.nc']
    assert m.grid_file == pathlib.Path('TPXO9_atlas_v5/grid_tpxo9_atlas_30_v5.nc')
    assert sorted(m.model_file['u']) == [pathlib.Path(f) for f in model_files]
    assert m.type == ['u', 'v']
    assert m.scale == 1.0/100.0
    assert m.compressed is False
    # test derived properties
    assert m.long_name['u'] == 'zonal_tidal_current'
    assert m.long_name['v'] == 'meridional_tidal_current'

# PURPOSE: test glob file functionality
def test_definition_TPXO9_glob():
    """Tests the reading of the TPXO9-atlas-v5 model definition file for currents
    with glob file searching
    """
    m = pyTMD.io.model().from_file(filepath.joinpath('model_TPXO9-atlas-v5_currents.def'))
    model_files = ['TPXO9_atlas_v5/u_2n2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_k1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_k2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_m2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_m4_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_mf_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_mm_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_mn4_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_ms4_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_n2_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_o1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_p1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_q1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_s1_tpxo9_atlas_30_v5.nc',
        'TPXO9_atlas_v5/u_s2_tpxo9_atlas_30_v5.nc']
    # create temporary files for testing glob functionality
    for model_file in model_files:
        local = filepath.joinpath(model_file)
        local.parent.mkdir(parents=True, exist_ok=True)
        local.touch(exist_ok=True)
    # create temporary grid file
    local = filepath.joinpath(m.grid_file)
    local.touch(exist_ok=True)
    # create model definition file
    fid = io.StringIO()
    attrs = ['name','format','compressed','type','scale']
    for attr in attrs:
        val = getattr(m,attr)
        if isinstance(val,list):
            fid.write('{0}\t{1}\n'.format(attr,','.join(val)))
        else:
            fid.write('{0}\t{1}\n'.format(attr,val))
    # append glob strings for model file
    glob_string = r'TPXO9_atlas_v5/u*.nc'
    fid.write('{0}\t{1}\n'.format('model_file',glob_string))
    fid.write('{0}\t{1}\n'.format('grid_file',m.grid_file))
    fid.seek(0)
    # use model definition file as input
    model = pyTMD.io.model(directory=filepath).from_file(fid)
    for attr in attrs:
        assert getattr(model,attr) == getattr(m,attr)
    # verify that the model files match
    for key,val in model.model_file.items():
        assert (len(val) == len(model_files))
        for f in model_files:
            assert pathlib.Path(filepath).joinpath(f) in model.model_file[key]
    # close the glob definition file
    fid.close()
    # clean up model
    shutil.rmtree(filepath.joinpath('TPXO9_atlas_v5'))

# parameterize model
@pytest.mark.parametrize("MODEL", pyTMD.io.model.FES())
def test_parse_FES_elevation(MODEL):
    """Tests the parsing of FES-type elevation model files
    """
    m = pyTMD.io.model(verify=False).elevation(MODEL)
    constituents = [pyTMD.io.model.parse_file(f) for f in m.model_file]
    assert (m.constituents == constituents)

# parameterize model
current_models = set(pyTMD.io.model.FES()) & set(pyTMD.io.model.ocean_current())
@pytest.mark.parametrize("MODEL", sorted(current_models))
def test_parse_FES_currents(MODEL):
    """Tests the parsing of FES-type current model files
    """
    # test ocean current constituents
    m = pyTMD.io.model(verify=False).current(MODEL)
    constituents = [pyTMD.io.model.parse_file(f) for f in m.model_file['u']]
    assert (m.constituents == constituents)
    constituents = [pyTMD.io.model.parse_file(f) for f in m.model_file['v']]
    assert (m.constituents == constituents)

# parameterize model
@pytest.mark.parametrize("MODEL", pyTMD.io.model.GOT())
def test_parse_GOT_elevation(MODEL):
    """Tests the parsing of GOT-type elevation model files
    """
    m = pyTMD.io.model(verify=False).elevation(MODEL)
    m.constituents = [pyTMD.io.model.parse_file(f) for f in m.model_file]
    constituents = ['q1','o1','p1','k1','n2','m2','s2','k2','s1','m4']
    assert all(c in m.constituents for c in constituents)

# parameterize model
@pytest.mark.parametrize("MODEL", pyTMD.io.model.ATLAS())
def test_parse_TPXO9_elevation(MODEL):
    """Tests the parsing of ATLAS-type elevation model files
    """
    m = pyTMD.io.model(verify=False).elevation(MODEL)
    m.constituents = [pyTMD.io.model.parse_file(f) for f in m.model_file]
    constituents = ['q1','o1','p1','k1','n2','m2','s2','k2','m4','ms4','mn4','2n2']
    assert all(c in m.constituents for c in constituents)
    # test additional constituents found in newer models
    if MODEL in ('TPXO9-atlas-v3','TPXO9-atlas-v4','TPXO9-atlas-v5'):
        assert all(c in m.constituents for c in ['mf','mm'])
    if MODEL in ('TPXO9-atlas-v5',):
        assert all(c in m.constituents for c in ['s1',])

# parameterize model
current_models = set(pyTMD.io.model.ATLAS()) & set(pyTMD.io.model.ocean_current())
@pytest.mark.parametrize("MODEL", sorted(current_models))
def test_parse_TPXO9_currents(MODEL):
    """Tests the parsing of ATLAS-type current model files
    """
    m = pyTMD.io.model(verify=False).current(MODEL)
    m.constituents = [pyTMD.io.model.parse_file(f) for f in m.model_file['u']]
    constituents = ['q1','o1','p1','k1','n2','m2','s2','k2','m4','ms4','mn4','2n2']
    assert all(c in m.constituents for c in constituents)
    # test additional constituents found in newer models
    if MODEL in ('TPXO9-atlas-v3','TPXO9-atlas-v4','TPXO9-atlas-v5'):
        assert all(c in m.constituents for c in ['mf','mm'])
    if MODEL in ('TPXO9-atlas-v5',):
        assert all(c in m.constituents for c in ['s1',])
