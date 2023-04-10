"""
test_model.py (04/2023)
Tests the reading of model definition files
"""
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
    assert sorted(m.model_file) == [pathlib.Path(f) for f in model_files]
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

def test_definition_FES_currents():
    """Tests the reading of the FES2014 model definition file
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
    assert sorted(m.model_file['u']) == [pathlib.Path(f) for f in model_files_u]
    assert sorted(m.model_file['v']) == [pathlib.Path(f) for f in model_files_v]
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
    assert sorted(m.model_file) == [pathlib.Path(f) for f in model_files]
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
        'TPXO9_atlas_v5/h_s2_tpxo9_atlas_30_v5']
    assert m.grid_file == pathlib.Path('TPXO9_atlas_v5/grid_tpxo9_atlas_30_v5.nc')
    assert sorted(m.model_file) == [pathlib.Path(f) for f in model_files]
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

def test_definition_TPXO9_currents():
    """Tests the reading of the TPXO9-atlas-v5 model definition file
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
        'TPXO9_atlas_v5/u_s2_tpxo9_atlas_30_v5']
    assert m.grid_file == pathlib.Path('TPXO9_atlas_v5/grid_tpxo9_atlas_30_v5.nc')
    assert sorted(m.model_file['u']) == [pathlib.Path(f) for f in model_files]
    assert m.type == ['u', 'v']
    assert m.scale == 1.0/100.0
    assert m.compressed is False
