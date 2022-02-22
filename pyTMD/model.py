#!/usr/bin/env python
u"""
model.py
Written by Tyler Sutterley (02/2022)
Retrieves tide model parameters for named tide models and
    from model definition files

UPDATE HISTORY:
    Updated 02/2022: added Arctic 2km model (Arc2kmTM) to list of models
    Updated 01/2022: added global Empirical Ocean Tide model (EOT20)
    Updated 12/2021: added TPXO9-atlas-v5 to list of available tide models
        added atl10 attributes for tidal elevation files
    Written 09/2021
"""
import os
import re
import io
import copy

class model:
    """Retrieves tide model parameters for named models or
    from a model definition file for use in the pyTMD tide
    prediction programs
    """
    def __init__(self, directory=os.getcwd(), **kwargs):
        # set default keyword arguments
        kwargs.setdefault('compressed',False)
        kwargs.setdefault('format','netcdf')
        kwargs.setdefault('verify',True)
        # set initial attributes
        self.atl03 = None
        self.atl06 = None
        self.atl07 = None
        self.atl10 = None
        self.atl11 = None
        self.atl12 = None
        self.compressed = copy.copy(kwargs['compressed'])
        self.constituents = None
        self.description = None
        self.directory = os.path.expanduser(directory)
        self.format = copy.copy(kwargs['format'])
        self.gla12 = None
        self.grid_file = None
        self.long_name = None
        self.model_file = None
        self.name = None
        self.projection = None
        self.reference = None
        self.scale = None
        self.type = None
        self.variable = None
        self.verify = copy.copy(kwargs['verify'])
        self.version = None

    def grid(self,m):
        """Create a model object from known tide grid files
        """
        # model name
        self.name = m
        # select between known tide models
        if (m == 'CATS0201'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'cats0201_tmd')
            self.grid_file = self.pathfinder('grid_CATS')
        elif (m == 'CATS2008'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'CATS2008')
            self.grid_file = self.pathfinder('grid_CATS2008')
        elif (m == 'CATS2008_load'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,
                'CATS2008a_SPOTL_Load')
            self.grid_file = self.pathfinder('grid_CATS2008a_opt')
        elif (m == 'TPXO9-atlas'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas')
            self.version = 'v1'
        elif (m == 'TPXO9-atlas-v2'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas_v2')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas_30_v2')
            self.version = 'v2'
        elif (m == 'TPXO9-atlas-v3'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas_v3')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas_30_v3')
            self.version = 'v3'
        elif (m == 'TPXO9-atlas-v4'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas_v4')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas_30_v4')
            self.version = 'v4'
        elif (m == 'TPXO9-atlas-v5'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas_v5')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas_30_v5')
            self.version = 'v5'
        elif (m == 'TPXO9.1'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'TPXO9.1','DATA')
            self.grid_file = self.pathfinder('grid_tpxo9')
            self.version = '9.1'
        elif (m == 'TPXO8-atlas'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'tpxo8_atlas')
            self.grid_file = self.pathfinder('grid_tpxo8atlas_30_v1')
            self.version = '8'
        elif (m == 'TPXO7.2'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'TPXO7.2_tmd')
            self.grid_file = self.pathfinder('grid_tpxo7.2')
            self.version = '7.2'
        elif (m == 'TPXO7.2_load'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'TPXO7.2_load')
            self.grid_file = self.pathfinder('grid_tpxo6.2')
            self.version = '7.2'
        elif (m == 'AODTM-5'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'aodtm5_tmd')
            self.grid_file = self.pathfinder('grid_Arc5km')
        elif (m == 'AOTIM-5'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'aotim5_tmd')
            self.grid_file = self.pathfinder('grid_Arc5km')
        elif (m == 'AOTIM-5-2018'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'Arc5km2018')
            self.grid_file = self.pathfinder('grid_Arc5km2018')
            self.version = '2018'
        elif (m == 'Arc2kmTM'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'Arc2kmTM')
            self.grid_file = self.pathfinder('grid_Arc2kmTM_v1')
        elif (m == 'Gr1km-v2'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'greenlandTMD_v2')
            self.grid_file = self.pathfinder('grid_Greenland8.v2')
            self.version = 'v2'
        else:
            raise Exception("Unlisted tide model")
        # return the model parameters
        return self

    def elevation(self,m):
        """Create a model object from known tidal elevation models
        """
        # model name
        self.name = m
        # model type
        self.type = 'z'
        # select between known tide models
        if (m == 'CATS0201'):
            self.model_directory = os.path.join(self.directory,'cats0201_tmd')
            self.grid_file = self.pathfinder('grid_CATS')
            self.model_file = self.pathfinder('h0_CATS02_01')
            self.format = 'OTIS'
            self.projection = '4326'
            # model description and references
            self.reference = ('https://mail.esr.org/polar_tide_models/'
                'Model_CATS0201.html')
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'CATS2008'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'CATS2008')
            self.grid_file = self.pathfinder('grid_CATS2008')
            self.model_file = self.pathfinder('hf.CATS2008.out')
            self.projection = 'CATS2008'
            # model description and references
            self.reference = ('https://www.esr.org/research/'
                'polar-tide-models/list-of-polar-tide-models/cats2008/')
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'CATS2008_load'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,
                'CATS2008a_SPOTL_Load')
            self.grid_file = self.pathfinder('grid_CATS2008a_opt')
            self.model_file = self.pathfinder('h_CATS2008a_SPOTL_load')
            self.projection = 'CATS2008'
            # model description and references
            self.reference = ('https://www.esr.org/research/'
                'polar-tide-models/list-of-polar-tide-models/cats2008/')
            self.atl03 = 'tide_load'
            self.atl06 = 'tide_load'
            self.atl07 = 'height_segment_load'
            self.atl10 = 'height_segment_load'
            self.atl11 = 'tide_load'
            self.atl12 = 'tide_load_seg'
            self.gla12 = 'd_ldElv'
            self.variable = 'tide_load'
            self.long_name = "Load Tide"
            self.description = ("Local displacement due to Ocean "
                "Loading (-6 to 0 cm)")
        elif (m == 'TPXO9-atlas'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas')
            model_files = ['h_q1_tpxo9_atlas_30','h_o1_tpxo9_atlas_30',
                'h_p1_tpxo9_atlas_30','h_k1_tpxo9_atlas_30',
                'h_n2_tpxo9_atlas_30','h_m2_tpxo9_atlas_30',
                'h_s2_tpxo9_atlas_30','h_k2_tpxo9_atlas_30',
                'h_m4_tpxo9_atlas_30','h_ms4_tpxo9_atlas_30',
                'h_mn4_tpxo9_atlas_30','h_2n2_tpxo9_atlas_30']
            self.model_file = self.pathfinder(model_files)
            self.projection = '4326'
            self.scale = 1.0/1000.0
            self.version = 'v1'
            # model description and references
            self.reference = ('http://volkov.oce.orst.edu/tides/'
                'tpxo9_atlas.html')
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'TPXO9-atlas-v2'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas_v2')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas_30_v2')
            model_files = ['h_q1_tpxo9_atlas_30_v2','h_o1_tpxo9_atlas_30_v2',
                'h_p1_tpxo9_atlas_30_v2','h_k1_tpxo9_atlas_30_v2',
                'h_n2_tpxo9_atlas_30_v2','h_m2_tpxo9_atlas_30_v2',
                'h_s2_tpxo9_atlas_30_v2','h_k2_tpxo9_atlas_30_v2',
                'h_m4_tpxo9_atlas_30_v2','h_ms4_tpxo9_atlas_30_v2',
                'h_mn4_tpxo9_atlas_30_v2','h_2n2_tpxo9_atlas_30_v2']
            self.model_file = self.pathfinder(model_files)
            self.projection = '4326'
            self.scale = 1.0/1000.0
            self.version = 'v2'
            # model description and references
            self.reference = 'https://www.tpxo.net/global/tpxo9-atlas'
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'TPXO9-atlas-v3'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas_v3')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas_30_v3')
            model_files = ['h_q1_tpxo9_atlas_30_v3','h_o1_tpxo9_atlas_30_v3',
                'h_p1_tpxo9_atlas_30_v3','h_k1_tpxo9_atlas_30_v3',
                'h_n2_tpxo9_atlas_30_v3','h_m2_tpxo9_atlas_30_v3',
                'h_s2_tpxo9_atlas_30_v3','h_k2_tpxo9_atlas_30_v3',
                'h_m4_tpxo9_atlas_30_v3','h_ms4_tpxo9_atlas_30_v3',
                'h_mn4_tpxo9_atlas_30_v3','h_2n2_tpxo9_atlas_30_v3',
                'h_mf_tpxo9_atlas_30_v3','h_mm_tpxo9_atlas_30_v3']
            self.model_file = self.pathfinder(model_files)
            self.projection = '4326'
            self.scale = 1.0/1000.0
            self.version = 'v3'
            # model description and references
            self.reference = 'https://www.tpxo.net/global/tpxo9-atlas'
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'TPXO9-atlas-v4'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas_v4')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas_30_v4')
            model_files = ['h_q1_tpxo9_atlas_30_v4','h_o1_tpxo9_atlas_30_v4',
                'h_p1_tpxo9_atlas_30_v4','h_k1_tpxo9_atlas_30_v4',
                'h_n2_tpxo9_atlas_30_v4','h_m2_tpxo9_atlas_30_v4',
                'h_s2_tpxo9_atlas_30_v4','h_k2_tpxo9_atlas_30_v4',
                'h_m4_tpxo9_atlas_30_v4','h_ms4_tpxo9_atlas_30_v4',
                'h_mn4_tpxo9_atlas_30_v4','h_2n2_tpxo9_atlas_30_v4',
                'h_mf_tpxo9_atlas_30_v4','h_mm_tpxo9_atlas_30_v4']
            self.model_file = self.pathfinder(model_files)
            self.projection = '4326'
            self.scale = 1.0/1000.0
            self.version = 'v4'
            # model description and references
            self.reference = 'https://www.tpxo.net/global/tpxo9-atlas'
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'TPXO9-atlas-v5'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas_v5')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas_30_v5')
            model_files = ['h_q1_tpxo9_atlas_30_v5','h_o1_tpxo9_atlas_30_v5',
                'h_p1_tpxo9_atlas_30_v5','h_k1_tpxo9_atlas_30_v5',
                'h_n2_tpxo9_atlas_30_v5','h_m2_tpxo9_atlas_30_v5',
                'h_s1_tpxo9_atlas_30_v5','h_s2_tpxo9_atlas_30_v5',
                'h_k2_tpxo9_atlas_30_v5','h_m4_tpxo9_atlas_30_v5',
                'h_ms4_tpxo9_atlas_30_v5','h_mn4_tpxo9_atlas_30_v5',
                'h_2n2_tpxo9_atlas_30_v5','h_mf_tpxo9_atlas_30_v5',
                'h_mm_tpxo9_atlas_30_v5']
            self.model_file = self.pathfinder(model_files)
            self.projection = '4326'
            self.scale = 1.0/1000.0
            self.version = 'v5'
            # model description and references
            self.reference = 'https://www.tpxo.net/global/tpxo9-atlas'
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'TPXO9.1'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'TPXO9.1','DATA')
            self.grid_file = self.pathfinder('grid_tpxo9')
            self.model_file = self.pathfinder('h_tpxo9.v1')
            self.projection = '4326'
            self.version = '9.1'
            # model description and references
            self.reference = ('http://volkov.oce.orst.edu/'
                'tides/global.html')
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'TPXO8-atlas'):
            self.format = 'ATLAS'
            self.model_directory = os.path.join(self.directory,'tpxo8_atlas')
            self.grid_file = self.pathfinder('grid_tpxo8atlas_30_v1')
            self.model_file = self.pathfinder('hf.tpxo8_atlas_30_v1')
            self.projection = '4326'
            self.version = '8'
            # model description and references
            self.reference = ('http://volkov.oce.orst.edu/'
                'tides/tpxo8_atlas.html')
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'TPXO7.2'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'TPXO7.2_tmd')
            self.grid_file = self.pathfinder('grid_tpxo7.2')
            self.model_file = self.pathfinder('h_tpxo7.2')
            self.projection = '4326'
            self.version = '7.2'
            # model description and references
            self.reference = ('http://volkov.oce.orst.edu/'
                'tides/global.html')
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'TPXO7.2_load'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'TPXO7.2_load')
            self.grid_file = self.pathfinder('grid_tpxo6.2')
            self.model_file = self.pathfinder('h_tpxo7.2_load')
            self.projection = '4326'
            self.version = '7.2'
            # model description and references
            self.reference = ('http://volkov.oce.orst.edu/'
                'tides/global.html')
            self.atl03 = 'tide_load'
            self.atl06 = 'tide_load'
            self.atl07 = 'height_segment_load'
            self.atl10 = 'height_segment_load'
            self.atl11 = 'tide_load'
            self.atl12 = 'tide_load_seg'
            self.gla12 = 'd_ldElv'
            self.variable = 'tide_load'
            self.long_name = "Load Tide"
            self.description = ("Local displacement due to Ocean "
                "Loading (-6 to 0 cm)")
        elif (m == 'AODTM-5'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'aodtm5_tmd')
            self.grid_file = self.pathfinder('grid_Arc5km')
            self.model_file = self.pathfinder('h0_Arc5km.oce')
            self.projection = 'PSNorth'
            # model description and references
            self.reference = ('https://www.esr.org/research/'
                'polar-tide-models/list-of-polar-tide-models/'
                'aodtm-5/')
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'AOTIM-5'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'aotim5_tmd')
            self.grid_file = self.pathfinder('grid_Arc5km')
            self.model_file = self.pathfinder('h_Arc5km.oce')
            self.projection = 'PSNorth'
            # model description and references
            self.reference = ('https://www.esr.org/research/'
                'polar-tide-models/list-of-polar-tide-models/'
                'aotim-5/')
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'AOTIM-5-2018'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'Arc5km2018')
            self.grid_file = self.pathfinder('grid_Arc5km2018')
            self.model_file = self.pathfinder('h_Arc5km2018')
            self.projection = 'PSNorth'
            self.version = '2018'
            # model description and references
            self.reference = ('https://www.esr.org/research/'
                'polar-tide-models/list-of-polar-tide-models/'
                'aotim-5/')
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'Arc2kmTM'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'Arc2kmTM')
            self.grid_file = self.pathfinder('grid_Arc2kmTM_v1')
            self.model_file = self.pathfinder('h_Arc2kmTM_v1')
            self.projection = '3413'
            self.version = 'v1'
            # model description and references
            self.reference = 'https://doi.org/10.18739/A2PV6B79W'
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'Gr1km-v2'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'greenlandTMD_v2')
            self.grid_file = self.pathfinder('grid_Greenland8.v2')
            self.model_file = self.pathfinder('h_Greenland8.v2')
            self.projection = '3413'
            self.version = 'v2'
            # model description and references
            self.reference = 'https://doi.org/10.1002/2016RG000546'
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'GOT4.7'):
            self.format = 'GOT'
            self.model_directory = os.path.join(self.directory,
                'GOT4.7','grids_oceantide')
            model_files = ['q1.d','o1.d','p1.d','k1.d','n2.d',
                'm2.d','s2.d','k2.d','s1.d','m4.d']
            self.model_file = self.pathfinder(model_files)
            self.scale = 1.0/100.0
            self.version = '4.7'
            # model description and references
            self.reference = ('https://denali.gsfc.nasa.gov/'
                'personal_pages/ray/MiscPubs/'
                '19990089548_1999150788.pdf')
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'GOT4.7_load'):
            self.format = 'GOT'
            self.model_directory = os.path.join(self.directory,
                'GOT4.7','grids_loadtide')
            model_files = ['q1load.d','o1load.d',
                'p1load.d','k1load.d','n2load.d',
                'm2load.d','s2load.d','k2load.d',
                's1load.d','m4load.d']
            self.model_file = self.pathfinder(model_files)
            self.scale = 1.0/1000.0
            self.version = '4.7'
            # model description and references
            self.reference = ('https://denali.gsfc.nasa.gov/'
                'personal_pages/ray/MiscPubs/'
                '19990089548_1999150788.pdf')
            self.atl03 = 'tide_load'
            self.atl06 = 'tide_load'
            self.atl07 = 'height_segment_load'
            self.atl10 = 'height_segment_load'
            self.atl11 = 'tide_load'
            self.atl12 = 'tide_load_seg'
            self.gla12 = 'd_ldElv'
            self.variable = 'tide_load'
            self.long_name = "Load Tide"
            self.description = ("Local displacement due to Ocean "
                "Loading (-6 to 0 cm)")
        elif (m == 'GOT4.8'):
            self.format = 'GOT'
            self.model_directory = os.path.join(self.directory,
                'got4.8','grids_oceantide')
            model_files = ['q1.d','o1.d','p1.d','k1.d','n2.d',
                'm2.d','s2.d','k2.d','s1.d','m4.d']
            self.model_file = self.pathfinder(model_files)
            self.scale = 1.0/100.0
            self.version = '4.8'
            # model description and references
            self.reference = ('https://denali.gsfc.nasa.gov/'
                'personal_pages/ray/MiscPubs/'
                '19990089548_1999150788.pdf')
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'GOT4.8_load'):
            self.format = 'GOT'
            self.model_directory = os.path.join(self.directory,
                'got4.8','grids_loadtide')
            model_files = ['q1load.d','o1load.d',
                'p1load.d','k1load.d','n2load.d',
                'm2load.d','s2load.d','k2load.d',
                's1load.d','m4load.d']
            self.model_file = self.pathfinder(model_files)
            self.scale = 1.0/1000.0
            self.version = '4.8'
            # model description and references
            self.reference = ('https://denali.gsfc.nasa.gov/'
                'personal_pages/ray/MiscPubs/'
                '19990089548_1999150788.pdf')
            self.atl03 = 'tide_load'
            self.atl06 = 'tide_load'
            self.atl07 = 'height_segment_load'
            self.atl10 = 'height_segment_load'
            self.atl11 = 'tide_load'
            self.atl12 = 'tide_load_seg'
            self.gla12 = 'd_ldElv'
            self.variable = 'tide_load'
            self.long_name = "Load Tide"
            self.description = ("Local displacement due to Ocean "
                "Loading (-6 to 0 cm)")
        elif (m == 'GOT4.10'):
            self.format = 'GOT'
            self.model_directory = os.path.join(self.directory,
                'GOT4.10c','grids_oceantide')
            model_files = ['q1.d','o1.d','p1.d','k1.d','n2.d',
                'm2.d','s2.d','k2.d','s1.d','m4.d']
            self.model_file = self.pathfinder(model_files)
            self.scale = 1.0/100.0
            self.version = '4.10'
            # model description and references
            self.reference = ('https://denali.gsfc.nasa.gov/'
                'personal_pages/ray/MiscPubs/'
                '19990089548_1999150788.pdf')
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'GOT4.10_load'):
            self.format = 'GOT'
            self.model_directory = os.path.join(self.directory,
                'GOT4.10c','grids_loadtide')
            model_files = ['q1load.d','o1load.d',
                'p1load.d','k1load.d','n2load.d',
                'm2load.d','s2load.d','k2load.d',
                's1load.d','m4load.d']
            self.model_file = self.pathfinder(model_files)
            self.scale = 1.0/1000.0
            self.version = '4.10'
            # model description and references
            self.reference = ('https://denali.gsfc.nasa.gov/'
                'personal_pages/ray/MiscPubs/'
                '19990089548_1999150788.pdf')
            self.atl03 = 'tide_load'
            self.atl06 = 'tide_load'
            self.atl07 = 'height_segment_load'
            self.atl10 = 'height_segment_load'
            self.atl11 = 'tide_load'
            self.atl12 = 'tide_load_seg'
            self.gla12 = 'd_ldElv'
            self.variable = 'tide_load'
            self.long_name = "Load Tide"
            self.description = ("Local displacement due to Ocean "
                "Loading (-6 to 0 cm)")
        elif (m == 'FES2014'):
            self.format = 'FES'
            self.model_directory = os.path.join(self.directory,
                'fes2014','ocean_tide')
            model_files = ['2n2.nc','eps2.nc','j1.nc','k1.nc',
                'k2.nc','l2.nc','la2.nc','m2.nc','m3.nc','m4.nc',
                'm6.nc','m8.nc','mf.nc','mks2.nc','mm.nc',
                'mn4.nc','ms4.nc','msf.nc','msqm.nc','mtm.nc',
                'mu2.nc','n2.nc','n4.nc','nu2.nc','o1.nc','p1.nc',
                'q1.nc','r2.nc','s1.nc','s2.nc','s4.nc','sa.nc',
                'ssa.nc','t2.nc']
            self.model_file = self.pathfinder(model_files)
            self.constituents = ['2n2','eps2','j1','k1','k2','l2',
                'lambda2','m2','m3','m4','m6','m8','mf','mks2','mm',
                'mn4','ms4','msf','msqm','mtm','mu2','n2','n4','nu2',
                'o1','p1','q1','r2','s1','s2','s4','sa','ssa','t2']
            self.scale = 1.0/100.0
            self.version = 'FES2014'
            # model description and references
            self.reference = ('https://www.aviso.altimetry.fr/'
                'en/data/products/auxiliary-products/'
                'global-tide-fes.html')
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'FES2014_load'):
            self.format = 'FES'
            self.model_directory = os.path.join(self.directory,
                'fes2014','load_tide')
            model_files = ['2n2.nc','eps2.nc','j1.nc','k1.nc',
                'k2.nc','l2.nc','la2.nc','m2.nc','m3.nc','m4.nc',
                'm6.nc','m8.nc','mf.nc','mks2.nc','mm.nc',
                'mn4.nc','ms4.nc','msf.nc','msqm.nc','mtm.nc',
                'mu2.nc','n2.nc','n4.nc','nu2.nc','o1.nc','p1.nc',
                'q1.nc','r2.nc','s1.nc','s2.nc','s4.nc','sa.nc',
                'ssa.nc','t2.nc']
            self.model_file = self.pathfinder(model_files)
            self.constituents = ['2n2','eps2','j1','k1','k2','l2',
                'lambda2','m2','m3','m4','m6','m8','mf','mks2','mm',
                'mn4','ms4','msf','msqm','mtm','mu2','n2','n4','nu2',
                'o1','p1','q1','r2','s1','s2','s4','sa','ssa','t2']
            self.scale = 1.0/100.0
            self.version = 'FES2014'
            # model description and references
            self.reference = ('https://www.aviso.altimetry.fr/'
                'en/data/products/auxiliary-products/'
                'global-tide-fes.html')
            self.atl03 = 'tide_load'
            self.atl06 = 'tide_load'
            self.atl07 = 'height_segment_load'
            self.atl10 = 'height_segment_load'
            self.atl11 = 'tide_load'
            self.atl12 = 'tide_load_seg'
            self.gla12 = 'd_ldElv'
            self.variable = 'tide_load'
            self.long_name = "Load Tide"
            self.description = ("Local displacement due to Ocean "
                "Loading (-6 to 0 cm)")
        elif (m == 'EOT20'):
            self.format = 'FES'
            self.model_directory = os.path.join(self.directory,
                'EOT20','ocean_tides')
            model_files = ['2N2_ocean_eot20.nc','J1_ocean_eot20.nc',
                'K1_ocean_eot20.nc','K2_ocean_eot20.nc',
                'M2_ocean_eot20.nc','M4_ocean_eot20.nc',
                'MF_ocean_eot20.nc','MM_ocean_eot20.nc',
                'N2_ocean_eot20.nc','O1_ocean_eot20.nc',
                'P1_ocean_eot20.nc','Q1_ocean_eot20.nc',
                'S1_ocean_eot20.nc','S2_ocean_eot20.nc',
                'SA_ocean_eot20.nc','SSA_ocean_eot20.nc',
                'T2_ocean_eot20.nc']
            self.model_file = self.pathfinder(model_files)
            self.constituents = ['2n2','j1','k1','k2','m2','m4',
                'mf','mm','n2','o1','p1','q1','s1','s2','sa',
                'ssa','t2']
            self.scale = 1.0/100.0
            self.version = 'EOT20'
            # model description and references
            self.reference = 'https://doi.org/10.17882/79489'
            self.atl03 = 'tide_ocean'
            self.atl06 = 'tide_ocean'
            self.atl07 = 'height_segment_ocean'
            self.atl10 = 'height_segment_ocean'
            self.atl11 = 'tide_ocean'
            self.atl12 = 'tide_ocean_seg'
            self.gla12 = 'd_ocElv'
            self.variable = 'tide_ocean'
            self.long_name = "Ocean Tide"
            self.description = ("Ocean Tides including diurnal and "
                "semi-diurnal (harmonic analysis), and longer period "
                "tides (dynamic and self-consistent equilibrium).")
        elif (m == 'EOT20_load'):
            self.format = 'FES'
            self.model_directory = os.path.join(self.directory,
                'EOT20','load_tides')
            model_files = ['2N2_load_eot20.nc','J1_load_eot20.nc',
                'K1_load_eot20.nc','K2_load_eot20.nc',
                'M2_load_eot20.nc','M4_load_eot20.nc',
                'MF_load_eot20.nc','MM_load_eot20.nc',
                'N2_load_eot20.nc','O1_load_eot20.nc',
                'P1_load_eot20.nc','Q1_load_eot20.nc',
                'S1_load_eot20.nc','S2_load_eot20.nc',
                'SA_load_eot20.nc','SSA_load_eot20.nc',
                'T2_load_eot20.nc']
            self.model_file = self.pathfinder(model_files)
            self.constituents = ['2n2','j1','k1','k2','m2','m4',
                'mf','mm','n2','o1','p1','q1','s1','s2','sa',
                'ssa','t2']
            self.scale = 1.0/100.0
            self.version = 'EOT20'
            # model description and references
            self.reference = 'https://doi.org/10.17882/79489'
            self.atl03 = 'tide_load'
            self.atl06 = 'tide_load'
            self.atl07 = 'height_segment_load'
            self.atl10 = 'height_segment_load'
            self.atl11 = 'tide_load'
            self.atl12 = 'tide_load_seg'
            self.gla12 = 'd_ldElv'
            self.variable = 'tide_load'
            self.long_name = "Load Tide"
            self.description = ("Local displacement due to Ocean "
                "Loading (-6 to 0 cm)")
        else:
            raise Exception("Unlisted tide model")
        # return the model parameters
        return self

    def current(self,m):
        """Create a model object from known tidal current models
        """
        # model name
        self.name = m
        # model type
        self.type = ['u','v']
        # select between tide models
        if (m == 'CATS0201'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'cats0201_tmd')
            self.grid_file = self.pathfinder('grid_CATS')
            self.model_file = dict(u=self.pathfinder('UV0_CATS02_01'))
            self.projection = '4326'
            # model description and references
            self.reference = ('https://mail.esr.org/polar_tide_models/'
                'Model_CATS0201.html')
        elif (m == 'CATS2008'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'CATS2008')
            self.grid_file = self.pathfinder('grid_CATS2008')
            self.model_file = dict(u=self.pathfinder('uv.CATS2008.out'))
            self.projection = 'CATS2008'
        elif (m == 'TPXO9-atlas'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas')
            model_files = {}
            model_files['u'] = ['u_q1_tpxo9_atlas_30','u_o1_tpxo9_atlas_30',
                'u_p1_tpxo9_atlas_30','u_k1_tpxo9_atlas_30',
                'u_n2_tpxo9_atlas_30','u_m2_tpxo9_atlas_30',
                'u_s2_tpxo9_atlas_30','u_k2_tpxo9_atlas_30',
                'u_m4_tpxo9_atlas_30','u_ms4_tpxo9_atlas_30',
                'u_mn4_tpxo9_atlas_30','u_2n2_tpxo9_atlas_30']
            model_files['v'] = ['v_q1_tpxo9_atlas_30','v_o1_tpxo9_atlas_30',
                'v_p1_tpxo9_atlas_30','v_k1_tpxo9_atlas_30',
                'v_n2_tpxo9_atlas_30','v_m2_tpxo9_atlas_30',
                'v_s2_tpxo9_atlas_30','v_k2_tpxo9_atlas_30',
                'v_m4_tpxo9_atlas_30','v_ms4_tpxo9_atlas_30',
                'v_mn4_tpxo9_atlas_30','v_2n2_tpxo9_atlas_30']
            self.model_file = {}
            for key,val in model_files.items():
                self.model_file[key] = self.pathfinder(val)
            self.projection = '4326'
            self.scale = 1.0/100.0
            self.version = 'v1'
            # model description and references
            self.reference = ('http://volkov.oce.orst.edu/tides/'
                'tpxo9_atlas.html')
        elif (m == 'TPXO9-atlas-v2'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas_v2')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas_30_v2')
            model_files = {}
            model_files['u'] = ['u_q1_tpxo9_atlas_30_v2','u_o1_tpxo9_atlas_30_v2',
                'u_p1_tpxo9_atlas_30_v2','u_k1_tpxo9_atlas_30_v2',
                'u_n2_tpxo9_atlas_30_v2','u_m2_tpxo9_atlas_30_v2',
                'u_s2_tpxo9_atlas_30_v2','u_k2_tpxo9_atlas_30_v2',
                'u_m4_tpxo9_atlas_30_v2','u_ms4_tpxo9_atlas_30_v2',
                'u_mn4_tpxo9_atlas_30_v2','u_2n2_tpxo9_atlas_30_v2']
            model_files['v'] = ['v_q1_tpxo9_atlas_30_v2','v_o1_tpxo9_atlas_30_v2',
                'v_p1_tpxo9_atlas_30_v2','v_k1_tpxo9_atlas_30_v2',
                'v_n2_tpxo9_atlas_30_v2','v_m2_tpxo9_atlas_30_v2',
                'v_s2_tpxo9_atlas_30_v2','v_k2_tpxo9_atlas_30_v2',
                'v_m4_tpxo9_atlas_30_v2','v_ms4_tpxo9_atlas_30_v2',
                'v_mn4_tpxo9_atlas_30_v2','v_2n2_tpxo9_atlas_30_v2']
            self.model_file = {}
            for key,val in model_files.items():
                self.model_file[key] = self.pathfinder(val)
            self.projection = '4326'
            self.scale = 1.0/100.0
            self.version = 'v2'
            # model description and references
            self.reference = 'https://www.tpxo.net/global/tpxo9-atlas'
        elif (m == 'TPXO9-atlas-v3'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas_v3')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas_30_v3')
            model_files = {}
            model_files['u'] = ['u_q1_tpxo9_atlas_30_v3','u_o1_tpxo9_atlas_30_v3',
                'u_p1_tpxo9_atlas_30_v3','u_k1_tpxo9_atlas_30_v3',
                'u_n2_tpxo9_atlas_30_v3','u_m2_tpxo9_atlas_30_v3',
                'u_s2_tpxo9_atlas_30_v3','u_k2_tpxo9_atlas_30_v3',
                'u_m4_tpxo9_atlas_30_v3','u_ms4_tpxo9_atlas_30_v3',
                'u_mn4_tpxo9_atlas_30_v3','u_2n2_tpxo9_atlas_30_v3']
            model_files['v'] = ['v_q1_tpxo9_atlas_30_v3','v_o1_tpxo9_atlas_30_v3',
                'v_p1_tpxo9_atlas_30_v3','v_k1_tpxo9_atlas_30_v3',
                'v_n2_tpxo9_atlas_30_v3','v_m2_tpxo9_atlas_30_v3',
                'v_s2_tpxo9_atlas_30_v3','v_k2_tpxo9_atlas_30_v3',
                'v_m4_tpxo9_atlas_30_v3','v_ms4_tpxo9_atlas_30_v3',
                'v_mn4_tpxo9_atlas_30_v3','v_2n2_tpxo9_atlas_30_v3']
            self.model_file = {}
            for key,val in model_files.items():
                self.model_file[key] = self.pathfinder(val)
            self.projection = '4326'
            self.scale = 1.0/100.0
            self.version = 'v3'
            # model description and references
            self.reference = 'https://www.tpxo.net/global/tpxo9-atlas'
        elif (m == 'TPXO9-atlas-v4'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas_v4')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas_30_v4')
            model_files = {}
            model_files['u'] = ['u_q1_tpxo9_atlas_30_v4','u_o1_tpxo9_atlas_30_v4',
                'u_p1_tpxo9_atlas_30_v4','u_k1_tpxo9_atlas_30_v4',
                'u_n2_tpxo9_atlas_30_v4','u_m2_tpxo9_atlas_30_v4',
                'u_s2_tpxo9_atlas_30_v4','u_k2_tpxo9_atlas_30_v4',
                'u_m4_tpxo9_atlas_30_v4','u_ms4_tpxo9_atlas_30_v4',
                'u_mn4_tpxo9_atlas_30_v4','u_2n2_tpxo9_atlas_30_v4']
            model_files['v'] = ['v_q1_tpxo9_atlas_30_v4','v_o1_tpxo9_atlas_30_v4',
                'v_p1_tpxo9_atlas_30_v4','v_k1_tpxo9_atlas_30_v4',
                'v_n2_tpxo9_atlas_30_v4','v_m2_tpxo9_atlas_30_v4',
                'v_s2_tpxo9_atlas_30_v4','v_k2_tpxo9_atlas_30_v4',
                'v_m4_tpxo9_atlas_30_v4','v_ms4_tpxo9_atlas_30_v4',
                'v_mn4_tpxo9_atlas_30_v4','v_2n2_tpxo9_atlas_30_v4']
            self.model_file = {}
            for key,val in model_files.items():
                self.model_file[key] = self.pathfinder(val)
            self.projection = '4326'
            self.scale = 1.0/100.0
            self.version = 'v4'
            # model description and references
            self.reference = 'https://www.tpxo.net/global/tpxo9-atlas'
        elif (m == 'TPXO9-atlas-v5'):
            self.model_directory = os.path.join(self.directory,'TPXO9_atlas_v5')
            self.grid_file = self.pathfinder('grid_tpxo9_atlas_30_v5')
            model_files = {}
            model_files['u'] = ['u_q1_tpxo9_atlas_30_v5','u_o1_tpxo9_atlas_30_v5',
                'u_p1_tpxo9_atlas_30_v5','u_k1_tpxo9_atlas_30_v5',
                'u_n2_tpxo9_atlas_30_v5','u_m2_tpxo9_atlas_30_v5',
                'u_s1_tpxo9_atlas_30_v5','u_s2_tpxo9_atlas_30_v5',
                'u_k2_tpxo9_atlas_30_v5','u_m4_tpxo9_atlas_30_v5',
                'u_ms4_tpxo9_atlas_30_v5','u_mn4_tpxo9_atlas_30_v5',
                'u_2n2_tpxo9_atlas_30_v5']
            model_files['v'] = ['v_q1_tpxo9_atlas_30_v5','v_o1_tpxo9_atlas_30_v5',
                'v_p1_tpxo9_atlas_30_v5','v_k1_tpxo9_atlas_30_v5',
                'v_n2_tpxo9_atlas_30_v5','v_m2_tpxo9_atlas_30_v5',
                'u_s1_tpxo9_atlas_30_v5','u_s2_tpxo9_atlas_30_v5',
                'u_k2_tpxo9_atlas_30_v5','u_m4_tpxo9_atlas_30_v5',
                'u_ms4_tpxo9_atlas_30_v5','u_mn4_tpxo9_atlas_30_v5',
                'u_2n2_tpxo9_atlas_30_v5']
            self.model_file = {}
            for key,val in model_files.items():
                self.model_file[key] = self.pathfinder(val)
            self.projection = '4326'
            self.scale = 1.0/100.0
            self.version = 'v5'
            # model description and references
            self.reference = 'https://www.tpxo.net/global/tpxo9-atlas'
        elif (m == 'TPXO9.1'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'TPXO9.1')
            self.grid_file = self.pathfinder('grid_tpxo9')
            self.model_file = dict(u=self.pathfinder('u_tpxo9.v1'))
            self.projection = '4326'
            self.version = '9.1'
            # model description and references
            self.reference = ('http://volkov.oce.orst.edu/tides/'
                'global.html')
        elif (m == 'TPXO8-atlas'):
            self.format = 'ATLAS'
            self.model_directory = os.path.join(self.directory,'tpxo8_atlas')
            self.grid_file = self.pathfinder('grid_tpxo8atlas_30_v1')
            self.model_file = dict(u=self.pathfinder('uv.tpxo8_atlas_30_v1'))
            self.projection = '4326'
            self.version = '8'
            # model description and references
            self.reference = ('http://volkov.oce.orst.edu/tides/'
                'tpxo8_atlas.html')
        elif (m == 'TPXO7.2'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'TPXO7.2_tmd')
            self.grid_file = self.pathfinder('grid_tpxo7.2')
            self.model_file = dict(u=self.pathfinder('u_tpxo7.2'))
            self.projection = '4326'
            self.version = '7.2'
            # model description and references
            self.reference = ('http://volkov.oce.orst.edu/tides/'
                'global.html')
        elif (m == 'AODTM-5'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'aodtm5_tmd')
            self.grid_file = self.pathfinder('grid_Arc5km')
            self.model_file = dict(u=self.pathfinder('UV0_Arc5km'))
            self.projection = 'PSNorth'
            # model description and references
            self.reference = ('https://www.esr.org/research/'
                'polar-tide-models/list-of-polar-tide-models/'
                'aodtm-5/')
        elif (m == 'AOTIM-5'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'aotim5_tmd')
            self.grid_file = self.pathfinder('grid_Arc5km')
            self.model_file = dict(u=self.pathfinder('UV_Arc5km'))
            self.projection = 'PSNorth'
            # model description and references
            self.reference = ('https://www.esr.org/research/'
                'polar-tide-models/list-of-polar-tide-models/'
                'aotim-5/')
        elif (m == 'AOTIM-5-2018'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'Arc5km2018')
            self.grid_file = self.pathfinder('grid_Arc5km2018')
            self.model_file = dict(u=self.pathfinder('UV_Arc5km2018'))
            self.projection = 'PSNorth'
            self.version = '2018'
            # model description and references
            self.reference = ('https://www.esr.org/research/'
                'polar-tide-models/list-of-polar-tide-models/'
                'aotim-5/')
        elif (m == 'Arc2kmTM'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'Arc2kmTM')
            self.grid_file = self.pathfinder('grid_Arc2kmTM_v1')
            self.model_file = dict(u=self.pathfinder('UV_Arc2kmTM_v1'))
            self.projection = '3413'
            self.version = 'v1'
            # model description and references
            self.reference = 'https://doi.org/10.18739/A2PV6B79W'
        elif (m == 'Gr1km-v2'):
            self.format = 'OTIS'
            self.model_directory = os.path.join(self.directory,'greenlandTMD_v2')
            self.grid_file = self.pathfinder('grid_Greenland8.v2')
            self.model_file = dict(u=self.pathfinder('u_Greenland8_rot.v2'))
            self.projection = '3413'
            self.version = 'v2'
            # model description and references
            self.reference = 'https://doi.org/10.1002/2016RG000546'
        elif (m == 'FES2014'):
            self.format = 'FES'
            model_directory = {}
            model_directory['u'] = os.path.join(self.directory,
                'fes2014','eastward_velocity')
            model_directory['v'] = os.path.join(self.directory,
                'fes2014','northward_velocity')
            model_files = ['2n2.nc','eps2.nc','j1.nc','k1.nc',
                'k2.nc','l2.nc','la2.nc','m2.nc','m3.nc','m4.nc',
                'm6.nc','m8.nc','mf.nc','mks2.nc','mm.nc',
                'mn4.nc','ms4.nc','msf.nc','msqm.nc','mtm.nc',
                'mu2.nc','n2.nc','n4.nc','nu2.nc','o1.nc','p1.nc',
                'q1.nc','r2.nc','s1.nc','s2.nc','s4.nc','sa.nc',
                'ssa.nc','t2.nc']
            self.model_file = {}
            for key,val in model_directory.items():
                self.model_directory = os.path.expanduser(model_directory)
                self.model_file[key] = self.pathfinder(val)
            self.constituents = ['2n2','eps2','j1','k1','k2','l2','lambda2',
                'm2','m3','m4','m6','m8','mf','mks2','mm','mn4','ms4','msf',
                'msqm','mtm','mu2','n2','n4','nu2','o1','p1','q1','r2','s1',
                's2','s4','sa','ssa','t2']
            self.scale = 1.0
            self.version = 'FES2014'
            # model description and references
            self.reference = ('https://www.aviso.altimetry.fr/en/data/products'
                'auxiliary-products/global-tide-fes.html')
        else:
            raise Exception("Unlisted tide model")
        # return the model parameters
        return self

    @property
    def gzip(self):
        """compression flag"""
        return '.gz' if self.compressed else ''

    @property
    def suffix(self):
        """format suffix flag"""
        return '.nc' if (self.format == 'netcdf') else ''

    def pathfinder(self,model_file):
        """Completes file paths and appends file and gzip suffixes
        """
        if isinstance(model_file,list):
            output_file = [os.path.join(self.model_directory,
                ''.join([f,self.suffix,self.gzip])) for f in model_file]
            valid = all([os.access(f, os.F_OK) for f in output_file])
        elif isinstance(model_file,str):
            output_file = os.path.join(self.model_directory,
                ''.join([model_file,self.suffix,self.gzip]))
            valid = os.access(output_file, os.F_OK)
        #-- check that (all) output files exist
        if self.verify and not valid:
            raise FileNotFoundError(output_file)
        #-- return the complete output path
        return output_file

    def from_file(self, definition_file):
        """Create a model object from an input definition file
        """
        # variable with parameter definitions
        parameters = {}
        # Opening definition file and assigning file ID number
        if isinstance(definition_file,io.IOBase):
            fid = copy.copy(definition_file)
        else:
            fid = open(os.path.expanduser(definition_file), 'r')
        # for each line in the file will extract the parameter (name and value)
        for fileline in fid:
            # Splitting the input line between parameter name and value
            part = fileline.rstrip().split(maxsplit=1)
            # filling the parameter definition variable
            parameters[part[0]] = part[1]
        # close the parameter file
        fid.close()
        # convert from dictionary to model variable
        temp = self.from_dict(parameters)
        # verify model name, format and type
        assert temp.name
        assert temp.format in ('OTIS','ATLAS','netcdf','GOT','FES')
        assert temp.type
        # verify necessary attributes are with model format
        assert temp.model_file
        # split model file into list if an ATLAS, GOT or FES file
        # model files can be comma, tab or space delimited
        # extract full path to tide model files
        if re.search(r'[\s\,]+', temp.model_file):
            temp.model_file = [os.path.expanduser(f) for f in
                re.split(r'[\s\,]+',temp.model_file)]
            temp.model_directory = os.path.dirname(temp.model_file[0])
        else:
            temp.model_file = os.path.expanduser(temp.model_file)
            temp.model_directory = os.path.dirname(temp.model_file)
        # extract full path to tide grid file
        if temp.format in ('OTIS','ATLAS','netcdf'):
            assert temp.grid_file
            temp.grid_file = os.path.expanduser(temp.grid_file)
        if temp.format in ('OTIS','ATLAS'):
            assert temp.projection
        # convert scale from string to float
        if temp.format in ('netcdf','GOT','FES'):
            assert temp.scale
            temp.scale = float(temp.scale)
        if temp.format in ('FES',):
            assert temp.version
        # split type into list if currents u,v
        if re.search(r'[\s\,]+', temp.type):
            temp.type = re.split(r'[\s\,]+',temp.type)
        # convert boolean strings
        if isinstance(temp.compressed,str):
            temp.compressed = self.to_bool(temp.compressed)
        # return the model parameters
        return temp

    def from_dict(self,d):
        """Create a model object from a python dictionary
        """
        for key,val in d.items():
            setattr(self,key,copy.copy(val))
        # return the model parameters
        return self

    def to_bool(self,val):
        """Converts strings of True/False to a boolean values
        """
        if val.lower() in ('y','yes','t','true','1'):
            return True
        elif val.lower() in ('n','no','f','false','0'):
            return False
        else:
            raise ValueError('Invalid boolean string {0}'.format(val))
