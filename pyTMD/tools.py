#!/usr/bin/env python
u"""
tools.py
Written by Tyler Sutterley (11/2022)
Jupyter notebook, user interface and plotting tools

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    ipywidgets: interactive HTML widgets for Jupyter notebooks and IPython
        https://ipywidgets.readthedocs.io/en/latest/
    ipyleaflet: Jupyter / Leaflet bridge enabling interactive maps
        https://github.com/jupyter-widgets/ipyleaflet
    matplotlib: Python 2D plotting library
        http://matplotlib.org/
        https://github.com/matplotlib/matplotlib

UPDATE HISTORY:
    Updated 11/2022: place more imports within try/except statements
    Updated 08/2022: place some imports behind try/except statements
    Updated 05/2022: include world copy jump in webmercator maps
    Updated 03/2022: add marker relocation routines from notebooks
    Updated 02/2022: add leaflet map projections
    Written 09/2021
"""
import io
import os
import copy
import base64
import datetime
import warnings
import numpy as np
import matplotlib
matplotlib.rcParams['axes.linewidth'] = 2.0
matplotlib.rcParams["animation.html"] = "jshtml"
import matplotlib.cm as cm
import matplotlib.colorbar
import matplotlib.animation
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pyTMD.model

# attempt imports
try:
    import IPython
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("IPython not available")
    warnings.warn("Some functions will throw an exception if called")
try:
    import ipyleaflet
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("ipyleaflet not available")
    warnings.warn("Some functions will throw an exception if called")
try:
    import ipywidgets
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("ipywidgets not available")
    warnings.warn("Some functions will throw an exception if called")
try:
    import pyproj
except (ImportError, ModuleNotFoundError) as e:
    warnings.filterwarnings("always")
    warnings.warn("pyproj not available")
    warnings.warn("Some functions will throw an exception if called")
# ignore warnings
warnings.filterwarnings("ignore")

class widgets:
    def __init__(self, **kwargs):
        # set default keyword options
        kwargs.setdefault('style', {})
        # set style
        self.style = copy.copy(kwargs['style'])

        # set the directory with tide models
        self.directory = ipywidgets.Text(
            value=os.getcwd(),
            description='Directory:',
            disabled=False
        )

        # dropdown menu for setting tide model
        model_list = sorted(pyTMD.model.ocean_elevation() +
            pyTMD.model.load_elevation())
        self.model = ipywidgets.Dropdown(
            options=model_list,
            value='GOT4.10',
            description='Model:',
            disabled=False,
            style=self.style,
        )

        # dropdown menu for setting ATLAS format model
        atlas_list = ['OTIS','netcdf']
        self.atlas = ipywidgets.Dropdown(
            options=atlas_list,
            value='netcdf',
            description='ATLAS:',
            disabled=False,
            style=self.style,
        )
        self.atlas.layout.display = 'none'

        # checkbox for setting if tide files are compressed
        self.compress = ipywidgets.Checkbox(
            value=True,
            description='Compressed?',
            disabled=False,
            style=self.style,
        )

        # date picker widget for setting time
        self.datepick = ipywidgets.DatePicker(
            description='Date:',
            value = datetime.date.today(),
            disabled=False,
            style=self.style,
        )

        # watch widgets for changes
        self.model.observe(self.set_atlas)

    # function for setting available map layers
    def set_atlas(self, sender):
        """function for updating ATLAS widget visibility
        """
        if (self.model.value in pyTMD.model.ATLAS()):
            self.atlas.layout.display = 'inline-flex'

# define projections for ipyleaflet tiles
projections = dict(
    # Alaska Polar Stereographic (WGS84)
    EPSG5936 = dict(
        Basemap = dict(
            name='EPSG:5936',
            custom=True,
            proj4def="""+proj=stere +lat_0=90 +lat_ts=90 +lon_0=-150 +k=0.994
                +x_0=2000000 +y_0=2000000 +datum=WGS84 +units=m +no_defs""",
            origin=[-2.8567784109255e+07, 3.2567784109255e+07],
            resolutions=[
                238810.813354,
                119405.406677,
                59702.7033384999,
                29851.3516692501,
                14925.675834625,
                7462.83791731252,
                3731.41895865639,
                1865.70947932806,
                932.854739664032,
                466.427369832148,
                233.213684916074,
                116.60684245803701,
                58.30342122888621,
                29.151710614575396,
                14.5758553072877,
                7.28792765351156,
                3.64396382688807,
                1.82198191331174,
                0.910990956788164,
                0.45549547826179,
                0.227747739130895,
                0.113873869697739,
                0.05693693484887,
                0.028468467424435
            ],
            bounds=[
                [-2623285.8808999992907047,-2623285.8808999992907047],
                [6623285.8803000003099442,6623285.8803000003099442]
            ]
        )
    ),
    # Polar Stereographic South (WGS84)
    EPSG3031 = dict(
        Basemap = dict(
            name='EPSG:3031',
            custom=True,
            proj4def="""+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1
                +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs""",
            origin=[-3.06361E7, 3.0636099999999993E7],
            resolutions=[
                67733.46880027094,
                33866.73440013547,
                16933.367200067736,
                8466.683600033868,
                4233.341800016934,
                2116.670900008467,
                1058.3354500042335,
                529.1677250021168,
                264.5838625010584,
            ],
            bounds=[
                [-4524583.19363305,-4524449.487765655],
                [4524449.4877656475,4524583.193633042]
            ]
        )
    )
)

# draw ipyleaflet map
class leaflet:
    def __init__(self, projection='Global', **kwargs):
        # set default keyword arguments
        kwargs.setdefault('map',None)
        kwargs.setdefault('attribution',True)
        kwargs.setdefault('zoom',1)
        kwargs.setdefault('zoom_control',False)
        kwargs.setdefault('scale_control',False)
        kwargs.setdefault('cursor_control',True)
        kwargs.setdefault('layer_control',True)
        kwargs.setdefault('center',(39,-108))
        # create basemap in projection
        if (projection == 'Global'):
            self.map = ipyleaflet.Map(center=kwargs['center'],
                zoom=kwargs['zoom'], max_zoom=15, world_copy_jump=True,
                attribution_control=kwargs['attribution'],
                basemap=ipyleaflet.basemaps.Esri.WorldTopoMap)
            self.crs = 'EPSG:3857'
        elif (projection == 'North'):
            self.map = ipyleaflet.Map(center=kwargs['center'],
                zoom=kwargs['zoom'], max_zoom=24,
                attribution_control=kwargs['attribution'],
                basemap=ipyleaflet.basemaps.Esri.ArcticOceanBase,
                crs=projections['EPSG5936']['Basemap'])
            self.map.add_layer(ipyleaflet.basemaps.Esri.ArcticOceanReference)
            self.crs = 'EPSG:5936'
        elif (projection == 'South'):
            self.map = ipyleaflet.Map(center=kwargs['center'],
                zoom=kwargs['zoom'], max_zoom=9,
                attribution_control=kwargs['attribution'],
                basemap=ipyleaflet.basemaps.Esri.AntarcticBasemap,
                crs=projections['EPSG3031']['Basemap'])
            self.crs = 'EPSG:3031'
        else:
            # use a predefined ipyleaflet map
            self.map = kwargs['map']
            self.crs = self.map.crs['name']
        # add control for layers
        if kwargs['layer_control']:
            self.layer_control = ipyleaflet.LayersControl(position='topleft')
            self.map.add_control(self.layer_control)
            self.layers = self.map.layers
        # add control for zoom
        if kwargs['zoom_control']:
            zoom_slider = ipywidgets.IntSlider(description='Zoom level:',
                min=self.map.min_zoom, max=self.map.max_zoom, value=self.map.zoom)
            ipywidgets.jslink((zoom_slider, 'value'), (self.map, 'zoom'))
            zoom_control = ipyleaflet.WidgetControl(widget=zoom_slider,
                position='topright')
            self.map.add_control(zoom_control)
        # add control for spatial scale bar
        if kwargs['scale_control']:
            scale_control = ipyleaflet.ScaleControl(position='topright')
            self.map.add_control(scale_control)
        # add control for cursor position
        if kwargs['cursor_control']:
            self.cursor = ipywidgets.Label()
            cursor_control = ipyleaflet.WidgetControl(widget=self.cursor,
                position='bottomleft')
            self.map.add_control(cursor_control)
            # keep track of cursor position
            self.map.on_interaction(self.handle_interaction)
        # add control for marker
        if kwargs['marker_control']:
            # add marker with default location
            self.marker = ipyleaflet.Marker(location=kwargs['center'],
                draggable=True)
            self.map.add_layer(self.marker)
            # add text with marker location
            self.marker_text = ipywidgets.Text(
                value='{0:0.8f},{1:0.8f}'.format(*kwargs['center']),
                description='Lat/Lon:',
                disabled=False)
            # watch marker widgets for changes
            self.marker.observe(self.set_marker_text)
            self.marker_text.observe(self.set_marker_location)
            self.map.observe(self.set_map_center)
            # add control for marker location
            marker_control = ipyleaflet.WidgetControl(
                widget=self.marker_text, position='bottomright')
            self.map.add_control(marker_control)

    # convert points to EPSG:4326
    def transform(self, x, y, proj4def):
        # convert geolocation variable to EPSG:4326
        crs1 = pyproj.CRS.from_string(proj4def)
        crs2 = pyproj.CRS.from_string('EPSG:4326')
        trans = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
        return trans.transform(x, y)

    # fix longitudes to be -180:180
    def wrap_longitudes(self, lon):
        phi = np.arctan2(np.sin(lon*np.pi/180.0),np.cos(lon*np.pi/180.0))
        # convert phi from radians to degrees
        return phi*180.0/np.pi

    # add function for setting marker text if location changed
    def set_marker_text(self, sender):
        LAT,LON = self.marker.location
        self.marker_text.value = '{0:0.8f},{1:0.8f}'.format(LAT,
            self.wrap_longitudes(LON))

    # add function for setting map center if location changed
    def set_map_center(self, sender):
        self.map.center = self.marker.location

    # add function for setting marker location if text changed
    def set_marker_location(self, sender):
        LAT,LON = [float(i) for i in self.marker_text.value.split(',')]
        self.marker.location = (LAT,LON)

    # handle cursor movements for label
    def handle_interaction(self, **kwargs):
        if (kwargs.get('type') == 'mousemove'):
            lat,lon = kwargs.get('coordinates')
            lon = self.wrap_longitudes(lon)
            self.cursor.value = u"""Latitude: {d[0]:8.4f}\u00B0,
                Longitude: {d[1]:8.4f}\u00B0""".format(d=[lat,lon])
