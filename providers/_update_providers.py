"""
_update_providers.py (08/2024)
Update the projection variables in the providers
"""
import json
import pyproj
import inspect
import pathlib
import argparse
import numpy as np
import pyTMD.utilities

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Update the projection variables in the providers"
            """,
        fromfile_prefix_chars="@"
    )
    # command line parameters
    parser.add_argument('--pretty', '-p',
        action='store_true',
        help='Pretty print the json file')
    parser.add_argument('--verbose', '-v',
        action='store_true',
        help='Verbose output')
    return parser

def get_projection(PROJ):
    # previously defined coordinate reference systems
    # python dictionary with named conversion functions
    CRS = {}
    CRS['3031'] = _EPSG3031
    CRS['3413'] = _EPSG3413
    CRS['CATS2008'] = _CATS2008
    CRS['3976'] = _EPSG3976
    CRS['AEDNorth'] = _AEDNorth
    CRS['4326'] = _EPSG4326
    # check that PROJ for conversion was entered correctly
    # run named conversion program and return values
    try:
        crs = CRS[PROJ].__call__()
    except KeyError as exc:
        pass
    else:
        # return the output variables
        return crs

def _EPSG3031():
    """
    CRS for models in EPSG:3031 (Antarctic Polar Stereographic)
    """
    # coordinate reference system information
    crs = pyproj.CRS.from_user_input({'proj':'stere',
        'lat_0':-90, 'lat_ts':-71, 'lon_0':0, 'x_0':0., 'y_0':0.,
        'ellps':'WGS84', 'datum':'WGS84', 'units':'km'})
    return crs

# wrapper function for models in EPSG 3413 (Sea Ice Polar Stereographic North)
def _EPSG3413():
    """
    CRS for models in EPSG:3413 (Sea Ice Polar Stereographic North)
    """
    # coordinate reference system information
    crs = pyproj.CRS.from_user_input({'proj':'stere',
        'lat_0':90, 'lat_ts':70, 'lon_0':-45, 'x_0':0., 'y_0':0.,
        'ellps':'WGS84', 'datum':'WGS84', 'units':'km'})
    return crs

# wrapper function for CATS2008 tide models
def _CATS2008():
    """
    CRS for Circum-Antarctic Tidal Simulation models
    """
    # coordinate reference system information
    crs = pyproj.CRS.from_user_input({'proj':'stere',
        'lat_0':-90, 'lat_ts':-71, 'lon_0':-70, 'x_0':0., 'y_0':0.,
        'ellps':'WGS84', 'datum':'WGS84', 'units':'km'})
    return crs

# wrapper function for models in EPSG 3976 (NSIDC Sea Ice Stereographic South)
def _EPSG3976():
    """
    CRS for models in EPSG:3976 (Sea Ice Polar Stereographic South)
    """
    # coordinate reference system information
    crs = pyproj.CRS.from_user_input({'proj':'stere',
        'lat_0':-90, 'lat_ts':-70, 'lon_0':0, 'x_0':0., 'y_0':0.,
        'ellps':'WGS84', 'datum':'WGS84', 'units':'km'})
    return crs

# function for models in (idealized) Azimuthal Equidistant projection
def _AEDNorth():
    """
    CRS for models in idealized Azimuthal Equidistant projections
    """
    # projections for converting to and from input EPSG
    R = 111700.0*180.0/np.pi
    crs = pyproj.CRS.from_user_input({'proj':'aeqd','lat_0':90,
        'lon_0':270,'x_0':0.,'y_0':0.,'ellps':'sphere',
        'R':R,'units':'km'})
    return crs

# wrapper function to pass lat/lon values or convert if EPSG
def _EPSG4326():
    """
    CRS for models in EPSG:4326 (WGS84 Latitude/Longitude)
    """
    crs = pyproj.CRS.from_epsg(4326)
    return crs

# wrapper function for using custom projections
def _custom(self, PROJ: int | str):
    """
    CRS for models in a custom projection

    Parameters
    ----------
    PROJ: int or str
        Spatial reference system code for coordinate transformations
    """
    # coordinate reference system information
    crs = self.from_input(PROJ)
    return crs

def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # find providers
    providers = [f for f in filepath.iterdir() if (f.suffix == '.json')]
    # for each provider
    for provider in providers:
        with provider.open('r', encoding='utf-8') as fid:
            p = json.load(fid)
            for key, val in p.items():
                # update projections
                for model, parameters in val.items():
                    if 'projection' in parameters:
                        crs = get_projection(parameters['projection'])
                        if crs.crs.to_epsg() is not None:
                            p[key][model]['projection'] = crs.crs.to_string()
                        else:
                            d = crs.crs.to_dict()
                            d.pop('no_defs')
                            p[key][model]['projection'] = d.copy()

        # writing model parameters to JSON database file
        with provider.open('w', encoding='utf-8') as fid:
            indent = 4 if args.pretty else None
            json.dump(p, fid, indent=indent, sort_keys=True)

if __name__ == '__main__':
    main()