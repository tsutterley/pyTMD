"""
model_to_json.py (08/2024)
Converts model definitions to a json file
"""
import re
import copy
import json
import inspect
import pathlib
import argparse
import pyTMD.io

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parents[1]

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Converts model definitions to a json file"
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

def serialize(d):
    """
    Formats paths to be JSON serializable
    """
    # all paths are relative to current
    directory = pathlib.Path().absolute()
    # iterate over keys
    for key, val in d.items():
        val = copy.copy(d[key])
        if isinstance(val, pathlib.Path):
            d[key] = str(val.relative_to(directory))
        elif isinstance(val, (list, tuple)) and isinstance(val[0], pathlib.Path):
            d[key] = [str(v.relative_to(directory)) for v in val]
        elif isinstance(val, dict):
            d[key] = serialize(val)
    # return the model dictionary
    return d

def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # output dictionary
    d = dict(elevation={}, current={})

    # ocean tide models
    for n in pyTMD.io.model.ocean_elevation():
        m = pyTMD.io.model(directory=None, verify=False, format='OTIS').elevation(n)
        d['elevation'][n] = serialize(m.to_dict())
    # ATLAS netcdf models
    for n in pyTMD.io.model.ATLAS():
        m = pyTMD.io.model(directory=None, verify=False).elevation(n)
        d['elevation'][f'{n}-nc'] = serialize(m.to_dict())
    # load tide models
    for n in pyTMD.io.model.load_elevation():
        m = pyTMD.io.model(directory=None, verify=False, format='OTIS').elevation(n)
        d['elevation'][n] = serialize(m.to_dict())
    # ocean current models
    for n in pyTMD.io.model.ocean_current():
        m = pyTMD.io.model(directory=None, verify=False, format='OTIS').current(n)
        d['current'][n] = serialize(m.to_dict())
    # ATLAS netcdf models
    for n in pyTMD.io.model.ATLAS():
        m = pyTMD.io.model(directory=None, verify=False).current(n)
        d['current'][f'{n}-nc'] = serialize(m.to_dict())

    # writing model parameters to JSON database file
    json_file = filepath.joinpath('pyTMD','data','database.json')
    print(f'\t{json_file}') if args.verbose else None
    with open(json_file, 'w') as fid:
        indent = 4 if args.pretty else None
        json.dump(d, fid, indent=indent, sort_keys=True)

if __name__ == '__main__':
    main()