"""
_providers_to_database.py (08/2024)
Compress providers to a single JSON database
"""
import re
import copy
import json
import inspect
import pathlib
import argparse
import pyTMD.utilities

# current file path
filename = inspect.getframeinfo(inspect.currentframe()).filename
filepath = pathlib.Path(filename).absolute().parent

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Compress providers to a single JSON database"
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

def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    # output dictionary
    d = dict(elevation={}, current={})
    # find providers
    providers = [f for f in filepath.iterdir() if (f.suffix == '.json')]
    # for each provider
    for provider in providers:
        with provider.open('r', encoding='utf-8') as fid:
            p = json.load(fid)
            for key, val in p.items():
                d[key].update(val)

    # writing model parameters to JSON database file
    json_file = pyTMD.utilities.get_data_path(['data','database.json'])
    print(f'\t{json_file}') if args.verbose else None
    with open(json_file, 'w') as fid:
        indent = 4 if args.pretty else None
        json.dump(d, fid, indent=indent, sort_keys=True)

if __name__ == '__main__':
    main()