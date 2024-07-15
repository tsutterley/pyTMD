"""
def_to_json.py (07/2024)
Converts a definition file to a json file
"""
import re
import json
import pathlib
import argparse

def read_definition_file(definition_file):
    parameters = {}
    fid = open(definition_file, 'r')
    for fileline in fid:
        # Splitting the input line between parameter name and value
        part = fileline.rstrip().split(maxsplit=1)
        # filling the parameter definition variable
        parameters[part[0]] = part[1]
    fid.close()
    return parameters

# PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Converts a definition file to a json file"
            """,
        fromfile_prefix_chars="@"
    )
    # command line parameters
    parser.add_argument('infile',
        type=pathlib.Path, nargs='+',
        help='Definition file to convert')
    parser.add_argument('--pretty', '-p',
        action='store_true',
        help='Pretty print the json file')
    parser.add_argument('--verbose', '-v',
        action='store_true',
        help='Verbose output')
    parser.add_argument('--cleanup', '-c',
        action='store_true',
        help='Remove original definition files')
    return parser

def main():
    # Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()
    # iterate over each input file
    for definition_file in args.infile:
        print(f'{definition_file} -->') if args.verbose else None
        # Reading each definition file
        parameters = read_definition_file(definition_file)
        if re.search(r';', parameters['model_file']):
            # split model into list of files for each direction
            model_file_u, model_file_v = parameters['model_file'].split(';')
            parameters['model_file'] = dict(
                u=re.split(r'[\s\,]+', model_file_u),
                v=re.split(r'[\s\,]+', model_file_v)
            )
        elif re.search(r',', parameters['model_file']):
            # split model into list of files
            parameters['model_file'] = re.split(r'[\s\,]+', parameters['model_file'])
        if 'constituents' in parameters and re.search(r',', parameters['constituents']):
            parameters['constituents'] = re.split(r'[\s\,]+', parameters['constituents'])
        if 'type' in parameters and re.search(r',', parameters['type']):
            parameters['type'] = re.split(r'[\s\,]+', parameters['type'])
        if 'compressed' in parameters:
            parameters['compressed'] = eval(parameters['compressed'])
        if 'scale' in parameters:
            parameters['scale'] = float(parameters['scale'])
        # Writing the parameters to a json file
        json_file = definition_file.with_suffix('.json')
        print(f'\t{json_file}') if args.verbose else None
        with open(json_file, 'w') as fid:
            indent = 4 if args.pretty else None
            json.dump(parameters, fid, indent=indent)
        # Removing the definition file
        if args.cleanup:
            definition_file.unlink()

if __name__ == '__main__':
    main()