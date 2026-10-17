#!/usr/bin/env python
"""
model.py
Written by Tyler Sutterley (10/2026)
Retrieves tide model parameters for named tide models and
    from model definition files

PYTHON DEPENDENCIES:
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
        https://pyproj4.github.io/pyproj/
    xarray: N-D labeled arrays and datasets in Python
        https://docs.xarray.dev/en/stable/

UPDATE HISTORY:
    Updated 10/2026: change parameter names from 'm' to 'name'
    Updated 09/2026: expose default xarray engines as possible readers
    Updated 06/2026: add validate argument to from_dict method
        split old parse json function into a series of validation functions
    Updated 04/2026: add __variables__ attribute containing model variables
    Updated 03/2026: add support for FES-native netCDF4 files
    Updated 02/2026: add HTML representation for model objects using xarray
        set tidal constituent units (if unset) in a loop
        check if units are compatible with known types before setting units
    Updated 11/2025: use default cache directory if directory is None
        added crs property for model coordinate reference system
        refactor to use new simpler (flattened) database format
    Updated 08/2025: use numpy degree to radian conversions
        update node equilibrium tide estimation
        add functions for converting a model to an xarray DataTree object
    Updated 06/2025: add function for reducing list of model files
        fix extra_databases to not overwrite the default database
        add capability to use a dictionary to expand the model database
    Updated 02/2025: fixed missing grid kwarg for reading from TMD3 models
    Updated 11/2024: use Love numbers for long-period tides in node equilibrium
    Updated 10/2024: add wrapper functions to read and interpolate constants
        add functions to append node tide equilibrium values to amplitudes
        remove redundant default keyword arguments to readers and interpolators
    Updated 09/2024: use JSON database for known model parameters
        drop support for the ascii definition file format
        add file_format and nodal correction attributes
        export database as a dataclass for easier access
        added variable name and descriptions for long period tides
    Updated 08/2024: added attribute for minor constituents to infer
        allow searching over iterable glob strings in definition files
        added option to try automatic detection of definition file format
        added new TPXO10-atlas-v2 to list of models
    Updated 07/2024: added new FES2022 and FES2022_load to list of models
        added JSON format for model definition files
        use parse function from constituents class to extract names
        renamed format for ATLAS to ATLAS-compact
        renamed format for netcdf to ATLAS-netcdf
        renamed format for FES to FES-netcdf and added FES-ascii
        renamed format for GOT to GOT-ascii and added GOT-netcdf
    Updated 05/2024: make subscriptable and allow item assignment
    Updated 04/2024: append v-components of velocity only to netcdf format
    Updated 11/2023: revert TPXO9-atlas currents changes to separate dicts
    Updated 09/2023: fix scale values for TPXO9-atlas currents
    Updated 08/2023: changed ESR netCDF4 format to TMD3 format
        updated filenames for CATS2008-v2023 to final version
    Updated 06/2023: remap FES2012 e2 constituent to eps2
    Updated 04/2023: added global HAMTIDE11 model
        made ICESat, ICESat-2 and output file attributes properties
        updated model definition read function for currents
        using pathlib to define and expand tide model paths
        add basic file searching with glob strings in definition files
        add long_name and description attributes for current variables
        added exceptions for files missing when using glob patterns
        simplify TPXO9-atlas currents dictionaries to single list
    Updated 03/2023: add basic variable typing to function inputs
    Updated 12/2022: moved to io and added deprecation warning to old
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 06/2022: added Greenland 1km model (Gr1kmTM) to list of models
        updated citation url for Goddard Ocean Tide (GOT) models
    Updated 05/2022: added ESR CATS2022 to list of models
        added attribute for flexure fields being available for model
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
        set default directory to None for documentation
    Updated 03/2022: added static decorators to define model lists
    Updated 02/2022: added Arctic 2km model (Arc2kmTM) to list of models
    Updated 01/2022: added global Empirical Ocean Tide model (EOT20)
    Updated 12/2021: added TPXO9-atlas-v5 to list of available tide models
        added atl10 attributes for tidal elevation files
    Written 09/2021
"""

from __future__ import annotations

import io
import copy
import json
import pyproj
import pathlib
import warnings
import xarray as xr
import pyTMD.utilities
from collections.abc import Iterable

# suppress warnings
warnings.filterwarnings("ignore", category=UserWarning)

__all__ = ["DataBase", "load_database", "model"]

# default working data directory for tide models
_default_directory = pyTMD.utilities.get_cache_path()
# available xarray engines
_available_engines = xr.backends.list_engines()


# allow model database to be subscriptable
# and have attribute access
class DataBase:
    """pyTMD model database and parameters"""

    def __init__(self, d: dict):
        self.__dict__ = copy.deepcopy(d)

    def update(self, d: dict):
        """Update the keys of the model database"""
        self.__dict__.update(d)

    def keys(self):
        """Returns the keys of the model database"""
        return self.__dict__.keys()

    def values(self):
        """Returns the values of the model database"""
        return self.__dict__.values()

    def items(self):
        """Returns the items of the model database"""
        return self.__dict__.items()

    def __str__(self):
        """String representation of the ``DataBase`` object"""
        return str(self.__dict__)

    def __repr__(self):
        """Representation of the ``DataBase`` object"""
        return self.__str__()

    def get(self, key, default=None):
        if not hasattr(self, key) or getattr(self, key) is None:
            return default
        else:
            return getattr(self, key, default)

    def __getitem__(self, key):
        return getattr(self, key)


# PURPOSE: load the JSON database of model files
def load_database(extra_databases: list = []):
    """
    Load the ``JSON`` database of model files

    Parameters
    ----------
    extra_databases: list, default []
        A list of additional databases to load, as either
        ``JSON`` file paths or dictionaries

    Returns
    -------
    parameters: dict
        Database of model parameters
    """
    # path to model database
    database = pyTMD.utilities.get_data_path(["data", "database.json"])
    # extract JSON data
    with database.open(mode="r", encoding="utf-8") as fid:
        parameters = json.load(fid)
    # verify that extra_databases is iterable
    if isinstance(extra_databases, (str, pathlib.Path, dict)):
        extra_databases = [extra_databases]
    # load any additional databases
    for db in extra_databases:
        # use database parameters directly if a dictionary
        if isinstance(db, dict):
            extra_database = copy.copy(db)
        # otherwise load parameters from JSON file path
        else:
            # verify that extra database file exists
            db = pyTMD.utilities.Path(db)
            if not db.exists():
                raise FileNotFoundError(db)
            # extract JSON data
            with db.open(mode="r", encoding="utf-8") as fid:
                extra_database = json.load(fid)
        # Add additional models to database
        parameters.update(extra_database)
    return DataBase(parameters)


class model:
    """Retrieves tide model parameters for named models or
    from a model definition file

    Attributes
    ----------
    compressed: bool
        Model files are gzip compressed
    directory: str, pathlib.Path or None, default None
        Working data directory for tide models
    extra_databases: list, default []
        Additional databases for model parameters
    verify: bool
        Verify that all model files exist
    """

    def __init__(
        self,
        directory: str | pathlib.Path | None = None,
        **kwargs,
    ):
        # set default keyword arguments
        kwargs.setdefault("compressed", False)
        kwargs.setdefault("verify", True)
        kwargs.setdefault("extra_databases", [])
        # set initial attributes
        self.compressed = copy.copy(kwargs["compressed"])
        self.constituents = None
        self.minor = None
        # set working data directory
        self.directory = None
        if directory is not None:
            self.directory = pyTMD.utilities.Path(directory)
        # set any extra databases
        self.extra_databases = copy.copy(kwargs["extra_databases"])
        self.format = kwargs.get("format", None)
        self.name = kwargs.get("name", None)
        self.verify = copy.copy(kwargs["verify"])
        self.__parameters__ = {}

    def from_database(
        self,
        name: str,
        group: tuple = ("z", "u", "v"),
    ):
        """
        Create a model object from database of known tidal models

        Parameters
        ----------
        name: str
            Model name
        group: tuple, default ('z', 'u', 'v')
            Model variable(s) to extract
        """
        # set working data directory if unset
        if self.directory is None:
            self.directory = pyTMD.utilities.Path(_default_directory)
        # select between known tide models
        parameters = load_database(extra_databases=self.extra_databases)
        # try to extract parameters for model
        try:
            self.from_dict(parameters[name], validate=False)
        except (ValueError, KeyError, AttributeError) as exc:
            raise ValueError(f"Unlisted tide model {name}") from exc
        # verify model types to extract
        if isinstance(group, str):
            group = (group,)
        # verify paths
        for g in group:
            # verify model group is valid
            g = g.lower()
            # skip if model group is unavailable
            if not hasattr(self, g):
                continue
            # validate paths: grid file for OTIS, ATLAS models
            if hasattr(self[g], "grid_file"):
                self[g].grid_file = self.pathfinder(self[g].grid_file)
            # validate paths: model constituent files
            self[g].model_file = self.pathfinder(self[g].model_file)
        # return the model parameters
        self.validate_name()
        self.validate_format()
        self.validate_projection()
        # set dictionary of parameters
        self.__parameters__ = self.to_dict(serialize=True)
        return self

    def from_file(
        self,
        definition_file: str | pathlib.Path | io.IOBase,
        **kwargs,
    ):
        """
        Create a model object from an input definition file

        Parameters
        ----------
        definition_file: str, pathlib.Path or io.IOBase
            Model definition file for creating model object
        """
        # load and parse definition file
        if isinstance(definition_file, io.IOBase):
            self.parse_json(definition_file)
        elif isinstance(definition_file, (str, pathlib.Path)):
            definition_file = pyTMD.utilities.Path(definition_file)
            with definition_file.open(mode="r", encoding="utf8") as fid:
                self.parse_json(fid)
        # return the model object
        return self

    def from_dict(self, d: dict, validate: bool = True):
        """
        Create a model object from a dictionary of parameters

        Parameters
        ----------
        d: dict
            Model object parameters
        validate: bool, default True
            Validate the model parameters after loading
        """
        # copy model parameters
        self.__parameters__ = copy.deepcopy(d)
        for key, val in d.items():
            if isinstance(val, dict) and key not in ("projection",):
                setattr(self, key, DataBase(val))
            else:
                setattr(self, key, copy.copy(val))
        # validate the model parameters
        if validate:
            self.validate_model()
        # return the model parameters
        return self

    def to_dict(self, **kwargs):
        """
        Create a dictionary from a model object

        Parameters
        ----------
        fields: list, default all
            List of model attributes to output
        serialize: bool, default False
            Serialize dictionary for ``JSON`` output
        """
        # default fields
        keys = ["name", "format", "projection", "reference", "z", "u", "v"]
        # set default keyword arguments
        kwargs.setdefault("fields", keys)
        kwargs.setdefault("serialize", False)
        # output dictionary
        d = {}
        # for each field
        for key in kwargs["fields"]:
            if hasattr(self, key) and getattr(self, key) is not None:
                d[key] = getattr(self, key)
        # serialize dictionary for JSON output
        if kwargs["serialize"]:
            d = self.serialize(d)
        # return the model dictionary
        return d

    @property
    def gzip(self) -> str:
        """Returns suffix for ``gzip`` compression"""
        return ".gz" if self.compressed else ""

    @property
    def corrections(self) -> str:
        """
        Returns the corrections group for the model
        """
        part1, _, part2 = self.format.partition("-")
        if self.format in ("GOT-ascii",):
            return "perth3"
        else:
            return part1

    @property
    def file_format(self) -> str:
        """
        Returns the file format for the model
        """
        part1, _, part2 = self.format.partition("-")
        if self.format in ("ATLAS-compact"):
            return part1
        elif "-" in self.format:
            return part2
        else:
            return self.format

    @property
    def multifile(self) -> bool:
        """Returns if the model uses individual files for constituents"""
        # try to find a valid mode group
        for g in ("z", "u", "v"):
            # verify case of model group
            g = g.lower()
            # skip if model group is unavailable
            if not hasattr(self, g):
                continue
            return isinstance(self[g].model_file, list)

    @property
    def crs(self):
        """Coordinate reference system of the model"""
        # default is EPSG:4326 (WGS84)
        CRS = self.get("projection", 4326)
        return pyproj.CRS.from_user_input(CRS)

    @staticmethod
    def known_formats(**kwargs) -> list:
        """
        Returns list of known model formats
        """
        # load the database of model parameters
        parameters = load_database(**kwargs)
        # extract all known formats
        format_list = []
        for model, val in parameters.items():
            format_list.append(val["format"])
        # append available xarray engines
        for key, val in _available_engines.items():
            format_list.append(key)
        # return unique list of formats
        return sorted(set(format_list))

    @staticmethod
    def ocean_elevation(**kwargs) -> list:
        """
        Returns list of ocean tide elevation models
        """
        # load the database of model parameters
        parameters = load_database(**kwargs)
        # extract all known ocean tide elevation models
        model_list = []
        for model, val in parameters.items():
            if ("z" in val) and (val["z"]["variable"] == "tide_ocean"):
                model_list.append(model)
            if ("z" in val) and (val["z"]["variable"] == "tide_lpe"):
                model_list.append(model)
        # return unique list of models
        return sorted(set(model_list))

    @staticmethod
    def load_elevation(**kwargs) -> list:
        """
        Returns list of load tide elevation models
        """
        # load the database of model parameters
        parameters = load_database(**kwargs)
        # extract all known load tide elevation models
        model_list = []
        for model, val in parameters.items():
            if ("z" in val) and (val["z"]["variable"] == "tide_load"):
                model_list.append(model)
        # return unique list of models
        return sorted(set(model_list))

    @staticmethod
    def ocean_current(**kwargs) -> list:
        """
        Returns list of tidal current models
        """
        # load the database of model parameters
        parameters = load_database(**kwargs)
        # extract all known ocean tide current models
        model_list = []
        for model, val in parameters.items():
            if ("u" in val) or ("v" in val):
                model_list.append(model)
        # return unique list of models
        return sorted(set(model_list))

    @staticmethod
    def OTIS(**kwargs) -> list:
        """
        Returns list of OTIS formatted models
        """
        # load the database of model parameters
        parameters = load_database(**kwargs)
        # extract all known OTIS models
        model_list = []
        for model, val in parameters.items():
            if val["format"] == "OTIS":
                model_list.append(model)
        # return unique list of models
        return sorted(set(model_list))

    @staticmethod
    def ATLAS_compact(**kwargs) -> list:
        """
        Returns list of ATLAS-compact formatted models
        """
        # load the database of model parameters
        parameters = load_database(**kwargs)
        # extract all known ATLAS-compact models
        model_list = []
        for model, val in parameters.items():
            if val["format"] == "ATLAS-compact":
                model_list.append(model)
        # return unique list of models
        return sorted(set(model_list))

    @staticmethod
    def TMD3(**kwargs) -> list:
        """
        Returns list of TMD3 formatted models
        """
        # load the database of model parameters
        parameters = load_database(**kwargs)
        # extract all known TMD3 models
        model_list = []
        for model, val in parameters.items():
            if val["format"] == "TMD3":
                model_list.append(model)
        # return unique list of models
        return sorted(set(model_list))

    @staticmethod
    def ATLAS(**kwargs) -> list:
        """
        Returns list of ATLAS-netcdf formatted models
        """
        # load the database of model parameters
        parameters = load_database(**kwargs)
        # extract all known TMD3 models
        model_list = []
        for model, val in parameters.items():
            if val["format"] == "ATLAS-netcdf":
                model_list.append(model)
        # return unique list of models
        return sorted(set(model_list))

    @staticmethod
    def GOT(**kwargs) -> list:
        """
        Returns list of GOT formatted models
        """
        # load the database of model parameters
        parameters = load_database(**kwargs)
        # extract all known GOT-ascii or GOT-netcdf models
        model_list = []
        for model, val in parameters.items():
            if val["format"] in ("GOT-ascii", "GOT-netcdf"):
                model_list.append(model)
        # return unique list of models
        return sorted(set(model_list))

    @staticmethod
    def FES(**kwargs) -> list:
        """
        Returns list of FES formatted models
        """
        # load the database of model parameters
        parameters = load_database(**kwargs)
        # extract all known FES-ascii or FES-netcdf models
        model_list = []
        for model, val in parameters.items():
            if val["format"] in ("FES-ascii", "FES-netcdf", "FES-native"):
                model_list.append(model)
        # return unique list of models
        return sorted(set(model_list))

    def pathfinder(
        self,
        model_file: str | pathlib.Path | list,
    ):
        """
        Completes file paths and appends ``gzip`` suffix

        Parameters
        ----------
        model_file: str, pathlib.Path or list
            Model file(s) to complete
        """
        # set working data directory if unset
        if self.directory is None:
            self.directory = pyTMD.utilities.Path(_default_directory)
        # complete model file paths
        if isinstance(model_file, list):
            output_file = [self.pathfinder(f) for f in model_file]
            valid = all([f.exists() for f in output_file])
        elif isinstance(model_file, str):
            output_file = self.directory.joinpath(
                "".join([model_file, self.gzip])
            )
            valid = output_file.exists()
        # check that (all) output files exist
        if self.verify and not valid and not self.compressed:
            # try seeing if there are compressed files
            self.compressed = True
            output_file = self.pathfinder(model_file)
        elif self.verify and not valid:
            raise FileNotFoundError(output_file)
        # return the complete output path
        return output_file

    def parse_json(self, fid: io.IOBase):
        """
        Load and parse ``JSON`` definition file

        Parameters
        ----------
        fid: io.IOBase
            Open definition file object
        """
        # attempt to read and parse a JSON file
        try:
            # load JSON file
            parameters = json.load(fid)
        except json.decoder.JSONDecodeError as exc:
            raise IOError("Cannot load model definition file") from exc
        # set model parameters from dictionary
        self.from_dict(parameters)

    def validate_model(self):
        """Validates the model parameters for a model object"""
        # validate model name
        self.validate_name()
        # validate model format
        self.validate_format()
        # validate model file paths
        self.validate_paths()
        # verify that projection attribute exists for projected models
        self.validate_projection()
        # set dictionary of parameters
        self.__parameters__ = self.to_dict(serialize=True)

    def validate_format(self):
        """Asserts that the model format is a known type"""
        # known remapped cases
        mapping = [
            ("ATLAS", "ATLAS-compact"),
            ("netcdf", "ATLAS-netcdf"),
            ("FES", "FES-netcdf"),
            ("GOT", "GOT-ascii"),
        ]
        # iterate over known remapped cases
        for m in mapping:
            # check if tide model is a remapped case
            if self.format == m[0]:
                self.format = m[1]
        # assert that tide model is a known format
        if self.format not in self.known_formats():
            raise ValueError(f"Unknown model format {self.format}")

    def validate_name(self):
        """Asserts that the model name is specified"""
        # verify model name is specified
        if not hasattr(self, "name") or self.name is None:
            raise ValueError("Model name is not specified")

    def validate_paths(self):
        """Resolves paths to model files for specific model formats"""
        # resolve grid files for projected (OTIS, ATLAS) model formats
        if self.format in ("OTIS", "ATLAS-compact", "ATLAS-netcdf"):
            for g in ("z", "u", "v"):
                # check that model group is available
                if not hasattr(self, g):
                    continue
                if (
                    not hasattr(self[g], "grid_file")
                    or self[g].grid_file is None
                ):
                    raise ValueError(
                        f"Grid file is not specified for group {g}"
                    )
                # check if grid file is relative or absolute
                if self.directory is not None:
                    self[g].grid_file = self.directory.joinpath(
                        self[g].grid_file
                    )
                else:
                    self[g].grid_file = pyTMD.utilities.Path(self[g].grid_file)
        # resolve model constituent files
        for g in ("z", "u", "v"):
            # check that model group is available
            if not hasattr(self, g):
                continue
            # get model files for model group
            if self.directory is not None:
                # use glob strings to find files in directory
                glob_string = copy.copy(self[g].model_file)
                # search singular glob string or iterable glob strings
                if isinstance(glob_string, str):
                    # singular glob string
                    self[g].model_file = list(self.directory.glob(glob_string))
                elif isinstance(glob_string, Iterable):
                    # iterable glob strings
                    self[g].model_file = []
                    for p in glob_string:
                        self[g].model_file.extend(self.directory.glob(p))
            elif isinstance(self[g].model_file, list):
                # resolve paths to model files
                self[g].model_file = [
                    pyTMD.utilities.Path(f) for f in self[g].model_file
                ]
            else:
                # fully defined single file case
                self[g].model_file = pyTMD.utilities.Path(self[g].model_file)

    def validate_projection(self):
        """Asserts that the projection attribute is specified for projected models"""
        if self.format in ("OTIS", "ATLAS-compact", "TMD3") and (
            not hasattr(self, "projection") or self.projection is None
        ):
            raise ValueError("Projection is not specified for model")

    def serialize(self, d: dict):
        """
        Encodes dictionary to be ``JSON`` serializable

        Parameters
        ----------
        d: dict
            Parameters to serialize
        """
        # iterate over keys
        for key, val in d.items():
            val = copy.copy(d[key])
            if isinstance(val, pathlib.Path):
                d[key] = str(val)
            elif isinstance(val, (list, tuple)) and isinstance(
                val[0], pathlib.Path
            ):
                d[key] = [str(v) for v in val]
            elif isinstance(val, dict):
                d[key] = self.serialize(val)
            elif isinstance(val, DataBase):
                d[key] = self.serialize(val.__dict__)
        # return the model dictionary
        return d

    def parse_constituents(
        self,
        group: str = "z",
        **kwargs,
    ) -> list:
        """
        Parses tide model files for a list of model constituents

        Parameters
        ----------
        group: str, default 'z'
            Model group
        """
        if isinstance(self[group].model_file, (str, pathlib.Path)):
            # single file case
            self.constituents = [
                self.parse_file(self[group].model_file, **kwargs)
            ]
        elif isinstance(self[group].model_file, list):
            # multiple file case
            self.constituents = [
                self.parse_file(f, **kwargs) for f in self[group].model_file
            ]
        # return the model parameters
        return self

    @staticmethod
    def parse_file(
        model_file: str | pathlib.Path,
        raise_error: bool = False,
    ):
        """
        Parses a model file for a tidal constituent name

        Parameters
        ----------
        model_file: str or pathlib.Path
            Tide model file to parse
        raise_error: bool, default False
            Raise ``ValueError`` if constituent is not found in file name

        Returns
        -------
        constituent: str or list
            Constituent name
        """
        # import constituents parser
        from pyTMD.constituents import _parse_name

        # convert to pathlib.Path
        model_file = pathlib.Path(model_file)
        # try to parse the constituent name from the file name
        try:
            return _parse_name(model_file.name)
        except ValueError:
            pass
        # if no constituent name is found
        if raise_error:
            raise ValueError(f"Constituent not found in file {model_file}")
        else:
            return None

    def reduce_constituents(
        self,
        constituents: str | list,
        group: str = "z",
    ):
        """
        Reduce model files to a subset of constituents

        Parameters
        ----------
        constituents: str or list
            List of constituents names
        group: str, default 'z'
            Model group
        """
        # if no constituents are specified, return self
        if constituents is None:
            return None
        # verify that constituents is a list
        if isinstance(constituents, str):
            constituents = [constituents]
        # parse constituents from model files
        try:
            self.parse_constituents(group=group, raise_error=True)
        except ValueError as exc:
            return None
        # only run for multiple files
        if isinstance(self[group].model_file, list):
            # multiple file case
            # filter model files to constituents
            self[group].model_file = [
                self[group].model_file[self.constituents.index(c)]
                for c in constituents
                if (c in self.constituents)
            ]
        # update list of constituents
        self.parse_constituents(group=group)
        # return self
        return self

    def open_dataset(self, **kwargs):
        """
        Open model files as an xarray Dataset

        Parameters
        ----------
        kwargs: dict
            Additional keyword arguments for opening model files

        Returns
        -------
        ds: xarray.Dataset
            Tide model data
        """
        # import tide model functions
        from pyTMD.io import OTIS, ATLAS, GOT, FES

        # set default keyword arguments
        kwargs.setdefault("group", "z")
        kwargs.setdefault("use_default_units", True)
        kwargs.setdefault("append_node", False)
        kwargs.setdefault("compressed", self.compressed)
        kwargs.setdefault("constituents", None)
        # model group
        group = kwargs["group"].lower()
        if group not in ("z", "u", "v"):
            raise ValueError(f"Invalid model group {group}")
        # reduce constituents if specified
        self.reduce_constituents(kwargs["constituents"])
        # extract model file
        model_file = self[group].get("model_file")
        if self.format in ("OTIS", "ATLAS-compact", "TMD3"):
            # open OTIS/TMD3/ATLAS-compact files as xarray Dataset
            ds = OTIS.open_dataset(
                model_file,
                grid_file=self[group].get("grid_file"),
                format=self.file_format,
                crs=self.crs,
                **kwargs,
            )
        elif self.format in ("ATLAS-netcdf",):
            # open ATLAS netCDF4 files as xarray Dataset
            ds = ATLAS.open_dataset(
                model_file,
                grid_file=self[group].get("grid_file"),
                format=self.file_format,
                **kwargs,
            )
        elif self.format in ("GOT-ascii", "GOT-netcdf"):
            # open GOT ASCII/netCDF4 files as xarray Dataset
            ds = GOT.open_mfdataset(
                model_file, format=self.file_format, **kwargs
            )
        elif self.format in ("FES-ascii", "FES-netcdf", "FES-native"):
            # open FES ASCII/netCDF4 files as xarray Dataset
            ds = FES.open_mfdataset(
                model_file, format=self.file_format, **kwargs
            )
        elif self.multifile and self.format in _available_engines.keys():
            # try to read a multi-file using default xarray engines
            ds = xr.open_mfdataset(model_file, engine=self.format, **kwargs)
        elif self.format in _available_engines.keys():
            # try to read using default xarray engines
            ds = xr.open_dataset(model_file, engine=self.format, **kwargs)
        else:
            raise ValueError(f"Unknown reader for format {self.format}")
        # append node equilibrium tide if not in constituents list
        if kwargs["append_node"] and ("node" not in ds.tmd.constituents):
            # calculate and append node equilibrium tide
            ds = ds.tmd.node_equilibrium()
        # add attributes
        ds.attrs["source"] = self.name
        # add coordinate reference system to Dataset
        ds.attrs["crs"] = self.crs.to_dict()
        # check if units attribute can be parsed and is a known type
        # if units cannot be parsed: use value defined in the model database
        for c in ds.tmd.constituents:
            if not ds[c].tmd._has_compatible_units:
                ds[c].attrs["units"] = self[group].units
        # convert to default units
        if kwargs["use_default_units"]:
            ds = ds.tmd.to_default_units()
        # return xarray dataset
        return ds

    def open_datatree(
        self,
        group: tuple = ("z", "u", "v"),
        **kwargs,
    ):
        """
        Open model files as an xarray DataTree

        Parameters
        ----------
        group: tuple, default ('z', 'u', 'v')
            List of model types to extract
        kwargs: dict
            Additional keyword arguments for opening model files

        Returns
        -------
        dtree: xr.DataTree
            Tide model data
        """
        # output dictionary of xarray Datasets
        ds = {}
        # try to read model files
        for g in group:
            # skip if model group is unavailable
            if not hasattr(self, g.lower()):
                continue
            # open xarray Dataset
            ds[g] = self.open_dataset(group=g, **kwargs)
        # create xarray DataTree from dictionary
        dtree = xr.DataTree.from_dict(ds)
        # return the model xarray DataTree
        return dtree

    def __str__(self):
        """String representation of the ``io.model`` object"""
        properties = ["pyTMD.io.model"]
        properties.append(f"    name: {self.name}")
        return "\n".join(properties)

    def __repr__(self):
        """Representation of the ``io.model`` object"""
        return self.__str__()

    def _repr_html_(self):
        """HTML representation of the ``io.model`` object"""
        header = "pyTMD.io.model"
        header_components = [f"<div class='xr-obj-type'>{header}</div>"]
        sections = []
        data_vars = self.__variables__.copy()
        parameters = {
            k: v for k, v in self.__parameters__.items() if k not in data_vars
        }
        sections.append(xr.core.formatting_html.attr_section(parameters))
        for v in data_vars:
            sections.append(
                xr.core.formatting_html._mapping_section(
                    mapping=self.__parameters__[v],
                    name=f"{v}-Attributes",
                    details_func=xr.core.formatting_html.summarize_attrs,
                    max_items_collapse=0,
                    expand_option_name="display_expand_attrs",
                )
            )
        return xr.core.formatting_html._obj_repr(
            self, header_components, sections
        )

    @property
    def __variables__(self):
        """List of model variables"""
        return [k for k in ("z", "u", "v") if k in self.__parameters__]

    def get(self, key, default=None):
        return getattr(self, key, default) or default

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)
