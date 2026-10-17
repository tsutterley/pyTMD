"""
Input/output functions for reading and writing tidal data
"""

import os
import pathlib
from . import ATLAS
from . import FES
from . import GOT
from . import OTIS
from . import IERS
from . import NOAA
from . import dataset
from .model import model, load_database

# set environmental variable for anonymous s3 access
os.environ["AWS_NO_SIGN_REQUEST"] = "YES"


def open_dataset(
    name: str,
    directory: str | pathlib.Path | None = None,
    group: str = "z",
    **kwargs,
):
    """
    Wrapper function for opening files as an xarray Dataset
    for a tide model in the database

    Parameters
    ----------
    name: str
        Model name
    directory: str, pathlib.Path or None, default None
        Working data directory for tide models
    group: str, default "z"
        Model variable to extract
    kwargs: dict
        Additional keyword arguments for opening model files
    """
    m = model(directory=directory).from_database(name, group=group)
    ds = m.open_dataset(group=group, **kwargs)
    return ds


def open_datatree(
    name: str,
    directory: str | pathlib.Path | None = None,
    group: tuple = ("z", "u", "v"),
    **kwargs,
):
    """
    Wrapper function for opening files as an xarray DataTree
    for a tide model in the database

    Parameters
    ----------
    name: str
        Model name
    directory: str, pathlib.Path or None, default None
        Working data directory for tide models
    group: tuple, default ('z', 'u', 'v')
        Model variable(s) to extract
    kwargs: dict
        Additional keyword arguments for opening model files
    """
    m = model(directory=directory).from_database(name, group=group)
    dtree = m.open_datatree(group=group, **kwargs)
    return dtree
