========
model.py
========

Class with parameters for named tide models

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/model.py

General Attributes and Methods
==============================

.. class:: model(object)

    .. attribute:: model.atl03

        HDF5 dataset string for output ATL03 tide heights

    .. attribute:: model.atl06

        HDF5 dataset string for output ATL06 tide heights

    .. attribute:: model.atl07

        HDF5 dataset string for output ATL07 tide heights

    .. attribute:: model.atl11

        HDF5 dataset string for output ATL11 tide heights

    .. attribute:: model.atl12

        HDF5 dataset string for output ATL12 tide heights

    .. attribute:: model.compressed

        Model files are gzip compressed

    .. attribute:: model.constituents

        Model constituents for ``FES`` models

    .. attribute:: model.description

        HDF5 ``description`` attribute string for output tide heights

    .. attribute:: model.directory

        Working data directory for tide models

    .. attribute:: model.format

        Model format (``OTIS``, ``ATLAS``, ``netcdf``, ``GOT``, ``FES``)

    .. attribute:: model.gla12

        HDF5 dataset string for output GLA12 tide heights

    .. attribute:: model.grid_file

        Model grid file for ``OTIS`` and ``ATLAS`` models

    .. attribute:: model.gzip

        Suffix if model is compressed

    .. attribute:: model.long_name

        HDF5 ``long_name`` attribute string for output tide heights

    .. attribute:: model.model_directory

        Full path to model directory

    .. attribute:: model.model_file

        Model constituent file or list of files

    .. attribute:: model.name

        Model name

    .. attribute:: model.projection

        Model projection for ``OTIS`` and ``ATLAS`` models

    .. attribute:: model.scale

        Model scaling factor for converting to output units

    .. attribute:: model.suffix

        Suffix if ATLAS model is ``'netcdf'`` format

    .. attribute:: model.type

        Model type (``z``, ``u``, ``v``)

    .. method:: model.pathfinder(model_file)

        Completes file paths and appends file and gzip suffixes

    .. method:: model.from_file(definition_file)

        Create a model object from an input definition file

    .. method:: model.from_dict(d)

        Create a model object from a python dictionary

    .. method:: model.to_bool(val)

        Converts strings of True/False to a boolean values
