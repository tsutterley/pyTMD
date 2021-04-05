======================
Setup and Installation
======================

Dependencies
############
``pyTMD`` is dependent on open source programs that can be installed using OS-specific package management systems,
`conda <https://anaconda.org/conda-forge/repo>`_ or from source:

- `GDAL <https://gdal.org/index.html>`_
- `GEOS <https://trac.osgeo.org/geos>`_
- `PROJ <https://proj.org/>`_
- `HDF5 <https://www.hdfgroup.org>`_
- `netCDF <https://www.unidata.ucar.edu/software/netcdf>`_
- `libxml2 <http://xmlsoft.org/>`_
- `libxslt <http://xmlsoft.org/XSLT/>`_

The version of GDAL used within ``pyTMD`` will match the version of the installed C program.
The path to the C program that will be used with ``pyTMD`` is given by:

.. code-block:: bash

    gdal-config --datadir

The ``pyTMD`` installation uses the ``gdal-config`` routines to set the GDAL package version.

Installation
############
Presently ``pyTMD`` is available for use as a `GitHub repository <https://github.com/tsutterley/pyTMD>`_ and
from the `Python Package Index (pypi) <https://pypi.org/project/pyTMD/>`_.
The contents of the repository can be download as a
`zipped file <https://github.com/tsutterley/pyTMD/archive/main.zip>`_  or cloned.

To use this repository, please fork into your own account and then clone onto your system:

.. code-block:: bash

    git clone https://github.com/tsutterley/pyTMD.git

Can then install using ``setuptools``:

.. code-block:: bash

    python3 setup.py install

or ``pip``

.. code-block:: bash

    python3 -m pip install --user .

Alternatively can install the ``pyTMD`` utilities directly from GitHub with ``pip``:

.. code-block:: bash

    python3 -m pip install --user git+https://github.com/tsutterley/pyTMD.git

Executable versions of this repository can also be tested using
`Binder <https://mybinder.org/v2/gh/tsutterley/pyTMD/main>`_ or
`Pangeo <https://binder.pangeo.io/v2/gh/tsutterley/pyTMD/main>`_.
