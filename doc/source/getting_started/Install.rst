======================
Setup and Installation
======================

Dependencies
############

``pyTMD`` is dependent on several open source programs that can be installed using
OS-specific package management systems (e.g. ``apt`` or ``homebrew``),
``conda`` or from source:

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

``pyTMD`` is available for download from the `GitHub repository <https://github.com/tsutterley/pyTMD>`_,
the `Python Package Index (pypi) <https://pypi.org/project/pyTMD/>`_,
and from `conda-forge <https://anaconda.org/conda-forge/pytmd>`_.


The simplest installation for most users will likely be using ``conda``:

.. code-block:: bash

    conda install -c conda-forge pytmd

``conda`` installed versions of ``pyTMD`` can be upgraded to the latest stable release:

.. code-block:: bash

    conda update pytmd

To use the development repository, please fork ``pyTMD`` into your own account and then clone onto your system:

.. code-block:: bash

    git clone https://github.com/tsutterley/pyTMD.git

``pyTMD`` can then be installed within the package directory using ``setuptools``:

.. code-block:: bash

    python3 setup.py install

or ``pip``

.. code-block:: bash

    python3 -m pip install --user .

The development version of ``pyTMD`` can also be installed directly from GitHub using ``pip``:

.. code-block:: bash

    python3 -m pip install --user git+https://github.com/tsutterley/pyTMD.git
