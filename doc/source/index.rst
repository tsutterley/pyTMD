===================
pyTMD Documentation
===================

Welcome to the documentation for ``pyTMD``, a Python-based tidal prediction software.
This documentation is intended to explain how to compute ocean, solid Earth, load and pole tide variations using the set of ``pyTMD`` programs.

Introduction
------------

.. grid:: 2 2 4 4
    :padding: 0

    .. grid-item-card::  Installation
      :text-align: center
      :link: ./getting_started/Install.html

      :material-outlined:`download;5em`

    .. grid-item-card::  Getting Started
      :text-align: center
      :link: ./getting_started/Getting-Started.html

      :material-outlined:`hiking;5em`

    .. grid-item-card::  Background
      :text-align: center
      :link: ./getting_started/Background.html

      :material-outlined:`library_books;5em`

    .. grid-item-card::  Examples
      :text-align: center
      :link: ./user_guide/Examples.html

      :material-outlined:`apps;5em`

Contribute
----------

.. grid:: 2 2 4 4
    :padding: 0

    .. grid-item-card::  Guidelines
      :text-align: center
      :link: ./getting_started/Contributing.html

      :material-outlined:`groups;5em`

    .. grid-item-card::  Code of Conduct
      :text-align: center
      :link: ./getting_started/Code-of-Conduct.html

      :material-outlined:`gavel;5em`

    .. grid-item-card::  Discussions
      :text-align: center
      :link: https://github.com/tsutterley/pyTMD/discussions

      :material-outlined:`forum;5em`

    .. grid-item-card::  Issues
      :text-align: center
      :link: https://github.com/tsutterley/pyTMD/issues

      :material-outlined:`bug_report;5em`

Project Details
---------------

.. grid:: 2 2 4 4
    :padding: 0

    .. grid-item-card::  Release Notes
      :text-align: center
      :link: ./release_notes/Release-Notes.html

      :material-outlined:`edit_note;5em`

    .. grid-item-card::  Contributors
      :text-align: center
      :link: ./project/Contributors.html

      :material-outlined:`diversity_1;5em`

    .. grid-item-card::  License
      :text-align: center
      :link: ./project/Licenses.html

      :material-outlined:`balance;5em`

    .. grid-item-card::  Citation Information
      :text-align: center
      :link: ./project/Citations.html

      :material-outlined:`alternate_email;5em`


.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Getting Started

    getting_started/Install.rst
    getting_started/Getting-Started.rst
    getting_started/Background.rst
    getting_started/Contributing.rst
    getting_started/Code-of-Conduct.rst
    getting_started/Resources.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: User Guide

    user_guide/Examples.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: API Reference

    api_reference/arguments.rst
    api_reference/astro.rst
    api_reference/compute.rst
    api_reference/crs.rst
    api_reference/ellipse.rst
    api_reference/interpolate.rst
    api_reference/io/io.rst
    api_reference/math.rst
    api_reference/predict.rst
    api_reference/solve/solve.rst
    api_reference/spatial.rst
    api_reference/utilities.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Utilities

    api_reference/arcticdata_tides.rst
    api_reference/aviso_fes_tides.rst
    api_reference/gsfc_got_tides.rst
    api_reference/reduce_OTIS_files.rst
    api_reference/usap_cats_tides.rst
    api_reference/verify_box_tpxo.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Use Cases

    api_reference/compute_LPET_elevations.rst
    api_reference/compute_LPT_displacements.rst
    api_reference/compute_OPT_displacements.rst
    api_reference/compute_SET_displacements.rst
    api_reference/compute_tidal_currents.rst
    api_reference/compute_tidal_elevations.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Project Details

    project/Contributors.rst
    project/Licenses.rst
    project/Citations.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Release Notes

    release_notes/Release-Notes.rst

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Bibliography

    project/Bibliography.rst
