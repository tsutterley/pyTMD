====================
output_otis_tides.py
====================

- Writes OTIS-format tide files

Calling Sequence
----------------

.. code-block:: python

    import pyTMD.output_otis_grid
    import pyTMD.output_otis_elevation
    import pyTMD.output_otis_transport
    output_otis_grid(grid_file,xlim,ylim,hz,mz,iob,dt)
    output_otis_elevation(elevation_file,h,xlim,ylim,constituents)
    output_otis_transport(transport_file,u,v,xlim,ylim,constituents)

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/output_otis_tides.py

.. autofunction:: pyTMD.output_otis_grid

.. autofunction:: pyTMD.output_otis_elevation

.. autofunction:: pyTMD.output_otis_transport
