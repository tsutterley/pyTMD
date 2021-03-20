output_otis_tides.py
====================

- Writes OTIS-format tide files

#### Calling Sequence
```python
from pyTMD.output_otis_tides import output_otis_grid, output_otis_elevation \
    output_otis_transport
output_otis_grid(grid_file,xlim,ylim,hz,mz,iob,dt)
output_otis_elevation(elevation_file,h,xlim,ylim,constituents)
output_otis_transport(transport_file,u,v,xlim,ylim,constituents)
```
[Source code](https://github.com/tsutterley/pyTMD/blob/main/pyTMD/output_otis_tides.py)

#### Arguments
 - `output_otis_grid()`
    1. `grid_file`: output OTIS grid file name
    2. `xlim`: x-coordinate grid-cell edges of output grid
    3. `ylim`: y-coordinate grid-cell edges of output grid
    4. `hz`: bathymetry
    5. `mz`: land/water mask
    6. `iob`: open boundary index
    7. `dt`: time step

 - `output_otis_elevation()`
    1. `elevation_file`: output OTIS elevation file name
    2. `h`: Eulerian form of tidal height oscillation
    3. `xlim`: x-coordinate grid-cell edges of output grid
    4. `ylim`: y-coordinate grid-cell edges of output grid
    5. `constituents`: tidal constituent IDs

 - `output_otis_transport()`
    1. `transport_file`: output OTIS transport file name
    2. `u`: Eulerian form of tidal zonal transport oscillation
    3. `v`: Eulerian form of tidal meridional transport oscillation
    4. `xlim`: x-coordinate grid-cell edges of output grid
    5. `ylim`: y-coordinate grid-cell edges of output grid
    6. `constituents`: tidal constituent IDs
