compute_tide_corrections.py
===========================

 - Calculates tidal elevations for correcting elevation data
 - Can use OTIS format tidal solutions provided by Ohio State University and ESR
 - Can use Global Tide Model (GOT) solutions provided by Richard Ray at GSFC

#### Calling Sequence
```
from gravity_toolkit.compute_tide_corrections import compute_tide_corrections
tide = compute_tide_corrections(x, y, delta_time, DIRECTORY=DIRECTORY,
    MODEL=MODEL,EPOCH=(2000,1,1,0,0,0), EPSG=3031)
```

#### Inputs
 1. `x`: x-coordinates in projection EPSG
 2. `y`: y-coordinates in projection EPSG
 3. `delta_time`: seconds since EPOCH

#### Options
 - `DIRECTORY`: working data directory for tide models
 - `MODEL`: Tide model to use in correction
 - `EPOCH`: time period for calculating delta times
     - default: J2000 (seconds since 2000-01-01T00:00:00)
 - `LEAPS`: need to compute leap seconds to convert to UTC
 - `EPSG`: input coordinate system
     - default: 3031 Polar Stereographic South, WGS84

#### Outputs
 - `tide`: tide height correction reconstructed using the nodal corrections
