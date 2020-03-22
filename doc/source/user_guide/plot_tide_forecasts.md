plot_tide_forecasts.py
======================

- Plots the daily tidal displacements for a given location
- Can use OTIS format tidal solutions provided by Ohio State University and ESR
- Can use Global Tide Model (GOT) solutions provided by Richard Ray at GSFC

#### Calling Sequence
```bash
python plot_tide_forecasts.py --directory=<path_to_directory> --tide=<model> \
    --coordinates=<coordinates> --date=<date_to_forecast>
```

#### Command Line Options
 - `-D X`, `--directory=X`: Working data directory
 - `-C X`, `--coordinates=X`: latitude and longitude of point to forecast
 - `--date=X`: date to forecast in ISO format (YYYY-MM-DD)
 - `-T X`, `--tide=X`: Tide model to use
     * CATS0201
     * CATS2008
     * CATS2008_load
     * TPXO9-atlas
     * TPXO9.1
     * TPXO8-atlas
     * TPXO7.2
     * TPXO7.2_load
     * AODTM-5
     * AOTIM-5
     * AOTIM-5-2018
     * GOT4.7
     * GOT4.7_load
     * GOT4.8
     * GOT4.8_load
     * GOT4.10
     * GOT4.10_load
