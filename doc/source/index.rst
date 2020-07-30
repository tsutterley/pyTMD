pyTMD
=====

Python-based tidal prediction software that reads OTIS, GOT and FES formatted
tidal solutions for predicting ocean and load tides and IERS conventions for
calculating radial pole tide displacements

.. graphviz::
    :caption: pyTMD Ocean and Load Tide Framework
    :align: center

    digraph {
        M [label="Tide Model" shape=box style="filled" color="darkorchid"]
        T [label="pyTMD" shape=box style="filled" color="gray"]
        P [label="Tide Predictions" shape=box style="filled" color="mediumseagreen"]
        H [label="Tide Heights" shape=box style="filled" color="mediumseagreen"]
        C [label="Average Tidal Currents" shape=box style="filled" color="mediumseagreen"]
        M -> T
        T -> P
        T -> H
        T -> C
    }

.. toctree::
    :maxdepth: 2
    :caption: Getting Started:

    getting_started/Install.md
    getting_started/Resources.md
    getting_started/Citations.md

.. toctree::
    :maxdepth: 1
    :caption: User Guide:

    user_guide/aviso_fes_tides.md
    user_guide/bilinear_interp.md
    user_guide/calc_astrol_longitudes.md
    user_guide/calc_delta_time.md
    user_guide/calc_iers_mean_pole.md
    user_guide/compute_LPT_displacements.md
    user_guide/compute_LPT_icebridge_data.md
    user_guide/compute_OPT_displacements.md
    user_guide/compute_OPT_icebridge_data.md
    user_guide/compute_tidal_elevations.md
    user_guide/compute_tide_corrections.md
    user_guide/compute_tides_icebridge_data.md
    user_guide/compute_tides_ICESat2_ATL03.md
    user_guide/compute_tides_ICESat2_ATL06.md
    user_guide/compute_tides_ICESat2_ATL07.md
    user_guide/compute_tides_ICESat2_ATL12.md
    user_guide/convert_calendar_decimal.md
    user_guide/convert_julian.md
    user_guide/convert_ll_xy.md
    user_guide/count_leap_seconds.md
    user_guide/iers_mean_pole.md
    user_guide/infer_minor_corrections.md
    user_guide/load_constituent.md
    user_guide/load_nodal_corrections.md
    user_guide/output_otis_tides.md
    user_guide/predict_tidal_ts.md
    user_guide/predict_tide_drift.md
    user_guide/predict_tide.md
    user_guide/read_FES_model.md
    user_guide/read_GOT_model.md
    user_guide/read_iers_EOP.md
    user_guide/read_netcdf_model.md
    user_guide/read_ocean_pole_tide.md
    user_guide/read_tide_model.md
    user_guide/reduce_OTIS_files.md
    user_guide/tidal_ellipse.md
