#!/usr/bin/env python
u"""
model.py
Written by Tyler Sutterley (11/2022)
Retrieves tide model parameters for named tide models and
    from model definition files

UPDATE HISTORY:
    Updated 11/2022: use f-strings for formatting verbose or ascii output
    Updated 06/2022: added Greenland 1km model (Gr1kmTM) to list of models
        updated citation url for Global Ocean Tide (GOT) models
    Updated 05/2022: added ESR CATS2022 to list of models
        added attribute for flexure fields being available for model
    Updated 04/2022: updated docstrings to numpy documentation format
        include utf-8 encoding in reads to be windows compliant
        set default directory to None for documentation
    Updated 03/2022: added static decorators to define model lists
    Updated 02/2022: added Arctic 2km model (Arc2kmTM) to list of models
    Updated 01/2022: added global Empirical Ocean Tide model (EOT20)
    Updated 12/2021: added TPXO9-atlas-v5 to list of available tide models
        added atl10 attributes for tidal elevation files
    Written 09/2021
"""
import warnings
import numpy as np
import pyTMD.io.model

class model(pyTMD.io.model):
    np.seterr(invalid='ignore')
    # inherit model class
    def __init__(self, *args, **kwargs):
        # raise warnings for deprecation of module
        warnings.filterwarnings("always")
        warnings.warn("Deprecated. Please use pyTMD.io instead", DeprecationWarning)
        super().__init__(*args, **kwargs)
