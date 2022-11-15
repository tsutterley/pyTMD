==================
aviso_fes_tides.py
==================

- Downloads the FES (Finite Element Solution) global tide model from `AVISO <https://www.aviso.altimetry.fr/en/data/products/auxiliary-products/global-tide-fes.html>`_
- FES outputs are licensed for scientific purposes only
- Decompresses the model tar files into the constituent files and auxiliary files
- Must have `data access to tide models from AVISO <https://www.aviso.altimetry.fr/en/data/data-access.html>`_

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/scripts/aviso_fes_tides.py

Calling Sequence
################

.. argparse::
    :filename: aviso_fes_tides.py
    :func: arguments
    :prog: aviso_fes_tides.py
    :nodescription:
    :nodefault:
