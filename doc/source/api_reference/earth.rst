===============
``pyTMD.earth``
===============

- Calculates Earth reference and deformation parameters
- Provides gravitational and ellipsoidal parameters from :cite:t:`HofmannWellenhof:2006hy,Petit:2010tp`

`Source code`__

.. __: https://github.com/pyTMD/pyTMD/blob/main/pyTMD/earth.py

General Methods
===============

.. autoclass:: pyTMD.earth.datum
   :members:

.. autofunction:: pyTMD.earth.love_numbers

.. autofunction:: pyTMD.earth.complex_love_numbers

.. autofunction:: pyTMD.earth.degree_love_numbers

.. autofunction:: pyTMD.earth.load_love_numbers

.. autofunction:: pyTMD.earth.adjust_load_love_numbers

.. autofunction:: pyTMD.earth.greens_function

.. autofunction:: pyTMD.earth._melchior_table_52

.. autofunction:: pyTMD.earth._infconv
