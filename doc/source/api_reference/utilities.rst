=========
utilities
=========

Download and management utilities for syncing time and auxiliary files

 - Can list a directory on a ftp host
 - Can download a file from a ftp or http host
 - Can download a file from CDDIS via https when NASA Earthdata credentials are supplied

`Source code`__

.. __: https://github.com/tsutterley/pyTMD/blob/main/pyTMD/utilities.py


General Methods
===============

.. autofunction:: pyTMD.utilities.get_data_path

.. autofunction:: pyTMD.utilities.import_dependency

.. autoclass:: pyTMD.utilities.reify
   :members:

.. autofunction:: pyTMD.utilities.file_opener

.. autofunction:: pyTMD.utilities.get_hash

.. autofunction:: pyTMD.utilities.get_git_revision_hash

.. autofunction:: pyTMD.utilities.get_git_status

.. autofunction:: pyTMD.utilities.url_split

.. autofunction:: pyTMD.utilities.roman_to_int

.. autofunction:: pyTMD.utilities.get_unix_time

.. autofunction:: pyTMD.utilities.isoformat

.. autofunction:: pyTMD.utilities.even

.. autofunction:: pyTMD.utilities.ceil

.. autofunction:: pyTMD.utilities.copy

.. autofunction:: pyTMD.utilities.check_ftp_connection

.. autofunction:: pyTMD.utilities.ftp_list

.. autofunction:: pyTMD.utilities.from_ftp

.. autofunction:: pyTMD.utilities._create_default_ssl_context

.. autofunction:: pyTMD.utilities._create_ssl_context_no_verify

.. autofunction:: pyTMD.utilities._set_ssl_context_options

.. autofunction:: pyTMD.utilities.check_connection

.. autofunction:: pyTMD.utilities.http_list

.. autofunction:: pyTMD.utilities.from_http

.. autofunction:: pyTMD.utilities.from_json

.. autofunction:: pyTMD.utilities.attempt_login

.. autofunction:: pyTMD.utilities.build_opener

.. autofunction:: pyTMD.utilities.get_token

.. autofunction:: pyTMD.utilities.list_tokens

.. autofunction:: pyTMD.utilities.revoke_token

.. autofunction:: pyTMD.utilities.check_credentials

.. autofunction:: pyTMD.utilities.cddis_list

.. autofunction:: pyTMD.utilities.from_cddis

.. autofunction:: pyTMD.utilities.iers_list

.. autofunction:: pyTMD.utilities.from_jpl_ssd

.. autofunction:: pyTMD.utilities.uhslc_list
