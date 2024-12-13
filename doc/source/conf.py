# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
import datetime
# sys.path.insert(0, os.path.abspath('.'))
import importlib.metadata


# -- Project information -----------------------------------------------------

# package metadata
metadata = importlib.metadata.metadata("pyTMD")
project = metadata["Name"]
year = datetime.date.today().year
copyright = f"2017\u2013{year}, Tyler C. Sutterley"
author = 'Tyler C. Sutterley'

# The full version, including alpha/beta/rc tags
# get semantic version from setuptools-scm
version = metadata["version"]
# append "v" before the version
release = f"v{version}"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_nb",
    "numpydoc",
    'sphinxcontrib.bibtex',
    "sphinx.ext.autodoc",
    "sphinx.ext.graphviz",
    "sphinx.ext.viewcode",
    "sphinx_design",
    "sphinxarg.ext"
]

# use myst for notebooks
source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
}
nb_execution_mode = "off"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['**.ipynb_checkpoints']

# location of master document (by default sphinx looks for contents.rst)
master_doc = 'index'

# -- Configuration options ---------------------------------------------------
autosummary_generate = True
autodoc_member_order = 'bysource'
numpydoc_show_class_members = False
pygments_style = 'native'
bibtex_bibfiles = ['_assets/references.bib']
bibtex_default_style = 'plain'

# -- Options for HTML output -------------------------------------------------

# html_title = metadata["Name"]
html_short_title = metadata["Name"]
html_show_sourcelink = False
html_show_sphinx = True
html_show_copyright = True

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    "logo_only": True,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_logo = "_assets/pyTMD_logo.png"
html_static_path = ['_static']
# fetch the project urls
project_urls = {}
for project_url in metadata.get_all('Project-URL'):
    name, _, url = project_url.partition(', ')
    project_urls[name.lower()] = url
# fetch the repository url
repository_url = project_urls.get('repository')
# add html context
html_context = {
    "menu_links": [
        (
            '<i class="fa fa-github fa-fw"></i> Source Code',
            repository_url,
        ),
        (
            '<i class="fa fa-book fa-fw"></i> License',
            f"{repository_url}/blob/main/LICENSE",
        ),
        (
            '<i class="fa fa-comment fa-fw"></i> Discussions',
            f"{repository_url}/discussions",
        ),
    ],
}

# Load the custom CSS files (needs sphinx >= 1.6 for this to work)
def setup(app):
    app.add_css_file("style.css")