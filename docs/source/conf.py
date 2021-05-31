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
# sys.path.insert(0, os.path.abspath('../../src/'))

from ipywidgets.embed import DEFAULT_EMBED_REQUIREJS_URL

# -- Project information -----------------------------------------------------

project = "CobraMod"
copyright = "2021, Stefano Camborda"
author = "Stefano Camborda"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "autoapi.extension",
    "sphinx.ext.autodoc",
    # Google docstring
    "sphinxcontrib.napoleon",
    # ipynb support
    "nbsphinx",
    "sphinx.ext.intersphinx",
]

# Extensions configuration
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "cobra": ("https://cobrapy.readthedocs.io/en/latest/", None),
    "escher": ("https://escher.readthedocs.io/en/latest/", None),
}
autoapi_dirs = ["../../src/cobramod"]
autoapi_root = "module"
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "imported-members",
]
autodoc_typehints = "description"
# Add any paths that contain templates here, relative to this directory.
# templates_path = []

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

# html_theme = "sphinx_rtd_theme"
# html_logo = "img/logo_2.png"
# html_title = (
#     "CobraMod: A pathway-centric curation tool for contraint-based "
#     + "metabolic models"
# )
# html_short_titel = "CobraMod"
# html_theme_options = {"logo_only": True}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = []

# Configuration for nbsphinx
nbsphinx_execute = "never"
html_js_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js",
    DEFAULT_EMBED_REQUIREJS_URL,
]
