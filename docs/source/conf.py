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
copyright = "2021, Stefano Camborda, Jan-Niklas Weder"
author = "Stefano Camborda, Jan-Niklas Weder, Nadine Töpfer"


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
    "sphinxcontrib.bibtex",
]

bibtex_bibfiles = ['biblio.bib']
bibtex_default_style = 'mystyle'
bibtex_reference_style = 'label'

html_css_files = [
    'custom.css',
]

import pybtex.plugin
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.style.template import node, FieldIsMissing, join, sentence  # ... and anything else needed

# Adapted from pybtex to use the customized functions
class MyStyle(UnsrtStyle):
    def format_names(self, role, as_sentence=True):
        formatted_names = names(role, sep=' and ', sep2 = ' and ', last_sep=' and ')
        if as_sentence:
            return sentence [formatted_names]
        else:
            return formatted_names


# Adapted from pybtex to get "," as separator instead of dots
@node
def sentence(children, data, capfirst=False, capitalize=False, add_period=True, sep=', '):

    text = join(sep)[children].format_data(data)
    if capfirst:
        text = text.capfirst()
    if capitalize:
        text = text.capitalize()
    if add_period:
        text = text + ","
    return text

# Adapted from pybtex to get a maximum of 3 authors or et al in the references.
@node
def names(children, context, role, **kwargs):
    """Return formatted names."""

    assert not children

    try:
        persons = context['entry'].persons[role]
    except KeyError:
        raise FieldIsMissing(role, context['entry'])

    style = context['style']
    print(style)

    if len(persons) <= 3:
        formatted_names = [style.format_name(person, style.abbreviate_names) for person in persons]
        print(persons)
    else:
        formatted_names = [style.format_name(persons[0], style.abbreviate_names)," et al."]
        try:
            kwargs["sep"] = " "
            kwargs["sep2"] = " "
        except KeyError:
            pass

    return join(**kwargs) [formatted_names].format_data(context)


pybtex.plugin.register_plugin('pybtex.style.formatting', 'mystyle', MyStyle)

# Extensions configuration
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "cobra": ("https://cobrapy.readthedocs.io/en/latest/", None),
    "escher": ("https://escher.readthedocs.io/en/latest/", None),
}
autoapi_dirs = ["../../src/cobramod"]
autoapi_root = "module"
autoapi_generate_api_docs = True
autoapi_add_toctree_entry = False
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "imported-members",
]
napoleon_custom_sections = [
    ("Arguments for reactions", "params_style"),
    ("Special arguments for databases", "params_style"),
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

html_theme = "sphinx_rtd_theme"
html_title = (
    "CobraMod: A pathway-centric curation tool for contraint-based "
    + "metabolic models"
)
html_short_titel = "CobraMod"
html_static_path = ["_static"]
html_logo = "logo_2.png"
html_theme_options = {"logo_only": True}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = []

# Configuration for nbsphinx
nbsphinx_execute = "never"
# This statement should avoid the creation of duplicates widgets
nbsphinx_widgets_path = ""
html_js_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js",
    DEFAULT_EMBED_REQUIREJS_URL,
]
