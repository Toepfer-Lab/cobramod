"""Test variables

This module creates the textbook models for testing. These are based on the
original model "e_coli_core" from BiGG database, which are also included in
COBRApy.
"""

from importlib import resources

import cobra
from cobra.io import read_sbml_model

from cobramod import data

with resources.as_file(
    resources.files(cobra).joinpath("data").joinpath("textbook.xml.gz")
) as textbook:
    textbook = read_sbml_model(textbook)

textbook_biocyc = read_sbml_model(
    resources.read_text(data, "textbook_biocyc.sbml")
)

textbook_biocyc_groups = read_sbml_model(
    resources.read_text(data, "textbook_biocyc_groups.sbml")
)

textbook_kegg = read_sbml_model(resources.read_text(data, "textbook_kegg.sbml"))
