"""Test variables

This module creates the textbook models for testing. These are based on the
original model "e_coli_core" from BiGG database, which are also included in
COBRApy.
"""
import pkg_resources
from cobra.io import read_sbml_model

textbook = read_sbml_model(
    pkg_resources.resource_filename("cobra", "data/textbook.xml.gz")
)

textbook_biocyc = read_sbml_model(
    pkg_resources.resource_filename("cobramod", "data/textbook_biocyc.sbml")
)
textbook_kegg = read_sbml_model(
    pkg_resources.resource_filename("cobramod", "data/textbook_kegg.sbml")
)
