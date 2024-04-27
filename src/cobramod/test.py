"""Test variables

This module creates the textbook models for testing. These are based on the
original model "e_coli_core" from BiGG database, which are also included in
COBRApy.
"""
from importlib import resources

import cobra
import cobramod
from cobra.io import read_sbml_model

textbook = read_sbml_model(
    resources.read_text(cobra, "data/textbook.xml.gz")
)

textbook_biocyc = read_sbml_model(
    resources.read_text(cobramod, "data/textbook_biocyc.sbml")
)

textbook_biocyc_groups = read_sbml_model(
    resources.read_text(cobramod, "data/textbook_biocyc_groups.sbml")
)

textbook_kegg = read_sbml_model(
    resources.read_text(cobramod, "data/textbook_kegg.sbml")
)
