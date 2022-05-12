#!/usr/bin/env python3
"""Test variables

This module creates the textbook models for testing. These are based on the
original model "e_coli_core" from BiGG database, which are also included in
COBRApy.
"""

from pathlib import Path

import cobra.data
from cobra import Model
from cobra.io import read_sbml_model
from importlib_resources import files, as_file

# This is the directory where the textbook models are located.
data_dir = Path(__file__).resolve().parent.joinpath("data")
# These are the three posibilities

textbook: Model

try:
    textbook
except NameError:
    textbook_raw = files(cobra.data).joinpath("textbook.xml.gz")
    with as_file(textbook_raw) as textbookXML:
        textbook = read_sbml_model(str(textbookXML))


textbook_biocyc = read_sbml_model(
    filename=str(data_dir.joinpath("textbook_biocyc.sbml"))
)
textbook_kegg = read_sbml_model(
    filename=str(data_dir.joinpath("textbook_kegg.sbml"))
)
