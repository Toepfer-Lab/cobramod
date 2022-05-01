#!/usr/bin/env python3
"""Test variables

This module creates the textbook models for testing. These are based on the
original model "e_coli_core" from BiGG database, which are also included in
COBRApy.
"""
import gzip

from importlib_resources import open_binary, files, as_file
from pathlib import Path

from cobra.io import read_sbml_model, validate_sbml_model
import cobra.data

# This is the directory where the textbook models are located.
data_dir = Path(__file__).resolve().parent.joinpath("data")
# These are the three posibilities

try:
    textbook
except NameError:
    textbook_raw = files(cobra.data).joinpath("textbook.xml")
    with as_file(textbook_raw) as textbookXML:
        textbook = read_sbml_model(str(textbookXML))


textbook_biocyc = read_sbml_model(
    filename=str(data_dir.joinpath("textbook_biocyc.sbml"))
)
textbook_kegg = read_sbml_model(
    filename=str(data_dir.joinpath("textbook_kegg.sbml"))
)
