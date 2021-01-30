#!/usr/bin/env python3
"""Test variables

This module creates the textbook models for testing. These are based on the
original model "e_coli_core" from BiGG database, which are also included in
COBRApy.
"""
from pathlib import Path

from cobra.io import read_sbml_model
from cobra.test import create_test_model

# This is the directory where the textbook models are located.
data_dir = Path(__file__).resolve().parent.joinpath("data")
# These are the three posibilities
textbook = create_test_model(model_name="textbook")
textbook_biocyc = read_sbml_model(
    filename=str(data_dir.joinpath("textbook_biocyc.sbml"))
)
textbook_kegg = read_sbml_model(
    filename=str(data_dir.joinpath("textbook_kegg.sbml"))
)
