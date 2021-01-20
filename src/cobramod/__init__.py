#!/usr/bin/env python3
"""CobraMod Package

It is a helpbox and extension package based on COBRApy . This package
facilitates the retrieval of biochemical data from multiple databases and help
users extend metabolic networks.

CobraMod ensures an easy-friendly API, customization and a record of all the
changes made to metabolic models.

Important functions from CobraMod that can be imported directly:

- get_data: Store and retrieve data from given database.
- create_object: Create a corresponding COBRApy object from given identifier.
- add_pathway: Extend given pathway into given model.
- test_result: Test given reaction to check for a feasible solution.
- translate: Check for cross-references.

The exclusive class :func:`cobramod.pathway.Pathway` is an extension of
:func:`cobra.core.group.Group`. CobraMod uses :func:`Escher` to visualize
the pathway and fluxes of that specific class. Some new methods added:

- visualize: Return an Builder for the representation of the pathway.
- solution: Filters solution and returns fluxes of only members of the class.

By default, CobraMod creates the file "debug.log", which displays the changes
occured when running a script. In order to modify the logging, import
"debug_log" from "cobramod.debug". The default logging level is defined as
INFO. Read the documentation of logging for more information.
"""
from cobramod.mod_parser import get_data, translate
from cobramod.creation import (
    create_object,
    add_meta_from_file,
    add_reactions_from_file,
    meta_string_to_model,
)
from cobramod.extension import add_pathway, test_result
from cobramod.pathway import Pathway
from cobramod.utils import model_convert

__all__ = [
    "get_data",
    "create_object",
    "add_meta_from_file",
    "add_reactions_from_file",
    "meta_string_to_model",
    "add_pathway",
    "test_result",
    "Pathway",
    "translate",
    "model_convert",
]

__version__ = "0.1.2"
