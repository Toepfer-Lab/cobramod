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
- add_metabolites: Add metabolites from different objects
- add_reactions: Add reactions from different objects
- add_pathway: Extend given pathway into given model.
- test_result: Test given reaction to check for a feasible solution.
- translate: Check for cross-references.

The exclusive class :func:`cobramod.pathway.Pathway` is an extension of
:class:`cobra.core.group.Group`. CobraMod uses :mod:`Escher` to visualize
the pathway and fluxes of that specific class. Some new methods added:

- visualize: Return an Builder for the representation of the pathway.
- solution: Filters solution and returns fluxes of only members of the class.

By default, CobraMod creates the file "debug.log", which displays the changes
occured when running a script. In order to modify the logging, import
"debug_log" from "cobramod.debug". The default logging level is defined as
INFO. Read the documentation of logging for more information.

For a list of databases, load variable :obj:`cobramod.available_databases`
"""
from cobramod.core.creation import (
    create_object,
    add_metabolites,
    add_reactions,
)
from cobramod.core.extension import add_pathway, test_result
from cobramod.core.pathway import Pathway, model_convert
from cobramod.core.retrieval import get_data, translate, available_databases

__all__ = [
    "get_data",
    "create_object",
    "add_reactions",
    "add_metabolites",
    "add_pathway",
    "test_result",
    "Pathway",
    "translate",
    "model_convert",
    "available_databases",
]

__version__ = "0.3.0"
