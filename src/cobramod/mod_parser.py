#!/usr/bin/env python3
from contextlib import suppress
from pathlib import Path

from cobramod.parsing.base import BaseParser
from cobramod.parsing.biocyc import BiocycParser
from cobramod.parsing.kegg import KeggParser
from cobramod.parsing.bigg import BiggParser

# FIXME: find a workaround to avoid linter problems
BiocycParser
KeggParser
BiggParser


def get_data(
    directory: Path,
    identifier: str,
    database: str,
    debug_level: int = 20,
    **kwargs
) -> dict:
    """
    Retrieves and tranform the data into a dictionary for given identifier
    in a specific database.

    Args:
        directory (Path): Directory to store and retrieve local data.
        identifier (str): original identifier
        database (str): Name of database. Options: "META", "ARA", "KEGG",
            "BIGG".
        debug_level (int, optional): Level of debugging. Read package logging
            for more info. Defaults to 20

    Keyword Arguments:
        model_id: Exclusive for BIGG. Original identifier of model to search.
            Some examples: "e_coli_core", "universal"
    Returns:
        dict: relevant data for given identifier
    """
    for parser in BaseParser.__subclasses__():
        with suppress(Warning):
            # This method is called to create Warning for the rest of
            # databases if database name is not the same
            parser._return_database(database=database)
            real_parser = parser
    return real_parser._retrieve_data(
        directory=directory,
        identifier=identifier,
        database=database,
        debug_level=debug_level,
        **kwargs
    )
