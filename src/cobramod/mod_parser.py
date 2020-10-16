#!/usr/bin/env python3
from pathlib import Path
from contextlib import suppress
from cobramod.parsing.base import BaseParser
from cobramod.parsing.biocyc import BiocycParser
from cobramod.parsing.kegg import KeggParser

# FIXME: find a workaround to avoid linter problems
BiocycParser
KeggParser


def get_data(
    directory: Path, identifier: str, database: str, debug_level: int = 20
) -> dict:
    """
    Retrieves and tranform the data into a dictionary for given identifier
    in a specific database.

    Args:
        directory (Path): Directory to store and retrieve local data.
        identifier (str): original identifier
        database (str): Name of database. Options: "META", "ARA", "KEGG"
        debug_level (int, optional): Level of debugging. Read package logging
            for more info. Defaults to 20

    Returns:
        dict: relevant data for given identifier
    """
    # for parser in [BiocycParseri]:
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
    )
