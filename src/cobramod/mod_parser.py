#!/usr/bin/env python3
from contextlib import suppress
from pathlib import Path
from typing import Type

from cobramod.error import WrongParserError, PatternNotFound
from cobramod.parsing.base import BaseParser
from cobramod.parsing.biocyc import BiocycParser
from cobramod.parsing.kegg import KeggParser
from cobramod.parsing.bigg import BiggParser
from cobramod.utils import _path_match, get_key_dict

parsers = [BiocycParser, KeggParser, BiggParser]
available_databases = ["META", "ARA", "META", "BIGG"]


def _get_parser(database: str) -> Type[BaseParser]:
    """
    Return the child instance of BaseParser that matches the name of the
    database.
    """
    for parser in BaseParser.__subclasses__():
        with suppress(WrongParserError):
            # This method will raise a WrongParserError. Skipping it will
            # return the the real parser
            parser._return_database(database=database)
            real_parser = parser
            break
    try:
        return real_parser
    except UnboundLocalError:
        raise WrongParserError("No parser found for that database")


def get_data(
    directory: Path,
    identifier: str,
    database: str,
    debug_level: int = 20,
    **kwargs,
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
    real_parser = _get_parser(database=database)
    return real_parser._retrieve_data(
        directory=directory,
        identifier=identifier,
        database=database,
        debug_level=debug_level,
        **kwargs,
    )


def _retrieve_dict(directory: Path, target: str) -> dict:
    """
    Search and return in given directory, specific target and return a
    dictionary with the parsed infomation.
    Args:
        directory (Path): Path to search. This includes subdirectories
        target (str): Pattern to search.

    Raises:
        FileNotFoundError: if target cannot be found
    """

    try:
        filename = _path_match(directory=directory, pattern=target)
    except StopIteration:
        raise FileNotFoundError(
            f"No file was found with the sub-string {target}"
        )
    for parser in BaseParser.__subclasses__():
        with suppress(WrongParserError, NotImplementedError):
            data_dict = parser._parse(
                root=parser._read_file(filename=filename)
            )["XREF"]
    try:
        return data_dict
    except UnboundLocalError:
        raise WrongParserError(
            "No parser could be identified. Please contact maintainers"
        )


def translate(directory: Path, target: str, database: str) -> str:
    """
    Return the identifier of crossref for given target. It can be a metabolite
    or a Reaction.

    Args:
        directory (Path): Path of stored data.
        target (str): Identifier to search for.
        database (str): Pattern for name of the cross-reference, e.g CAS, BIGG.

    Returns
        str: corresponding identifier for cross-reference.

    Raises:
        FileNotFoundError: If no target can be found
        WrongParserError: If target cannot be properly identified
    """
    data_dict = _retrieve_dict(directory=directory, target=target)
    try:
        key = get_key_dict(dictionary=data_dict, pattern=database)
        return data_dict[key]
    except PatternNotFound:
        raise PatternNotFound(
            "No parser could be identified. Please contact maintainers"
        )
