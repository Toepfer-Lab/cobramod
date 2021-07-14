#!/usr/bin/env python3
"""Data retrieval

This module creates the function :func:`cobramod.core.retrieval.get_data`,
which gets the data from a local directory or from different databases.
The available databases for data retrieval can be found in the variable
:obj:`cobramod.core.retrieval.available_databases`
"""
from contextlib import suppress
from collections import UserList
from pathlib import Path
from typing import Type

from cobramod.error import WrongParserError, PatternNotFound
from cobramod.parsing.base import BaseParser
from cobramod.parsing.biocyc import BiocycParser
from cobramod.parsing.kegg import KeggParser
from cobramod.parsing.bigg import BiggParser
from cobramod.parsing.plantcyc import PlantCycParser
from cobramod.utils import _path_match, get_key_dict

parsers = [BiocycParser, PlantCycParser, KeggParser, BiggParser]


class ListCobramod(UserList):
    """
    Simple list that prints out a message about the enormous size of Biocyc.
    """

    def __init__(self, initlist=[]):
        super().__init__(initlist=initlist)
        self.msg = (
            "Biocyc includes around 18.000 sub-databases. The complete list "
            + "can be found in 'https://biocyc.org/biocyc-pgdb-list.shtml'. "
            + "Please use the corresponding object identifier. e.g: 'ARA', "
            + "'GCF_000963925'"
        )

    def __repr__(self):
        print(self.msg)
        return super().__repr__()

    def __str__(self):
        return self.msg + "\n" + super().__str__()


available_databases = ListCobramod(["META", "PLANT", "KEGG", "BIGG"])


def _get_parser(database: str) -> Type[BaseParser]:
    """
    Return the child instance of BaseParser that matches the name of the
    database.
    """
    for parser in BaseParser.__subclasses__():
        with suppress(WrongParserError):
            # This method will raise a WrongParserError. Skipping it will
            # return the the real parser
            parser._check_database(database=database)
            real_parser = parser
            break
    try:
        return real_parser
    except UnboundLocalError:
        raise WrongParserError(
            "No parser found for that database. Please check "
            + "cobramod.retrieval.available_databases"
        )


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
        database (str): Name of database. Check
            :obj:`cobramod.available_databases` for a list of names.
        debug_level (int, optional): Level of debugging. Read package logging
            for more info. Defaults to 20

    Keyword Arguments:
        model_id: Exclusive for BIGG. Original identifier of model to search.
            Some examples: "e_coli_core", "universal"
        genome (str, optional): Exclusive for KEGG. Abbreviation for the
            specie involved. Genes will be obtained from this specie.
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
    dictionary with the parsed information.
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
            try:
                data_dict = parser._parse(
                    root=parser._read_file(filename=filename),
                    directory=directory,
                )["XREF"]
            except TypeError:
                data_dict = parser._parse(  # type: ignore
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
        PatternNotFound: If target cannot be properly identified
    """
    # Return parsed information
    data_dict = _retrieve_dict(directory=directory, target=target)
    try:
        # Search for the name of the database as a pattern
        key = get_key_dict(dictionary=data_dict, pattern=database)
        return data_dict[key]
    except PatternNotFound:
        raise PatternNotFound(
            "No could be identified. Probably the target does not include the "
            + "given database"
        )
