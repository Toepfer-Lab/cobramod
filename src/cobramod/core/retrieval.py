#!/usr/bin/env python3
"""Data retrieval

This module creates the function :func:`cobramod.core.retrieval.get_data`,
which gets the data from a local directory or from different databases.
The available databases for data retrieval can be found in the variable
:obj:`cobramod.core.retrieval.available_databases`
"""
from contextlib import suppress
from pathlib import Path
from typing import Type, Optional

from requests import HTTPError

from cobramod.error import WrongParserError, PatternNotFound
from cobramod.parsing.base import BaseParser
from cobramod.parsing.biocyc import BiocycParser
from cobramod.parsing.kegg import KeggParser
from cobramod.parsing.bigg import BiggParser
from cobramod.parsing.plantcyc import PlantCycParser
from cobramod.utils import _path_match, get_key_dict

parsers = [BiocycParser, PlantCycParser, KeggParser, BiggParser]


class Databases(object):
    """
    Simple object that shows the information about the database
    """

    def __init__(self):
        self.msg = (
            "CobraMod supports Biocyc, the Plant Metabolic Network (PMN), KEGG"
            " and BiGG Models repository. Biocyc includes around 18.000 "
            "sub-databases and a complete list for BioCyc can be found at "
            "'https://biocyc.org/biocyc-pgdb-list.shtml'. "
            "The database-specific identifiers can be found in the URL of the "
            'respective data. For instance, for "diphosphate" this is:'
        )

    def __repr__(self):
        """
        Returns a string about the definition of "identifier" and where
        to obtain them.
        """
        self.msg2 = """\n
BioCyc, sub-database ECOLI -> https://biocyc.org/compound?orgid=ECOLI&id=PPI ->
PPI\nPlant Metabolic Network, sub-database CORN ->
https://pmn.plantcyc.org/compound?orgid=CORN&id=PPI -> PPI\nKEGG ->
https://www.genome.jp/entry/C00013 -> C00013\nBiGG Models Repository,
universal model -> http://bigg.ucsd.edu/universal/metabolites/ppi -> ppi\n
CobraMod uses abbreviations to represent the databases or sub-databases:\n
Database -> Abbreviation\nBioCyc -> META or identifier of sub-database e.g:
ECOLI, ARA, GCF_000010885\nPlant Metabolic Network -> pmn:PLANT or identifier
of sub-database e.g pmn:ARA, pmn:CORN\nKEGG -> KEGG\nBiGG Models Repository ->
BIGG
        """
        return self.msg + self.msg2

    def _repr_html_(self):
        """
        Returns a HTML string about the definition of "identifier" and where
        to obtain them.
        """
        return f"""<p>{self.msg}</p>
<table style="width: 100%; border-collapse: collapse; float: left;" border="1">
<tbody> <tr> <td style="width: 50%;"> <h3>Database</h3> </td> <td style="width:
50%;"> <h3>URL with identifier (bold)</h3> </td> </tr> <tr> <td style="width:
50%;">BioCyc, sub-database ECOLI</td> <td style="width:
50%;">https://biocyc.org/compound?orgid=ECOLI&amp;id=<strong>PPI</strong></td>
</tr> <tr> <td style="width: 50%;">Plant Metabolic Network, sub-database
CORN</td> <td style="width:
50%;">https://pmn.plantcyc.org/compound?orgid=CORN&amp;id=<strong>PPI</strong></td>
</tr> <tr> <td style="width: 50%;">KEGG</td> <td style="width:
50%;">https://www.genome.jp/entry/<strong>C00013</strong></td> </tr> <tr> <td
style="width: 50%;">BiGG Models Repository, universal model</td> <td
style="width:
50%;">http://bigg.ucsd.edu/universal/metabolites/<strong>ppi</strong></td>
</tr> </tbody> </table> <p>CobraMod uses abbreviations to represent the
databases or sub-databases:</p> <table style="width: 100%; border-collapse:
collapse; float: left;" border="1"> <tbody> <tr> <td style="width: 50%;">
<h3>Database</h3> </td> <td style="width: 50%;"> <h3>Abbreviation</h3> </td>
</tr> <tr> <td style="width: 50%;">BioCyc</td> <td style="width: 50%;">META or
identifier of sub-database e.g: ECOLI, ARA, GCF_000010885</td> </tr> <tr> <td
style="width: 50%;">Plant Metabolic Network</td> <td style="width: 50%;">Prefix
"pmn:" with the sub-database identifier, e.g pmn:PLANT, pmn:ARA, pmn:CORN</td>
</tr> <tr> <td style="width: 50%;">KEGG</td> <td style="width: 50%;">KEGG</td>
</tr> <tr> <td style="width: 50%;">BiGG Models Repository</td> <td
style="width: 50%;">BIGG</td> </tr> </tbody> </table>
"""


available_databases = Databases()


def _get_parser(directory: Path, database: str) -> Type[BaseParser]:
    """
    Return the child instance of BaseParser that matches the name of the
    database.
    """
    for parser in BaseParser.__subclasses__():
        with suppress(WrongParserError):
            # This method will raise a WrongParserError. Skipping it will
            # return the the real parser
            parser._check_database(directory=directory, database=database)
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
    Retrieves and transforms the data into a dictionary for the given
    identifier from a specific database.

    Args:
        directory (Path): Directory to store and retrieve local data.
        identifier (str): Original identifier.
        database (str): Name of database. Check
            :obj:`cobramod.available_databases` for a list of names.
        debug_level (int, optional): Level of debugging. Read package logging
            for more info. Defaults to 20.

    Keyword Arguments:
        model_id: Exclusive for BIGG. Original identifier of the model to be
            searched. Some examples: "e_coli_core", "universal".
        genome (str, optional): Exclusive for KEGG. Abbreviation for the
            species involved. Genes will be obtained for this species.
    Returns:
        dict: Relevant data for the given identifier.
    """
    real_parser = _get_parser(directory=directory, database=database)
    return real_parser._retrieve_data(
        directory=directory,
        identifier=identifier,
        database=database,
        debug_level=debug_level,
        **kwargs,
    )


def _get_correct_data(
    replacement: dict,
    directory: Path,
    database: str,
    identifier: str,
    model_id: str,
    genome: Optional[str],
):
    # TODO: docstrings
    try:
        replacement[identifier]
        data_dict = get_data(
            directory=directory,
            database=database,
            identifier=replacement[identifier],
            model_id=model_id,
            genome=genome,
        )
    except KeyError:
        data_dict = get_data(
            directory=directory,
            database=database,
            identifier=identifier,
            model_id=model_id,
            genome=genome,
        )
    except HTTPError:
        data_dict = get_data(
            directory=directory,
            database=database,
            identifier=identifier,
            model_id=model_id,
            genome=genome,
        )
        data_dict["ENTRY"] = replacement[identifier]
    return data_dict


def _retrieve_dict(directory: Path, target: str) -> dict:
    """
    Search and return in given directory, specific target and return a
    dictionary with the parsed information.
    Args:
        directory (Path): Path to search. This includes subdirectories
        target (str): Pattern to search.

    Raises:
        FileNotFoundError: If target cannot be found.
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
    Return the identifier of crossref for the given target. It can be a
    metabolite or a Reaction.

    Args:
        directory (Path): Path of stored data.
        target (str): Identifier to search for.
        database (str): Pattern for the name of the cross-reference,
            e.g CAS, BIGG.

    Returns
        str: Corresponding identifier for cross-reference.

    Raises:
        PatternNotFound: If the target cannot be properly identified.
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
