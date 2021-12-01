#!/usr/bin/env python3
"""Data parsing for BiGG

This module handles the retrieval of data from BiGG into a local directory.
The possible type of data that can be download:

- Metabolites: Normally, simple names.
- Reactions: Mostly abbreviations.
- Genes: It is included in the Reactions. Names and gene-reaction-rule is also
acquired

They change identifiers depending on the model given. BiGG have multiple
models.
Contact maintainers if other types should be added.

Important class of the module:
- BiggParser: Child of the abstract class
:class:`cobramod.parsing.base.BaseParser`.
"""
from contextlib import suppress
from json import loads, JSONDecodeError
from pathlib import Path
from typing import Any, Tuple
from warnings import warn

from requests import get, HTTPError, Response

from cobramod.debug import debug_log
from cobramod.error import WrongParserError
from cobramod.parsing.base import BaseParser


def _find_url(model_id: str, identifier: str) -> Tuple[Response, str]:
    """
    Tries to find a valid URL for the API of BIGG. It will return a
    :class:`requests.Response` if URL is valid.
    object

    Args:
        model_id (str): Name for the specific BIGG model identifier.
        identifier (str): Identifier for the item. Can be a reaction or
            compound

    Returns:
        (Response,str): A valid request Response followed by the database
            Version.

    Raises:
        HTTPError: If identifier is not found in BIGG database.
    """
    for object_type in ("reactions", "metabolites"):
        with suppress(HTTPError):
            # Check that status is available
            debug_log.debug(
                f'{object_type.capitalize()[:-1]} "{identifier}" not found in '
                f'directory "BIGG", subdirectory "{model_id}".'
            )
            # Retrieve from URL
            url_text = (
                f"http://bigg.ucsd.edu/api/v2/models/{model_id}/{object_type}"
                f"/{identifier}"
            )
            debug_log.debug(f"Searching {url_text} for biochemical data.")
            # Get and check for errors
            response = get(url_text)
            response.raise_for_status()

            info_response = get("http://bigg.ucsd.edu/api/v2/database_version")
            info_response.raise_for_status()
            db_version = info_response.json()["bigg_models_version"]

            return response, db_version
    # Otherwise
    raise HTTPError(f"Identifier '{identifier}' not found in BIGG.")


def retrieve_data(directory: Path, identifier: str, model_id: str) -> dict:
    """
    Searchs in given parent directory if data is located in their respective
    model_id directory. Else, data will be retrieved from the corresponding
    model identifier. Returns dictionary from JSON.

    Args:
        directory (Path): Path to directory where data is located.
        identifier (str): Identifier for given database.
        model_id (str): Identifier of model for BiGG. Examples: *e_coli_core*,
            *universal*. Check http://bigg.ucsd.edu/models for a complete
            list of models.

    Raises:
        HTTPError: If object is not available in given database
        NotADirectoryError: If parent directory is not found

    Returns:
        dict: directory transformed from a JSON.
    """
    if directory.exists():
        # Create subdirectory data
        data_dir = directory.joinpath("BIGG").joinpath(model_id)
        filename = data_dir.joinpath(f"{identifier}.json")
        debug_log.debug(
            f'Searching "{identifier}" in directory BIGG, '
            f'sub-directory "{model_id}".'
        )
        try:
            return BiggParser._read_file(filename=filename)
        except FileNotFoundError:
            # It will raise an error by itself if needed
            response, db_version = _find_url(
                model_id=model_id, identifier=identifier
            )

            BaseParser.check_database_version(directory, "bigg", db_version)

            debug_log.info(
                f'Object "{identifier}" found. Saving in '
                f'directory "BIGG", subdirectory "{model_id}"'
            )
            # r.text is the raw text
            data_dir.mkdir(parents=True, exist_ok=True)
            with open(file=filename, mode="w+") as file:
                file.write(response.text)
            return loads(s=response.text)
    else:
        msg = "Directory not found. Please create the given directory."
        debug_log.critical(msg)
        raise NotADirectoryError(msg)


def _build_reference(json_data: dict) -> dict:
    """
    Return a dictionary of cross-references, where the keys are the
    cross-references and values the their identifiers.
    """
    references = json_data["database_links"]
    try:
        return {item: references[item][0]["id"] for item in references}
    except KeyError:
        raise KeyError(
            "Problem with given json dictionary. Inform maintainers"
        )


def _p_compound(json_data: dict) -> dict:
    """
    Returns a dictionary with the most important attributes of the compound.
    """
    # For generic compounds, a list comes, otherwise it is a regular string
    try:
        formula = json_data["formulae"][0]
        charge = json_data["charges"][0]
    except KeyError:
        formula = json_data["formula"]
        charge = json_data["charge"]
    return {
        "TYPE": "Compound",
        "ENTRY": json_data["bigg_id"],
        "NAME": json_data["name"],
        "FORMULA": formula,
        "CHARGE": charge,
        "DATABASE": "BIGG",
        "XREF": _build_reference(json_data=json_data),
    }


def _get_metabolites(json_data: dict) -> dict:
    """
    Returns a dictionary with participant-metabolites
    """
    meta_dict = dict()
    for item in json_data["metabolites"]:
        # left
        if item["stoichiometry"] < 0:
            meta_dict[f'l_{item["bigg_id"]}'] = item["stoichiometry"]
        # right
        elif item["stoichiometry"] > 0:
            meta_dict[f'r_{item["bigg_id"]}'] = item["stoichiometry"]
        else:
            raise Warning("Error in the stoichiometry. Verify original JSON")
    return meta_dict


def _p_genes(json_data: dict) -> dict:
    """
    Return a dictionary that includes the genes and gene-reaction-rule for
    given reaction. The dictionary include the key "genes", which is dictionary
    with the identifier and name of the gene; and the key "rule" for the
    COBRApy representation of the gene-reaction rule.
    """
    # results comes in a single list
    genes = dict()
    rule = str()
    # Try to obtain and create a dictionary with the gene identifier and their
    # corresponding names
    with suppress(KeyError):
        json_genes = json_data["results"][0]["genes"]
        for single in json_genes:
            genes[single["bigg_id"]] = single["name"]
    with suppress(KeyError):
        rule = json_data["results"][0]["gene_reaction_rule"]
    return {"genes": genes, "rule": rule}


def _p_reaction(json_data: dict) -> dict:
    """
    Returns a dictionary with the most important attributes of the reaction.
    """
    temp_dict = {
        "TYPE": "Reaction",
        "ENTRY": json_data["bigg_id"],
        "NAME": json_data["name"],
        "EQUATION": _get_metabolites(json_data=json_data),
        "BOUNDS": (-1000, 1000),
        "TRANSPORT": BaseParser._check_transport(
            data_dict=_get_metabolites(json_data=json_data)
        ),
        "DATABASE": "BIGG",
        "XREF": _build_reference(json_data=json_data),
        "GENES": _p_genes(json_data=json_data),
    }
    msg = (
        f'Reaction "{json_data["bigg_id"]}" from BIGG set to reversible. '
        "Please modify if necessary."
    )
    debug_log.warning(msg=msg)
    warn(message=msg, category=UserWarning)
    return temp_dict


class BiggParser(BaseParser):
    @staticmethod
    def _retrieve_data(
        directory: Path,
        identifier: str,
        database: str,
        debug_level: int,
        **kwargs,
    ) -> dict:
        """
        Retrieves data from given model in BiGG and parses the most important
        attributes into a dictionary.

        Args:
            directory (Path): Directory to store and retrieve local data.
            identifier (str): original identifier
            debug_level(int): Level of debugging. Read package logging
            for more info.

        Keyword Arguments:
            model_id (str): Identifier of the model. Some examples:
                "e_coli_core", "universal"

        Returns:
            dict: relevant data for given identifier
        """
        BiggParser._check_database(directory=directory, database=database)
        try:
            model_id = kwargs["model_id"]
        except KeyError:
            raise KeyError(
                'Argument "model_id" is missing. Please specify it in one of'
                + " the main functions."
            )
        # It will raise an error for HTTPError
        json_data = retrieve_data(
            directory=directory, identifier=identifier, model_id=model_id
        )
        debug_log.log(
            level=debug_level,
            msg=f'Data for "{identifier}" was retrieved ' + "from BIGG.",
        )
        return BiggParser._parse(root=json_data)

    @staticmethod
    def _parse(root: Any) -> dict:  # type: ignore
        """
        Parses the JSON dictionary and returns a dictionary with the most
        important attributes depending of the type of the identifier.
        """
        try:
            for method in (_p_reaction, _p_compound):
                with suppress(KeyError):
                    bigg_dict = method(json_data=root)
            return bigg_dict
        except TypeError:
            raise WrongParserError

    @staticmethod
    def _check_database(directory: Path, database: str):
        """
        Returns name of the database. It will raise a Error if name is
        incorrect.
        """
        if database != "BIGG":
            raise WrongParserError

    @staticmethod
    def _read_file(filename: Path) -> dict:
        """
        Reads the given file a returns a JSON dictionary with most important
        information from it.
        """
        try:
            with open(file=filename, mode="r") as file:
                unformatted_data = file.read()
            return loads(s=unformatted_data)
        except JSONDecodeError:
            raise WrongParserError("Wrong filetype. Please use a JSON file.")
