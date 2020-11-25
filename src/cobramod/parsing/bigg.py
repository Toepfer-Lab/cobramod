#!/usr/bin/env python3
from contextlib import suppress
from json import loads, JSONDecodeError
from pathlib import Path
from typing import Any

from requests import get, HTTPError

from cobramod.debug import debug_log
from cobramod.error import WrongParserError
from cobramod.utils import get_key_dict
from cobramod.parsing.base import BaseParser


def _get_json_bigg(
    directory: Path, identifier: str, model_id: str, object_type: str, **kwargs
) -> dict:
    """
    Searchs in given parent directory if data is located in their respective
    database directory. If not, data will be retrieved from the corresponding
    database. Returns dictionary from JSON.

    Args:
        directory (Path): Path to directory where data is located.
        identifier (str): Identifier for given database.
        model_id (str): Identifier of model for BiGG. Examples: *e_coli_core*,
            *universal*
        object_type (str): Either 'reaction' or 'metabolite'.

    Raises:
        Warning: If object is not available in given database
        NotADirectoryError: If parent directory is not found

    Returns:
        dict: directory transformed from a JSON.
    """
    if directory.exists():
        database = "BIGG"
        data_dir = directory.joinpath(database).joinpath(model_id)
        data_dir.mkdir(parents=True, exist_ok=True)
        filename = data_dir.joinpath(f"{identifier}.json")
        debug_log.debug(
            f'Searching {object_type} "{identifier}" in directory "BIGG"'
        )
        try:
            debug_log.debug(f"Identifier '{identifier}' found.")
            return BiggParser._read_file(filename=filename)
        except FileNotFoundError:
            debug_log.debug(
                f'{object_type.capitalize()} "{identifier}" not found in '
                f'directory "BIGG", subdirectory "{model_id}".'
            )
            # Retrieve from URL
            url_text = (
                f"http://bigg.ucsd.edu/api/v2/models/{model_id}/{object_type}s"
                f"/{identifier}"
            )
            debug_log.debug(f"Searching in {url_text}")
            r = get(url_text)
            if r.status_code >= 400:
                msg = (
                    f'"{identifier}" could not be retrieved. Status code: '
                    f"{r.status_code}"
                )
                debug_log.error(msg)
                raise HTTPError(msg)
            else:
                unformatted_data = r.text
                debug_log.info(
                    f'Object "{identifier}" found. Saving in '
                    f'directory "BIGG", subdirectory "{model_id}"'
                )
                with open(file=filename, mode="w+") as f:
                    f.write(unformatted_data)
                return loads(s=unformatted_data)
    else:
        msg = "Directory not found"
        debug_log.critical(msg)
        raise NotADirectoryError(msg)


def _build_reference(json_data: dict) -> dict:
    """
    From the json dictionary, return a dictionary with keys as cross-references
    and values as their identifiers.
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
    Parse the information of the JSON dictionary and returns a dictionary with
    the most important attributes of the compound.
    """
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


def _p_metabolites(json_data: dict) -> dict:
    """
    Checks the JSON data and returns a dictionary with participant-metabolites
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


def _p_reaction(json_data: dict) -> dict:
    """
    Parses the data of the JSON dictionary a returns a dictionary with the most
    important attributes of the reaction.
    """
    temp_dict = {
        "TYPE": "Reaction",
        "ENTRY": json_data["bigg_id"],
        "NAME": json_data["name"],
        "EQUATION": _p_metabolites(json_data=json_data),
        # FIXME: assuming bounds. Create proper method
        "BOUNDS": (-1000, 1000),
        "TRANSPORT": BaseParser._check_transport(
            data_dict=_p_metabolites(json_data=json_data)
        ),
        "DATABASE": "BIGG",
        "XREF": _build_reference(json_data=json_data),
    }
    return temp_dict


def _get_db_names(data_dict: dict, pattern: str) -> str:
    """
    From the original data of a BiGG JSON, retrieve the corresponding
    identifier for its crossref database.
    """
    try:
        database = get_key_dict(
            dictionary=data_dict["database_links"], pattern=pattern
        )
        return data_dict["database_links"][database][0]["id"]
    except Warning:
        msg = "Given database is not found in the crossref section"
        debug_log.error(msg)
        raise KeyError(msg)


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
            "e_coli_core", "universal"
            debug_level(int): Level of debugging. Read package logging
            for more info.

        Keyword Arguments:
            model_id (str): Identifier of the model. Some examples:

        Returns:
            dict: relevant data for given identifier
        """
        try:
            json_data = _get_json_bigg(
                directory=directory,
                database=database,
                identifier=identifier,
                object_type="metabolite",
                **kwargs,
            )
            msg_extra = "Identified as a metabolite"
        except KeyError:
            json_data = _get_json_bigg(
                directory=directory,
                database=database,
                identifier=identifier,
                object_type="reaction",
                **kwargs,
            )
            msg_extra = "Identified as a reaction"
        debug_log.log(
            level=debug_level,
            msg=f'Data for "{identifier}" retrieved. ' + msg_extra,
        )
        return BiggParser._parse(root=json_data)

    @staticmethod
    def _parse(root: Any) -> dict:
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
    def _return_database(database: str) -> str:
        """
        Returns name of the database. It will raise a Error if name is
        incorrect.
        """
        if database == "BIGG":
            return database
        else:
            raise WrongParserError

    @staticmethod
    def _read_file(filename: Path) -> dict:
        """
        Reads the given file a returns a JSON dictionary with most important
        information from it.
        """
        try:
            with open(file=filename, mode="r") as f:
                unformatted_data = f.read()
            return loads(s=unformatted_data)
        except JSONDecodeError:
            # TODO find exception type
            raise WrongParserError("Wrong filetype")
