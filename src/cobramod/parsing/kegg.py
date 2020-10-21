#!/usr/bin/env python3
from typing import Generator, Union, Dict, NamedTuple
from itertools import chain
from cobramod.debug import debug_log
from pathlib import Path
import requests
from cobramod.parsing.base import BaseParser


class MetaboliteTuple(NamedTuple):
    identifier: str
    coefficient: float


def _get_keys(raw: str) -> Generator:
    """
    Returns as Generator the keys for KEGGs data. Input is raw text
    """
    lines = (line for line in raw.split(sep="\n"))
    for line in lines:
        segment = line.split(" ")
        if segment[0] == "///":
            continue
        if segment[0] != "":
            yield segment[0]


def _find_key(line: str, keys: set) -> Union[str, None]:
    """
    Return key if found in given line. Else, None is returned.
    """
    for key in keys:
        if key in line:
            return key
    return None


def _create_dict(raw: str) -> dict:
    """
    Formats most of the keys for Keggs data.
    """
    lines = (line for line in raw.split("\n"))
    keys = set(_get_keys(raw=raw))
    actual_key = str()
    kegg_dict: Dict[str, list] = {"": []}
    for line in lines:
        key = _find_key(line=line, keys=keys)
        if key is None:
            # Append to older key.
            kegg_dict[actual_key].append(line.strip().rstrip())
        else:
            actual_key = key
            line = line[line.find(key) + len(key) :].strip().rstrip()
            kegg_dict[actual_key] = [line]
    del kegg_dict[""]
    return kegg_dict


def _get_coefficient(line: str, SIDE: int, prefix: str) -> MetaboliteTuple:
    """
    Returns a NamedTuple with the metabolite identifier and coefficient
    that appears in given line.

    Args:
        line (str): string with information
        SIDE (int): Constant to multiply depending on the side
        prefix (str): prefix for the metabolite

    Returns:
        NamedTuple: tupple with identifier and coefficient
    """
    # " 4 H+ "
    line_list = line.strip().rstrip().split(" ")
    try:
        return MetaboliteTuple(
            identifier=prefix + line_list[1],
            coefficient=float(line_list[0]) * SIDE,
        )
    except IndexError:
        return MetaboliteTuple(
            identifier=prefix + line_list[0], coefficient=1.0 * SIDE
        )


def _give_metabolites(line: str) -> dict:
    """
    Returns dictionary with metabolites identifiers and corresponding
    coeffient
    """
    middle = line.find("=")
    temp_dict = dict()
    reactants = line[: middle - 1]
    products = line[middle + 2 :]
    for item in reactants.split(" + "):
        MetaTuple = _get_coefficient(line=item, SIDE=-1, prefix="r_")
        temp_dict[MetaTuple.identifier] = MetaTuple.coefficient
    for item in products.split(" + "):
        MetaTuple = _get_coefficient(line=item, SIDE=1, prefix="l_")
        temp_dict[MetaTuple.identifier] = MetaTuple.coefficient
    return temp_dict


def _get_reversibility(line: str) -> tuple:
    """
    Returns the bounds for reaction depending of the string
    """
    # FIXME: Direction depends also from extra keys
    middle = line.find("=")
    line = line[middle - 1 : middle + 2]
    if line == "<=>":
        bounds = (-1000, 1000)
    elif line == "==>":
        bounds = (0, 1000)
    elif line == "<==":
        bounds = (-1000, 0)
    else:
        raise Warning(f'"Reversibility for "{line}" could not be found.')
    return bounds


class KeggParser(BaseParser):
    @staticmethod
    def _parse(raw: str) -> dict:
        """
        Parses raw text from KEGG database and transforms it in a dictionary
        that can be used later as commom data type for cobramod. It identifies
        the type of object depending of it information.

        Args:
            raw (str): raw string of text from KEGGs

        Returns:
            dict: data from KEGG
        """
        kegg_dict = _create_dict(raw=raw)
        # Will always have an ENTRY KEY
        # FIXME:  Molecular
        kegg_dict["TYPE"] = list(
            chain.from_iterable(item.split() for item in kegg_dict["ENTRY"])
        )[-1]
        kegg_dict["ENTRY"] = list(
            chain.from_iterable(item.split() for item in kegg_dict["ENTRY"])
        )[-2]
        # Only first name
        kegg_dict["NAME"] = kegg_dict["NAME"][0]
        kegg_dict["DATABASE"] = "KEGG"
        if kegg_dict["TYPE"] == "Compound":
            kegg_dict["FORMULA"] = (
                "".join(kegg_dict["FORMULA"]).strip().rstrip()
            )
            # FIXME: search information about charges
            kegg_dict["CHARGE"] = 0
        # For reactions only
        elif kegg_dict["TYPE"] == "Reaction":
            kegg_dict["BOUNDS"] = _get_reversibility(
                line=kegg_dict["EQUATION"][0]
            )
            kegg_dict["EQUATION"] = _give_metabolites(
                line=kegg_dict["EQUATION"][0]
            )
            kegg_dict["TRANSPORT"] = BaseParser._check_transport(
                data_dict=kegg_dict["EQUATION"]
            )
        else:
            raise NotImplementedError(
                "Could not parse given root. Please inform maintainers."
            )
        return kegg_dict

    @staticmethod
    def _retrieve_data(
        directory: Path, identifier: str, database: str, debug_level: int
    ) -> dict:
        # FIXME: temporal solution for database
        """
        Retrieves data from KEGG database and parses the most important
        attributes into a dictionary

        Args:
            directory (Path): Directory to store and retrieve local data.
            identifier (str): original identifier
            database (str): name of database.
        debug_level (int): Level of debugging. Read package logging
            for more info.

        Returns:
            dict: relevant data for given identifier
        """
        raw = _get_unformatted_kegg(directory=directory, identifier=identifier)
        debug_log.log(
            level=debug_level, msg=f'Data for "{identifier}" retrieved.'
        )
        return KeggParser._parse(raw=raw)

    @staticmethod
    def _return_database(database: str) -> str:
        """
        Returns the name of the database. This method is used to compare with
        given database name. It will raise a warning if both names are not
        equal or belong to the list of proper names.
        """
        if database == "KEGG":
            return database
        else:
            raise Warning(f'Given database "{database}" does not exist')


def _get_unformatted_kegg(directory: Path, identifier: str) -> str:
    """
    Retrieves and stores the data of given identifier for the KEGG database.

    Args:
        directory (Path): Path for data storage.
        identifier (str): official name for given object in KEGGs.

    Raises:
        NotImplementedError: If identifier contains a semicolon (:)
        Warning: If given identifier is not found in KEGGs.
        NotADirectoryError:  If parent directory is not found.

    Returns:
        str: raw string from given KEGGs identifier
    """
    # NOTE: As KEGG only support xml for pathways, Metabolites and Reactions
    # have to be parsed differently. Only specific keys are necesary
    database = "KEGG"
    if ":" in identifier:
        raise NotImplementedError(
            "Identifiers with semicolons cannot be parse yet"
        )
    if directory.exists():
        data_dir = KeggParser._define_base_dir(
            directory=directory, database=database
        )
        filename = data_dir.joinpath(f"{identifier}.txt")
        debug_log.debug(f'Searching "{identifier}" in directory "{database}".')
        try:
            with open(file=filename, mode="r") as f:
                unformatted_data = f.read()
                debug_log.debug(f"Identifier '{identifier}' found.")
            return unformatted_data
        except FileNotFoundError:
            debug_log.debug(
                f'"{identifier}" not found in directory "{database}".'
            )
            # Retrieve from URL
            url_text = f"http://rest.kegg.jp/get/{identifier}/"
            debug_log.debug(f"Searching in {url_text}")
            r = requests.get(url_text)
            if r.status_code == 404:
                msg = f'"{identifier}" not available in "{database}".'
                debug_log.error(msg)
                raise Warning(msg)
            else:
                unformatted_data = r.text
                debug_log.info(
                    f'Object "{identifier}" found in database. Saving in '
                    f'directory "{database}".'
                )
                with open(file=filename, mode="w+") as f:
                    f.write(unformatted_data)
                return unformatted_data
    else:
        msg = "Directory not found"
        debug_log.critical(msg)
        raise NotADirectoryError(msg)
