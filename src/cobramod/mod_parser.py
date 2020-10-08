#!/usr/bin/env python3
from pathlib import Path
from cobramod.debug import debug_log
from cobramod.parsing import biocyc, kegg
from typing import Any
import requests
from itertools import chain
import xml.etree.ElementTree as ET


def _define_base_dir(directory: Path, database: str) -> Path:
    """
    Returns Path object for given database. If directory does not exist.
    It will be created.

    Args:
        directory (Path): Parent directory.
        database (str): Name of database. Options: "META", "ARA". "KEGG"

    Returns:
        Path: Path object for database.
    """
    if directory.joinpath(database).exists():
        return directory.joinpath(database)
    else:
        directory.joinpath(database).mkdir()
        return directory.joinpath(database)


def _get_unformatted_kegg(directory: Path, identifier: str) -> str:
    """
    Retrieves and stores the data of given identifier for the KEGG database

    Args:
        directory (Path): Path for data storage.
        identifier (str): official name for given object in KEGGs

    Raises:
        NotImplementedError: If identifier contains a semicolon (:)
        Warning: If given identifier is not found in KEGGs
        NotADirectoryError:  If parent directory is not found

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
        data_dir = _define_base_dir(
            directory=directory, database=database)
        filename = data_dir.joinpath(f'{identifier}.txt')
        debug_log.debug(f'Searching "{identifier}" in directory "{database}"')
        try:
            with open(file=filename, mode="r") as f:
                unformatted_data = f.read()
                debug_log.debug('Found')
            return unformatted_data
        except FileNotFoundError:
            debug_log.warning(
                f'"{identifier}" not found in directory "{database}".')
            # Retrieve from URL
            url_text = (
                f'http://rest.kegg.jp/get/{identifier}/')
            debug_log.debug(f'Searching in {url_text}')
            r = requests.get(url_text)
            if r.status_code == 404:
                msg = f'"{identifier}" not available in "{database}"'
                debug_log.error(msg)
                raise Warning(msg)
            else:
                unformatted_data = r.text
                debug_log.debug(
                    f'Object found and saved in directory "{database}".')
                with open(file=filename, mode="w+") as f:
                    f.write(unformatted_data)
                return unformatted_data
    else:
        msg = "Directory not found"
        debug_log.error(msg)
        raise NotADirectoryError(msg)


def _parse_kegg(raw: str) -> dict:
    """
    Parses raw text from KEGG database and transforms it in a dictionary that
    can be used later as commom data type for cobramod. It idenfies the type
    of object depending of it information

    Args:
        raw (str): raw string of text from KEGGs

    Returns:
        dict: data from KEGG
    """
    kegg_dict = kegg._create_dict(raw=raw)
    # Will always have an ENTRY KEY
    # FIXME:  Molecular
    kegg_dict["TYPE"] = list(
        chain.from_iterable(item.split() for item in kegg_dict["ENTRY"]))[-1]
    kegg_dict["ENTRY"] = list(
        chain.from_iterable(item.split() for item in kegg_dict["ENTRY"]))[-2]
    # Only first name
    kegg_dict["NAME"] = kegg_dict["NAME"][0]
    kegg_dict["DATABASE"] = "KEGG"
    if kegg_dict["TYPE"] == "Compound":
        kegg_dict["FORMULA"] = "".join(kegg_dict["FORMULA"]).strip().rstrip()
        # FIXME: search information about charges
        kegg_dict["CHARGE"] = 0
    # For reactions only
    elif kegg_dict["TYPE"] == "Reaction":
        kegg_dict["BOUNDS"] = kegg._get_reversibility(
            line=kegg_dict["EQUATION"][0])
        kegg_dict["EQUATION"] = kegg._give_metabolites(
            line=kegg_dict["EQUATION"][0])
    else:
        raise NotImplementedError(
            'Could not parse given root. Please inform maintainers.')
    return kegg_dict


def _get_xml_from_biocyc(
        directory: Path, identifier: str, database: str) -> ET.Element:
    """
    Searchs in given parent directory if data is located in their respective
    database directory. If not, data will be retrievied from the corresponding
    database. Returns root of given identifier.

    Args:
        directory (Path): Path to directory where data is located.
        identifier (str): identifier for given database.
        database (str): Name of database. Options: "META", "ARA".

    Raises:
        Warning: If object is not available in given database
        NotADirectoryError: If parent directory is not found

    Returns:
        ET.Element: root of XML file
    """
    if directory.exists():
        data_dir = _define_base_dir(
            directory=directory, database=database)
        filename = data_dir.joinpath(f'{identifier}.xml')
        debug_log.debug(f'Searching "{identifier}" in directory "{database}"')
        try:
            root = ET.parse(str(filename)).getroot()
            debug_log.debug('Found')
            return root
        except FileNotFoundError:
            debug_log.warning(
                f'"{identifier}" not found in directory "{database}".')
            # Retrieve from URL
            url_text = (
                f'https://websvc.biocyc.org/getxml?{database}:{identifier}')
            debug_log.debug(f'Searching in {url_text}')
            r = requests.get(url_text)
            if r.status_code == 404:
                msg = f'"{identifier}" not available in "{database}"'
                debug_log.error(msg)
                raise Warning(msg)
            else:
                root = ET.fromstring(r.text)  # defining root
                tree = ET.ElementTree(root)
                tree.write(str(filename))
                debug_log.debug(
                    f'Object found and saved in directory "{database}".')
                return root
    else:
        msg = "Directory not found"
        debug_log.error(msg)
        raise NotADirectoryError(msg)


def _parse_biocyc(root: Any) -> dict:
    """
    Parses the root and returns a dictionary with the dictionary of the root,
    with the most important information depending on its object type.
    """
    if root.find("Compound") or root.find("Protein"):
        biocyc_dict = biocyc._p_compound(root=root)
    elif root.find("Reaction"):
        biocyc_dict = biocyc._p_reaction(root=root)
    elif root.find("Pathway"):
        biocyc_dict = biocyc._p_pathway(root=root)
    else:
        raise NotImplementedError(
            'Could not parse given root. Please inform maintainers.')
    biocyc_dict["DATABASE"] = root.find("*[@frameid]").attrib["orgid"]
    return biocyc_dict


def get_data(
        directory: Path,
        identifier: str,
        database: str) -> dict:
    """
    Retrieves and tranform the data into a dictionary for given identifier
    in a specific database.

    Args:
        directory (Path): Directory to store and retrieve local data.
        identifier (str): original identifier
        database (str): Name of database. Options: "META", "ARA", "KEGG"

    Returns:
        dict: relevant data for given identifier
    """
    def _get_kegg(directory: Path, identifier: str) -> dict:
        raw = _get_unformatted_kegg(directory=directory, identifier=identifier)
        return _parse_kegg(raw=raw)

    def _get_biocyc(
            directory: Path, identifier: str, database: str) -> dict:
        root = _get_xml_from_biocyc(
            directory=directory, identifier=identifier, database=database)
        return _parse_biocyc(root=root)

    biocyc_names = ["META", "ARA"]
    if database in biocyc_names:
        return _get_biocyc(
            directory=directory, identifier=identifier, database=database)
    elif database == "KEGG":
        return _get_kegg(directory=directory,  identifier=identifier)
    else:
        raise Warning("Give database does not exits.") from None
