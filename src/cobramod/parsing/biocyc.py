#!/usr/bin/env python3
from cobramod.debug import debug_log
from cobramod.parsing.base import BaseParser
from typing import Any, Dict
import xml.etree.ElementTree as ET
from pathlib import Path
import requests


def _p_compound(root: Any) -> dict:
    """
    Parses given xml root into a dictionary for metabolite from biocyc
    """
    identifier = root.find("*[@frameid]").attrib["frameid"]
    try:
        formula = (
            root.find("./*/cml/*/formula").attrib["concise"].replace(" ", "")
        )
        charge = int(root.find("./*/cml/molecule").attrib["formalCharge"])
    # Must be Protein
    except AttributeError:
        formula = "X"
        charge = 0
        debug_log.warning(f'Biocyc ID "{identifier}" treated as Protein')
    # For names only
    try:
        name = root.find("./*/cml/*").attrib["title"]
    except AttributeError:
        name = identifier
    # Temporal fixup
    if formula == "":
        formula = "X"
    temp_dict = {
        "TYPE": root.find("*[@frameid]").tag,
        "ENTRY": identifier,
        "NAME": name,
        "FORMULA": formula,
        "CHARGE": charge,
    }
    return temp_dict


def _p_metabolites(root: Any) -> dict:
    """
    Parses the data of given root and returns a dictionary with the
    information.
    """
    meta_dict: Dict[str, float] = dict()
    left_metabolites = root.findall("./Reaction/left")
    right_metabolites = root.findall("./Reaction/right")
    for meta in left_metabolites:
        try:
            coef = float(meta.find("coefficient").text) * -1
        except AttributeError:
            coef = -1  # default
        try:
            meta_identifier = (
                meta.find("*/[@frameid]").attrib["frameid"].strip().rstrip()
            )
        except AttributeError:
            raise AttributeError("Reaction cannot find participants")
        meta_dict[f"l_{meta_identifier}"] = coef
    for meta in right_metabolites:
        try:
            coef = float(meta.find("coefficient").text)
        except AttributeError:
            coef = 1  # default
        try:
            meta_identifier = (
                meta.find("*/[@frameid]").attrib["frameid"].strip().rstrip()
            )
        except AttributeError:
            raise AttributeError("Reaction cannot find participants")
        meta_dict[f"r_{meta_identifier}"] = coef
    return meta_dict


def _check_direction(root: Any) -> tuple:
    """
    Veifies that the direction of the reactions is the same as stated in
    the root file.
    """
    # Reversible <->
    text = root.find("*/reaction-direction").text
    if "REVERSIBLE" in text:
        bounds = (-1000, 1000)
    elif "RIGHT-TO-LEFT" in text:
        bounds = (-1000, 0)
    elif "LEFT-TO-RIGHT" in text:
        bounds = (0, 1000)
    else:
        raise Warning(
            f"Reversibility for "
            f'"{root.find("*[@frameid]").attrib["frameid"]}" could not '
            f"be found"
        )
    return bounds


def _p_reaction(root: Any) -> dict:
    """
    Parses the root file and return the data in a dictionary
    """
    identifier = root.find("*[@frameid]").attrib["frameid"]
    try:
        name = root.find("*[@ID]/enzymatic-reaction/*/common-name").text
    except AttributeError:
        name = identifier
    temp_dict = {
        "TYPE": root.find("*[@frameid]").tag,
        "ENTRY": identifier,
        "NAME": name,
        "EQUATION": _p_metabolites(root=root),
        "BOUNDS": _check_direction(root=root),
        "TRANSPORT": BaseParser._check_transport(
            data_dict=_p_metabolites(root=root)
        ),
    }
    return temp_dict


def _get_unsorted_pathway(root: Any) -> tuple:
    """
    Returns a dictionary with sequences (edges) for a graph; and a set
    with all participants (vertex).
    """
    reaction_dict = dict()
    reaction_set = set()
    for rxn_line in root.findall("*/reaction-ordering"):
        current = rxn_line.find("*/[@frameid]").attrib["frameid"]
        prior = rxn_line.find("predecessor-reactions/").attrib["frameid"]
        # If the direction of keys and values changes, then
        # many reactions would get lost. This way, even with
        # multiple compounds, pathways remain
        # NOTE: check for behaviour
        # Replacing values produces cuts with are needed to avoid cyclic
        # No reactions are missed
        reaction_dict[current] = prior
        # TODO: add information
        reaction_set.add(current)
        reaction_set.add(prior)
    name = root.find("Pathway").attrib["frameid"]
    # If dictinary has only one element, raise an error.
    if len(reaction_dict) == 0:
        raise NotImplementedError("Path has only a reaction. Add separately")
    debug_log.debug(f'Dictionary for pathway "{name}" succesfully created')
    return reaction_dict, reaction_set


def _p_pathway(root: Any) -> dict:
    """
    Return dictionary with data for given xml root.
    """
    identifier = root.find("*[@frameid]").attrib["frameid"]
    reaction_dict, reaction_set = _get_unsorted_pathway(root=root)
    temp_dict = {
        "TYPE": root.find("*[@frameid]").tag,
        "NAME": root.find("*[@frameid]").attrib["frameid"],
        "ENTRY": identifier,
        "PATHWAY": reaction_dict,
        "SET": reaction_set,
    }
    return temp_dict


class BiocycParser(BaseParser):
    @staticmethod
    def _parse(root: Any) -> dict:
        """
        Parses the root and returns a dictionary with the dictionary of the
        root, with the most important information depending on its object type.
        """
        if root.find("Compound") or root.find("Protein"):
            biocyc_dict = _p_compound(root=root)
        elif root.find("Reaction"):
            biocyc_dict = _p_reaction(root=root)
        elif root.find("Pathway"):
            biocyc_dict = _p_pathway(root=root)
        else:
            raise NotImplementedError(
                "Could not parse given root. Please inform maintainers."
            )
        biocyc_dict["DATABASE"] = root.find("*[@frameid]").attrib["orgid"]
        return biocyc_dict

    @staticmethod
    def _retrieve_data(
        directory: Path, identifier: str, database: str, debug_level: int
    ) -> dict:
        """
        Retrieves data from biocyc and parses the most important attributes
        into a dictionary.

        Args:
            directory (Path): Directory to store and retrieve local data.
            identifier (str): original identifier
            database (str): Name of the database. Some options: "META", "ARA"
        debug_level (int): Level of debugging. Read package logging
            for more info.

        Returns:
            dict: relevant data for given identifier
        """
        root = _get_xml_from_biocyc(
            directory=directory, identifier=identifier, database=database
        )
        debug_log.log(
            level=debug_level, msg=f'Data for "{identifier}" retrieved.'
        )
        return BiocycParser._parse(root=root)

    @staticmethod
    def _return_database(database: str) -> str:
        """
        Returns the name of the database. This method is used to compare with
        given database name. It will raise a warning if both names are not
        equal or belong to the list of proper names.
        """
        names = ["META", "ARA"]
        if database in names:
            return database
        else:
            raise Warning(f'Given database "{database}" does not exist')


def _get_xml_from_biocyc(
    directory: Path, identifier: str, database: str
) -> ET.Element:
    """
    Searchs in given parent directory if data is located in their respective
    database directory. If not, data will be retrievied from the corresponding
    database. Returns root of given identifier.

    Args:
        directory (Path): Path to directory where data is located.
        identifier (str): identifier for given database.
        database (str): Name of database. Options: "META", "ARA".

    Raes:
        Warning: If object is not available in given database
        NotADirectoryError: If parent directory is not found

    Returns:
        ET.Element: root of XML file
    """
    if directory.exists():
        data_dir = BiocycParser._define_base_dir(
            directory=directory, database=database
        )
        filename = data_dir.joinpath(f"{identifier}.xml")
        debug_log.debug(f'Searching "{identifier}" in directory "{database}"')
        try:
            root = ET.parse(str(filename)).getroot()
            debug_log.debug(f"Identifier '{identifier}' found.")
            return root
        except FileNotFoundError:
            debug_log.debug(
                f'"{identifier}" not found in directory "{database}".'
            )
            # Retrieve from URL
            url_text = (
                f"https://websvc.biocyc.org/getxml?{database}:{identifier}"
            )
            debug_log.debug(f"Searching in {url_text}")
            r = requests.get(url_text)
            if r.status_code == 404:
                msg = f'"{identifier}" not available in "{database}"'
                debug_log.error(msg)
                raise Warning(msg)
            else:
                root = ET.fromstring(r.text)  # defining root
                tree = ET.ElementTree(root)
                tree.write(str(filename))
                debug_log.info(
                    f'Object "{identifier}" found in database. Saving in '
                    f'directory "{database}".'
                )
                return root
    else:
        msg = "Directory not found"
        debug_log.critical(msg)
        raise NotADirectoryError(msg)
