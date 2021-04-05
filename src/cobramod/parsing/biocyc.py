#!/usr/bin/env python3
"""Data parsing for Biocyc

This module handles the retrieval of data from BioCyc into a local directory.
The posible type of data that can be download:

- Metabolites: Normally have an abbreviation or short name.
- Reactions: Can have the words "RXN" in the identifier. Enzymes can sometimes
be used instead.
- Pathways

Contact maintainers if other types should be added.

Important class of the module:
- BiocycParser: Child of the abstract class
:class:`cobramod.parsing.base.BaseParser`.
"""
from contextlib import suppress
from xml.etree.ElementTree import (
    fromstring,
    Element,
    ElementTree,
    parse,
    ParseError,
)
from pathlib import Path
from typing import Any, Dict

from requests import get, HTTPError

from cobramod.debug import debug_log
from cobramod.error import WrongParserError
from cobramod.parsing.base import BaseParser


def _build_reference(root: Any) -> dict:
    """
    From the original root object, returns a dictionary with the corresponding
    crossref and its identifiers.
    """
    references = dict()
    try:
        for child in root.findall("*/dblink"):
            references[child.find("dblink-db").text] = child.find(
                "dblink-oid"
            ).text
        return references
    except IndexError:
        raise IndexError("Current root could not be find references.")


def _p_compound(root: Any) -> dict:
    """
    Parses given xml root into a dictionary for metabolite from biocyc
    """
    identifier = root.find("*[@frameid]").attrib["frameid"]
    root_type = root.find("*[@frameid]").tag
    try:
        formula = (
            root.find("./*/cml/*/formula").attrib["concise"].replace(" ", "")
        )
        charge = int(root.find("./*/cml/molecule").attrib["formalCharge"])
    # Must be Protein
    except AttributeError:
        formula = "X"
        charge = 0
        debug_log.warning(
            f'Biocyc ID "{identifier}" could not find a chemical formula. '
            f'Formula set to "X" and charge to 0'
        )
    # For names only
    try:
        name = root.find("./*/cml/*").attrib["title"]
    except AttributeError:
        name = identifier
    # Temporal fixup
    if formula == "":
        formula = "X"
    return {
        "TYPE": root_type,
        "ENTRY": identifier,
        "NAME": name,
        "FORMULA": formula,
        "CHARGE": charge,
        "DATABASE": root.find("*[@frameid]").attrib["orgid"],
        "XREF": _build_reference(root=root),
    }


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
        except (AttributeError, ValueError):
            # Value errors are for the strings like 'n'
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
        except (AttributeError, ValueError):
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


def _p_genes(root: Any, identifier: str):
    database = "META"
    url_text = (
        f"https://websvc.biocyc.org/apixml?fn=genes-of-reaction&id={database}:"
        f"{identifier}&detail=full"
    )
    r = get(url_text)
    if r.status_code >= 400:
        rule = None
        genes = None
    else:
        root = fromstring(r.text)  # defining root
        tree: Any = ElementTree(root)
        genes = []
        for gene in tree.findall("Gene"):
            try:
                name = gene.find("common-name").text
            except AttributeError:
                name = identifier
            genes.append({"identifier": gene.attrib["frameid"], "name": name})
        # FIXME: Temporal OR rule
        rule = " or ".join([test["identifier"] for test in genes])
    return {"genes": genes, "rule": rule}


def _p_reaction(root: Any) -> dict:
    """
    Parses the root file and return the data in a dictionary
    """
    identifier = root.find("*[@frameid]").attrib["frameid"]
    try:
        name = root.find("*[@ID]/enzymatic-reaction/*/common-name").text
    except AttributeError:
        name = identifier
    return {
        "TYPE": root.find("*[@frameid]").tag,
        "ENTRY": identifier,
        "NAME": name,
        "EQUATION": _p_metabolites(root=root),
        "BOUNDS": _check_direction(root=root),
        "TRANSPORT": BaseParser._check_transport(
            data_dict=_p_metabolites(root=root)
        ),
        "DATABASE": root.find("*[@frameid]").attrib["orgid"],
        "XREF": _build_reference(root=root),
        "GENES": _p_genes(root=root, identifier=identifier),
    }


def get_graph(root: Any) -> dict:
    """
    Return a directed graph from given XML Element. A key represent the
    parent reaction and the value is the child. A parent can have multiple
    childs as tuples. End-reactions have a NoneType as child.

    Args:
        root (Element): An XML element, which represent a XML file with the
            information of a pathway.

    Returns:
        Dict: Directed graph from root.

    Raises:
        WrongParserError: If given root does not represent a pathway.
    """
    # Parent: child
    graph = dict()  # type: ignore
    reactions = set()
    # line: Element
    for line in root.findall("*/reaction-ordering"):
        # Get reactions and add them to set
        child = line.find("*/[@frameid]").attrib["frameid"]
        parent = line.find("predecessor-reactions/").attrib["frameid"]
        reactions.add(parent)
        reactions.add(child)
        # Expand keys
        try:
            # If at least 2
            if isinstance(graph[parent], tuple):
                graph[parent] += child
                continue
            # Else transform single to tuple
            graph[parent] = (graph[parent], child)
        except KeyError:
            # Create new one
            graph[parent] = child
    # Find missing elements and insert a None
    for reaction in reactions:
        try:
            graph[reaction]
        except KeyError:
            graph[reaction] = None
    if not graph:
        # It could be a single-reaction pathway. Give a None since it is just
        # one
        try:
            name = root.find("*/reaction-layout/Reaction").attrib["frameid"]
            graph[name] = None
        except AttributeError:
            raise WrongParserError("Given root does not belong to a Pathway")
    return graph


def _p_pathway(root: Any) -> dict:
    """
    Return dictionary with data for given xml root.
    """
    identifier = root.find("*[@frameid]").attrib["frameid"]
    reaction_dict = get_graph(root=root)
    temp_dict = {
        "TYPE": root.find("*[@frameid]").tag,
        "NAME": root.find("*[@frameid]").attrib["frameid"],
        "ENTRY": identifier,
        "PATHWAY": reaction_dict,
        "DATABASE": root.find("*[@frameid]").attrib["orgid"],
        "XREF": _build_reference(root=root),
    }
    return temp_dict


class BiocycParser(BaseParser):
    @staticmethod
    def _parse(root: Any) -> dict:
        """
        Parses the root and returns a dictionary with the dictionary of the
        root, with the most important information depending on its object type.
        """
        # try to define dictionary
        for method in (_p_reaction, _p_pathway, _p_compound):
            with suppress(WrongParserError, AttributeError):
                biocyc_dict = method(root=root)
                biocyc_dict["DATABASE"] = root.find("*[@frameid]").attrib[
                    "orgid"
                ]
                return biocyc_dict
        raise NotImplementedError(
            "Given root could not be parsed properly. Contact maintainers"
        )

    @staticmethod
    def _retrieve_data(
        directory: Path,
        identifier: str,
        database: str,
        debug_level: int,
        **kwargs,
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
            raise WrongParserError

    @staticmethod
    def _read_file(filename: Path) -> Element:
        """
        Reads and return given filename as a Element. I will raise a error if
        file is not an valid xml.
        """
        try:
            return parse(source=str(filename)).getroot()
        except ParseError:
            raise WrongParserError("Wrong filetype")


def _get_xml_from_biocyc(
    directory: Path, identifier: str, database: str
) -> Element:
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
            return BiocycParser._read_file(filename=filename)
        except FileNotFoundError:
            debug_log.debug(
                f'"{identifier}" not found in directory "{database}".'
            )
            # Retrieve from URL
            url_text = (
                f"https://websvc.biocyc.org/getxml?{database}:{identifier}"
            )
            debug_log.debug(f"Searching in {url_text}")
            r = get(url_text)
            if r.status_code >= 400:
                msg = f'"{identifier}" is not available in "{database}"'
                debug_log.error(msg)
                # Try with META
                if database != "META":
                    # This will raise an error if not found in META
                    root = _get_xml_from_biocyc(
                        directory=directory,
                        identifier=identifier,
                        database="META",
                    )
                    return root
                raise HTTPError(msg)
            else:
                root = fromstring(r.text)  # defining root
                tree = ElementTree(root)
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
