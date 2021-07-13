#!/usr/bin/env python3
"""Data parsing for PlantCyc

This module handles the retrieval of data from PlantCyc into a local directory.
The possible type of data that can be download:

- Metabolites: Normally have an abbreviation or short name.
- Reactions: Can have the words "RXN" in the identifier. Enzymes can sometimes
be used instead. The gene information for the reactions is include if found.
- Pathways

Contact maintainers if other types should be added.

Important class of the module:
- PlantCycParser: Child of the abstract class
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
from warnings import warn

from requests import get, HTTPError

from cobramod.debug import debug_log
from cobramod.error import (
    WrongParserError,
    NoGeneInformation,
    SuperpathwayWarning,
)
from cobramod.parsing.base import BaseParser


def _build_reference(root: Any) -> dict:
    """
    Returns a dictionary with the corresponding cross-references and their
    identifiers.
    """
    references = dict()
    try:
        for child in root.findall("*/dblink"):
            references[child.find("dblink-db").text] = child.find(
                "dblink-oid"
            ).text
        return references
    except IndexError:
        raise IndexError("No references can be found in given XML.")


def _p_compound(root: Any) -> dict:
    """
    Returns a dictionary for the representation of a metabolite
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
            # ToDo hier ändern
            f'Chemical formula for the PlantCyc metabolite "{identifier}" '
            'could not be found. Formula set to "X" and charge to 0. '
            "Please modify it if necessary."
        )
    # In case formula is a empty string
    if not formula:
        formula = "X"
        debug_log.warning(
            # ToDo hier ändern
            f'Chemical formula for the PlantCyc metabolite "{identifier}" '
            'could not be found. Formula set to "X". '
            "Please modify it if necessary."
        )
    # obtain name
    try:
        name = root.find("./*/cml/*").attrib["title"]
    except AttributeError:
        name = identifier
    return {
        "TYPE": root_type,
        "ENTRY": identifier,
        "NAME": name,
        "FORMULA": formula,
        "CHARGE": charge,
        "DATABASE": root.find("*[@frameid]").attrib["orgid"],
        "XREF": _build_reference(root=root),
    }


def _get_metabolites(root: Any) -> dict:
    """
    Returns a dictionary with the metabolites that participates in the reaction
    of given XML.
    """
    metabolites: Dict[str, float] = dict()
    reactants = root.findall("./Reaction/left")
    products = root.findall("./Reaction/right")
    for side in (reactants, products):
        # Define factor
        FACTOR = 1
        PREFIX = "r"
        if side is reactants:
            FACTOR = -1
            PREFIX = "l"
        # Add metabolites to the dictionary
        for metabolite in side:
            try:
                coefficient = (
                        float(metabolite.find("coefficient").text) * FACTOR
                )
            except (AttributeError, ValueError):
                # Value errors are for the strings like 'n'
                coefficient = 1 * FACTOR  # default
            try:
                identifier = (
                    metabolite.find("*/[@frameid]")
                        .attrib["frameid"]
                        .strip()
                        .rstrip()
                )
            except AttributeError:
                raise AttributeError("Reaction cannot find participants")
            metabolites[f"{PREFIX}_{identifier}"] = coefficient
    return metabolites


def _check_direction(root: Any) -> tuple:
    """
    Verifies that the direction of the reactions is the same as stated in
    the root file.
    """
    text = root.find("*/reaction-direction").text
    reversibility = {
        "REVERSIBLE": (-1000, 1000),
        "RIGHT-TO-LEFT": (-1000, 0),
        "PHYSIOL-RIGHT-TO-LEFT": (-1000, 0),
        "LEFT-TO-RIGHT": (0, 1000),
        "PHYSIOL-LEFT-TO-RIGHT": (0, 1000),
    }
    try:
        return reversibility[text]
    except KeyError:
        raise KeyError(
            "Reversibility for "
            + f'"{root.find("*[@frameid]").attrib["frameid"]}" could not '
            + f"be found. The string shows: \n {text}"
        )


def _p_genes(root: Any, identifier: str, directory: Path) -> dict:
    """
    Returns a dictionary with the corresponding genes and gene-reaction rules.
    This function will try to read a file for given identifier. If nothing is
    found, the dictionary will have empty entries.
    """
    rule = str()
    genes = dict()
    # Get the information and check if Genes can be found to be parsed
    with suppress(FileNotFoundError):
        tree: Any = PlantCycParser._read_file(
            filename=directory.joinpath("GENES").joinpath(
                f"{identifier}_genes.xml"
            )
        )

        for gene in tree.findall("Gene"):
            try:
                # Sometimes the genes have a name
                name = gene.find("common-name").text
            except AttributeError:
                name = identifier
            genes[gene.attrib["frameid"]] = name
        # Assuming rule
        rule = " or ".join(genes.keys())
        debug_log.warning(
            f'Gene-reaction rule for reaction "{identifier}" assumed'
            + ' to be "OR". Please modify it if necessary.'
        )
    return {"genes": genes, "rule": rule}


def _p_reaction(root: Any, directory: Path) -> dict:
    """
    Return the dictionary for the representation of a reaction
    """

    identifier = root.find("*[@frameid]").attrib["frameid"]
    try:
        name = root.find("*[@ID]/enzymatic-reaction/*/common-name").text
    except AttributeError:
        name = identifier
    database = root.find("*[@frameid]").attrib["orgid"]
    return {
        "TYPE": root.find("*[@frameid]").tag,
        "ENTRY": identifier,
        "NAME": name,
        "EQUATION": _get_metabolites(root=root),
        "BOUNDS": _check_direction(root=root),
        "TRANSPORT": BaseParser._check_transport(
            data_dict=_get_metabolites(root=root)
        ),
        "DATABASE": database,
        "XREF": _build_reference(root=root),
        "GENES": _p_genes(
            root=root,
            identifier=identifier,
            directory=directory.joinpath(database),
        ),
    }


def get_graph(root: Any) -> dict:
    """
    Return a directed graph from given XML Element. A key represent the
    parent reaction and the value is the child. A parent can have multiple
    children as tuples. End-reactions have a NoneType as child.

    Args:
        root (Element): An XML element, which represent a XML file with the
            information of a pathway.

    Returns:
        Dict: Directed graph from root.

    Raises:
        WrongParserError: If given root does not represent a pathway.
    """
    # Verify non-superpathway
    with suppress(TypeError):
        identifier = root.find("*[@frameid]").attrib["frameid"]
        for parent in root.findall("Pathway/parent/*"):
            if "Super" in parent.get("frameid"):
                msg = (
                        f'Pathway "{identifier}" was identified as a '
                        + "superpathway. This type of pathway does not normally "
                        + "included all reactions. Please add the corresponding "
                        + "sub-pathways singlely!"
                )
                debug_log.warning(msg=msg)
                warn(message=msg, category=SuperpathwayWarning)
                break
    # Parent: child
    graph: Dict[str, Any] = dict()
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
                graph[parent] = graph[parent] + (child,)
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
    # Verify if Superpathway
    return graph


def _p_pathway(root: Any) -> dict:
    """
    Return dictionary for the representation of a pathway.
    """
    identifier = root.find("*[@frameid]").attrib["frameid"]
    reaction_dict = get_graph(root=root)
    pathway = {
        "TYPE": root.find("*[@frameid]").tag,
        "NAME": root.find("*[@frameid]").attrib["frameid"],
        "ENTRY": identifier,
        "PATHWAY": reaction_dict,
        "DATABASE": root.find("*[@frameid]").attrib["orgid"],
        "XREF": _build_reference(root=root),
    }
    return pathway


# ToDo hier ändern
class PlantCycParser(BaseParser):
    @staticmethod
    def _parse(root: Any, directory: Path) -> dict:
        """
        Parses the root and returns a dictionary with the dictionary of the
        root, with the most important information depending on its object type.
        """
        # try to define dictionary
        for method in (_p_reaction, _p_pathway, _p_compound):
            with suppress(WrongParserError, AttributeError):
                try:
                    plantcyc_dict = method(  # type: ignore
                        root=root, directory=directory
                    )
                except TypeError:
                    plantcyc_dict = method(root=root)  # type: ignore
                # This might change due to sub-database
                plantcyc_dict["DATABASE"] = root.find("*[@frameid]").attrib[
                    "orgid"
                ]
                return plantcyc_dict
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
        Retrieves data from PlantCyc and parses the most important attributes
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

        root = retrieve_data(
            directory=directory, identifier=identifier, database=database
        )
        debug_log.log(
            level=debug_level, msg=f'Data for "{identifier}" retrieved.'
        )
        return PlantCycParser._parse(root=root, directory=directory)

    @staticmethod
    def _check_database(database: str):
        """
        Returns name of the database. It will raise a Error if name is
        incorrect.
        """
        if not isinstance(database, str) or not database.startswith("pmn:"):
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


def _get_gene_xml(directory: Path, identifier: str, database: str):
    """
    Retrieve the gene information for given reaction in a specific database

    Args:
        directory (Path): Directory to store and retrieve local data.
        identifier (str): original identifier
        database (str): Name of the database. Some options: "META", "ARA"

    Raises:
        HTTPError: If given reaction does not have gene information available
    """
    # GENES directory will depend from the sub-database
    directory = directory.joinpath("GENES")
    if not directory.exists():
        directory.mkdir()
    # Retrieval of the Gene information
    filename = directory.joinpath(f"{identifier}_genes.xml")
    if not filename.exists():
        # This URL will not necessarily raise exception
        url_text = (
            f"https://pmn.plantcyc.org/apixml?fn=genes-of-reaction&id="
            f"{database}:{identifier}&detail=full"
        )
        response = get(url_text)
        try:
            response.raise_for_status()
            root = fromstring(response.text)
            # Check for results in root
            tree: Any = ElementTree(root)

            if int(tree.find("*/num_results").text) == 0:
                raise NoGeneInformation
            if database == "META":
                msg = (
                    f'Object "{identifier}" comes from "META". Please use '
                    "another sub-database from PlantCyc to add proper genes. "
                    "Skipping retrieval of gene information."
                )
                debug_log.warning(msg)
                warn(message=msg, category=UserWarning)
                raise NoGeneInformation(msg)
            tree.write(str(filename))
            debug_log.info(
                f'Object "{identifier}_gene.xml" saved in '
                f'directory "{database}/GENES".'
            )
        except HTTPError:
            # Warning
            raise NoGeneInformation(
                f"Object {identifier} does not have gene information"
            )


def retrieve_data(directory: Path, identifier: str, database: str) -> Element:
    """
    Searchs in given parent directory if data is located in their respective
    database directory. If not, data will be retrieved from the corresponding
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
    # remove pmn: from database str
    database = database[4:]

    if directory.exists():
        data_dir = directory.joinpath("PMN", database)
        data_dir.mkdir(parents=True, exist_ok=True)
        filename = data_dir.joinpath(f"{identifier}.xml")
        debug_log.debug(f'Searching "{identifier}" in directory "{database}"')
        # Search for the file on directory. Otherwise retrive from database
        try:
            return PlantCycParser._read_file(filename=filename)
        except FileNotFoundError:
            debug_log.debug(
                f'"{identifier}" not found in directory "{database}".'
            )
            # Retrieve from URL
            url_text = (
                f"https://pmn.plantcyc.org/getxml?{database}:{identifier}"
            )
            debug_log.debug(f"Searching in {url_text}")
            response = get(url_text)
            try:
                response.raise_for_status()
                # defining root
                root = fromstring(response.text)
                tree = ElementTree(root)
                tree.write(str(filename))
                debug_log.info(
                    f'Object "{identifier}" found in database. Saving in '
                    f'directory "{database}".'
                )
                # Obtain genes if possible. This should be only call one time
                # If information is already available, the genes should be
                # available if found
                with suppress(NoGeneInformation):
                    # This will include Paths
                    _get_gene_xml(
                        directory=data_dir,
                        identifier=identifier,
                        database=database,
                    )
                return root
            except HTTPError:
                msg = f'"{identifier}" is not available in "{database}"'
                debug_log.critical(msg)
                # Try with META
                if database != "PLANT":
                    # This will raise an error if not found in META
                    with suppress(HTTPError):
                        root = retrieve_data(
                            directory=directory,
                            identifier=identifier,
                            database="pmn:PLANT",
                        )
                        return root
                # In case there is nothing in META
                raise HTTPError(msg)
    else:
        msg = "Directory not found. Please create the given directory."
        debug_log.critical(msg)
        raise NotADirectoryError(msg)
