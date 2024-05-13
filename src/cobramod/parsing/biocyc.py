"""Data parsing for Biocyc

This module handles the retrieval of data from BioCyc into a local directory.
The possible type of data that can be downloaded:

- Metabolites: Normally have an abbreviation or short name.
- Reactions: Can have the words "RXN" in the identifier. Enzymes can sometimes
be used instead. The gene information for the reactions is included if found.
"""

from __future__ import annotations

import urllib.parse
import xml.etree.ElementTree as et
from contextlib import suppress
from pathlib import Path
from typing import Any

import requests

import cobramod.error as cmod_error
import cobramod.utils as cmod_utils
from cobramod.debug import debug_log


def build_cross_references_xml(root: Any) -> dict[str, str]:
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


def reaction_str_from_xml(root: et.Element) -> str:
    reactants = root.findall("./Reaction/left")
    products = root.findall("./Reaction/right")

    if not reactants or not products:
        raise AttributeError(
            "The current et.Element does not show any properties of being a"
            "reaction."
        )
    reaction_str: str = ""

    for side in (reactants, products):
        prefix = "r"

        if side is reactants:
            prefix = "l"

        for i, item in enumerate(side):
            element = item.find("coefficient")
            if not element:
                coefficient: float = 1.0

            else:
                value = element.text
                if not value:
                    raise AttributeError("Given Elements does not include text")

                coefficient = float(value)

            element = item.find("*/[@frameid]")

            if element is None:
                raise AttributeError("Metabolite was not found in the XML root")
            else:
                identifier = element.attrib["frameid"].strip().rstrip()

            reaction_str += f"{coefficient} {prefix}_{identifier}"

            if i < len(side) - 1:
                reaction_str += " + "

            if side is reactants and i == len(side) - 1:
                reaction_str += f" {get_direction_from_xml(root)} "

    return reaction_str


def get_direction_from_xml(root: et.Element) -> str:
    """
    Returns the direction of the chemical reaction
    """
    element = root.find("*/reaction-direction")
    if element is None:
        raise AttributeError(
            "Attribute 'reaction-direction' was not found in the given Element"
        )
    text = element.text
    if not text:
        raise AttributeError(
            "Given element does not include text and thus, reversibility "
            "cannot be infered"
        )

    reversibility = {
        "REVERSIBLE": "<->",
        "RIGHT-TO-LEFT": "<--",
        "PHYSIOL-RIGHT-TO-LEFT": "<--",
        "LEFT-TO-RIGHT": "-->",
        "PHYSIOL-LEFT-TO-RIGHT": "-->",
    }
    return reversibility[text]


def parse_reaction_attributes(
    root: et.Element, entry: str, gene_directory: Path
) -> dict[str, Any]:
    name_element = root.find("*[@ID]/enzymatic-reaction/*/common-name")
    if name_element is not None and name_element.text:
        name = name_element.text
    else:
        name = entry
    equation = reaction_str_from_xml(root)

    attributes = {
        "name": name,
        "equation": equation,
        "xref": build_cross_references_xml(root),
        "replacements": {},
        "transport": cmod_utils.is_transport(equation),
        "genes": parse_genes(entry, gene_directory),
    }
    return attributes


def parse_metabolite_attributes(root: et.Element, entry: str):
    try:
        element = root.find("./*/cml/*")

        if element is None:
            raise AttributeError("No compound was found in given XML root")
        name = element.attrib["title"]

    except AttributeError:
        name = entry

    try:
        element = root.find("./*/cml/*/formula")
        if element is None:
            raise AttributeError("No formula was found in given XML root")
        formula = element.attrib["concise"].replace(" ", "")

        element = root.find("./*/cml/molecule")
        if element is None:
            raise AttributeError("No formula was found in given XML root")
        charge = float(element.attrib["formalCharge"])

    # Must be Protein
    except AttributeError:
        formula = "X"
        charge = 0

    attributes: dict[str, Any] = {
        "name": name,
        "formula": formula,
        "charge": charge,
        "xref": build_cross_references_xml(root),
    }
    return attributes


def get_graph(root: et.Element) -> dict:
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
    element = root.find("*[@frameid]")
    if element is None:
        raise TypeError("Not a valid Element")

    identifier = element.get("frameid")

    for pathway in element.findall("*/Pathway"):
        frameid = pathway.get("frameid")

        if frameid is None:
            raise TypeError("Not a valid Element")

        if "Super" in frameid:
            msg = (
                f'Pathway "{identifier}" was identified as a '
                + "superpathway. This type of pathway does not normally "
                + "included all reactions. Please add the corresponding "
                + "sub-pathways singlely!. "
                # NOTE: ToDo
                + "For the moment, CobraMod does not support them"
            )
            debug_log.error(msg=msg)
            raise cmod_error.SuperpathwayException(msg)

    # identifiers = [
    #     a.attrib["frameid"]
    #     for a in root.findall("*/reaction-ordering/Reaction")
    # ]
    graph: dict[str, Any] = dict()
    reactions = set()

    for line in root.findall("*/reaction-ordering"):
        # Get reactions and add them to set
        element = line.find("*/[@frameid]")
        if element is None:
            raise TypeError("Not a valid Element")
        child = element.attrib["frameid"]

        element = line.find("predecessor-reactions/")
        if element is None:
            raise TypeError("Not a valid Element")
        parent = element.attrib["frameid"]

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
            element = root.find("*/reaction-layout/Reaction")
            if element is None:
                raise TypeError("Not a valid Element")
            name = element.attrib["frameid"]
            graph[name] = None
        except AttributeError:
            raise cmod_error.WrongParserError(
                "Given root does not belong to a Pathway"
            )

    return graph


def parse_pathway_attributes(root: et.Element, entry: str) -> dict[str, Any]:
    graph = get_graph(root)
    element = root.find("*[@frameid]")

    if element is not None:
        name = element.attrib["frameid"]
    else:
        name = entry

    attributes = {
        "name": name,
        "xref": build_cross_references_xml(root),
        "pathway": graph,
    }
    return attributes


def retrieve_gene_information(directory: Path, identifier: str, database: str):
    """
    Retrieve the gene information for given reaction in a specific database
    Args:
        directory (Path): Directory to store and retrieve local data.
        identifier (str): original identifier
        database (str): Name of the database. Some options: "META", "ARA".
            Full list in: "https://biocyc.org/biocyc-pgdb-list.shtml"
    Raises:
        HTTPError: If given reaction does not have gene information available
    """
    # GENES directory will depend from the sub-database
    directory = directory.joinpath(database, "GENES")

    if not directory.exists():
        directory.mkdir()

    # Retrieval of the Gene information
    filename = directory.joinpath(f"{identifier}_genes.xml")
    if database == "META":
        msg = (
            f'Object "{identifier}" retrieved from "META". Please use '
            "a sub-database to add proper genes. "
            "Skipping retrieval of gene information."
        )
        debug_log.warning(msg)

        return

    if not filename.exists():
        # This URL will not necessarily raise exception
        encoded_id = urllib.parse.quote(identifier, safe="")

        url_text = (
            f"https://websvc.biocyc.org/apixml?fn=genes-of-reaction&id="
            f"{database}:{encoded_id}&detail=full"
        )

        # FIXME: find a better way
        user, pwd = cmod_utils.get_credentials(
            Path.cwd().joinpath("credentials.txt")
        )
        s = requests.Session()
        s.post(
            "https://websvc.biocyc.org/credentials/login/",
            data={"email": user, "password": pwd},
        )
        response = s.get(url_text)
        s.close()
        try:
            response.raise_for_status()
            root = et.fromstring(response.text)
            # Check for results in root
            tree = et.ElementTree(root)

            element = tree.find("*/num_results")

            if not isinstance(element, et.Element):
                raise AttributeError("No gene information available")

            # Ignore metabolites
            if element.text == "0":
                return

            tree.write(str(filename))
            debug_log.info(
                f'Object "{identifier}_gene.xml" saved in '
                f'directory "{database}/GENES".'
            )
        except requests.HTTPError:
            raise AttributeError(
                f"Object {identifier} does not have gene information"
            )


def parse_genes(identifier: str, directory: Path) -> dict:
    """
    Returns a dictionary with the corresponding genes and gene-reaction rules.
    This function will try to read a file for given identifier. If nothing is
    found, the dictionary will have empty entries.
    """
    rule = str()
    genes = dict()
    # Get the information and check if Genes can be found to be parsed
    with suppress(FileNotFoundError):
        tree = et.parse(directory.joinpath(f"{identifier}_genes.xml")).getroot()

        if not isinstance(tree, et.Element):
            raise TypeError("Given root is not a valid Element object")

        for gene in tree.findall("Gene"):
            name = gene.attrib.get("frameid")

            if name is None:
                raise TypeError("Given root is not a valid Element object")

            genes[name] = name

        # Assuming rule
        rule = " or ".join(genes.keys())
        debug_log.warning(
            f'Gene-reaction rule for reaction "{identifier}" set'
            + ' to "OR". Please modify it if necessary.'
        )
    return {"genes": genes, "rule": rule}
