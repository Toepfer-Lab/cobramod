"""Data parsing for Biocyc

This module handles the retrieval of data from BioCyc into a local directory.
The possible type of data that can be downloaded:

- Metabolites: Normally have an abbreviation or short name.
- Reactions: Can have the words "RXN" in the identifier. Enzymes can sometimes
be used instead. The gene information for the reactions is included if found.
- Pathways

Contact maintainers if other types should be added.

Important class of the module:
- BiocycParser: Child of the abstract class
:class:`cobramod.parsing.base.BaseParser`.
"""
from __future__ import annotations

import xml.etree.ElementTree as et
from contextlib import suppress
from pathlib import Path
from typing import Any
from warnings import warn

import cobramod.utils as cmod_utils
from cobramod.debug import debug_log
from cobramod.error import (
    SuperpathwayWarning,
    WrongParserError,
)


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


def parse_reaction_attributes(root: et.Element, entry: str) -> dict[str, Any]:
    try:
        name = root.find("*[@ID]/enzymatic-reaction/*/common-name")
        if name is not None and name.text:
            name = name.text
    except AttributeError:
        name = entry
    equation = reaction_str_from_xml(root)
    attributes = {
        "name": name,
        "equation": equation,
        "xref": build_cross_references_xml(root),
        "replacements": {},
        "transport": cmod_utils.is_transport(equation),
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
    graph: dict[str, Any] = dict()
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
        tree = et.parse(
            directory.joinpath("GENES").joinpath(f"{identifier}_genes.xml")
        ).getroot()

        if not isinstance(tree, et.Element):
            raise TypeError("Given root is not a valid Element object")

        for gene in tree.findall("Gene"):
            element = gene.find("common-name")

            if not isinstance(element, et.Element):
                raise TypeError("Given root is not a valid Element object")

            name = element.text

            if not name:
                name = identifier

            genes[gene.attrib["frameid"]] = name

        # Assuming rule
        rule = " or ".join(genes.keys())
        debug_log.warning(
            f'Gene-reaction rule for reaction "{identifier}" set'
            + ' to "OR". Please modify it if necessary.'
        )
    return {"genes": genes, "rule": rule}
