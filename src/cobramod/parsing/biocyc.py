#!/usr/bin/env python3
"""Data parsing for Biocyc

This module handles the retrieval of data from BioCyc into a local directory.
The posible type of data that can be download:

- Metabolites: Normally have an abbreviation or short name.
- Reactions: Can have the words "RXN" in the identifier. Enzymes can sometimes
be used instead. The gene information for the reactions is include if found.
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
from cobramod.error import WrongParserError, NoGeneInformation
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


def _p_genes(root: Any, identifier: str, directory: Path):
    """
    Parses given root and returns a dictionary with the corresponding genes and
    gene-reaction rules. This function will try to read a file for given
    identifier. If nothing is found, the dictionary will have empty entries.
    """
    rule = str()
    genes = dict()
    # Get the information and check if Genes can be found to be parsed
    with suppress(FileNotFoundError):
        tree: Any = BiocycParser._read_file(
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
        # FIXME: Temporal OR rule
        rule = " or ".join(genes.keys())
    return {"genes": genes, "rule": rule}


def _p_reaction(root: Any, directory: Path) -> dict:
    """
    Parses the root file and return the data in a dictionary
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
        "EQUATION": _p_metabolites(root=root),
        "BOUNDS": _check_direction(root=root),
        "TRANSPORT": BaseParser._check_transport(
            data_dict=_p_metabolites(root=root)
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
    def _parse(root: Any, directory: Path) -> dict:
        """
        Parses the root and returns a dictionary with the dictionary of the
        root, with the most important information depending on its object type.
        """
        # try to define dictionary
        for method in (_p_reaction, _p_pathway, _p_compound):
            with suppress(WrongParserError, AttributeError):
                try:
                    biocyc_dict = method(  # type: ignore
                        root=root, directory=directory
                    )
                except TypeError:
                    biocyc_dict = method(root=root)  # type: ignore
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
        return BiocycParser._parse(root=root, directory=directory)

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


def _get_gene_xml(directory: Path, identifier: str, database: str):
    # GENES directory will depend from the subdatabase
    directory = directory.joinpath("GENES")
    if not directory.exists():
        directory.mkdir()
    # Retrieval of the Gene information
    filename = directory.joinpath(f"{identifier}_genes.xml")
    if not filename.exists():
        # This URL will not necessarily raise exceptio
        url_text = (
            f"https://websvc.biocyc.org/apixml?fn=genes-of-reaction&id="
            f"{database}:{identifier}&detail=full"
        )
        r = get(url_text)
        try:
            r.raise_for_status()
            root = fromstring(r.text)
            # Check for results in root
            tree: Any = ElementTree(root)
            if int(tree.find("*/num_results").text) == 0:
                raise HTTPError
            tree.write(str(filename))
            debug_log.info(
                f'Object "{identifier}_gene.xml" saved in '
                f'directory "{database}/GENES".'
            )
        except HTTPError:
            # Warning
            NoGeneInformation(
                f"Object {identifier} does not have gene information"
            )


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
        # Search for the file on directory. Otherwise retrive from database
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
            try:
                r.raise_for_status()
                # defining root
                root = fromstring(r.text)
                tree = ElementTree(root)
                tree.write(str(filename))
                debug_log.info(
                    f'Object "{identifier}" found in database. Saving in '
                    f'directory "{database}".'
                )
                # Obtain genes if possible. This should be only call one time
                # If information is already available, the genes should be
                # available if found
                with suppress(HTTPError):
                    # This will include Paths
                    _get_gene_xml(
                        directory=data_dir,
                        identifier=identifier,
                        database=database,
                    )
                return root
            except HTTPError:
                msg = f'"{identifier}" is not available in "{database}"'
                debug_log.error(msg)
                # Try with META
                if database != "META":
                    # This will raise an error if not found in META
                    with suppress(HTTPError):
                        root = _get_xml_from_biocyc(
                            directory=directory,
                            identifier=identifier,
                            database="META",
                        )
                        return root
                # In case there is nothing in META
                raise HTTPError(msg)
    else:
        msg = "Directory not found. Please create the given directory."
        debug_log.critical(msg)
        raise NotADirectoryError(msg)
