"""Data parsing for KEGG

This module handles the retrieval of data from KEGG into a local directory.
The possible type of data that can be downloaded:

- Metabolite: Identifiers that start with the letter C, e.g C00001
- Reactions: Identifiers that start with the letter R, e.g R00001. The
gene information for reactions is also included if the specie is specified
- Module Pathways: Identifiers that start with the letter M, e.g M00001
"""

from __future__ import annotations

from contextlib import suppress
from itertools import chain
from pathlib import Path
from typing import Any, Generator, Optional

import requests

import cobramod.utils as cmod_utils
from cobramod.debug import debug_log
from cobramod.error import WrongParserError

KO_LINK = "http://rest.kegg.jp/link/ko/"


def parse_metabolite_attributes(
    data: dict[str, list[str]], entry: str
) -> dict[str, Any]:
    formula = data.get("FORMULA", ["X"]).pop()

    # NOTE: replace e.g (HPO3)n

    if formula.find("("):
        formula = formula.replace("(", "")
        formula = formula.replace(")", "")
        formula = formula.replace("n", "")

    attributes: dict[str, Any] = {
        "name": data.get("NAME", [entry]).pop(),
        "formula": formula,
        # NOTE: BiGG compounds does not show charges at a specific
        # ph-Value
        "charge": 0,
        "xref": build_references(data),
    }
    return attributes


def parse_reaction_attributes(
    data: dict[str, list[str]], entry: str, genome: str, gene_directory: Path
) -> dict[str, Any]:
    reaction_str = data.get("EQUATION", []).pop()

    if not reaction_str:
        raise AttributeError(f"There is no reaction found for '{entry}'")

    position = cmod_utils.get_arrow_position(reaction_str)

    if position == -1:
        raise AttributeError(
            f"Reaction string for {entry} did not include an arrow in the "
            "equation"
        )
    reactants = [
        item.rstrip().strip()
        for item in reaction_str[:position].rstrip().strip().split("+")
    ]
    products = [
        item.rstrip().strip()
        for item in reaction_str[position + 3 :].rstrip().strip().split("+")
    ]

    arrow = reaction_str[position : position + 3]

    reaction_str = str()

    for side in (reactants, products):
        prefix = "r"

        if side is reactants:
            prefix = "l"

        for i, item in enumerate(side):
            try:
                coefficient, identifier = [
                    item.strip().rstrip() for item in item.split(" ")
                ]
            except ValueError:
                # In case no coefficient is given, it must be 1
                coefficient = "1"
                identifier = item

            reaction_str += f"{coefficient} {prefix}_{identifier}"

            if i < len(side) - 1:
                reaction_str += " + "

            if side is reactants and i == len(side) - 1:
                reaction_str += f" {arrow} "

    attributes = {
        "name": data.get("NAME", [entry]).pop(),
        "equation": reaction_str,
        "xref": build_references(data),
        "replacements": {},
        "transport": cmod_utils.is_transport(reaction_str),
        "genes": parse_genes(
            identifier=entry,
            directory=gene_directory,
            genome=genome,
        ),
    }
    return attributes


def data_from_string(raw: str) -> dict[str, list[str]]:
    """
    Formats most of the keys for KEGG data and returns a dictionary.
    """
    lines = [line for line in raw.split(sep="\n")]
    data: dict[str, list[str]] = dict()

    for line in lines:
        segment = line.split(" ")

        if segment[0] != "///" and segment[0]:
            data[segment[0]] = []

    actual_key: str = ""

    for line in lines:
        if "///" in line:
            break

        # Create a list for each new key
        key = None
        for entry in data.keys():
            if entry in line:
                key = entry
                break

        if not key:
            data[actual_key].append(line.strip().rstrip())
        else:
            actual_key = key
            line = line[line.find(key) + len(key) :].strip().rstrip()
            data[actual_key] = [line]
    return data


def build_references(data_dict: dict[str, list[str]]) -> dict[str, str]:
    """
    Return a dictionary, where the keys are the names of cross-references
    and keys their identifiers. If nothing is found, it will return None
    """
    references = dict()

    for line in data_dict.get("DBLINKS", []):
        items = [item.rstrip().strip() for item in line.split(":")]
        references[items[0]] = items[1]

    return references


def ko_generator(identifier: str) -> Generator[str, None, None]:
    """
    Returns a list with the corresponding KO-entries for given identifier.
    Otherwise it will raise a HTTPError
    """
    url_text = f"http://rest.kegg.jp/link/ko/{identifier}"
    try:
        response = requests.get(url_text)
        response.raise_for_status()
        return (
            single
            for single in response.text.split()
            # This statement removes identifier from the split. e.g
            # [rn:R00001, ko:k00001]
            if identifier not in single
        )
    except requests.HTTPError:
        raise requests.HTTPError(
            f'No K-Identifiers could be retrieved for "{identifier}".'
        )


def parse_ko_to_genes(
    string: str, reaction: str, genome: Optional[str]
) -> list[str]:
    """
    Returns with a list with the corresponding genes for given genome. String
    is the raw text that include the gene information. If Abbreviation or no
    genome is present then, an empty list is returned.
    """
    # Retrive the KO-identifiers
    text = string.replace("\t", " ").splitlines()

    genes: dict[str, list] = dict()
    # Separate the lines and append the gene identifier to the corresponding
    # genome
    for line in text:
        _, single = line.split()
        specie, identifier = single.split(":")
        try:
            genes[specie].append(identifier)
        except KeyError:
            genes[specie] = [identifier]
    # Return the corresponding gene or empty list
    if genome:
        try:
            return genes[genome]
        except KeyError:
            msg = (
                f'Reaction "{reaction}" does not have a "{genome}" '
                + "abbreviation as a specie. No genes will be added."
            )
            debug_log.warning(msg=msg)
            return []
    msg = (
        f'Nothing was specified in argument "genome". Reaction "{reaction}"'
        " will not include genes. Please modify if necessary."
    )
    debug_log.warning(msg=msg)
    return []


def retrieve_kegg_genes(directory: Path, identifier: str):
    """
    Stores the genes for given reaction in given directory. Function will call
    a HTTPError if nothing is found.
    """
    # Ignore pathways and compounds
    if identifier.startswith("M") or identifier.startswith("C"):
        return

    directory = directory.joinpath("KEGG", "GENES")

    if not directory.exists():
        directory.mkdir()

    # Retrieval of the Gene information
    filename = directory.joinpath(f"{identifier}_genes.txt")

    if not filename.exists():
        try:
            string = "+".join(ko_generator(identifier))
            url_text = f"http://rest.kegg.jp/link/genes/{string}"
            response = requests.get(url_text)
            response.raise_for_status()

            with open(file=filename, mode="w") as file:
                file.write(response.text)

        except requests.ConnectionError:
            # NOTE: Using Generator otherwise, server breaks connection
            generator = ko_generator(identifier)
            while True:
                try:
                    ko = next(generator)
                    url_text = f"http://rest.kegg.jp/link/genes/{ko}"
                    response = requests.get(url_text)
                    response.raise_for_status()

                    with open(file=filename, mode="a+") as file:
                        file.write(response.text)
                except requests.HTTPError:
                    raise requests.HTTPError(
                        f'Gene information for "{identifier}" unavailable'
                    )
                except StopIteration:
                    break


def parse_genes(
    directory: Path, identifier: str, genome: str
) -> dict[str, Any]:
    """
    From given KEGG dictionary returns a dictionary with the key
    "genes" which include a dictionary with the identifier and name of the
    gene; and the key "rule" for the COBRApy representation of the
    gene-reaction rule
    """
    genes: dict[str, str] = dict()
    rule = str()

    filename = directory.joinpath(f"{identifier}_genes.txt")
    with suppress(FileNotFoundError):
        with open(file=filename, mode="r") as file:
            genes_list = parse_ko_to_genes(
                string=str(file.read()), reaction=identifier, genome=genome
            )
            for gene in genes_list:
                genes[gene] = ""

            rule = " or ".join(genes.keys())

    if genes:
        rule = " or ".join(genes.keys())
        debug_log.warning(
            f'Gene-reaction rule for reaction "{identifier}" from KEGG '
            'set to "OR".'
        )
    return {"genes": genes, "rule": rule}


def get_graph(kegg_dict: dict) -> dict:
    """
    Returns dictionary with sequences for a graph, where the key the prior
    reaction is and the value, its successor; and a set with all participant
    reaction (vertex)
    """
    # Obtain sequences as a set and regular reactions
    sequences = [item.split(sep=" ")[0] for item in kegg_dict["REACTION"]]
    whole_set = set(
        chain.from_iterable(line.split(sep=",") for line in sequences)
    )
    reactions = [line.split(sep=",") for line in sequences]
    graph: dict[str, Any] = dict()

    line: list
    for index_j, line in enumerate(iterable=reactions):
        reaction: str
        for reaction in line:
            # If index not found, then ignore it
            with suppress(IndexError):
                # In KEGG the order of the reactions is defined by its order in
                # the sequence.
                parent = reaction
                children = reactions[index_j + 1]
                if len(children) == 1:
                    # Take as a string
                    children = children[0]
                else:
                    # All elements as a tuple
                    children = tuple(children)
                # Expand keys
                try:
                    # If at least 2
                    if isinstance(graph[parent], tuple):
                        graph[parent] += children
                        continue
                    # Else transform single string to tuple
                    graph[parent] = (graph[parent], children)
                except KeyError:
                    # Create new one
                    graph[parent] = children
    # Check that the graph includes all reactions
    for reaction in whole_set:
        try:
            graph[reaction]
        except KeyError:
            graph[reaction] = None
    if not graph:
        raise WrongParserError("Given root does not belong to a Pathway")
    return graph


def parse_pathway_attributes(
    data: dict[str, list[str]], entry: str
) -> dict[str, Any]:
    graph = get_graph(data)
    name = data.get("NAME", [entry])

    attributes = {
        "name": name.pop(),
        "xref": build_references(data),
        "pathway": graph,
    }
    return attributes
