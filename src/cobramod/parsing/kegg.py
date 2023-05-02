"""Data parsing for KEGG

This module handles the retrieval of data from KEGG into a local directory.
The possible type of data that can be downloaded:

- Metabolite: Identifiers that start with the letter C, e.g C00001
- Reactions: Identifiers that start with the letter R, e.g R00001. The
gene information for reactions is also included if the specie is specified
- Module Pathways: Identifiers that start with the letter M, e.g M00001

Contact maintainers if other types should be added.

Important class of the module:
- KeggParser: Child of the abstract class
:class:`cobramod.parsing.base.BaseParser`.
"""
from __future__ import annotations

import io
from contextlib import suppress
from itertools import chain
from pathlib import Path
from typing import Any, Generator, Optional, Union
from warnings import warn

from requests import HTTPError, get

import cobramod.utils as cmod_utils
from cobramod.debug import debug_log
from cobramod.error import AbbreviationWarning, WrongParserError

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
    data: dict[str, list[str]], entry: str
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
            directory=Path(),
            genome="",
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


def _kegg_info_to_version(info: str) -> str:
    with io.StringIO(info) as lines:
        for num, line in enumerate(lines, 1):
            if num == 2:
                return line[line.find("Release") + 8 :].rstrip()

    msg = (
        "Error determining the kegg version. "
        'Instead, "Undefined" is used as version.'
    )
    warn(
        message=msg,
        category=UserWarning,
    )
    return "Undefined"


def _get_keys(raw: str) -> Generator:
    """
    Returns as Generator the keys for KEGG data. Input is raw text
    """
    lines = (line for line in raw.split(sep="\n"))
    for line in lines:
        segment = line.split(" ")
        if segment[0] != "///" and segment[0]:
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
    Formats most of the keys for KEGG data and returns a dictionary.
    """
    lines = (line for line in raw.split("\n"))
    keys = set(_get_keys(raw=raw))
    actual_key = str()
    kegg_dict: dict[str, list] = {}
    # Create dictionary
    for line in lines:
        if "///" in line:
            break
        key = _find_key(line=line, keys=keys)
        if not key:
            # Append to older key.
            kegg_dict[actual_key].append(line.strip().rstrip())
        else:
            actual_key = key
            line = line[line.find(key) + len(key) :].strip().rstrip()
            kegg_dict[actual_key] = [line]
    return kegg_dict


# def _get_metabolites(line: str) -> tuple:
#     """
#     Returns the bounds of the reaction as a tuple and a dictionary with the
#     participant metabolites of the reaction
#
#     Args:
#         line (str): Equation in string form
#
#     Returns:
#         Tuple: bounds (tuple), dictionary with metabolites
#
#     Raises:
#         UnboundLocalError: If reversibility cannot be found for given reaction
#     """
#     # Obtaining bounds
#     arrows = {"<=>": (-1000, 1000), "==>": (0, 1000), "<==": (-1000, 1000)}
#     for item in arrows.keys():
#         if item in line:
#             position = line.find(item)
#             arrow = line[position : position + 3]
#     try:
#         arrow
#     except UnboundLocalError:
#         raise UnboundLocalError(
#             f'Reversibility for "{line}" could not be found.'
#         )
#     # Getting metabolites
#     metabolites = dict()
#     reactants, products = [
#         side.rstrip().strip().split("+") for side in line.split(arrow)
#     ]
#     for side in (reactants, products):
#         FACTOR = 1
#         prefix = "r"
#         if side is reactants:
#             FACTOR = -1
#             prefix = "l"
#         for metabolite in side:
#             try:
#                 coefficient, identifier = [
#                     item.strip().rstrip()
#                     for item in metabolite.rstrip().strip().split(" ")
#                 ]
#             except ValueError:
#                 # In case no coefficient is given, it must be 1
#                 identifier = metabolite.strip().rstrip()
#                 coefficient = "1"
#             metabolites[f"{prefix}_{identifier}"] = float(coefficient) * FACTOR
#     return arrows[arrow], metabolites


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


def _get_ko(identifier: str) -> list:
    """
    Returns a list with the corresponding KO-entries for given identifier.
    Otherwise it will raise a HTTPError
    """
    url_text = f"http://rest.kegg.jp/link/ko/{identifier}"
    try:
        response = get(url_text)
        response.raise_for_status()
        return [
            single
            for single in response.text.split()
            # This statement removes identifier from the split. e.g
            # [rn:R00001, ko:k00001]
            if identifier not in single
        ]
    except HTTPError:
        raise HTTPError(
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
        ko, single = line.split()
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
            warn(message=msg, category=AbbreviationWarning)
            debug_log.warning(msg=msg)
            return []
    msg = (
        f'Nothing was specified in argument "genome". Reaction "{reaction}"'
        " will not include genes. Please modify if necessary."
    )
    warn(message=msg, category=UserWarning)
    debug_log.warning(msg=msg)
    return []


def retrieve_kegg_genes(directory: Path, identifier: str):
    """
    Stores the genes for given reaction in given directory. Function will call
    a HTTPError if nothing is found.
    """
    # GENES directory will depend from the subdatabase
    directory = directory.joinpath("KEGG").joinpath("GENES")
    if not directory.exists():
        directory.mkdir()
    # Retrieval of the Gene information
    filename = directory.joinpath(f"{identifier}_genes.txt")
    if not filename.exists():
        # Obtain list of KO
        # This URL will not necessarily raise exception
        string = "+".join(_get_ko(identifier=identifier))
        url_text = f"http://rest.kegg.jp/link/genes/{string}"
        try:
            response = get(url_text)
            response.raise_for_status()
            with open(file=filename, mode="w+") as file:
                file.write(response.text)
        except HTTPError:
            raise HTTPError(f'Gene information for "{identifier}" unavailable')


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

    filename = (
        directory.joinpath("KEGG")
        .joinpath("GENES")
        .joinpath(f"{identifier}_genes.txt")
    )
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
            + ' set to "OR".'
        )
    return {"genes": genes, "rule": rule}


# def _p_reaction(kegg_dict: dict, directory: Path, genome: str) -> dict:
#     """
#     Returns a dictionary that is the representation of a reaction
#     """
#     try:
#         object_type = list(
#             chain.from_iterable(item.split() for item in kegg_dict["ENTRY"])
#         )[-1]
#         if object_type != "Reaction":
#             raise WrongParserError(
#                 "Given dictionary does not belong to a reaction."
#             )
#         identifier = list(
#             chain.from_iterable(item.split() for item in kegg_dict["ENTRY"])
#         )[-2]
#         # Obtain bounds
#         bounds, metabolites = _get_metabolites(line=kegg_dict["EQUATION"][0])
#         with suppress(KeyError):
#             if "irreversible" in ",".join(kegg_dict["COMMENT"]):
#                 bounds = (0, 1000)
#         return {
#             "TYPE": object_type,
#             "ENTRY": identifier,
#             "NAME": kegg_dict["NAME"][0],
#             "DATABASE": "KEGG",
#             "EQUATION": metabolites,
#             "BOUNDS": bounds,
#             "TRANSPORT": BaseParser._check_transport(data_dict=metabolites),
#             "XREF": build_references(data_dict=kegg_dict),
#             "GENES": parse_genes(
#                 kegg_dict=kegg_dict,
#                 identifier=identifier,
#                 directory=directory,
#                 genome=genome,
#             ),
#         }
#     except KeyError:
#         raise WrongParserError
#
#
# def _p_compound(kegg_dict: dict) -> dict:
#     """
#     Returns a dictionary that is the representation of a metabolite
#     """
#     try:
#         object_type = list(
#             chain.from_iterable(item.split() for item in kegg_dict["ENTRY"])
#         )[-1]
#         if object_type != "Compound":
#             raise WrongParserError(
#                 "Given dictionary does not belong to a metabolite."
#             )
#         identifier = list(
#             chain.from_iterable(item.split() for item in kegg_dict["ENTRY"])
#         )[-2]
#         try:
#             # When formulas includes 'n' e.g 'H4P2O7(HPO3)n'
#             if "(" in kegg_dict["FORMULA"][0]:
#                 formula = "X"
#                 msg = (
#                     f'Sum formula for metabolite "{identifier}" from KEGG '
#                     'could not be found. Formula set to "X" and charge to 0. '
#                     "Please modify it if necessary."
#                 )
#                 debug_log.warning(msg=msg)
#                 warn(msg)
#             else:
#                 formula = "".join(kegg_dict["FORMULA"]).strip().rstrip()
#         except KeyError:
#             # Some general compounds do not come with formulas
#             formula = "X"
#             msg = (
#                 f'Sum formula for metabolite "{identifier}" from KEGG '
#                 'could not be found. Formula set to "X" and charge to 0. '
#                 "Please modify it if necessary."
#             )
#             debug_log.warning(msg=msg)
#             warn(msg)
#         return {
#             "TYPE": list(
#                 chain.from_iterable(item.split() for item in kegg_dict["ENTRY"])
#             )[-1],
#             "ENTRY": identifier,
#             "NAME": kegg_dict["NAME"][0],
#             "DATABASE": "KEGG",
#             "FORMULA": formula,
#             # FIXME: search information about charges
#             "CHARGE": 0,
#             "XREF": build_references(data_dict=kegg_dict),
#         }
#     except KeyError:
#         raise WrongParserError


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


def _p_pathway(kegg_dict: dict) -> dict:
    """
    Returns a dictionary, which is the representation of a pathway.
    """
    try:
        if "Pathway" not in kegg_dict["ENTRY"][0]:
            raise WrongParserError(
                "Given dictionary does not belong to a pathway."
            )
        graph = get_graph(kegg_dict=kegg_dict)
        return {
            "TYPE": "Pathway",
            "ENTRY": list(
                chain.from_iterable(item.split() for item in kegg_dict["ENTRY"])
            )[0],
            "NAME": kegg_dict["NAME"][0],
            "DATABASE": "KEGG",
            "PATHWAY": graph,
            "XREF": build_references(data_dict=kegg_dict),
        }
    except KeyError:
        raise WrongParserError


def _p_enzyme(kegg_dict: dict):
    """
    If enzyme information is found, raise NotImplementedError.
    """
    try:
        if "Enzyme" in kegg_dict["ENTRY"][0]:
            raise NotImplementedError(
                "Enzymes are currently not implemented. Please contact "
                + "maintainers."
            )
        else:
            raise WrongParserError(
                "Given dictionary does not belong to an Enzyme."
            )
    except KeyError:
        raise WrongParserError
