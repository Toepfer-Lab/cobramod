"""Data parsing for BiGG

This module handles the retrieval of data from BiGG into a local directory.
The possible type of data that can be downloaded:

- Metabolites: Normally, simple names.
- Reactions: Mostly abbreviations.
- Genes: It is included in the Reactions. Names and gene-reaction-rule is also
acquired

They change identifiers depending on the model given. BiGG have multiple models
"""

from __future__ import annotations

from contextlib import suppress
from typing import Any

import requests

import cobramod.utils as cmod_utils
from cobramod.debug import debug_log


def find_url(model_id: str, query: str) -> tuple[requests.Response, str]:
    """
    Tries to find a valid URL for the API of BIGG. It will return a
    :class:`requests.Response` if URL is valid.
    object

    Args:
        model_id (str): Name for the specific BIGG model identifier.
        identifier (str): Identifier for the item. Can be a reaction or
            compound

    Returns:
        (Response,str): A valid request Response followed by the database
            Version.

    Raises:
        HTTPError: If identifier is not found in BIGG database.
    """
    for object_type in ("reactions", "metabolites"):
        with suppress(requests.HTTPError):
            # Check that status is available
            debug_log.debug(
                f'{object_type.capitalize()[:-1]} "{query}" not found in '
                f'directory "BIGG", subdirectory "{model_id}".'
            )
            # FIXME: BIGG combines compartments in their internal identifier.
            # To overcome this, use universal
            # if object_type == "metabolites":
            #     model_id = "universal"

            # Retrieve from URL
            url_text = (
                f"http://bigg.ucsd.edu/api/v2/models/{model_id}/{object_type}"
                f"/{query}"
            )
            debug_log.debug(f"Searching {url_text} for biochemical data.")
            # Get and check for errors
            response = requests.get(url_text)
            response.raise_for_status()

            info_response = requests.get(
                "http://bigg.ucsd.edu/api/v2/database_version"
            )
            info_response.raise_for_status()
            db_version = info_response.json()["bigg_models_version"]

            return response, db_version
    # Otherwise
    raise requests.HTTPError(
        f"Identifier '{query}' not found in BIGG (model: '{model_id}')."
    )


def build_reference(json_data: dict[str, Any]) -> dict[str, str]:
    """
    Return a dictionary of cross-references, where the keys are the
    cross-references and values the their identifiers.
    """
    references = json_data["database_links"]
    try:
        return {item: references[item][0]["id"] for item in references}
    except KeyError:
        raise KeyError("Problem with given json dictionary. Inform maintainers")


def parse_reaction_attributes(
    data: dict[str, Any], entry: str
) -> dict[str, Any]:
    results = data.get("results", None)

    if results:
        databases = data.get("database_links", None)
        data = results.pop()
        data["database_links"] = databases

    reaction_str = data.get("reaction_string", "")

    reaction_str = reaction_str.replace("&#8652;", "<->")
    reaction_str = reaction_str.replace("&#x2192;", "-->")
    reaction_str = reaction_str.replace("&rarr;", "-->")
    reaction_str = reaction_str.replace("&larr;", "<--")

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

            reaction_str += f"{coefficient} {prefix}_{identifier[:-2]}"

            if i < len(side) - 1:
                reaction_str += " + "

            if side is reactants and i == len(side) - 1:
                reaction_str += f" {arrow} "

    genes = parse_genes(
        data.get("genes", []), data.get("gene_reaction_rule", "")
    )
    attributes: dict[str, Any] = {
        "name": data.get("name", data.get("exported_reaction_id", "")),
        "equation": reaction_str,
        "xref": build_reference(data),
        "replacements": {},
        "transport": cmod_utils.is_transport(reaction_str),
        "genes": genes,
    }
    return attributes


def parse_metabolite_attributes(data: dict[str, Any]) -> dict[str, Any]:
    # For generic compounds, a list comes, otherwise it is a regular
    # string
    try:
        formula: str = data.get("formulae", data.get("formula", [])).pop()
        charge = int(data["charges"].pop())

    except AttributeError:
        # NOTE: add message
        formula = data.get("formulae", data.get("formula", "X"))
        charge = int(data["charge"])

    attributes = {
        "name": data.get("name", ""),
        "formula": formula,
        "charge": charge,
        "xref": build_reference(data),
    }
    return attributes


def parse_genes(data: list[dict[str, Any]], rule: str):
    genes: dict[str, str] = {}
    for gene in data:
        genes[gene["bigg_id"]] = gene.get("name", "")

    return {"genes": genes, "rule": rule}
