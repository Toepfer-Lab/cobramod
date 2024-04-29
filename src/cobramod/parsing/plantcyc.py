"""Gene parsing for Plantcyc

The main functions to parse Metabolites and Reactions are found in the module
Biocyc
"""

import urllib.parse
import xml.etree.ElementTree as et
from pathlib import Path
from typing import Any

import requests

from cobramod.debug import debug_log


def retrieve_gene_information(directory: Path, identifier: str, database: str):
    """
    Retrieve the gene information for given reaction in a specific database

    Args:
        directory (Path): Directory to store and retrieve local data.
        identifier (str): original identifier
        database (str): Name of the database. Some options: "META", "ARA".
            Full list in: "https://plantcyc.org/list-of-pgdbs"

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
            f"https://pmn.plantcyc.org/apixml?fn=genes-of-reaction&id="
            f"{database}:{encoded_id}&detail=full"
        )
        response = requests.get(url_text)
        try:
            response.raise_for_status()
            root = et.fromstring(response.text)
            # Check for results in root
            tree: Any = et.ElementTree(root)

            element = tree.find("*/num_results")
            if not isinstance(element, et.Element):
                raise AttributeError("No gene information available")

            tree.write(str(filename))
            debug_log.info(
                f'Object "{identifier}_gene.xml" saved in '
                f'directory "{database}/GENES".'
            )
        except requests.HTTPError:
            raise AttributeError(
                f"Object {identifier} does not have gene information"
            )
