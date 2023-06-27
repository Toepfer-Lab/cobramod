"""Data parsing for SolanaCyc

This module handles the retrieval of data from SolanaCyc into a local
directory. The possible type of data that can be downloaded:

- Metabolites: Normally have an abbreviation or short name.
- Reactions: Can have the words "RXN" in the identifier. Enzymes can sometimes
be used instead. The gene information for the reactions is included if found.
- Pathways

Contact maintainers if other types should be added.

Important class of the module:
- SolCycParser: Child of the abstract class
:class:`cobramod.parsing.base.BaseParser`.
"""
import urllib.parse
import xml.etree.ElementTree as et
from pathlib import Path
from typing import Any
from warnings import warn

import requests

from cobramod.debug import debug_log


def retrieve_gene_information(directory: Path, identifier: str, database: str):
    """
    Retrieve the gene information for given reaction in a specific database

    Args:
        directory (Path): Directory to store and retrieve local data.
        identifier (str): original identifier
        database (str): Name of the database. Some options: "META", "LYCO".
            Full list in: "https://solcyc.sgn.cornell.edu/xmlquery?dbs"

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
        warn(message=msg, category=UserWarning)

        return

    if not filename.exists():
        # This URL will not necessarily raise exception
        encoded_id = urllib.parse.quote(identifier, safe="")

        url_text = (
            f"https://solcyc.sgn.cornell.edu/apixml?fn=genes-of-reaction&id="
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
