import requests
import logging
from typing import Optional

from cobramod.dbwalker.dataclasses import GenerellIdentifiers

logger = logging.getLogger("cobramod.DBWalker.SwissLipids")
logger.propagate = True

# Ensure console output
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def get_swisslipids_id(
    identifier: str, identifier_type: str = "auto"
) -> Optional[str]:
    """
    Get SwissLipids ID from InChI Key or SMILES string.

    Args:
        identifier: InChI Key or SMILES string
        identifier_type: Type of identifier ("inchikey" or "smiles")

    Returns:
        SwissLipids ID if found, None otherwise
    """
    base_url = "https://www.swisslipids.org/api/index.php/"
    try:
        identifier_type = identifier_type.lower()

        if identifier_type == "inchikey":
            url = base_url + f"advancedSearch?InChIkey={identifier}"
        elif identifier_type == "smiles":
            # seems to only work for simple SMILES strings that do not
            # need encoding
            # Example for SMILES also does not work (06.08.2025):
            # https://www.swisslipids.org/#/advanced
            escapedID = identifier.replace("(", "\(")
            escapedID = escapedID.replace(")", "\)")
            escapedID = escapedID.replace("=", "%3D")

            url = base_url + f"advancedSearch?SMILES={escapedID}"
        else:
            logger.error(f"Unsupported identifier type: {identifier_type}")
            return None

        logger.info(
            f"Requesting SwissLipids ID for {identifier_type}: {identifier}"
        )

        response = requests.get(url, timeout=30)

        if (
            response.status_code == 404
            and response.json()["message"] == "Sorry, no result found"
        ):
            logger.error(
                f"No result available in database for: {identifier_type}={identifier}"
            )
            return None

        print(response.status_code)
        print(response.json())
        response.raise_for_status()

        data = response.json()

        # Extract SwissLipids ID from response
        if isinstance(data, list) and len(data) > 0:
            # For search results, take the first match
            logger.info(
                f"Found {len(data)} results for {identifier_type}: {identifier}"
            )
            return data[0].get("entity_id")
        elif isinstance(data, dict):
            # For direct lookups
            swissID = data.get("entity_id")
            logger.info(f"Found {swissID} for {identifier_type}:{identifier}")
            return swissID

        return None

    except (requests.RequestException, KeyError, ValueError) as e:
        logger.error(f"Error fetching SwissLipids ID: {e}")
        return None


def get_swisslipids_data(swisslipids_id: str) -> GenerellIdentifiers:
    """
    Get chemical identifiers from SwissLipids ID.

    Args:
        swisslipids_id: SwissLipids database identifier

    Returns:
        GenerellIdentifiers containing inchi, inchi_key, and smiles if available
    """
    try:
        # SwissLipids API endpoint for entity details
        url = (
            f"https://www.swisslipids.org/api/index.php/entity/{swisslipids_id}"
        )

        logger.info(
            f"Requesting data from SwissLipids for ID: {swisslipids_id}"
        )

        response = requests.get(url, timeout=30)

        if (
            response.status_code == 500
            and response.json()["message"] == "ERROR: entity is unknown"
        ):
            logger.warning(f"No entry in SwissLipids for: {swisslipids_id}")
            return GenerellIdentifiers()

        response.raise_for_status()

        data = response.json()

        # Extract chemical identifiers from response
        inchi = data["structures"].get("inchi")
        if inchi is not None and inchi.startswith("InChI="):
            inchi = inchi[6:]
            inchi = inchi.replace("\\/", "/")

        inchi_key = data["structures"].get("inchikey")
        smiles = data["structures"].get("smiles")

        return GenerellIdentifiers(
            inchi=inchi, inchi_key=inchi_key, smiles=smiles
        )

    except (requests.RequestException, KeyError, ValueError) as e:
        logger.error(f"Error fetching data from SwissLipids: {e}")
        return GenerellIdentifiers()


def swisslipids_search_by_name(name: str) -> Optional[str]:
    """
    Search SwissLipids by compound name and return the first matching ID.

    Args:
        name: Compound name to search for

    Returns:
        SwissLipids ID if found, None otherwise
    """
    try:
        url = "https://www.swisslipids.org/api/entity/search"
        params = {"query": name, "format": "json"}

        logger.info(f"Searching SwissLipids for name: {name}")

        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()

        data = response.json()

        if isinstance(data, list) and len(data) > 0:
            return data[0].get("id") or data[0].get("swisslipids_id")

        return None

    except (requests.RequestException, KeyError, ValueError) as e:
        logger.error(f"Error searching SwissLipids by name: {e}")
        return None
