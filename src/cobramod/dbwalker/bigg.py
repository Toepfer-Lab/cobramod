import requests
import logging
from typing import Union, Tuple

from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import GenerellIdentifiers

logger = logging.getLogger("cobramod.DBWalker.Bigg")
logger.propagate = True

# Ensure console output
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


class Bigg(Database):
    def getDBIdentifierFromSmiles(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> str:
        raise NotImplementedError

    def getDBIdentifierFromInchi(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> str:
        raise NotImplementedError

    def getDBIdentifierFromInchiKey(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> str:
        raise NotImplementedError

    def validateGenerellIdentifiers(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> Tuple[GenerellIdentifiers, GenerellIdentifiers]:
        raise NotImplementedError

    def __init__(self):
        super().__init__()

    def getGenerellIdentifier(self, dbIdentifier: str, **kwargs) -> GenerellIdentifiers:
        """
        Convert a BiGG metabolite ID to available general identifiers.

        Args:
            bigg_id: BiGG database metabolite identifier

        Returns:
            Dictionary containing inchi, inchi_key, and smiles if available
        """
        foundIDs = GenerellIdentifiers(inchi=None, inchi_key=None, smiles=None)

        try:
            url = f"http://bigg.ucsd.edu/api/v2/universal/metabolites/{dbIdentifier}"

            logger.info(
                f"Requesting data from Bigg for metabolite with id {dbIdentifier}"
            )
            logger.debug(f"Request URL: {url}")

            response = requests.get(url)
            response.raise_for_status()
        except requests.RequestException as e:
            logger.error(f"Request to BiGG API failed: {e}")
            return foundIDs

        data = response.json()

        try:
            data = data["database_links"]
        except KeyError:
            logger.warning("No database links found in BiGG response")
            return foundIDs

        try:
            foundIDs.inchi = data["InChI"][0]["id"]
            logger.debug(f"Found InChI: {foundIDs.inchi}")
        except (KeyError, IndexError):
            logger.debug("No InChI found in BiGG response")

        try:
            foundIDs.inchi_key = data["InChI Key"][0]["id"]
            logger.debug(f"Found InChI Key: {foundIDs.inchi_key}")
        except (KeyError, IndexError):
            logger.debug("No InChI Key found in BiGG response")

        try:
            foundIDs.smiles = data["SMILES"][0]["id"]
            logger.debug(f"Found SMILES: {foundIDs.smiles}")
        except (KeyError, IndexError):
            logger.debug("No SMILES found in BiGG response")

        return foundIDs
