import logging

from typing import Optional, Union, Tuple

import requests

from cobramod import Settings
from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import GenerellIdentifiers

logger = logging.getLogger("cobramod.DBWalker.Chebi")
logger.propagate = True

# Ensure console output
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


class Chebi(Database):

    chebi_ftp = "ftp://ftp.ebi.ac.uk/pub/databases/chebi/flat_files"
    settings = Settings()

    def __init__(self):
        super().__init__()

        self.cached_structure = self.__getStructureFile()

    def __getStructureFile(self):
        """
        Check if a local structure file is available and of not get it from chebis ftp server.

        """

        logger.debug("Checking if local structure file")

        # structure file
        url = self.chebi_ftp + "/chebi-structure.tsv"

        cached_file = self.settings.cacheDir / "chebi" / "chebi-structure.tsv"

        if cached_file.exists():
            logger.debug("Found local structure file")
            return cached_file

        logger.debug("ChEBI structure file not found. Downloading...")

        response = requests.get(url)

        if response.status_code == 200:
            with open(cached_file, 'wb') as file:
                file.write(response.content)
            logger.debug(
                f"ChEBI structure file donwloaded and saved at {str(cached_file)}."
            )
        else:
            logger.debug(
                f"Download failed with status code {response.status_code}"
            )
            response.raise_for_status()

        return cached_file

    def getDBIdentifier(self, identifier: GenerellIdentifiers) -> Optional[str]:
        pass

    def getGenerellIdentifier(self, dbIdentifier: str) -> GenerellIdentifiers:
        raise NotImplementedError

    def getDBIdentifierFromSmiles(self, smiles: Union[str, GenerellIdentifiers]) -> Optional[str]:
        raise NotImplementedError

    def getDBIdentifierFromInchi(self, inchi: Union[str, GenerellIdentifiers]) -> Optional[str]:
        raise NotImplementedError

    def getDBIdentifierFromInchiKey(self, inchikey: Union[str, GenerellIdentifiers]) -> str:
        raise NotImplementedError

    def validateGenerellIdentifiers(self, smiles: Union[str, GenerellIdentifiers]) -> Tuple[
        GenerellIdentifiers, GenerellIdentifiers]:
        raise NotImplementedError

