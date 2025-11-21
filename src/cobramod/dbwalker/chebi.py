import logging

from typing import Optional, Union, Tuple

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

