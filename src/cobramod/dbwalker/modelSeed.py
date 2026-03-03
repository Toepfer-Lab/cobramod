import logging
from typing import Union, Optional

from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import GenerellIdentifiers

logger = logging.getLogger("cobramod.DBWalker.ModelSeed")
logger.propagate = True


class ModelSeed(Database):
    def __init__(self):
        super().__init__()

        self.__getStructureFile()
        self.__load_structure_file()

    def __getCompoundFlatFile(self):
        pass

    def getGenerellIdentifier(self, dbIdentifier: str) -> GenerellIdentifiers:
        pass

    @property
    def name(self) -> str:
        return "ModelSeed"

    def save_cache(self):
        pass

    @property
    def AnnotationPrefix(self) -> str:
        return "seed.compound"

    def getDBIdentifierFromSmiles(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> Optional[str]:
        pass

    def getDBIdentifierFromInchi(
        self, inchi: Union[str, GenerellIdentifiers]
    ) -> Optional[str]:
        pass

    def getDBIdentifierFromInchiKey(
        self, inchikey: Union[str, GenerellIdentifiers]
    ) -> str:
        pass
