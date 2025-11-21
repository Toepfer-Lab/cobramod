from abc import ABC, abstractmethod
from typing import Union, Optional

from typing_extensions import Tuple

from cobramod.dbwalker.dataclasses import GenerellIdentifiers


class Database(ABC):
    """
    A Class defining the interface for interacting with databsed to aquire databse specific identifiers, SMILES, InChI and InChIKey.
    """

    @abstractmethod
    def getGenerellIdentifier(
        self, dbIdentifier: str
    ) -> GenerellIdentifiers:
        """
        This methods queries SMILES, InChI and InChIKey from a database. If any is not available it will be set to None.
        Returns: A GenerellIdentifiers object containing containing SMILES, InChI and InChIKey if available.
        """

        ...

    @abstractmethod
    def getDBIdentifier(self, identifier:GenerellIdentifiers) -> Optional[str]:
        """
        This methods queries the database identifier from a given GenerellIdentifiers.
        Therefore using all defined Generellidentifier and validating that they result in the same database identifier.
        """

        ...

    @abstractmethod
    def getDBIdentifierFromSmiles(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> Optional[str]:
        """
        This method queries the database identifier from a given SMILES string.
        """

        ...

    @abstractmethod
    def getDBIdentifierFromInchi(
        self, inchi: Union[str, GenerellIdentifiers]
    ) -> Optional[str]:
        """
        This method queries the database identifier from a given Inchi.
        """

        ...

    @abstractmethod
    def getDBIdentifierFromInchiKey(
        self, inchikey: Union[str, GenerellIdentifiers]
    ) -> str:
        """
        This method queries the database identifier from a given InchIKey.
        """

        ...
