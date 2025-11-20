
from abc import ABC, abstractmethod
from typing import Union

from typing_extensions import Tuple

from cobramod.dbwalker.dataclasses import GenerellIdentifiers


class Database(ABC):
    """
    A Class defining the interface for interacting with databsed to aquire databse specific identifiers, SMILES, InChI and InChIKey.
    """

    @abstractmethod
    def getGenerellIdentifier(self, dbIdentifier:str, **kwargs) -> GenerellIdentifiers:
        """
        This methods queries SMILES, InChI and InChIKey from a database. If any is not available it will be set to None.
        Returns: A GenerellIdentifiers object containing containing SMILES, InChI and InChIKey if available.
        """

        ...

    @abstractmethod
    def getDBIdentifierFromSmiles(self, smiles: Union[str,GenerellIdentifiers]) -> str:
        """
        This method queries the database identifier from a given SMILES string.
        """

    @abstractmethod
    def getDBIdentifierFromInchi(self, smiles: Union[str,GenerellIdentifiers]) -> str:
        """
        This method queries the database identifier from a given Inchi.
        """

        ...

    @abstractmethod
    def getDBIdentifierFromInchiKey(self, smiles: Union[str,GenerellIdentifiers]) -> str:
        """
        This method queries the database identifier from a given InchIKey.
        """

        ...

    @abstractmethod
    def validateGenerellIdentifiers(self, smiles: Union[str,GenerellIdentifiers]) -> Tuple[GenerellIdentifiers,GenerellIdentifiers]:
        """
        This method queries available generell identifiers based on a databse identifier and validates them via checking if a generell identifier maps back to the original database identifier. If this is the case it will be added to the GenerellIdentifier object containing validated identifiers. If not it will be added to the GenerellIdentifiers object containing non verified identifiers..
        """
        ...