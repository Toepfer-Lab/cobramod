import logging
from abc import ABC, abstractmethod
from typing import Union, Optional

from typing_extensions import Tuple

from cobramod.dbwalker.dataclasses import GenerellIdentifiers


class Database(ABC):
    """
    A Class defining the interface for interacting with databsed to aquire databse specific identifiers, SMILES, InChI and InChIKey.
    """
    logger = logging.getLogger("cobramod.DBWalker.DataBase")

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

    def _validateGeneralIdentifiersWithDBIDs(self, generelID:GenerellIdentifiers, identifier:str):
        """
        Internal method that maps the entries of a GenerellIdenfiers object to a database identifier.
         With this it validates that they are equal to the passed identifier. This method is internally used to validate
         that general identifers map back to the DBIdentifer and are therefore within this Database coherent.
         Missmatches will be logged and removed from the GenerellIdentifiers object.
        Returns:

        """

        if generelID.smiles is not None:
            smiles_databaseID = self.getDBIdentifierFromSmiles(generelID.smiles)
            if identifier != smiles_databaseID:
                self.logger.error(
                    f"The SMILES string ({generelID.smiles}) does not map back to the original BioCyc ID ({identifier}), instead it points to {smiles_databaseID}. Ignoring SMILES"
                )
                generelID.smiles = None

        if generelID.inchi is not None:
            inchi_databaseID = self.getDBIdentifierFromInchi(generelID.inchi)
            if identifier != inchi_databaseID:
                self.logger.error(
                    f"The InChi ({generelID.inchi}) does not map back to the original Biocyc ID ({identifier}), instead it points to {inchi_databaseID}. Ignoring the InChi"
                )

        if generelID.inchi_key is not None:
            inchikey_databaseID = self.getDBIdentifierFromInchiKey(generelID.inchi_key)
            if identifier != inchikey_databaseID:
                self.logger.error(
                    f"The InChiKey ({generelID.inchi_key}) does not map back to the original Biocyc ID ({identifier}), instead it points to {inchikey_databaseID}. Ignoring the InChi"
                )



    @abstractmethod
    def getDBIdentifierFromInchiKey(
        self, inchikey: Union[str, GenerellIdentifiers]
    ) -> str:
        """
        This method queries the database identifier from a given InchIKey.
        """

        ...
