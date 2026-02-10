import logging
from abc import ABC, abstractmethod
from typing import Union, Optional

import numpy
import pandas as pd
from typing_extensions import Tuple

from cobramod import Settings
from cobramod.dbwalker.dataclasses import GenerellIdentifiers


class Database(ABC):
    """
    A Class defining the interface for interacting with databsed to aquire databse specific identifiers, SMILES, InChI and InChIKey.
    """

    logger = logging.getLogger("cobramod.DBWalker.DataBase")

    @abstractmethod
    def getGenerellIdentifier(self, dbIdentifier: str) -> GenerellIdentifiers:
        """
        This methods queries SMILES, InChI and InChIKey from a database. If any is not available it will be set to None.
        Returns: A GenerellIdentifiers object containing containing SMILES, InChI and InChIKey if available.
        """

        ...

    @property
    @abstractmethod
    def name(self) -> str: ...

    @abstractmethod
    def save_cache(self): ...

    @property
    @abstractmethod
    def AnnotationPrefix(self) -> str: ...

    def getDBIdentifier(self, identifier: GenerellIdentifiers) -> Optional[str]:
        """
        This methods queries the database identifier from a given GenerellIdentifiers.
        Therefore using all defined Generellidentifier and validating that they result in the same database identifier.
        """

        smilesBasedID = None
        inchiBasedID = None
        inchikeyBasedID = None

        self.logger.debug(
            "Querying all available GenerellIdentifier, to check whether they point to the same database ID."
        )
        possible_IDs = []

        if identifier.smiles is not None:
            smilesBasedID = self.getDBIdentifierFromSmiles(identifier)
            possible_IDs.append(smilesBasedID)

        if identifier.inchi is not None:
            inchiBasedID = self.getDBIdentifierFromInchi(identifier)
            possible_IDs.append(inchiBasedID)

        if identifier.inchi_key is not None:
            inchikeyBasedID = self.getDBIdentifierFromInchiKey(identifier)
            possible_IDs.append(inchikeyBasedID)

        self.logger.debug(
            "Queried all available GenerellIdentifier. Checking if they point to the same database ID."
        )

        if len(possible_IDs) == 0:
            self.logger.warning(
                "Did not find any DB IDs any of the GenerellIdentifiers."
            )
            return None

        allEqual = all(entry == possible_IDs[0] for entry in possible_IDs)

        if not allEqual:
            self.logger.error(
                f"Generell Identifier for supposedly the same object result in different DB IDs in {self.name}."
                f"\n\t\t{'DB ID':<10}{str(inchiBasedID):<10}<-{'InChi':<10}{identifier.inchi:<10}"
                f"\n\t\t{'DB ID':<10}{str(inchikeyBasedID):<10}<-{'InChiKey':<10}{identifier.inchi_key}"
                f"\n\t\t{'DB ID':<10}{str(smilesBasedID):<10}<-{'Smiles':<10}{identifier.smiles}"
            )
            return None
        else:
            self.logger.debug(
                f"All available Identifier point towards the same DB ID in {self.name}:"
                f"\n\t\t{'DB ID':<10}{str(inchiBasedID):<10}<-{'InChi':<10}{identifier.inchi:<10}"
                f"\n\t\t{'DB ID':<10}{str(inchikeyBasedID):<10}<-{'InChiKey':<10}{identifier.inchi_key}"
                f"\n\t\t{'DB ID':<10}{str(smilesBasedID):<10}<-{'Smiles':<10}{identifier.smiles}"
            )

            if not isinstance(possible_IDs[0], str):
                return str(possible_IDs[0])

            return possible_IDs[0]

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

    def _validateGeneralIdentifiersWithDBIDs(
        self, generelID: GenerellIdentifiers, identifier: str
    ):
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
            inchikey_databaseID = self.getDBIdentifierFromInchiKey(
                generelID.inchi_key
            )
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


class DependentDatabase(ABC):
    """
    A dependent Database is a database that maps its metabolites to another Database instead of having associations between their IDs and general IDs like smiles, inchi and InChiKey.
    """

    dependentOn: Database

    @abstractmethod
    def getGenerellIdentifier(self, dbIdentifier: str) -> GenerellIdentifiers:
        """
        This methods queries SMILES, InChI and InChIKey from a database. If any is not available it will be set to None.
        Returns: A GenerellIdentifiers object containing containing SMILES, InChI and InChIKey if available.
        """

        ...

    @abstractmethod
    def getDBIdentifier(self, identifier: GenerellIdentifiers) -> Optional[str]:
        """
        This methods queries the database identifier from a given GenerellIdentifiers.
        Therefore using all defined Generellidentifier and validating that they result in the same database identifier.
        """

        ...
