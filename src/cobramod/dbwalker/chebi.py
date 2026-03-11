import logging
from pathlib import Path
from typing import Optional, Tuple, Union

import pandas as pd
import requests

from cobramod import Settings
from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import (
    GenerellIdentifiers,
    Unavailable,
    UnavailableType,
    Uncertain,
)

logger = logging.getLogger("cobramod.DBWalker.Chebi")
logger.propagate = True


class Chebi(Database):
    chebi_ftp = "https://ftp.ebi.ac.uk/pub/databases/chebi"
    settings = Settings()

    def __init__(self):
        super().__init__()

        self.__getStructureFile()
        self.__load_structure_file()

    @property
    def AnnotationPrefix(self) -> str:
        return "CHEBI"

    @property
    def name(self) -> str:
        return "Chebi"

    def save_cache(self):
        pass

    def __getStructureFile(self):
        """
        Check if a local structure file is available and of not get it from chebis ftp server.

        """

        logger.debug("Checking if local structure file exists")

        # structure file
        url = self.chebi_ftp + "/flat_files/structures.tsv.gz"
        cached_file = Path(
            self.settings.cacheDir / "chebi" / "chebi-structure.tsv.gz"
        )
        cached_file.parent.mkdir(parents=True, exist_ok=True)

        if cached_file.exists():
            logger.debug("Found local structure file")
            return

        logger.debug("ChEBI structure file not found. Downloading...")

        response = requests.get(url)

        if response.status_code == 200:
            with open(cached_file, "wb") as file:
                file.write(response.content)
            logger.debug(
                f"ChEBI structure file downloaded and saved at {str(cached_file)}."
            )
        else:
            logger.error(
                f"Download failed with status code {response.status_code}"
            )
            response.raise_for_status()

        return

    def __load_structure_file(self):
        structure_file = (
                self.settings.cacheDir / "chebi" / "chebi-structure.tsv.gz"
        )

        self.structure_file = pd.read_csv(
            structure_file, sep="\t", compression="gzip"
        )

    def getGenerellIdentifier(self, dbIdentifier: str) -> GenerellIdentifiers:

        if isinstance(dbIdentifier, str):
            dbIdentifier = int(dbIdentifier)

        row = self.structure_file.loc[self.structure_file["compound_id"] == dbIdentifier]

        if len(row) == 0:
            return Unavailable

        if len(row) > 1:
            raise ValueError("The CHEBI flat file contains one ID two times. Please check that the flat file is valid and create an issue at our GitHub repository: https://github.com/Toepfer-Lab/cobramod")

        smiles = row["smiles"].iloc[0]
        inchi = row["standard_inchi"].iloc[0]
        inchikey = row["standard_inchi_key"].iloc[0]

        gID = GenerellIdentifiers(
            smiles=smiles if pd.notna(smiles) else Unavailable,
            inchi=inchi if pd.notna(inchi) else Unavailable,
            inchi_key=inchikey if pd.notna(inchikey) else Unavailable,
        )

        return gID

    def getDBIdentifierFromSmiles(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> Optional[str]:
        if isinstance(smiles, GenerellIdentifiers):
            smiles = smiles.smiles

        result = self.structure_file.loc[
            self.structure_file["smiles"] == smiles, "compound_id"
        ]

        if len(result) == 1:
            return str(result.iloc[0])

        elif len(result) > 1:
            logger.warning(
                f"Found multiple ({len(result)}) entries for SMILES ({smiles}) in Chebi"
            )

            return Uncertain(
                type="DBID",
                possibilities=result.astype(str).tolist(),
            )

        else:
            logger.warning(f"Found no match for SMILES ({smiles}) in Chebi)")

            return Unavailable

    def getDBIdentifierFromInchi(
        self, inchi: Union[str, GenerellIdentifiers]
    ) -> Optional[str]:
        if isinstance(inchi, GenerellIdentifiers):
            inchi = inchi.inchi

        result = self.structure_file.loc[
            self.structure_file["standard_inchi"] == inchi, "compound_id"
        ]

        if len(result) == 1:
            return str(result.iloc[0])

        elif len(result) > 1:
            logger.warning(
                f"Found multiple ({len(result)}) entries for InChI ({inchi}) in Chebi"
            )

            return Uncertain(
                type="DBID",
                possibilities=result,
            )
        else:
            logger.warning(f"Found no match for InChI ({inchi}) in Chebi)")

            return Unavailable

    def getDBIdentifierFromInchiKey(
        self, inchikey: Union[str, GenerellIdentifiers]
    ) -> Union[str, UnavailableType]:
        if isinstance(inchikey, GenerellIdentifiers):
            inchikey = inchikey.inchi_key

        result = self.structure_file.loc[
            self.structure_file["standard_inchi_key"] == inchikey, "compound_id"
        ]

        if len(result) == 1:
            return str(result.iloc[0])

        elif len(result) > 1:
            logger.warning(
                f"Found multiple ({len(result)}) entries for InChIKey ({inchikey}) in Chebi"
            )

            return Uncertain(
                type="DBID",
                possibilities=result,
            )
        else:
            logger.warning(
                f"Found no match for InChIKey ({inchikey}) in Chebi)"
            )

            return Unavailable

    def validateGenerellIdentifiers(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> Tuple[GenerellIdentifiers, GenerellIdentifiers]:
        raise NotImplementedError
