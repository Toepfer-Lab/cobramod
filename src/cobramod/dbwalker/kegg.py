import logging
from typing import Union, Tuple, Optional, overload

import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry

from cobramod import Settings
from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import GenerellIdentifiers
from cobramod.dbwalker.pubchem import PubChem

logger = logging.getLogger("cobramod.DBWalker.Kegg")
logger.propagate = True

# Ensure console output
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


class Kegg(Database):

    def __init__(self):
        super().__init__()
        self.__pubchem = PubChem()
        self.__session = requests.Session()
        self.__session.mount('https://', HTTPAdapter(max_retries=Retry(total=5,
                backoff_factor=0.5)
        ))

        self.__settings = Settings()
        self.__cache_file = self.__settings.cacheDir / f"{self.name}.csv"
        self._cache_added = 0

        if not self.__cache_file.exists():
            self.cache = pd.DataFrame(columns=["DB-ID", "CID"])
            self.cache.to_csv(self.__cache_file, index=False)
        else:
            self.cache = pd.read_csv(self.__cache_file)


    @property
    def name(self) -> str:
        return "Kegg"

    @property
    def AnnotationPrefix(self) -> str:
        return "kegg.compound"

    def getGenerellIdentifier(
        self, dbIdentifier: str, **kwargs
    ) -> GenerellIdentifiers:
        cid = self.get_cid_from_kegg_id(dbIdentifier=dbIdentifier)

        # Kegg itself does not provide identifiers instead it maps them to PubChem

        generellIdentifier = self.__pubchem.getGenerellIdentifier(dbIdentifier=cid)

        return generellIdentifier

    def getDBIdentifierFromSmiles(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> str:
        cid = self.__pubchem.getDBIdentifierFromSmiles(smiles=smiles)
        kegg_id = self.get_kegg_id_from_cid(cid= cid)

        return kegg_id

    def getDBIdentifierFromInchi(
        self, inchi: Union[str, GenerellIdentifiers]
    ) -> str:

        cid = self.__pubchem.getDBIdentifierFromInchi(inchi = inchi)
        kegg_id = self.get_kegg_id_from_cid(cid=cid)

        return kegg_id

    def getDBIdentifierFromInchiKey(
        self, inchikey: Union[str, GenerellIdentifiers]
    ) -> str:

        cid = self.__pubchem.getDBIdentifierFromInchiKey(inchikey = inchikey)
        kegg_id = self.get_kegg_id_from_cid(cid=cid)

        return kegg_id

    def get_kegg_id_from_cid(self, cid):
        """
        Given a PubChem CID, return the corresponding KEGG compound ID (e.g., 'C00022').
        """

        try:
            cached_results = self.cache.loc[self.cache["CID"] == cid]["DB-ID"]
            if len(cached_results) > 0:
                return cached_results.iloc[0]

        except KeyError:
            pass

        self.__settings.limiter.try_acquire("kegg")
        url = f"https://rest.kegg.jp/conv/compound/pubchem:{cid}"
        response = self.__session.get(url)

        if response.status_code != 200:
            logger.error(
                f"Error ({response.status_code}) getting KEGG compound ID from {cid}"
            )
            return None

        lines = response.text.strip().split("\n")
        if not lines:
            return None

        # Parse the first matching line
        try:
            kegg_id = lines[0].split("\t")[1]
        except IndexError:
            return None

        if ":" in kegg_id:
            kegg_id = kegg_id.split(":")[1]

        #kegg_id = f"kegg.compound:{kegg_id}"

        self.cache.loc[len(self.cache.index)] = [kegg_id,cid]
        self._cache_added += 1
        if self._cache_added % 10 == 0:
            self.save_cache()

        return kegg_id

    def get_cid_from_kegg_id(self, dbIdentifier):
        """
        Given a KEGG compound ID (e.g., 'C00022' or 'kegg.compound:C00022'),
        return the corresponding PubChem CID.
        """
        # Extract the compound ID if it's in full format
        if ":" in dbIdentifier:
            compound_id = dbIdentifier.split(":")[1]
        else:
            compound_id = dbIdentifier

        try:
            cached_results = self.cache.loc[self.cache["DB-ID"] == dbIdentifier]["CID"]
            if len(cached_results) > 0:
                return cached_results.iloc[0]

        except KeyError:
            pass

        self.__settings.limiter.try_acquire("kegg")
        url = f"https://rest.kegg.jp/conv/pubchem/cpd:{compound_id}"
        response = self.__session.get(url)

        if response.status_code != 200:
            logger.error(
                f"Error ({response.status_code}) while getting CID for {compound_id}"
            )
            return None

        lines = response.text.strip().split("\n")
        if not lines:
            logger.info(f"No entry found for KEGG ID {dbIdentifier}")
            return None

        # Parse the first matching line
        pubchem_entry = lines[0].split("\t")[1]

        if ":" in pubchem_entry:
            cid = pubchem_entry.split(":")[1]
        else:
            cid = pubchem_entry

        self.cache.loc[len(self.cache.index)] = [dbIdentifier, cid]
        self._cache_added += 1
        if self._cache_added % 10 == 0:
            self.save_cache()

        return cid

    def save_cache(self):
        self.cache.to_csv(
            self.__cache_file,
            index=False,
        )

    def __del__(self):
        self.save_cache()
        self.__session.close()