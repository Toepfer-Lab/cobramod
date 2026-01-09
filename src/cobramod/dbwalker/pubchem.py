import logging
from typing import Union, Tuple, List

import pandas as pd
import requests
from typing_extensions import overload

from cobramod import Settings
from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import GenerellIdentifiers

logger = logging.getLogger("cobramod.DBWalker.PubChem")
logger.propagate = True


class PubChem(Database):

    def __init__(self):
        super().__init__()
        self.__settings = Settings()
        self.__cache_file = self.__settings.cacheDir / f"{self.name}.csv"
        self._cache_added = 0
        self.__session = requests.Session()

        if not self.__cache_file.exists():
            self.cache = pd.DataFrame(columns=["DB-ID", "SMILES", "InChI", "InChIKey"])
            self.cache.to_csv(self.__cache_file, index=False)
        else:
            self.cache = pd.read_csv(self.__cache_file)

    @property
    def name(self) -> str:
        return "PubChem"

    @property
    def AnnotationPrefix(self) -> str:
        return "pubchem.compound"

    def getGenerellIdentifier(
        self, dbIdentifier: Union[str,int], **kwargs
    ) -> GenerellIdentifiers:
        if isinstance(dbIdentifier, str) and "pubchem.compound:" in dbIdentifier:
            dbIdentifier = int(dbIdentifier.replace("pubchem.compound:", ""))
        elif isinstance(dbIdentifier, str) and dbIdentifier.isdigit():
            dbIdentifier = int(dbIdentifier)

        try:
            cached_results = self.cache.loc[self.cache["DB-ID"] == dbIdentifier]
            gID = GenerellIdentifiers()
            gID.smiles = cached_results["SMILES"]
            gID.inchi = cached_results["InChI"]
            gID.inchi_key = cached_results["InChIKey"]
            return gID

        except KeyError:
            pass

        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{dbIdentifier}/property/MolecularFormula,InChIKey,InChI,SMILES/json"

        logger.debug(
            f"Getting identifiers from PubChem using the following url:\n{url}"
        )

        response = self.__session.get(url=url, timeout=30)
        response.raise_for_status()

        value = response.json()

        generell_identifiers = GenerellIdentifiers()

        for entry in value["PropertyTable"]["Properties"]:
            inchi = entry.get("InChI")
            if generell_identifiers.inchi is not None:
                assert generell_identifiers.inchi == inchi
            else:
                generell_identifiers.inchi = inchi

            inchikey = entry.get("InChIKey")
            if generell_identifiers.inchi_key is not None:
                if not generell_identifiers.inchi_key == inchikey:
                    logger.error(
                        "Results from PubChem for InChiKeys dont match. Got at least two diffrent InChIKeys:"
                        f"{generell_identifiers.inchi_key} and {inchikey}."
                    )
                    raise ValueError()
            else:
                generell_identifiers.inchi_key = inchikey

            smiles = entry.get("SMILES")
            if generell_identifiers.smiles is not None:
                assert generell_identifiers.smiles == smiles
            else:
                generell_identifiers.smiles = smiles

        self.cache.loc[len(self.cache.index)] = [dbIdentifier, generell_identifiers.smiles, generell_identifiers.inchi, generell_identifiers.inchi_key]

        self._cache_added += 1

        if self._cache_added % 10 == 0:
            self.save_cache()

        return generell_identifiers


    def getDBIdentifierFromSmiles(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(smiles, GenerellIdentifiers):
            assert smiles.smiles is not None
            smiles = smiles.smiles

        try:
            cached_results = self.cache.loc[self.cache["SMILES"] == smiles]["DB-ID"]
            if len(cached_results)> 0:
                return cached_results.iloc[0]

        except KeyError:
            pass

        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/TXT"
        data = {
            "smiles": smiles,
        }

        self.logger.debug(
            f"Getting identifiers from PubChem using the following url:\n{url}\n & Post data: {data}"
        )
        response = self.__session.post(url, files=data, timeout=30)

        response.raise_for_status()

        value = response.text.rstrip()

        return value

    def getDBIdentifierFromInchi(
        self, inchi: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(inchi, GenerellIdentifiers):
            assert inchi.inchi is not None
            inchi = inchi.inchi

        try:
            cached_results = self.cache.loc[self.cache["InChI"] == inchi]["DB-ID"]
            if len(cached_results)> 0:
                return cached_results.iloc[0]

        except KeyError:
            pass

        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/TXT"

        data = {
            "inchi": inchi,
        }

        self.logger.debug(
            f"Getting identifiers from PubChem using the following url:\n{url}\n & Post data: {data}"
        )
        response = self.__session.post(url, files=data, timeout=30)

        response.raise_for_status()

        value = response.text.rstrip()

        return value

    def getDBIdentifierFromInchiKey(
        self, inchikey: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(inchikey, GenerellIdentifiers):
            assert inchikey.inchi_key is not None
            inchikey = inchikey.inchi_key

        try:
            cached_results = self.cache.loc[self.cache["InChIKey"] == inchikey]["DB-ID"]
            if len(cached_results)> 0:
                return cached_results.iloc[0]

        except KeyError:
            pass

        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/cids/TXT"

        data = {
            "inchikey": inchikey,
        }

        response = self.__session.post(url, files=data, timeout=30)

        response.raise_for_status()

        value = response.text.rstrip()

        return value

    def validateGenerellIdentifiers(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> Tuple[GenerellIdentifiers, GenerellIdentifiers]:
        raise NotImplementedError

    def getChemicalFormular(self, cid: str) -> str:
        """
        Query PubChem by CID to retrieve the chemical formula.

        Args:
            cid (int or str): PubChem Compound ID

        Returns:
            str: Chemical formula or None if not found
        """

        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{cid}/property/MolecularFormula/txt"
        )

        response = self.__session.get(url, timeout=30)
        response.raise_for_status()

        value = response.text.rstrip()

        return value

    def save_cache(self):
        self.cache.to_csv(
            self.__cache_file,
            index=False,
        )

    def __del__(self):
        self.save_cache()
        self.__session.close()