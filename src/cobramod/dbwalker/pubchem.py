import logging
from typing import Union, Tuple, List

import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from typing_extensions import overload

from cobramod import Settings
from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.cache import Cache
from cobramod.dbwalker.dataclasses import GenerellIdentifiers, Unavailable

class PubChem(Database):

    def __init__(self):
        super().__init__()
        self.__settings = Settings()
        self.__cache_dir = self.__settings.cacheDir / self.name
        self._cache = Cache(cache_dir= self.__cache_dir)
        self.__session = requests.Session()
        self.__session.mount('https://', HTTPAdapter(max_retries=Retry(total=5,
                backoff_factor=0.5,allowed_methods=frozenset(['GET', 'POST']))
        ))

        self.logger = logging.getLogger("cobramod.DBWalker.PubChem")
        self.logger.propagate = True

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

        cached = self._cache.getByID(dbIdentifier)

        if cached is not None and not cached.anyNoneEntries():
            return cached


        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{dbIdentifier}/property/MolecularFormula,InChIKey,InChI,SMILES/json"

        self.logger.debug(
            f"Getting identifiers from PubChem using the following url: {url}"
        )

        self.__settings.limiter.try_acquire("pubchem")
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
                    self.logger.error(
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

        self._cache.addGenerellIdentifiers(generell_identifiers, dbID= dbIdentifier)
        return generell_identifiers


    def getDBIdentifierFromSmiles(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(smiles, GenerellIdentifiers):
            assert smiles.smiles is not None
            smiles = smiles.smiles

        cached = self._cache.getBySmiles(smiles=smiles)
        if isinstance(cached, Unavailable):
            return None

        elif isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/TXT"
        data = {
            "smiles": smiles,
        }

        self.logger.debug(
            f"Getting identifiers from PubChem using the following url: {url} & Post data: {data}"
        )

        self.__settings.limiter.try_acquire("pubchem")
        response = self.__session.post(url, files=data, timeout=30)

        response.raise_for_status()

        value = response.text.rstrip()


        self._cache.addSmiles(smiles=smiles, dbID= value)
        return value

    def getDBIdentifierFromInchi(
        self, inchi: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(inchi, GenerellIdentifiers):
            assert inchi.inchi is not None
            inchi = inchi.inchi

        cached = self._cache.getByInchi(inchi= inchi)
        if isinstance(cached, Unavailable):
            return cached

        elif isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/TXT"

        data = {
            "inchi": inchi,
        }

        self.logger.debug(
            f"Getting identifiers from PubChem using the following url: {url} & Post data: {data}"
        )

        self.__settings.limiter.try_acquire("pubchem")
        response = self.__session.post(url, files=data, timeout=30)

        response.raise_for_status()

        value = response.text.rstrip()

        self._cache.addInchi(inchi= inchi, dbID= value)
        return value

    def getDBIdentifierFromInchiKey(
        self, inchikey: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(inchikey, GenerellIdentifiers):
            assert inchikey.inchi_key is not None
            inchikey = inchikey.inchi_key

        cached = self._cache.getByInchiKey(inchikey = inchikey)
        if isinstance(cached, Unavailable):
            return None

        elif isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/cids/TXT"

        data = {
            "inchikey": inchikey,
        }

        self.__settings.limiter.try_acquire("pubchem")
        response = self.__session.post(url, files=data, timeout=30)

        response.raise_for_status()

        value = response.text.rstrip()

        self._cache.addInchiKey(inchikey=inchikey, dbID= value)
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

        self.__settings.limiter.try_acquire("pubchem")
        response = self.__session.get(url, timeout=30)
        response.raise_for_status()

        value = response.text.rstrip()

        return value

    def getCIDfromSID(self, sid: str) -> str:

        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/cids/TXT"
        )

        self.__settings.limiter.try_acquire("pubchem")
        response = self.__session.get(url, timeout=30)
        response.raise_for_status()

        value = response.text.rstrip()

        return value

    def getSIDfromCID(self, cid: str) -> str:

        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/cid/{cid}/sids/TXT"
        )

        self.__settings.limiter.try_acquire("pubchem")
        response = self.__session.get(url, timeout=30)
        response.raise_for_status()

        value = response.text.rstrip()

        return value

    def save_cache(self):
        self._cache.save_cache()

    def __del__(self):
        self.save_cache()
        self.__session.close()