import atexit
import logging
from typing import Optional, Tuple, Union, overload

import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry

from cobramod import Settings
from cobramod.dbwalker.cache import Cache
from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import GenerellIdentifiers, Unavailable
from cobramod.dbwalker.pubchem import PubChem

logger = logging.getLogger("cobramod.DBWalker.Kegg")
logger.propagate = True


class Kegg(Database):
    def __init__(self):
        super().__init__()
        self.__pubchem = PubChem()
        self.__session = requests.Session()
        self.__session.mount(
            "https://",
            HTTPAdapter(max_retries=Retry(total=5, backoff_factor=0.5)),
        )

        self.__settings = Settings()
        self._cache_folder = self.__settings.cacheDir / self.name
        self._cache = Cache(cache_dir=self._cache_folder)

        atexit.register(self.save_cache)

    @property
    def name(self) -> str:
        return "Kegg"

    @property
    def AnnotationPrefix(self) -> str:
        return "kegg.compound"

    def getGenerellIdentifier(
        self, dbIdentifier: str, **kwargs
    ) -> GenerellIdentifiers:
        cached = self._cache.getByID(dbIdentifier)

        if cached is not None and not cached.anyNoneEntries():
            return cached

        cid = self.get_cid_from_kegg_id(dbIdentifier=dbIdentifier)
        if isinstance(cid, Unavailable):
            return Unavailable()

        # Kegg itself does not provide identifiers instead it maps them to PubChem

        generellIdentifier = self.__pubchem.getGenerellIdentifier(
            dbIdentifier=cid
        )

        self._cache.addGenerellIdentifiers(
            generellIdentifier, dbID=dbIdentifier
        )
        return generellIdentifier

    def getDBIdentifierFromSmiles(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(smiles, GenerellIdentifiers):
            smiles = smiles.smiles

        cached = self._cache.getBySmiles(smiles=smiles)
        if isinstance(cached, Unavailable):
            return cached

        elif isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        cid = self.__pubchem.getDBIdentifierFromSmiles(smiles=smiles)
        sid = self.__pubchem.getSIDsFromCIDs(cid=cid)

        kegg_id = self.get_kegg_id_from_sid(sid=sid)

        self._cache.addSmiles(smiles=smiles, dbID=kegg_id)
        return kegg_id

    def getDBIdentifierFromInchi(
        self, inchi: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(inchi, GenerellIdentifiers):
            inchi = inchi.inchi

        cached = self._cache.getByInchi(inchi=inchi)
        if isinstance(cached, Unavailable):
            return cached

        elif isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        cid = self.__pubchem.getDBIdentifierFromInchi(inchi=inchi)

        if isinstance(cid, Unavailable):
            self._cache.addInchi(inchi=inchi, dbID=Unavailable())
            return Unavailable()

        sid = self.__pubchem.getSIDsFromCIDs(cid=cid)

        kegg_id = self.get_kegg_id_from_sid(sid=sid)

        self._cache.addInchi(inchi=inchi, dbID=kegg_id)
        return kegg_id

    def getDBIdentifierFromInchiKey(
        self, inchikey: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(inchikey, GenerellIdentifiers):
            inchikey = inchikey.inchi_key

        cached = self._cache.getByInchiKey(inchikey=inchikey)
        if isinstance(cached, Unavailable):
            return cached

        elif isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        cid = self.__pubchem.getDBIdentifierFromInchiKey(inchikey=inchikey)
        sid = self.__pubchem.getSIDsFromCIDs(cid=cid)

        kegg_id = self.get_kegg_id_from_sid(sid=sid)

        self._cache.addInchiKey(inchikey=inchikey, dbID=kegg_id)
        return kegg_id

    def get_kegg_id_from_sid(self, sid):
        """
        Given a PubChem CID, return the corresponding KEGG compound ID (e.g., 'C00022').
        """

        self.__settings.limiter.try_acquire("kegg")
        url = f"https://rest.kegg.jp/conv/compound/pubchem:{sid}"
        response = self.__session.get(url)

        if response.status_code != 200:
            logger.error(
                f"Error ({response.status_code}) getting KEGG compound ID from {sid}"
            )
            response.raise_for_status()

        lines = response.text.strip().split("\n")
        if not lines:
            return Unavailable()

        # Parse the first matching line
        try:
            kegg_id = lines[0].split("\t")[1]
        except IndexError:
            return Unavailable()

        if ":" in kegg_id:
            kegg_id = kegg_id.split(":")[1]

        # kegg_id = f"kegg.compound:{kegg_id}"

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

        self.__settings.limiter.try_acquire("kegg")
        url = f"https://rest.kegg.jp/conv/pubchem/cpd:{compound_id}"
        response = self.__session.get(url)

        if response.status_code != 200:
            logger.error(
                f"Error ({response.status_code}) while getting CID for {compound_id}"
            )
            response.raise_for_status()

        lines = response.text.strip().split("\n")
        if not lines:
            logger.info(f"No entry found for KEGG ID {dbIdentifier}")
            return Unavailable()

        # Parse the first matching line
        pubchem_entry = lines[0].split("\t")[1]

        if ":" in pubchem_entry:
            sid = pubchem_entry.split(":")[1]
        else:
            sid = pubchem_entry

        cid = self.__pubchem.getCIDsFromSIDs(sid=sid)

        if isinstance(cid, Unavailable):
            return cid
        elif len(cid) == 1:
            cid = cid[0]
        else:
            return Unavailable()

        return cid

    def save_cache(self):
        self._cache.save_cache()
