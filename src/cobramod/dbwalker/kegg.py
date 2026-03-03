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
        if cid is Unavailable:
            return Unavailable

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
        if cached is Unavailable:
            return Unavailable

        elif isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        cid = self.__pubchem.getDBIdentifierFromSmiles(smiles=smiles)
        if cid is Unavailable:
            logger.info(
                f"Did not find a hit for SMILES ({smiles}) in PubChem."
            )
            self._cache.addSmiles(smiles=smiles, dbID=Unavailable)
            return Unavailable

        sid = self.__pubchem.getSIDsFromCIDs(cid=cid)
        if sid is Unavailable:
            logger.info(
                f"No SIDs found for CID ({cid})"
            )
            self._cache.addSmiles(smiles=smiles, dbID=Unavailable)
            return Unavailable

        kegg_specific_sid = self.__pubchem.getKeggSpecificSIDs(sid)
        if kegg_specific_sid is Unavailable:
            logger.info(
                f"No Kegg specific SID found in SIDs ({sid})"
            )
            self._cache.addSmiles(smiles=smiles, dbID=Unavailable)
            return Unavailable

        kegg_id = self.get_kegg_id_from_sid(sid=kegg_specific_sid)
        logger.info(f"Got Kegg ID ({kegg_id}) for Kegg specific SID ({sid})")

        self._cache.addSmiles(smiles=smiles, dbID=kegg_id)
        return kegg_id

    def getDBIdentifierFromInchi(
        self, inchi: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(inchi, GenerellIdentifiers):
            inchi = inchi.inchi

        cached = self._cache.getByInchi(inchi=inchi)
        if cached is Unavailable:
            return cached

        elif isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        cid = self.__pubchem.getDBIdentifierFromInchi(inchi=inchi)

        if cid is Unavailable:
            logger.info(
                f"No CID found for InChI ({inchi})"
            )
            self._cache.addInchi(inchi=inchi, dbID=Unavailable)
            return Unavailable

        sid = self.__pubchem.getSIDsFromCIDs(cid=cid)
        if sid is Unavailable:
            logger.info(
                f"No SIDs found for CID ({cid})"
            )
            self._cache.addInchi(inchi=inchi, dbID=Unavailable)
            return Unavailable

        kegg_specific_sid = self.__pubchem.getKeggSpecificSIDs(sid)
        if kegg_specific_sid is Unavailable:
            logger.info(
                f"No Kegg specific SIDs found in SIDs ({sid})")
            self._cache.addInchi(inchi=inchi, dbID=Unavailable)
            return Unavailable

        kegg_id = self.get_kegg_id_from_sid(sid=kegg_specific_sid)
        logger.info(f"Found Kegg ID ({kegg_id}) for Kegg specific SID ({sid})")

        self._cache.addInchi(inchi=inchi, dbID=kegg_id)
        return kegg_id

    def getDBIdentifierFromInchiKey(
        self, inchikey: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(inchikey, GenerellIdentifiers):
            inchikey = inchikey.inchi_key

        cached = self._cache.getByInchiKey(inchikey=inchikey)
        if cached is Unavailable:
            return Unavailable

        elif isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        cid = self.__pubchem.getDBIdentifierFromInchiKey(inchikey=inchikey)
        if cid is Unavailable:
            logger.info(
                f"No CID found for InChIKey ({inchikey})"
            )
            self._cache.addInchiKey(inchikey=inchikey, dbID=Unavailable)
            return Unavailable

        sid = self.__pubchem.getSIDsFromCIDs(cid=cid)
        if sid is Unavailable:
            logger.info(
                f"No SIDs found for CID ({cid})"
            )
            self._cache.addInchiKey(inchikey=inchikey, dbID=Unavailable)
            return Unavailable

        kegg_specific_sid = self.__pubchem.getKeggSpecificSIDs(sid)
        if kegg_specific_sid is Unavailable:
            logger.info(
                f"No Kegg specific SID found in SIDs ({sid})"
            )

        kegg_id = self.get_kegg_id_from_sid(sid=kegg_specific_sid)

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
            logger.info(f"No KEGG compound ID found for SID({sid})")
            return Unavailable

        # Parse the first matching line
        try:
            kegg_id = lines[0].split("\t")[1]
        except IndexError:
            logger.info(f"No KEGG compound ID found for SID({sid})")
            return Unavailable

        if ":" in kegg_id:
            kegg_id = kegg_id.split(":")[1]

        # kegg_id = f"kegg.compound:{kegg_id}"

        logger.info(f"Got Kegg ID ({kegg_id}) for SID ({sid})")

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
            return Unavailable

        # Parse the first matching line
        pubchem_entry = lines[0].split("\t")[1]

        if ":" in pubchem_entry:
            sid = pubchem_entry.split(":")[1]
        else:
            sid = pubchem_entry

        logger.info(f"Found SID ({sid}) for KEGG ID {dbIdentifier}")
        cid = self.__pubchem.getCIDsFromSIDs(sid=sid)

        if cid is Unavailable:
            return Unavailable
        elif len(cid) == 1:
            cid = cid[0]
        else:
            return Unavailable  # ToDo Uncertain?

        return cid

    def save_cache(self):
        self._cache.save_cache()
