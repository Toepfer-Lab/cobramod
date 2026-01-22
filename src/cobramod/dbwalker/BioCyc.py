import urllib

import pandas as pd
import requests
from typing import Optional, Union, Tuple, overload, List

import xml.etree.ElementTree as ET
import logging

from rdkit.Chem import inchi

from cobramod.dbwalker.cache import Cache
from cobramod.settings import Settings
from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import GenerellIdentifiers, Unavailable
import cobramod


class BioCyc(Database):

    def __init__(self):
        super().__init__()
        self.__settings = Settings()
        self.__cache_base_dir = self.__settings.cacheDir / self.name
        self._caches = {}

        self.logger = logging.getLogger("cobramod.DBWalker.BioCyc")
        self.logger.propagate = True

    @property
    def name(self) -> str:
        return "BioCyc"

    @property
    def AnnotationPrefix(self) -> str:
        return "biocyc"

    def _get_cache(self, subDB):
        try:
            return self._caches[subDB]
        except KeyError:
            cache = Cache(cache_dir= self.__cache_base_dir / subDB)
            self._caches[subDB] = cache
            return cache

    @overload
    def getGenerellIdentifier(
        self, dbIdentifier: str
    ) -> GenerellIdentifiers: ...

    @overload
    def getGenerellIdentifier(
        self, dbIdentifier: str, BioCycSubDB: str
    ) -> GenerellIdentifiers: ...

    def getGenerellIdentifier(
        self, dbIdentifier: str, BioCycSubDB: Optional[str] = None
    ) -> GenerellIdentifiers:
        """
        Get InChI, InChI Key, and SMILES for a compound using BioCyc ID.

        Args:
            dbIdentifier: BioCyc compound ID (e.g., 'CPD-123')
            BioCycSubDB: BioCyc database name (default: 'META' for MetaCyc)

        Returns:
            Dictionary with 'inchi', 'inchi_key', and 'smiles' keys
        """

        if BioCycSubDB is not None:
            db, biocyc_id = BioCycSubDB, dbIdentifier
            pass
        elif ":" in dbIdentifier:
            db, biocyc_id = dbIdentifier.split(":")
        else:
            db, biocyc_id = "META", dbIdentifier


        cached = self._get_cache(db).getByID(dbIdentifier)

        if cached is not None and not cached.anyNoneEntries():
            return cached


        self.__settings.limiter.try_acquire("biocyc")
        url = "https://websvc.biocyc.org/getxml"

        params = {"id": f"{db}:{biocyc_id}"}

        # BioCyc cuts the SMILES string short otherwise
        headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/xml"}

        self.logger.debug(
            f"Querying BioCyc for compound {biocyc_id} in database {db}"
        )
        self.logger.debug(f"Request URL: {url}")
        self.logger.debug(f"Request parameters: {params}")

        result = GenerellIdentifiers()
        try:
            session = self.__settings._biocycSession
            response = session.get(
                url, params=params, headers=headers, timeout=30
            )
            response.raise_for_status()

            self.logger.debug(f"Response status code: {response.status_code}")

            # Check if BioCyc is asking for account creation
            if (
                response.text
                and isinstance(response.text, str)
                and "<title>Create Account</title>" in response.text
            ):
                self.logger.warning("BioCyc is requesting account creation.")
                if not self.__settings.BioCycLoggedIn:
                    self.logger.error(
                        "BioCyc account information are not provided can't login."
                    )
                return result

            # Parse XML response
            root = ET.fromstring(response.content)

            # Extract SMILES from XML
            smiles_element = root.find(".//string[@title='smiles']")
            if smiles_element is not None and smiles_element.text is not None:
                result.smiles = smiles_element.text.strip()
            else:
                result.smiles = Unavailable()

            # Extract InChI from XML
            inchi_element = root.find(".//inchi[@datatype='string']")
            if inchi_element is not None and inchi_element.text is not None:
                inchi = inchi_element.text.strip()
                result.inchi = inchi

            else:
                result.inchi = Unavailable()

            # Extract InChI Key from XML
            inchikey_element = root.find(".//inchi-key[@datatype='string']")
            if (
                inchikey_element is not None
                and inchikey_element.text is not None
            ):
                inchikey = inchikey_element.text.strip()

                if inchikey.startswith("InChIKey="):
                    result.inchi_key = inchikey[9:]
                else:
                    result.inchi_key = inchikey
            else:
                result.inchi_key = Unavailable()

        except requests.RequestException as e:
            self.logger.error(f"Error fetching data from BioCyc: {e}")
        except ET.ParseError as e:
            self.logger.error(f"Error parsing XML response: {e}")

        self.logger.debug(f"Final result: {result}")

        self.logger.debug(
            "Checking that the general identifiers map back to the original ID"
        )

        # ToDo validate that the moethod works on the object and not on a copy due to the namespace
        self._validateGeneralIdentifiersWithDBIDs(
            generelID=result, identifier=dbIdentifier
        )

        self._get_cache(db).addGenerellIdentifiers(result, dbID=(db, dbIdentifier))

        return result

    def getDBIdentifierFromSmiles(
        self,
        smiles: Union[str, GenerellIdentifiers],
        BioCycSubDB: Optional[str] = None,
    ) -> Optional[Union[str, List[str]]]:
        """
        Convert SMILES to BioCyc identifiers.

        Args:
            smiles: SMILES string of the compound

        Returns:
            GenerellIdentifiers with BioCyc ID if found, otherwise empty
        """

        if BioCycSubDB is None:
            BioCycSubDB = "META"

        if isinstance(smiles, GenerellIdentifiers):
            assert smiles.smiles is not None
            smiles = smiles.smiles

        assert isinstance(smiles, str)

        cached = self._get_cache(BioCycSubDB).getBySmiles(smiles=smiles)

        if isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        try:
            self.logger.info(f"Looking up BioCyc ID for SMILES: {smiles}")

            session = self.__settings._biocycSession

            self.__settings.limiter.try_acquire("biocyc")

            response = session.post(
                f"https://websvc.biocyc.org/{BioCycSubDB}/smiles-search",
                data={"smiles": smiles, "exact": "T", "fmt": "json"},
            )

            response.raise_for_status()

            data = response.json()

            results = data["RESULTS"]

            if results is None:
                self.logger.debug(
                    f"BioCyc did not report any matches for SMILES ({smiles}). "
                )

                self._get_cache(BioCycSubDB).addSmiles(smiles= smiles, dbID= Unavailable())
                return Unavailable()

            if len(results) > 1:
                self.logger.debug(
                    "BioCyc reported more than one exact match for SMILES. Now checking for string equality."
                )

            hits = 0
            matches: List[str] = []
            printable: List[str, str] = []

            for result in results:
                if result["SMILES"] == smiles:
                    hits = hits + 1
                    matches.append(result["OBJECT-ID"])
                    printable.append(
                        (result["COMMON-NAME"], result["OBJECT-ID"])
                    )

            if hits > 1:
                self.logger.error(
                    f"Found more than one entry for SMILES ({smiles}) in BioCyc ({BioCycSubDB}). Wont add anything due to being uncertain which ones is correct. The following IDs were found for: {printable})"
                )
                # ToDo raise error for uncertain match
                biocycID = Unavailable()
            elif hits == 1:
                self.logger.debug(
                    "Only one of the matches was exactly equal to the query. Proceeding normally."
                )
                biocycID = matches[0]
            else:
                self.logger.debug(
                    "No matches were found. Proceeding without match."
                )
                biocycID = Unavailable()

            self._get_cache(BioCycSubDB).addSmiles(smiles=smiles, dbID= Unavailable())
            return biocycID

        except requests.RequestException as e:
            self.logger.error(f"Error fetching data from BioCyc: {e}")
            return None

    def getDBIdentifierFromInchi(
        self,
        inchi: Union[str, GenerellIdentifiers],
        BioCycSubDB: Optional[str] = None,
    ) -> Optional[str]:
        """
        Convert InChI to BioCyc identifiers.

        Args:
            inchi: InChI string of the compound

        Returns:
            GenerellIdentifiers with BioCyc ID if found, otherwise empty
        """
        if BioCycSubDB is None:
            BioCycSubDB = "META"

        if isinstance(inchi, GenerellIdentifiers):
            assert inchi.inchi is not None
            inchi = inchi.inchi

        cached = self._get_cache(BioCycSubDB).getByInchi(inchi=inchi)

        if isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        url = f"https://websvc.biocyc.org/{BioCycSubDB}/inchi-search?inchi={inchi}&exact=T&fmt=json"

        self.logger.info(f"Requesting BioCyc for InChI: {inchi} using {url}")

        session = self.__settings._biocycSession

        self.__settings.limiter.try_acquire("biocyc")
        response = session.get(url, timeout=30)

        try:
            response.raise_for_status()

            data = response.json()
            if data["RESULTS"] is None:
                self.logger.warning(
                    f"Could not find a Entry in BioCyc({BioCycSubDB}) for InChI: {inchi}"
                )
                self._get_cache(BioCycSubDB).addInchi(inchi=inchi, dbID= Unavailable())
                return Unavailable()

            biocycID = data["RESULTS"][0]["OBJECT-ID"]

            if biocycID is None:
                self._get_cache(BioCycSubDB).addInchi(inchi=inchi, dbID= Unavailable())
                return Unavailable()

            return biocycID

        except requests.RequestException as e:
            self.logger.error(f"Error fetching data from BioCyc: {e}")
            self._get_cache(BioCycSubDB).addInchi(inchi=inchi, dbID=Unavailable())
            return Unavailable()

    def getDBIdentifierFromInchiKey(
        self,
        inchikey: Union[str, GenerellIdentifiers],
        BioCycSubDB: Optional[str] = None,
    ) -> Optional[str]:
        if BioCycSubDB is None:
            BioCycSubDB = "META"

        if isinstance(inchikey, GenerellIdentifiers):
            assert inchikey.inchi is not None
            inchikey = inchikey.inchi

        cached = self._get_cache(BioCycSubDB).getByInchiKey(inchikey = inchikey)

        if isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        url = f"https://biocyc.org/{BioCycSubDB}/search-query?type=COMPOUND&inchikey={inchikey}&exact=T&fmt=json"
        self.logger.info(
            f"Requesting BioCyc for InChIKey: {inchikey} using: {url}"
        )

        session = self.__settings._biocycSession

        self.__settings.limiter.try_acquire("biocyc")

        response = session.get(url, timeout=30)

        try:
            response.raise_for_status()

            data = response.json()
            biocycID = data["RESULTS"][0]["OBJECT-ID"]

            if biocycID is None:
                self._get_cache(BioCycSubDB).addInchiKey(inchikey=inchikey, dbID=Unavailable())
                return Unavailable()

            self._get_cache(BioCycSubDB).addInchiKey(inchikey=inchikey, dbID=biocycID)
            return biocycID

        except requests.RequestException as e:
            self.logger.error(f"Error fetching data from BioCyc: {e}")

            self._get_cache(BioCycSubDB).addInchiKey(inchikey=inchikey, dbID=Unavailable())
            return Unavailable()

    def __del__(self):
        for cache in self._caches.values():
            cache.save_cache()
