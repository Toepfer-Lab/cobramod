import urllib

import pandas as pd
import requests
from typing import Optional, Union, Tuple, overload, List

import xml.etree.ElementTree as ET
import logging

from rdkit.Chem import inchi

from cobramod.settings import Settings
from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import GenerellIdentifiers
import cobramod


class BioCyc(Database):

    def __init__(self):
        super().__init__()
        self.__settings = Settings()
        self._cacheDir = self.__settings.cacheDir / "biocyc"

        if not self._cacheDir.exists():
            self._cacheDir.mkdir()


        self.__cache_file = self._cacheDir / f"{self.name}.csv"
        self._cache_added = 0

        self._cache_smiles_not_found = set()
        self._cache_inchi_not_found = set()
        self._cache_inchikey_notfound = set()

        if not self.__cache_file.exists():
            index = pd.MultiIndex(levels=[[], []], names=[u'SUBDB', u'BioCycID'], codes=[[], []])
            self.cache = pd.DataFrame(index= index, columns=["SMILES", "InChI", "InChIKey"])
            self.cache.to_csv(self.__cache_file)
        else:
            self.cache = pd.read_csv(self.__cache_file, index_col = ["SUBDB", "BioCycID"])

        self.logger = logging.getLogger("cobramod.DBWalker.BioCyc")
        self.logger.propagate = True




        # Ensure console output
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)

    @property
    def name(self) -> str:
        return "BioCyc"

    @property
    def AnnotationPrefix(self) -> str:
        return "biocyc"

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

        try:
            cached_results = self.cache.loc[(db, biocyc_id)]
            gID = GenerellIdentifiers()
            gID.smiles = cached_results["SMILES"]
            gID.inchi = cached_results["InChI"]
            gID.inchi_key = cached_results["InChIKey"]
            return gID

        except KeyError:
            pass

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

            # Extract InChI from XML
            inchi_element = root.find(".//inchi[@datatype='string']")
            if inchi_element is not None and inchi_element.text is not None:
                inchi = inchi_element.text.strip()
                if inchi.startswith("InChI="):
                    result.inchi = inchi[6:]
                else:
                    result.inchi = inchi

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
        self.cache.loc[(db,biocyc_id),:] = [result.smiles, result.inchi, result.inchi_key]

        self._cache_added += 1
        if self._cache_added % 10 == 0:
            self.save_cache()

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

        try:

            cached_results = self.cache.loc[self.cache["SMILES"] == smiles].index.to_list()

            if len(cached_results) == 1:
                return cached_results
            elif len(cached_results) > 1:
                raise LookupError

        except KeyError:
            pass

        if smiles in self._cache_smiles_not_found:
            return None

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

                self._cache_smiles_not_found.add(smiles)
                return None

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
                    f"Found more than one entry for SMILES ({smiles}) in BioCyc ({BioCycSubDB}). Wont add anything due to being uncertain which ones is correct. The following IDs were found for:\n {printable})"
                )
                # ToDo raise error for uncertain match
                biocycID = None
            elif hits == 1:
                self.logger.debug(
                    "Only one of the matches was exactly equal to the query. Proceeding normally."
                )
                biocycID = matches[0]
            else:
                self.logger.debug(
                    "No matches were found. Proceeding without match."
                )
                self._cache_smiles_not_found.add(smiles)
                biocycID = None

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

        if inchi in self._cache_inchi_not_found:
            return None

        try:
            cached_results = self.cache.loc[self.cache["InChI"] == inchi]["DB-ID"]
            return cached_results.to_list()
        except KeyError:
            pass

        url = f"https://websvc.biocyc.org/{BioCycSubDB}/inchi-search?inchi={inchi}&exact=T&fmt=json"

        self.logger.info(f"Requesting BioCyc for InChI: {inchi}\nusing {url}")

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
                self._cache_inchi_not_found.add(inchi)
                return None

            biocycID = data["RESULTS"][0]["OBJECT-ID"]

            if biocycID is None:
                self._cache_inchi_not_found.add(inchi)

            return biocycID

        except requests.RequestException as e:
            self.logger.error(f"Error fetching data from BioCyc: {e}")
            self._cache_inchi_not_found.add(inchi)
            return None

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

        if inchikey in self._cache_inchikey_notfound:
            return None

        try:
            cached_results = self.cache.loc[self.cache["InChIKey"] == inchikey]["DB-ID"]
            return cached_results.to_list()
        except KeyError:
            pass

        url = f"https://biocyc.org/{BioCycSubDB}/search-query?type=COMPOUND&inchikey={inchikey}&exact=T&fmt=json"
        self.logger.info(
            f"Requesting BioCyc for InChIKey: {inchikey}\nusing: {url}"
        )

        session = self.__settings._biocycSession

        self.__settings.limiter.try_acquire("biocyc")

        response = session.get(url, timeout=30)

        try:
            response.raise_for_status()

            data = response.json()
            biocycID = data["RESULTS"][0]["OBJECT-ID"]

            if biocycID is None:
                self._cache_inchikey_notfound.add(inchikey)

            return biocycID

        except requests.RequestException as e:
            self.logger.error(f"Error fetching data from BioCyc: {e}")

            self._cache_inchikey_notfound.add(inchikey)
            return None

    def save_cache(self):
        self.cache.to_csv(
            self.__cache_file,
            index=False,
        )

    def __del__(self):
        self.save_cache()