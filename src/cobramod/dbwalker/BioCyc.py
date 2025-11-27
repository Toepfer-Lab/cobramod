import urllib

import requests
from typing import Optional, Union, Tuple, overload, List

import xml.etree.ElementTree as ET
import logging

from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import GenerellIdentifiers
import cobramod


class BioCyc(Database):
    logger = logging.getLogger("cobramod.DBWalker.BioCyc")
    logger.propagate = True

    # Ensure console output
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

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

        url = "https://websvc.biocyc.org/getxml"

        params = {"id": f"{db}:{biocyc_id}"}

        # BioCyc cuts the SMILES string short otherwise
        headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/xml"}

        BioCyc.logger.debug(
            f"Querying BioCyc for compound {biocyc_id} in database {db}"
        )
        self.logger.debug(f"Request URL: {url}")
        self.logger.debug(f"Request parameters: {params}")

        result = GenerellIdentifiers()
        try:
            cSettings = cobramod.Settings()
            session = cSettings._biocycSession
            response = session.get(
                url, params=params, headers=headers, timeout=30
            )
            response.raise_for_status()

            if cSettings.autoOpenCloseBioCycSession:
                cSettings._closeBiocycSession()

            self.logger.debug(f"Response status code: {response.status_code}")

            # Check if BioCyc is asking for account creation
            if (
                response.text
                and isinstance(response.text, str)
                and "<title>Create Account</title>" in response.text
            ):
                self.logger.warning("BioCyc is requesting account creation.")
                if not cSettings.BioCycLoggedIn:
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

        return result

    def getDBIdentifierFromSmiles(
        self,
        smiles: Union[str, GenerellIdentifiers],
        BioCycSubDB: Optional[str] = None,
    ) -> Optional[str]:
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

        try:
            self.logger.info(f"Looking up BioCyc ID for SMILES: {smiles}")

            cSettings = cobramod.Settings()
            session = cSettings._biocycSession

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

        url = f"https://websvc.biocyc.org/{BioCycSubDB}/inchi-search?inchi={inchi}&exact=T&fmt=json"

        self.logger.info(f"Requesting BioCyc for InChI: {inchi}\nusing {url}")

        cSettings = cobramod.Settings()
        session = cSettings._biocycSession

        response = session.get(url, timeout=30)

        try:
            response.raise_for_status()

            data = response.json()

            biocycID = data["RESULTS"][0]["OBJECT-ID"]
            return biocycID

        except requests.RequestException as e:
            self.logger.error(f"Error fetching data from BioCyc: {e}")
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

        url = f"https://biocyc.org/{BioCycSubDB}/search-query?type=COMPOUND&inchikey={inchikey}&exact=T&fmt=json"
        self.logger.info(
            f"Requesting BioCyc for InChIKey: {inchikey}\nusing: {url}"
        )

        cSettings = cobramod.Settings()
        session = cSettings._biocycSession

        response = session.get(url, timeout=30)
        print(response.json())

        try:
            response.raise_for_status()

            data = response.json()
            print(data)

            biocycID = data["RESULTS"][0]["OBJECT-ID"]
            return biocycID

        except requests.RequestException as e:
            self.logger.error(f"Error fetching data from BioCyc: {e}")

            return None

    def getDBIdentifier(self, identifier: GenerellIdentifiers) -> Optional[str]:
        smilesBasedID = None
        inchiBasedID = None
        inchikeyBasedID = None

        self.logger.debug(
            "Querying all available GenerellIdentifier, to check whether they point to the same database ID."
        )

        if identifier.smiles is not None:
            smilesBasedID = self.getDBIdentifierFromSmiles(identifier)

        if identifier.inchi is not None:
            inchiBasedID = self.getDBIdentifierFromInchi(identifier)

        if identifier.inchi_key is not None:
            inchikeyBasedID = self.getDBIdentifierFromInchiKey(identifier)

        self.logger.debug(
            "Queried all available GenerellIdentifier. Checking if they point to the same database ID."
        )

        missmatch = False

        # ugly implementation but easiest to understand and backwards compatible
        if smilesBasedID is None:
            if inchiBasedID is None:
                if inchikeyBasedID is None:
                    return None

                return inchikeyBasedID

            else:
                # inchibasedID is defined
                if inchiBasedID == inchikeyBasedID:
                    return inchiBasedID

                else:
                    missmatch = True

        else:
            # smilesbasedOD is defined
            if inchiBasedID is None:
                if inchikeyBasedID is None:
                    return smilesBasedID
                else:
                    if smilesBasedID == inchikeyBasedID:
                        return smilesBasedID
            else:
                if inchikeyBasedID is None:
                    if smilesBasedID == inchiBasedID:
                        return smilesBasedID
                    else:
                        missmatch = True
                else:
                    if smilesBasedID == inchiBasedID == inchikeyBasedID:
                        return smilesBasedID
                    else:
                        missmatch = True

        if missmatch:
            self.logger.error(
                "Generell Identifier for supposedly the same object result in different DB IDs."
                f"\n InChi: {identifier.inchi} -> DB ID {inchiBasedID}"
                f"\n InChiKey: {identifier.inchi_key} -> DB ID {inchikeyBasedID}"
                f"\n Smiles: {smilesBasedID} -> DB ID {smilesBasedID}"
            )
        else:
            self.logger.debug(
                "All available Identifier point towards the same DB ID:"
                f"\n InChi: {identifier.inchi} -> DB ID {inchiBasedID}"
                f"\n InChiKey: {identifier.inchi_key} -> DB ID {inchikeyBasedID}"
                f"\n Smiles: {smilesBasedID} -> DB ID {smilesBasedID}"
            )

        # ToDo aise error due to uncertainty
        return None
