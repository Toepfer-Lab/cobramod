import atexit
import logging
from pathlib import Path
from typing import List, Tuple, Union

import requests
from requests.adapters import HTTPAdapter, Retry

from cobramod import Settings
from cobramod.dbwalker.cache import Cache
from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import (
    GenerellIdentifiers,
    Unavailable,
    UnavailableType,
    Uncertain,
)


class PubChem(Database):
    def __init__(self, cachedir: Union[Path, str, None] = None):
        super().__init__()
        self.__settings = Settings()

        if cachedir:
            if isinstance(cachedir, str):
                cachedir = Path(cachedir)
            self.__cache_dir = cachedir
        else:
            self.__cache_dir = self.__settings.cacheDir / self.name

        self.__keggSIDs = None
        self._cache = Cache(cache_dir=self.__cache_dir)
        self.session = requests.Session()
        self.session.mount(
            "https://",
            HTTPAdapter(
                max_retries=Retry(
                    total=5,
                    backoff_factor=0.5,
                    allowed_methods=frozenset(["GET", "POST"]),
                )
            ),
        )

        self.logger = logging.getLogger("cobramod.DBWalker.PubChem")
        self.logger.propagate = True

        atexit.register(self.save_cache)

    @property
    def name(self) -> str:
        return "PubChem"

    @property
    def AnnotationPrefix(self) -> str:
        return "pubchem.compound"

    @property
    def KEGGprovidedSIDs(self) -> List[str]:
        if self.__keggSIDs is not None:
            return self.__keggSIDs

        file = self.__cache_dir / "KEGGprovidedSIDs.txt"

        if file.exists():
            with file.open() as f:
                keggSIDs = f.readlines()
                self.__keggSIDs = keggSIDs
                return keggSIDs

        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceall/KEGG/sids/TXT"

        self.__settings.limiter.try_acquire("pubchem")
        response = self.session.get(url=url, timeout=30)
        response.raise_for_status()

        keggSIDs = response.text.splitlines()
        self.__keggSIDs = keggSIDs

        with open(file, "w") as file_writer:
            data_to_write = "\n".join(keggSIDs)
            file_writer.write(data_to_write)

        return keggSIDs

    def getSIDsFromCIDs(self, cid: str) -> List[str]:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/sids/TXT"

        self.__settings.limiter.try_acquire("pubchem")
        response = self.session.get(url=url, timeout=30)
        response.raise_for_status()

        SIDs = response.text.splitlines()

        return SIDs

    def getKeggSpecificSIDs(self, cid) -> List[str]:
        all_sids = self.getSIDsFromCIDs(cid=cid)
        sidsFromKegg = self.KEGGprovidedSIDs

        keggSIDs = [x for x in all_sids if x in sidsFromKegg]

        return keggSIDs

    def getCIDsFromSIDs(self, sid: str) -> Union[List[str], UnavailableType]:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/cids/TXT"

        self.__settings.limiter.try_acquire("pubchem")
        response = self.session.get(url=url, timeout=30)

        not_found_text = """Status: 404
Code: PUGREST.NotFound
Message: No CIDs found for the given SID(s)"""

        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as err:
            if (
                err.response.status_code == 404
                and err.response.text.strip() == not_found_text.strip()
            ):
                return Unavailable
            else:
                raise err

        cids = response.text.splitlines()
        logging.info(f"Got CID ({cids}) for SID({sid})")

        return cids

    def getGenerellIdentifier(
        self, dbIdentifier: Union[str, int], **kwargs
    ) -> GenerellIdentifiers:
        if (
            isinstance(dbIdentifier, str)
            and "pubchem.compound:" in dbIdentifier
        ):
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
        response = self.session.get(url=url, timeout=30)
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

        self._cache.addGenerellIdentifiers(
            generell_identifiers, dbID=dbIdentifier
        )
        return generell_identifiers

    def getDBIdentifierFromSmiles(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(smiles, GenerellIdentifiers):
            assert smiles.smiles is not None
            smiles = smiles.smiles

        cached = self._cache.getBySmiles(smiles=smiles)

        if cached is Unavailable:
            return Unavailable

        elif isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/TXT"
        )
        data = {
            "smiles": smiles,
        }

        self.logger.debug(
            f"Getting identifiers from PubChem using the following url: {url} & Post data: {data}"
        )

        self.__settings.limiter.try_acquire("pubchem")
        response = self.session.post(url, files=data, timeout=30)

        response.raise_for_status()

        value = response.text.rstrip()
        if (
            value == "0"
        ):  # PubChem returns 0 for SMILES with a valid structure but no entry
            return Unavailable

        self._cache.addSmiles(smiles=smiles, dbID=value)
        return value

    def getDBIdentifierFromInchi(
        self, inchi: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(inchi, GenerellIdentifiers):
            assert inchi.inchi is not None
            inchi = inchi.inchi

        cached = self._cache.getByInchi(inchi=inchi)
        if cached is Unavailable:
            return Unavailable

        elif isinstance(cached, set):
            if len(cached) == 1:
                return list(cached)[0]

            else:
                return cached

        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/TXT"
        )

        data = {
            "inchi": inchi,
        }

        self.logger.debug(
            f"Getting identifiers from PubChem using the following url: {url} & Post data: {data}"
        )

        self.__settings.limiter.try_acquire("pubchem")
        response = self.session.post(url, files=data, timeout=30)

        response.raise_for_status()

        value = response.text.rstrip()

        if value == 0 or value == "0":
            value = Unavailable

        self._cache.addInchi(inchi=inchi, dbID=value)
        return value

    def getDBIdentifierFromInchiKey(
        self, inchikey: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(inchikey, GenerellIdentifiers):
            assert inchikey.inchi_key is not None
            inchikey = inchikey.inchi_key

        cached = self._cache.getByInchiKey(inchikey=inchikey)
        if cached is Unavailable:
            return Unavailable

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
        response = self.session.post(url, files=data, timeout=30)

        response.raise_for_status()

        value = response.text.rstrip()

        if "\n" in value:
            return Uncertain(possibilities=value.splitlines(), type="DBID")

        if value == 0:
            value = Unavailable

        self._cache.addInchiKey(inchikey=inchikey, dbID=value)
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
        response = self.session.get(url, timeout=30)
        response.raise_for_status()

        value = response.text.rstrip()

        return value

    def save_cache(self):
        self._cache.save_cache()
