import logging
from typing import Union, Tuple, List

import requests

from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import GenerellIdentifiers

logger = logging.getLogger("cobramod.DBWalker.PubChem")
logger.propagate = True


class PubChem(Database):
    def getGenerellIdentifier(self, dbIdentifier: str, **kwargs) -> GenerellIdentifiers:

        if "pubchem.compound:" in dbIdentifier:
            dbIdentifier = dbIdentifier.replace("pubchem.compound:", "")

        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{dbIdentifier}/property/MolecularFormula,InChIKey,InChI,SMILES/json"

        logger.debug(f"Getting identifiers from PubChem using the following url:\n{url}")

        response = requests.get(url=url)
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

        return generell_identifiers

        pass

    def getDBIdentifierFromSmiles(
        self, smiles: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(smiles, GenerellIdentifiers):
            assert smiles.smiles is not None
            smiles = smiles.smiles

        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/TXT"

        response = requests.get(url)
        response.raise_for_status()

        value = response.text.rstrip()

        return value

    def getDBIdentifierFromInchi(
        self, inchi: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(inchi, GenerellIdentifiers):
            assert inchi.inchi is not None
            inchi = inchi.inchi

        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{inchi}/cids/TXT"

        response = requests.get(url)
        response.raise_for_status()

        value = response.text.rstrip()

        return value

    def getDBIdentifierFromInchiKey(
        self, inchikey: Union[str, GenerellIdentifiers]
    ) -> str:
        if isinstance(inchikey, GenerellIdentifiers):
            assert inchikey.inchi_key is not None
            inchikey = inchikey.inchi_key

        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{inchikey}/cids/TXT"

        response = requests.get(url)
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

        response = requests.get(url)
        response.raise_for_status()

        value = response.text.rstrip()

        return value
