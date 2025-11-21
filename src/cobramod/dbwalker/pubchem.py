from typing import Union, Tuple

import requests

from cobramod.dbwalker.DataBase import Database
from cobramod.dbwalker.dataclasses import GenerellIdentifiers


class PubChem(Database):
    def getGenerellIdentifier(self, dbIdentifier: str, **kwargs) -> GenerellIdentifiers:
        generell_identifiers = GenerellIdentifiers()

        for property_name in ["InChi", "InChIKey"]:
            url = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
                f"{dbIdentifier}/property/{property_name}/txt"
            )

            response = requests.get(url=url)
            response.raise_for_status()

            value = response.text.rstrip()

            # PubChem returns one value per line if multiple valid
            if "\n" in value:
                value = value.split("\n")

            if property_name == "InChI":
                generell_identifiers.inchi = value
            else:
                generell_identifiers.inchi_key = value

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
