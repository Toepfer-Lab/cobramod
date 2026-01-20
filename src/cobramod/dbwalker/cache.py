from typing import Union, Dict

import pandas as pd
from cobramod.dbwalker.dataclasses import GenerellIdentifiers


class MissmatchError:
    pass


class Cache:
    def __init__(self):
        super().__init__()

        self.id_dict: Dict[str, GenerellIdentifiers] = {}
        self.smiles_dict: Dict[str, set] = {}
        self.inchi_dict: Dict[str, set] = {}
        self.inchi_key_dict: Dict[str, set] = {}

        self._cache_smiles_not_found = set()
        self._cache_inchi_not_found = set()
        self._cache_inchikey_not_found = set()


    def addSmiles(self, smiles, dbID):

        if smiles not in self.smiles_dict:
            self.smiles_dict[smiles] = set()

        self.smiles_dict[smiles].add(dbID)

        if dbID not in self.id_dict:
            self.id_dict[dbID] = GenerellIdentifiers(smiles = smiles)

        else:
            entry = self.id_dict[dbID]
            if entry.smiles != smiles:
                raise MissmatchError

    def addInchi(self, inchi, dbID):

        if inchi not in self.inchi_dict:
            self.inchi_dict[inchi] = set()

        self.inchi_dict[inchi].add(dbID)

        if dbID not in self.id_dict:
            self.id_dict[dbID] = GenerellIdentifiers(inchi = inchi)

        else:
            entry = self.id_dict[dbID]
            if entry.inchi != inchi:
                raise MissmatchError

    def addInchiKey(self, inchikey, dbID):

        if inchikey not in self.inchi_key_dict:
            self.inchi_key_dict[inchikey] = set()

        self.inchi_key_dict[inchikey].add(dbID)

        if dbID not in self.id_dict:
            self.id_dict[dbID] = GenerellIdentifiers(inchi_key = inchikey)

        else:
            entry = self.id_dict[dbID]
            if entry.inchi_key != inchikey:
                raise MissmatchError

    def addGenerellIdentifiers(self, gID:GenerellIdentifiers, dbID):

        if dbID not in self.id_dict:
            self.id_dict[dbID] = gID

        else:
            entry = self.id_dict[dbID]
            if not entry.weakEq(gID):
                raise MissmatchError

            entry += gID
            self.id_dict[dbID] = entry

        if gID.smiles is not None:
            if gID.smiles not in self.smiles_dict:
                self.smiles_dict[gID.smiles] = set()

            self.smiles_dict[gID.smiles].add(dbID)

        if gID.inchi is not None:
            if gID.inchi not in self.inchi_dict:
                self.inchi_dict[gID.inchi] = set()

            self.inchi_dict[gID.inchi].add(dbID)

        if gID.inchi_key is not None:
            if gID.inchi_key not in self.inchi_key_dict:
                self.inchi_key_dict[gID.inchi_key] = set()

            self.inchi_key_dict[gID.inchi_key].add(dbID)


    def getBySmiles(self, smiles):
        return self.smiles_dict[smiles]

    def getByInchi(self, inchi):
        return self.inchi_dict[inchi]

    def getByInchiKey(self, inchikey):
        return self.inchi_key_dict[inchikey]

    def getByID(self, dbIdentifier):
        return self.id_dict[dbIdentifier]