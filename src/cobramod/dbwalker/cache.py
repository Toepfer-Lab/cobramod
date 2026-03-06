import json
import logging
from pathlib import Path
from typing import Dict, Set, Union

import pandas as pd
from notebook.terminal import initialize

from cobramod.dbwalker.dataclasses import GenerellIdentifiers, Unavailable


class MissmatchError(BaseException):
    pass


logger = logging.getLogger("cobramod.DBWalker.cache")
logger.propagate = True


class Cache:
    def __init__(self, cache_dir: Union[Path, str, None] = None):
        super().__init__()
        self.initialized = False

        self.id_dict: Dict[str, GenerellIdentifiers] = {}
        self.smiles_dict: Dict[str, set] = {}
        self.inchi_dict: Dict[str, set] = {}
        self.inchi_key_dict: Dict[str, set] = {}

        self._cache_smiles_not_found: Set[str] = set()
        self._cache_inchi_not_found: Set[str] = set()
        self._cache_inchikey_not_found: Set[str] = set()

        if isinstance(cache_dir, str):
            cache_dir = Path(cache_dir)

        self._added: int = 0
        self._cache_folder = cache_dir

        if cache_dir is not None and (cache_dir / "cache.xml").exists():
            self.load_cache()

        self.initialized = True

    def clear_cache(self):
        self.id_dict: Dict[str, GenerellIdentifiers] = {}
        self.smiles_dict: Dict[str, set] = {}
        self.inchi_dict: Dict[str, set] = {}
        self.inchi_key_dict: Dict[str, set] = {}

        self._cache_smiles_not_found: Set[str] = set()
        self._cache_inchi_not_found: Set[str] = set()
        self._cache_inchikey_not_found: Set[str] = set()

    def save_cache(self):
        self._cache_folder.mkdir(exist_ok=True, parents=True)

        if self._cache_folder is None:
            raise KeyError

        with open(self._cache_folder / "cache.xml", "w") as file:
            dict_format = {}
            for key, value in self.id_dict.items():
                dict_format[key] = value.to_dict()

            file.write(json.dumps(dict_format, indent=4))

        with open(self._cache_folder / "smilesNotFound.txt", "w") as file:
            file.writelines(
                line + "\n" for line in self._cache_smiles_not_found
            )

        with open(self._cache_folder / "inchiNotFound.txt", "w") as file:
            file.writelines(line + "\n" for line in self._cache_inchi_not_found)

        with open(self._cache_folder / "inchikeyNotFound.txt", "w") as file:
            file.writelines(
                line + "\n" for line in self._cache_inchikey_not_found
            )

    def load_cache(self):
        if self._cache_folder is None:
            raise KeyError

        with open(self._cache_folder / "cache.xml", "r") as file:
            load_data = json.load(file)

            for key, value in load_data.items():
                gID = GenerellIdentifiers.from_dict(value)

                self.addGenerellIdentifiers(gID, key)

        with open(self._cache_folder / "smilesNotFound.txt", "r") as file:
            for line in file:
                self._cache_smiles_not_found.add(line.strip())

        with open(self._cache_folder / "inchiNotFound.txt", "r") as file:
            for line in file:
                self._cache_inchi_not_found.add(line.strip())

        with open(self._cache_folder / "inchikeyNotFound.txt", "r") as file:
            for line in file:
                self._cache_inchikey_not_found.add(line.strip())

    def __added(self):
        if not self.initialized:
            return

        self._added += 1
        if self._added % 10 == 0:
            self.save_cache()

    def addSmiles(self, smiles, dbID):
        if not isinstance(smiles, str):
            raise ValueError

        self.__added()
        logger.debug(
            f"Adding entry (SMILES:{smiles}, DBID:{dbID}) to the cache."
        )
        if dbID is Unavailable:
            self._cache_smiles_not_found.add(smiles)
            return

        if smiles not in self.smiles_dict:
            self.smiles_dict[smiles] = set()

        self.smiles_dict[smiles].add(dbID)

        if dbID not in self.id_dict:
            self.id_dict[dbID] = GenerellIdentifiers(smiles=smiles)

        else:
            entry = self.id_dict[dbID]
            if entry.smiles != smiles and entry.smiles is not None:
                raise MissmatchError
            entry.smiles = smiles

    def addInchi(self, inchi, dbID):
        if not isinstance(inchi, str):
            raise ValueError

        self.__added()

        logger.debug(f"Adding entry (InChI:{inchi}, DBID:{dbID}) to the cache.")

        if dbID is Unavailable:
            self._cache_inchi_not_found.add(inchi)
            return

        if inchi not in self.inchi_dict:
            self.inchi_dict[inchi] = set()

        self.inchi_dict[inchi].add(dbID)

        if dbID not in self.id_dict:
            self.id_dict[dbID] = GenerellIdentifiers(inchi=inchi)

        else:
            entry = self.id_dict[dbID]
            if entry.inchi != inchi and entry.inchi is not None:
                raise MissmatchError(
                    f"Found {entry.inchi}, while trying to add {inchi}, for DB-ID {dbID}"
                )
            
            entry.inchi = inchi

    def addInchiKey(self, inchikey, dbID):
        if not isinstance(inchikey, str):
            raise ValueError

        self.__added()

        logger.debug(
            f"Adding entry (InChIKey:{inchikey}, DBID:{dbID}) to the cache."
        )

        if dbID is Unavailable:
            self._cache_inchikey_not_found.add(inchikey)
            return

        if inchikey not in self.inchi_key_dict:
            self.inchi_key_dict[inchikey] = set()

        self.inchi_key_dict[inchikey].add(dbID)

        if dbID not in self.id_dict:
            self.id_dict[dbID] = GenerellIdentifiers(inchi_key=inchikey)

        else:
            entry = self.id_dict[dbID]
            if entry.inchi_key != inchikey and entry.inchi_key is not None:
                raise MissmatchError
            entry.inchi_key = inchikey

    def addGenerellIdentifiers(self, gID: GenerellIdentifiers, dbID):
        self.__added()

        logger.debug(f"Adding for ID ({dbID}) GID ({gID}) to the cache.")

        if dbID not in self.id_dict:
            self.id_dict[dbID] = gID

        else:
            entry = self.id_dict[dbID]
            if not entry.weakEq(gID):
                raise MissmatchError

            entry += gID
            self.id_dict[dbID] = entry

        if isinstance(gID.smiles, str):
            if gID.smiles not in self.smiles_dict:
                self.smiles_dict[gID.smiles] = set()

            self.smiles_dict[gID.smiles].add(dbID)

        if isinstance(gID.inchi, str):
            if gID.inchi not in self.inchi_dict:
                self.inchi_dict[gID.inchi] = set()

            self.inchi_dict[gID.inchi].add(dbID)

        if isinstance(gID.inchi_key, str):
            if gID.inchi_key not in self.inchi_key_dict:
                self.inchi_key_dict[gID.inchi_key] = set()

            self.inchi_key_dict[gID.inchi_key].add(dbID)

    def getBySmiles(self, smiles):
        if smiles in self._cache_smiles_not_found:
            return Unavailable

        return self.smiles_dict.get(smiles, None)

    def getByInchi(self, inchi):
        if inchi in self._cache_inchi_not_found:
            return Unavailable

        return self.inchi_dict.get(inchi, None)

    def getByInchiKey(self, inchikey):
        if inchikey in self._cache_inchikey_not_found:
            return Unavailable

        return self.inchi_key_dict.get(inchikey, None)

    def getByID(self, dbIdentifier):
        if dbIdentifier in self.id_dict:
            return self.id_dict[dbIdentifier]

        else:
            return None
