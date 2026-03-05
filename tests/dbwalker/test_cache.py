import tempfile
from pathlib import Path
from typing import List
from unittest import TestCase
from unittest.mock import MagicMock, patch

import pandas as pd
from rdkit import Chem

from cobramod.dbwalker.cache import Cache
from cobramod.dbwalker.dataclasses import GenerellIdentifiers, Unavailable


class TestCache(TestCase):

    def test_addSmiles(self):
        cache = Cache()

        cache.addSmiles("[H]O[H]", "CHEBI:15377")

        self.assertEqual(cache.getBySmiles("[H]O[H]"), set(["CHEBI:15377"]))

        expected = GenerellIdentifiers(smiles="[H]O[H]")
        self.assertEqual(cache.getByID("CHEBI:15377"),expected)
        self.assertEqual(len(cache.smiles_dict), 1)
        self.assertEqual(len(cache.id_dict), 1)

        self.assertEqual(len(cache.inchi_dict), 0)
        self.assertEqual(len(cache.inchi_key_dict), 0)
        self.assertEqual(len(cache._cache_smiles_not_found), 0)
        self.assertEqual(len(cache._cache_inchi_not_found), 0)
        self.assertEqual(len(cache._cache_inchikey_not_found), 0)

        cache.addSmiles("[H]O[H]", Unavailable)

        hit = cache.getBySmiles("[H]O[H]")
        self.assertIsInstance(hit, Unavailable)

    def test_addSmilesNotFound(self):
        cache = Cache()
        cache.addSmiles("[H]", "CHEBI:15377")

    def test_addInchi(self):
        cache = Cache()

        cache.addInchi("InChI=1S/H2O/h1H2", "CHEBI:15377")

        self.assertEqual(cache.getByInchi("InChI=1S/H2O/h1H2"), set(["CHEBI:15377"]))

        expected = GenerellIdentifiers(inchi="InChI=1S/H2O/h1H2")
        self.assertEqual(cache.getByID("CHEBI:15377"),expected)


    def test_addInchiKey(self):
        cache = Cache()

        cache.addInchiKey("XLYOFNOQVPJJNP-UHFFFAOYSA-N", "CHEBI:15377")

        self.assertEqual(cache.getByInchiKey("XLYOFNOQVPJJNP-UHFFFAOYSA-N"), set(["CHEBI:15377"]))

        expected = GenerellIdentifiers(inchi_key="XLYOFNOQVPJJNP-UHFFFAOYSA-N")
        self.assertEqual(cache.getByID("CHEBI:15377"),expected)

    def test_addGenerellIdentifiers(self):
        cache = Cache()

        generellIdentifiers = GenerellIdentifiers(
            smiles="[H]O[H]",
            inchi="InChI=1S/H2O/h1H2",
            inchi_key="XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        )
        cache.addGenerellIdentifiers(generellIdentifiers,"CHEBI:15377")

        self.assertEqual(cache.getByID("CHEBI:15377"), generellIdentifiers)


    def test_save_cache(self):
        cache = Cache()


        cache.addInchiKey("XLYOFNOQVPJJNP-UHFFFAOYSA-N", "CHEBI:15377")
        cache.addInchiKey(
            "InChI=1S/H2O/h1H2", Unavailable)

        #with tempfile.TemporaryDirectory() as tmpdirname:
        cache._cache_folder = Path("./") #= tmpdirname
        cache.save_cache()

        cache2 = Cache()
        cache2._cache_folder = Path("./") #= tmpdirname

        print(cache.id_dict)
        print(cache.inchi_key_dict)
        print(cache._cache_inchikey_not_found)

        cache2.load_cache()
        print(cache.id_dict)
        print(cache.inchi_key_dict)
        print(cache._cache_inchikey_not_found)

