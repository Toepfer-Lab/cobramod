from typing import List
from unittest import TestCase
from unittest.mock import patch, MagicMock
import pandas as pd
from rdkit import Chem

from cobramod.dbwalker.cache import Cache
from cobramod.dbwalker.dataclasses import GenerellIdentifiers


class TestCache(TestCase):

    def test_addSmiles(self):
        cache = Cache()

        cache.addSmiles("[H]O[H]", "CHEBI:15377")

        self.assertEqual(cache.getBySmiles("[H]O[H]"), set(["CHEBI:15377"]))

        expected = GenerellIdentifiers(smiles="[H]O[H]")
        self.assertEqual(cache.getByID("CHEBI:15377"),expected)


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