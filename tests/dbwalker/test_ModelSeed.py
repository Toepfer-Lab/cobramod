from unittest import TestCase

import pytest
import requests
from cobramod.dbwalker.modelSeed import (
    get_compound_info_by_modelseed_id,
    smiles2ModelSeed,
)
from cobramod.dbwalker.dataclasses import GenerellIdentifiers


class TestGetCompoundInfoByModelSeedId(TestCase):
    def test_live_api_glucose(self):
        """Test with glucose (cpd00027) - a well-known compound."""
        result = get_compound_info_by_modelseed_id("cpd00027")

        self.assertIsInstance(result, GenerellIdentifiers)
        # Glucose should have chemical identifiers
        has_identifier = any([result.smiles, result.inchi, result.inchi_key])
        self.assertTrue(
            has_identifier, "Glucose should have chemical identifiers"
        )
        expected = GenerellIdentifiers(
            smiles="OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
            inchi=None,
            inchi_key="WQZGKKKJIJFFOK-GASJEMHNSA-N",
        )

        self.assertEqual(result.smiles, expected.smiles)
        self.assertEqual(result.inchi, expected.inchi)
        self.assertEqual(result.inchi_key, expected.inchi_key)

    def test_live_api_empty_string(self):
        """Test with empty string as compound ID."""
        result = get_compound_info_by_modelseed_id("")

        self.assertIsInstance(result, GenerellIdentifiers)
        # Should return empty GenerellIdentifiers for empty string
        self.assertIsNone(result.smiles)
        self.assertIsNone(result.inchi)
        self.assertIsNone(result.inchi_key)


class TestSmiles2ModelSeed(TestCase):
    def test_live_api_water_thioglycol(self):
        """Test conversion of water SMILES to ModelSEED ID."""
        result = smiles2ModelSeed("OCCS")

        if result is not None:
            self.assertIsInstance(result, str)
            self.assertGreater(len(result), 0)

        expected = "cpd00688"

        self.assertEqual(result, expected)
