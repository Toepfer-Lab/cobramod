from unittest import TestCase
from unittest.mock import patch, Mock
import requests
from requests import Session

from cobramod.dbwalker.dataclasses import GenerellIdentifiers

from cobramod.dbwalker.BioCyc import BioCyc

import logging

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)

logger = logging.getLogger("cobramod")
logger.propagate = True
logger.setLevel(logging.DEBUG)


class TestGetCompoundInfoByBiocycId(TestCase):
    def setUp(self):
        self.biocyc = BioCyc()

    def test_get_BioCycIDviaSmiles(self):

        #simple smiles but not unique
        biocycID = self.biocyc.getDBIdentifierFromSmiles("O")
        self.assertIsNone(biocycID)

        # identifier not present
        smiles = "CC(=O)OCOC(=O)C"
        biocycID = self.biocyc.getDBIdentifierFromSmiles(smiles)
        self.assertIsNone(biocycID)

        #complicated id
        smiles = "C(O)[C@H]2(O[C@@H](O)[C@H](O[C@@H]1([C@H](O)[C@H]([C@H](O)[C@@H](CO)O1)O))[C@@H](O)[C@H](O)2)"
        biocycID = self.biocyc.getDBIdentifierFromSmiles(smiles)
        self.assertEqual(first=biocycID, second="BETA-KOJIBIOSE")

