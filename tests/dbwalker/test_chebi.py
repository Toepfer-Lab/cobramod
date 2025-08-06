from unittest import TestCase
import pytest
from unittest.mock import patch, Mock
import requests

from cobramod.dbwalker.chebi import get_chebi_compound_info, getChebiID
from cobramod.debug import change_to_debug

from cobramod.dbwalker.bigg import biggID2Identifier

import logging

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()]
)

# Ensure the logger in bigg.py propagates logs
logger = logging.getLogger("cobramod")
logger.propagate = True
logger.setLevel(logging.DEBUG)


class TestIdentifier2ChebiID(TestCase):

    def test_live_api_call(self):
        """Test against the livebiggID2Identifier BiGG API (may fail if API is unavailable)."""
        result = getChebiID(identifier='Cn1cnc2n(C)c(=O)n(C)c(=O)c12', identifier_type="SMILES")
        self.assertEqual("CHEBI:27732", result)

class TestChebiID2Identifier(TestCase):
    def test_live_api_call(self):
        """Test against the livebiggID2Identifier BiGG API (may fail if API is unavailable)."""

        result = get_chebi_compound_info(chebi_id="CHEBI:27732")
        self.assertEqual("1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3", result.inchi)
        self.assertEqual("RYYVLZVUVIJVGH-UHFFFAOYSA-N", result.inchi_key)
        self.assertEqual("Cn1cnc2n(C)c(=O)n(C)c(=O)c12", result.smiles)

        result = get_chebi_compound_info(chebi_id="27732")
        self.assertEqual("1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3", result.inchi)
        self.assertEqual("RYYVLZVUVIJVGH-UHFFFAOYSA-N", result.inchi_key)
        self.assertEqual("Cn1cnc2n(C)c(=O)n(C)c(=O)c12", result.smiles)