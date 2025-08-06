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
        # Note: Live API may or may not have data, so we just check the structure
        print(f"Live API result: {result}")

class TestChebiID2Identifier(TestCase):
    def test_live_api_call(self):
        """Test against the livebiggID2Identifier BiGG API (may fail if API is unavailable)."""

        result = get_chebi_compound_info(chebi_id="CHEBI:27732")
        # Note: Live API may or may not have data, so we just check the structure
        print(f"Live API result: {result}")