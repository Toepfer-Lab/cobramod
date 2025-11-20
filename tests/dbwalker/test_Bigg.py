from unittest import TestCase
import pytest
from unittest.mock import patch, Mock
import requests

from cobramod.dbwalker.bigg import Bigg
from cobramod.debug import change_to_debug

from cobramod.dbwalker.dataclasses import GenerellIdentifiers

import logging

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()],
)

# Ensure the logger in bigg.py propagates logs
logger = logging.getLogger("cobramod.DBWalker.Bigg")
logger.propagate = True
logger.setLevel(logging.DEBUG)


class TestBiggID2Identifier(TestCase):
    def setUp(self):
        self.connector = Bigg()

    def test_inchikey(self):
        """Test successful API call with all chemical identifiers present."""
        mock_response_data = {
            "name": "D-Glucose",
            "old_identifiers": ["glc__D"],
            "formulae": ["C6H12O6"],
            "compartments_in_models": [],
            "bigg_id": "glc__D",
            "database_links": {
                "InChI Key": [
                    {
                        "id": "WQZGKKKJIJFFOK-GASJEMHNSA-N",
                        "link": "https://identifiers.org/inchikey/WQZGKKKJIJFFOK-GASJEMHNSA-N",
                    }
                ],
                "BioCyc": [
                    {
                        "id": "META:Glucopyranose",
                        "link": "http://identifiers.org/biocyc/META:Glucopyranose",
                    }
                ],
            },
        }

        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result = self.connector.getGenerellIdentifier("glc__D")

            self.assertIsInstance(result, GenerellIdentifiers)
            # assert result.inchi == 'InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1'
            self.assertEqual(result.inchi_key, "WQZGKKKJIJFFOK-GASJEMHNSA-N")
            # assert result.smiles == 'C(C1C(C(C(C(O1)O)O)O)O)O'

    def test_request_exception_handling(self):
        """Test handling of requests exceptions."""
        with patch("requests.get") as mock_get:
            mock_get.side_effect = requests.RequestException("Connection error")

            result = self.connector.getGenerellIdentifier("test_id")

            assert isinstance(result, GenerellIdentifiers)
            assert result.inchi is None
            assert result.inchi_key is None
            assert result.smiles is None

    def test_empty_api_response(self):
        """Test handling of empty API response."""
        mock_response_data = {}

        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result = self.connector.getGenerellIdentifier("nonexistent_id")

            assert isinstance(result, GenerellIdentifiers)
            assert result.inchi is None
            assert result.inchi_key is None
            assert result.smiles is None

    def test_api_url_construction(self):
        """Test that correct API URL is constructed."""
        mock_response_data = {"inchi": "test_inchi"}

        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            self.connector.getGenerellIdentifier("atp")

            mock_get.assert_called_once_with(
                "http://bigg.ucsd.edu/api/v2/universal/metabolites/atp"
            )

    def test_live_api_call(self):
        """Test against the livebiggID2Identifier BiGG API (may fail if API is unavailable)."""
        result = self.connector.getGenerellIdentifier("glc__D")
        self.assertTrue(isinstance(result, GenerellIdentifiers))

        # assert result.inchi == 'InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1'
        self.assertEqual(result.inchi_key, "WQZGKKKJIJFFOK-GASJEMHNSA-N")
        # assert result.smiles == 'C(C1C(C(C(C(O1)O)O)O)O)O'
