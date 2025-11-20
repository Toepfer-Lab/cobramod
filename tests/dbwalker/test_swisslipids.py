import logging
from unittest import TestCase

import pytest
from unittest.mock import patch, Mock
import requests
from cobramod.dbwalker.swisslipids import (
    get_swisslipids_id,
    get_swisslipids_data,
    swisslipids_search_by_name,
)
from cobramod.dbwalker.dataclasses import GenerellIdentifiers

logger = logging.getLogger("cobramod")
logger.propagate = True
logger.setLevel(logging.DEBUG)


class TestGetSwissLipidsId(TestCase):
    def test_get_id_from_inchikey_success(self):
        """Test successful SwissLipids ID retrieval from InChI Key."""
        mock_response_data = [
            {
                "entity_id": "SLM:000000651",
                "entity_name": "1-hexadecanoyl-2-(9Z-octadecenoyl)-sn-glycero-3-phosphocholine",
                "external_id": "",
                "entity_type": "metabolite",
                "nb_exp_annot": 0,
                "classification_level": "Isomeric subspecies",
                "stats": {
                    "entity_id": "SLM:000000651",
                    "nb_cheminfo": "0",
                    "nb_classification": "0",
                    "nb_metabolism": "52",
                    "nb_infered_metabolism": "0",
                    "nb_functions": "0",
                    "nb_go": "0",
                    "nb_interactions": "0",
                    "nb_locations": "0",
                    "nb_genetics": "0",
                    "nb_total": "52",
                },
            }
        ]

        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result = get_swisslipids_id(
                "WQZGKKKJIJFFOK-GASJEMHNSA-N", "inchikey"
            )

            self.assertEqual("SLM:000000651", result)
            mock_get.assert_called_once_with(
                "https://www.swisslipids.org/api/index.php/advancedSearch?InChIkey=WQZGKKKJIJFFOK-GASJEMHNSA-N",
                timeout=30,
            )

    def test_get_id_from_smiles_success(self):
        """Test successful SwissLipids ID retrieval from SMILES."""
        mock_response_data = [
            {
                "entity_id": "SLM:000597962",
                "entity_name": "propyl propanoate",
                "external_id": "",
                "entity_type": "metabolite",
                "nb_exp_annot": 0,
                "classification_level": "not classified",
                "stats": {
                    "entity_id": "SLM:000597962",
                    "nb_cheminfo": "1",
                    "nb_classification": "1",
                    "nb_metabolism": "1",
                    "nb_infered_metabolism": "0",
                    "nb_functions": "0",
                    "nb_go": "0",
                    "nb_interactions": "0",
                    "nb_locations": "0",
                    "nb_genetics": "0",
                    "nb_total": "1",
                },
            }
        ]

        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result = get_swisslipids_id("CCCOC(CC)=O", "smiles")

            self.assertEqual("SLM:000597962", result)
            mock_get.assert_called_once_with(
                "https://www.swisslipids.org/api/index.php/advancedSearch?SMILES=CCCOC\(CC\)%3DO",
                timeout=30,
            )

    def test_unsupported_identifier_type(self):
        """Test handling of unsupported identifier types."""
        result = get_swisslipids_id("test", "invalid_type")
        self.assertIsNone(result)

    def test_api_error_handling(self):
        """Test handling of API errors."""
        with patch("requests.get") as mock_get:
            mock_get.side_effect = requests.RequestException("API Error")

            result = get_swisslipids_id("TEST-KEY", "inchikey")
            self.assertIsNone(result)

    def test_non_existing_id(self):
        """Test handling of non-existing IDs."""
        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = {
                "status": 404,
                "message": "Sorry, no result found",
                "title": "ERROR",
            }
            mock_response.status_code = 404
            mock_response.raise_for_status.return_value = 404
            mock_get.return_value = mock_response

            with self.assertLogs(
                "cobramod.DBWalker.SwissLipids", level="ERROR"
            ) as cm:
                result = get_swisslipids_id("NONEXISTENT", "SMILES")

                self.assertIsNone(result)
                self.assertEqual(
                    cm.output,
                    [
                        "ERROR:cobramod.DBWalker.SwissLipids:No result available in database for: smiles=NONEXISTENT"
                    ],
                )

    def test_response_without_id_field(self):
        """Test handling of response without ID field."""
        mock_response_data = {"name": "Test without ID"}

        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result = get_swisslipids_id("TEST-KEY", "inchikey")
            assert result is None

    def test_live_api_integration(self):
        """Integration test against live SwissLipids API."""
        result = get_swisslipids_data("SLM:000000273")

        self.assertIsInstance(result, GenerellIdentifiers)
        self.assertEqual(
            result.inchi,
            "1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1",
        )
        self.assertEqual(
            result.smiles, "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
        )
        self.assertEqual(result.inchi_key, "WQZGKKKJIJFFOK-GASJEMHNSA-N")


class TestGetSwissLipidsData(TestCase):
    def test_get_data_success_all_identifiers(self):
        """Test successful data retrieval with all chemical identifiers."""
        mock_response_data = {
            "structures": {
                "structure": "",
                "inchikey": "WQZGKKKJIJFFOK-GASJEMHNSA-N",
                "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
                "inchi": "InChI=1S\/C6H12O6\/c7-1-2-3(8)4(9)5(10)6(11)12-2\/h2-11H,1H2\/t2-,3-,4+,5-,6?\/m1\/s1",
            },
            "chemical_data": {
                "formula": "C6H12O6",
                "mass": 180.15588,
                "charge": 0,
                "mz": [],
            },
            "entity_id": "SLM:000000273",
            "entity_name": "D-glucose",
            "entity_type": "metabolite",
            "creation_date": "2013-06-24 20:32:23",
            "last_update": "2013-08-15",
            "tax_id": "null",
            "taxon": "null",
            "nb_exp_annot": 0,
            "synonyms": [
                {"name": "D-glucose", "type": "name", "source": "SwissLipids"}
            ],
            "xrefs": [
                {"source": "Rhea", "url": "", "id": "CHEBI:4167"},
                {
                    "source": "HMDB",
                    "url": "http:\/\/www.hmdb.ca\/metabolites\/HMDB00122",
                    "id": "HMDB00122",
                },
                {
                    "source": "ChEBI",
                    "url": "http:\/\/www.ebi.ac.uk\/chebi\/searchId.do?chebiId=4167",
                    "id": "4167",
                },
            ],
        }

        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result = get_swisslipids_data("SLM:000000273")

            self.assertIsInstance(result, GenerellIdentifiers)
            self.assertEqual(
                "1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1",
                result.inchi,
            )
            self.assertEqual("WQZGKKKJIJFFOK-GASJEMHNSA-N", result.inchi_key)
            self.assertEqual(
                "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", result.smiles
            )
            mock_get.assert_called_once_with(
                "https://www.swisslipids.org/api/index.php/entity/SLM:000000273",
                timeout=30,
            )

    def test_get_data_partial_identifiers(self):
        """Test data retrieval with only some identifiers present."""
        mock_response_data = {
            "structures": {
                "structure": "",
                "inchikey": "WQZGKKKJIJFFOK-GASJEMHNSA-N",
            },
            "chemical_data": {
                "formula": "C6H12O6",
                "mass": 180.15588,
                "charge": 0,
                "mz": [],
            },
            "entity_id": "SLM:000000273",
            "entity_name": "D-glucose",
            "entity_type": "metabolite",
            "creation_date": "2013-06-24 20:32:23",
            "last_update": "2013-08-15",
            "tax_id": "null",
            "taxon": "null",
            "nb_exp_annot": 0,
            "synonyms": [
                {"name": "D-glucose", "type": "name", "source": "SwissLipids"}
            ],
            "xrefs": [
                {"source": "Rhea", "url": "", "id": "CHEBI:4167"},
                {
                    "source": "HMDB",
                    "url": "http:\/\/www.hmdb.ca\/metabolites\/HMDB00122",
                    "id": "HMDB00122",
                },
                {
                    "source": "ChEBI",
                    "url": "http:\/\/www.ebi.ac.uk\/chebi\/searchId.do?chebiId=4167",
                    "id": "4167",
                },
            ],
        }

        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result = get_swisslipids_data("SLM:000000004")

            self.assertIsInstance(result, GenerellIdentifiers)
            self.assertIsNone(result.inchi)
            self.assertEqual("WQZGKKKJIJFFOK-GASJEMHNSA-N", result.inchi_key)
            self.assertIsNone(result.smiles)
            mock_get.assert_called_once_with(
                "https://www.swisslipids.org/api/index.php/entity/SLM:000000004",
                timeout=30,
            )

    def test_get_data_api_error(self):
        """Test handling of API errors."""
        with patch("requests.get") as mock_get:
            mock_get.side_effect = requests.RequestException("API Error")

            with self.assertLogs(
                "cobramod.DBWalker.SwissLipids", level="ERROR"
            ) as cm:
                result = get_swisslipids_data("SLM:000000999")
                expected = GenerellIdentifiers()

                self.assertEqual(result, expected)
                self.assertIsNone(result.inchi)
                self.assertIsNone(result.inchi_key)
                self.assertIsNone(result.smiles)

                self.assertEqual(
                    cm.output,
                    [
                        "ERROR:cobramod.DBWalker.SwissLipids:Error fetching data from SwissLipids: API Error"
                    ],
                )

    def test_get_data_empty_response(self):
        """Test handling of empty API response."""
        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = {}
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result = get_swisslipids_data("SLM:000000000")

            assert isinstance(result, GenerellIdentifiers)
            assert result.inchi is None
            assert result.inchi_key is None
            assert result.smiles is None


class TestSwissLipidsSearchByName:
    def test_search_by_name_success(self):
        """Test successful search by compound name."""
        mock_response_data = [
            {
                "id": "SLM:000000001",
                "name": "Palmitic acid",
                "synonyms": ["Hexadecanoic acid"],
            }
        ]

        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result = swisslipids_search_by_name("Palmitic acid")

            assert result == "SLM:000000001"
            mock_get.assert_called_once_with(
                "https://www.swisslipids.org/api/entity/search",
                params={"query": "Palmitic acid", "format": "json"},
                timeout=30,
            )

    def test_search_by_name_alternative_id_field(self):
        """Test search with alternative ID field name."""
        mock_response_data = [
            {"swisslipids_id": "SLM:000000002", "name": "Oleic acid"}
        ]

        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result = swisslipids_search_by_name("Oleic acid")

            assert result == "SLM:000000002"

    def test_search_by_name_no_results(self):
        """Test search with no results."""
        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = []
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result = swisslipids_search_by_name("Nonexistent compound")

            assert result is None

    def test_search_by_name_api_error(self):
        """Test handling of API errors during search."""
        with patch("requests.get") as mock_get:
            mock_get.side_effect = requests.RequestException("Network error")

            result = swisslipids_search_by_name("Test compound")

            assert result is None

    def test_search_by_name_malformed_response(self):
        """Test handling of malformed API response."""
        mock_response_data = [
            {
                "name": "Test compound"
                # Missing ID field
            }
        ]

        with patch("requests.get") as mock_get:
            mock_response = Mock()
            mock_response.json.return_value = mock_response_data
            mock_response.raise_for_status.return_value = None
            mock_get.return_value = mock_response

            result = swisslipids_search_by_name("Test compound")

            assert result is None

    def test_live_api_integration(self):
        """Integration test against live SwissLipids API."""
        # Test with a well-known lipid (may fail if API is unavailable)
        result = swisslipids_search_by_name("palmitic acid")

        # Just check that we get a reasonable response format if API is available
        if result is not None:
            assert isinstance(result, str)
            assert len(result) > 0
            # SwissLipids IDs typically start with 'SLM:'
            assert "SLM:" in result or result.startswith("SL")

        print(f"Live API search result for 'palmitic acid': {result}")
