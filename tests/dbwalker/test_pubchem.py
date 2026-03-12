import shutil
import tempfile
from unittest import TestCase
from unittest.mock import MagicMock, patch

import requests

from cobramod.dbwalker.dataclasses import (
    GenerellIdentifiers,
    Unavailable,
    Uncertain,
)
from cobramod.dbwalker.pubchem import PubChem


class TestGetPubChem(TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.pubchem = PubChem(cachedir=self.test_dir)

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_name(self):
        self.assertEqual("PubChem", self.pubchem.name)

    def test_AnnotationPrefix(self):
        self.assertEqual("pubchem.compound", self.pubchem.AnnotationPrefix)

    def test_getGenerellIdentifier(self):
        mock_session = MagicMock()
        mock_response = MagicMock()
        mock_response.status_code = 200

        mock_response.json.return_value = {
            "PropertyTable": {
                "Properties": [
                    {
                        "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                        "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                        "SMILES": "CCO",
                        "MolecularFormula": "C2H6O",
                    }
                ]
            }
        }
        mock_session.get.return_value = mock_response

        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            gID = self.pubchem.getGenerellIdentifier("702")

            self.assertEqual("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", gID.inchi)
            self.assertEqual("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", gID.inchi_key)
            self.assertEqual("CCO", gID.smiles)

            mock_session.get.assert_called_once()

        # Including prefix
        mock_session.reset_mock()
        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            gID = self.pubchem.getGenerellIdentifier("pubchem.compound:702")

            self.assertEqual("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", gID.inchi)
            self.assertEqual("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", gID.inchi_key)
            self.assertEqual("CCO", gID.smiles)
            mock_session.get.assert_called_once()

        # Cache
        mock_session.reset_mock()
        with patch.object(self.pubchem, "session", mock_session):
            gID_cached = self.pubchem.getGenerellIdentifier(702)

            self.assertEqual(
                "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", gID_cached.inchi
            )
            mock_session.get.assert_not_called()

        # Multiple potential GIDs
        mock_response.json.return_value = {
            "PropertyTable": {
                "Properties": [
                    {
                        "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                        "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                        "SMILES": "CCO",
                        "MolecularFormula": "C2H6O",
                    },
                    {
                        "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                        "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                        "SMILES": "CCO",
                        "MolecularFormula": "C2H6O",
                    },
                ]
            }
        }

        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            gID = self.pubchem.getGenerellIdentifier(703)

            self.assertEqual("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", gID.inchi)
            self.assertEqual("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", gID.inchi_key)

        # Test case 6: Edge case - mismatch in InChIKeys should raise ValueError
        mock_response.json.return_value = {
            "PropertyTable": {
                "Properties": [
                    {
                        "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                        "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                        "SMILES": "CCO",
                        "MolecularFormula": "C2H6O",
                    },
                    {
                        "InChI": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                        "InChIKey": "DIFFERENT-KEY-UHFFFAOYSA-N",
                        "SMILES": "CCO",
                        "MolecularFormula": "C2H6O",
                    },
                ]
            }
        }

        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            with self.assertRaises(ValueError):
                gid = self.pubchem.getGenerellIdentifier(704)

        # Test case 7: Edge case - None values in response
        mock_response.json.return_value = {
            "PropertyTable": {
                "Properties": [
                    {
                        "InChI": None,
                        "InChIKey": None,
                        "SMILES": "CCO",
                        "MolecularFormula": "C2H6O",
                    }
                ]
            }
        }

        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            gID = self.pubchem.getGenerellIdentifier(705)

            self.assertIsNone(gID.inchi)
            self.assertIsNone(gID.inchi_key)
            self.assertEqual("CCO", gID.smiles)

        # Test case 8: Live test
        self.pubchem._cache.clear_cache()
        gID = self.pubchem.getGenerellIdentifier(702)

        self.assertIsNotNone(gID.inchi)
        self.assertIsNotNone(gID.inchi_key)
        self.assertIsNotNone(gID.smiles)

    def test_getDBIdentifierFromSmiles(self):
        mock_session = MagicMock()
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "702\n"
        mock_session.post.return_value = mock_response

        # normal case
        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            pubchemID = self.pubchem.getDBIdentifierFromSmiles("CCO")
            self.assertEqual("702", pubchemID)
            mock_session.post.assert_called_once()

        # using cache from previous query
        mock_session.reset_mock()
        with patch.object(self.pubchem, "session", mock_session):
            pubchemID = self.pubchem.getDBIdentifierFromSmiles("CCO")
            self.assertEqual("702", pubchemID)
            mock_session.post.assert_not_called()

        # identifier not present
        mock_session.reset_mock()
        mock_response.text = "0\n"
        mock_session.post.return_value = mock_response

        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            pubchemID = self.pubchem.getDBIdentifierFromSmiles(
                "C1[C@H]([C@@H]([C@@H]([C@H](C(O1)OO)O)O)O)O"
            )
            self.assertEqual(Unavailable, pubchemID)

        # with GenerellIdentifiers object
        mock_session.reset_mock()
        mock_response.text = "702\n"
        mock_session.post.return_value = mock_response

        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            gID = GenerellIdentifiers(smiles="CCO")
            pubchemID = self.pubchem.getDBIdentifierFromSmiles(gID)
            self.assertEqual("702", pubchemID)

        # Live test
        self.pubchem._cache.clear_cache()
        pubchemID = self.pubchem.getDBIdentifierFromSmiles("CCO")
        self.assertIsNotNone(pubchemID)
        self.assertTrue(pubchemID.isdigit())

        pubchemID = self.pubchem.getDBIdentifierFromSmiles(
            "C1[C@H]([C@@H]([C@@H]([C@H](C(O1)OO)O)O)O)O"
        )
        self.assertEqual(Unavailable, pubchemID)

    def test_getDBIdentifierFromInchi(self):
        mock_session = MagicMock()
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "702\n"
        mock_session.post.return_value = mock_response

        # normal case
        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            pubchemID = self.pubchem.getDBIdentifierFromInchi(
                "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
            )
            self.assertEqual("702", pubchemID)
            mock_session.post.assert_called_once()

        # using cache from previous query
        mock_session.reset_mock()
        with patch.object(self.pubchem, "session", mock_session):
            pubchemID = self.pubchem.getDBIdentifierFromInchi(
                "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
            )
            self.assertEqual("702", pubchemID)
            mock_session.post.assert_not_called()

        # identifier not present (empty response)
        mock_session.reset_mock()
        mock_response.text = ""
        mock_session.post.return_value = mock_response

        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            pubchemID = self.pubchem.getDBIdentifierFromInchi(
                "InChI=1S/INVALID"
            )
            self.assertEqual(Unavailable, pubchemID)

        # with GenerellIdentifiers object
        mock_session.reset_mock()
        mock_response.text = "702\n"
        mock_session.post.return_value = mock_response

        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            gID = GenerellIdentifiers(inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
            pubchemID = self.pubchem.getDBIdentifierFromInchi(gID)
            self.assertEqual("702", pubchemID)

        # Live test
        self.pubchem._cache.clear_cache()
        pubchemID = self.pubchem.getDBIdentifierFromInchi(
            "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        )
        self.assertIsNotNone(pubchemID)
        self.assertTrue(pubchemID.isdigit())

    # def test_getDBIdentifierFromInchi(self):

    def test_getDBIdentifierFromInchiKey(self):
        mock_session = MagicMock()

        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "945872\n"
        mock_session.post.return_value = mock_response

        # normal case
        with patch.object(self.pubchem, "session", mock_session):
            self.pubchem._cache.clear_cache()
            pubchemID = self.pubchem.getDBIdentifierFromInchiKey(
                "NELACKPYIURSBY-MRXNPFEDSA-N"
            )
            self.assertEqual("945872", pubchemID)

        # using cache from previous query (should not call API)
        mock_response.text = "1264\n102212881\n"
        mock_session.post.return_value = mock_response

        with patch.object(self.pubchem, "session", mock_session):
            mock_session.reset_mock()
            pubchemID = self.pubchem.getDBIdentifierFromInchiKey(
                "NELACKPYIURSBY-MRXNPFEDSA-N"
            )
            self.assertEqual("945872", pubchemID)
            mock_session.post.assert_not_called()

            # multiple results returning uncertain
            self.pubchem._cache.clear_cache()
            pubchemID = self.pubchem.getDBIdentifierFromInchiKey(
                "PFRKDKQPQSBYQX-UHFFFAOYSA-O"
            )
            expected = Uncertain(
                possibilities=["1264", "102212881"], type="DBID"
            )
            self.assertEqual(expected, pubchemID)

        # live test
        self.pubchem._cache.clear_cache()
        pubchemID = self.pubchem.getDBIdentifierFromInchiKey(
            "JDSRHVWSAMTSSN-UHFFFAOYSA-N "
        )
        self.assertEqual("1426", pubchemID)

    def test_KEGGprovidedSIDs(self):
        # recreate pubchem object to ensure it does not reuse cache files
        test_dir = tempfile.mkdtemp()
        pubchem = PubChem(cachedir=test_dir)

        mock_session = MagicMock()
        mock_response = MagicMock()
        mock_response.status_code = 200

        # original source is file with one SID per line
        mock_response.text = "\n".join(str(i) for i in range(3303, 3365))
        mock_session.get.return_value = mock_response

        with patch.object(pubchem, "session", mock_session):
            pubchem._PubChem__keggSIDs = None
            keggSIDs = pubchem.KEGGprovidedSIDs

            mock_session.get.assert_called_with(
                url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceall/KEGG/sids/TXT",
                timeout=30,
            )

            expected = [str(i) for i in range(3303, 3365)]
            self.assertEqual(keggSIDs, expected)

    def test_getSIDsFromCIDs(self):
        mock_session = MagicMock()
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "12345\n67890\n11111\n"
        mock_session.get.return_value = mock_response

        with patch.object(self.pubchem, "session", mock_session):
            sids = self.pubchem.getSIDsFromCIDs("702")

            self.assertEqual(["12345", "67890", "11111"], sids)
            mock_session.get.assert_called_once()

    def test_getKeggSpecificSIDs(self):
        mock_session = MagicMock()
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "3303\n3304\n3305\n"
        mock_session.get.return_value = mock_response

        with patch.object(self.pubchem, "session", mock_session):
            # Force fresh download of KEGG SIDs
            self.pubchem._PubChem__keggSIDs = None
            self.pubchem._PubChem__cache_dir.mkdir(parents=True, exist_ok=True)

            # Input SIDs where only some overlap with KEGG list
            input_sids = ["3303", "99999", "3305", "00000"]
            kegg_sids = self.pubchem.getKeggSpecificSIDs(input_sids)

            self.assertEqual(["3303", "3305"], kegg_sids)

    def test_getCIDsFromSIDs(self):
        mock_session = MagicMock()
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "702\n5793\n"
        mock_session.get.return_value = mock_response

        # Normal case
        with patch.object(self.pubchem, "session", mock_session):
            cids = self.pubchem.getCIDsFromSIDs("12345")
            self.assertEqual(["702", "5793"], cids)
            mock_session.get.assert_called_once()

        # 404 not found returns Unavailable
        mock_session.reset_mock()
        mock_404_response = MagicMock()
        mock_404_response.status_code = 404
        mock_404_response.text = (
            "Status: 404\n"
            "Code: PUGREST.NotFound\n"
            "Message: No CIDs found for the given SID(s)"
        )
        http_error = requests.exceptions.HTTPError(response=mock_404_response)
        mock_404_response.raise_for_status.side_effect = http_error
        mock_session.get.return_value = mock_404_response

        with patch.object(self.pubchem, "session", mock_session):
            result = self.pubchem.getCIDsFromSIDs("99999")
            self.assertIs(Unavailable, result)

        # Other HTTP errors are re-raised
        mock_session.reset_mock()
        mock_500_response = MagicMock()
        mock_500_response.status_code = 500
        mock_500_response.text = "Server Error"
        http_error = requests.exceptions.HTTPError(response=mock_500_response)
        mock_500_response.raise_for_status.side_effect = http_error
        mock_session.get.return_value = mock_500_response

        with patch.object(self.pubchem, "session", mock_session):
            with self.assertRaises(requests.exceptions.HTTPError):
                self.pubchem.getCIDsFromSIDs("88888")

    def test_getChemicalFormular(self):
        mock_session = MagicMock()
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "C2H6O\n"
        mock_session.get.return_value = mock_response

        with patch.object(self.pubchem, "session", mock_session):
            formula = self.pubchem.getChemicalFormular("702")
            self.assertEqual("C2H6O", formula)
            mock_session.get.assert_called_once()

    def test_save_cache(self):
        with patch.object(self.pubchem._cache, "save_cache") as mock_save:
            self.pubchem.save_cache()
            mock_save.assert_called_once()

    def test_getDBIdentifier(self):
        # All three identifiers agree on the same DB ID
        with patch.object(self.pubchem, "getDBIdentifierFromSmiles", return_value="702"), \
            patch.object(self.pubchem, "getDBIdentifierFromInchi", return_value="702"), \
            patch.object(self.pubchem, "getDBIdentifierFromInchiKey", return_value="702"):

            gid = GenerellIdentifiers(
                smiles="CCO",
                inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                inchi_key="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            )
            result = self.pubchem.getDBIdentifier(gid)
            self.assertEqual("702", result)

        # Mismatch between identifiers returns Unavailable
        with patch.object(self.pubchem, "getDBIdentifierFromSmiles", return_value="702"), \
            patch.object(self.pubchem, "getDBIdentifierFromInchi", return_value="999"), \
            patch.object(self.pubchem, "getDBIdentifierFromInchiKey", return_value="702"):

            gid = GenerellIdentifiers(
                smiles="CCO",
                inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                inchi_key="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            )
            result = self.pubchem.getDBIdentifier(gid)
            self.assertIs(Unavailable, result)

        # No identifiers at all returns Unavailable
        gid = GenerellIdentifiers()
        result = self.pubchem.getDBIdentifier(gid)
        self.assertIs(Unavailable, result)

        # Some identifiers Unavailable, others agree
        with patch.object(self.pubchem, "getDBIdentifierFromSmiles", return_value="702"), \
            patch.object(self.pubchem, "getDBIdentifierFromInchi", return_value=Unavailable), \
            patch.object(self.pubchem, "getDBIdentifierFromInchiKey", return_value="702"):

            gid = GenerellIdentifiers(
                smiles="CCO",
                inchi="InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                inchi_key="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            )
            result = self.pubchem.getDBIdentifier(gid)
            self.assertEqual("702", result)