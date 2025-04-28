import tempfile
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

from cobra import Metabolite

from cobramod.core.crossreferences import (
    inchikey2pubchem_cid,
    add_crossreferences,
)


class TestCrossReferencesPubChem(TestCase):
    @patch("requests.get")
    def test_inchikey2pubchem_cid(self, mocked_post):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mocked_post.return_value.status_code = 200
            mocked_post.return_value.text = "154"

            value = inchikey2pubchem_cid("Test_ID", directory)
            self.assertEqual("154", value)

            # check that the local cache is used
            mocked_post.return_value.text = "155"
            value = inchikey2pubchem_cid("Test_ID", directory)
            self.assertEqual(value, "154")

            value = inchikey2pubchem_cid(["Test_ID2", "Test_ID3"], directory)
            self.assertEqual(["155", "155"], value)

            # test with more than one cid
            # PubChem may return multiple with one per line
            mocked_post.return_value.text = "154\n155\n"

            Path.unlink(directory / "XRef" / str("pubchem" + ".feather"))
            value = inchikey2pubchem_cid("Test_ID", directory)
            self.assertEqual(["154", "155"], value)

    @patch("requests.get")
    @patch("requests.post")
    def test_pubchem(self, mocked_post, mock_get):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)

            mocked_post.return_value.status_code = 200
            mocked_post.return_value.json.return_value = {}

            mock_get.return_value.status_code = 200
            mock_get.return_value.text = "9547068\n46891690\n"

            metabolite = Metabolite()
            metabolite.annotation = {"inchikey": "BBYWOYAFBUOUFP-JOCHJYFZSA-N"}
            add_crossreferences(metabolite, directory)

            expected_annotations = {
                "inchikey": "BBYWOYAFBUOUFP-JOCHJYFZSA-N",
                "pubchem.compound": ["9547068", "46891690"],
            }

            # check keys
            self.assertCountEqual(expected_annotations, metabolite.annotation)

            for key, value in expected_annotations.items():
                self.assertCountEqual(value, metabolite.annotation[key])
